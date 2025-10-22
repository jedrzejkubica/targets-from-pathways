import os
import sys
import logging
import argparse
import pathlib
import urllib.request
import zipfile
import shutil

import networkx

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)

# Known sources for Reactome files
REACTOME_MAPPING_URL = "https://download.reactome.org/94/Ensembl2Reactome_PE_All_Levels.txt"
REACTOME_INTERACTIONS_ZIP_URL = (
    "http://cpws.reactome.org/caBigR3WebApp2025/FIsInGene_04142025_with_annotations.txt.zip"
)
REACTOME_INTERACTIONS_TXT = "FIsInGene_04142025_with_annotations.txt"


def _ensure_reactome_mapping(path: pathlib.Path) -> None:
    """Download the Reactome Ensembl2Reactome mapping if missing and filename matches expected.

    Only auto-downloads when the basename is 'Ensembl2Reactome_PE_All_Levels.txt'.
    """
    if path.exists():
        return
    if path.name != "Ensembl2Reactome_PE_All_Levels.txt":
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Downloading Reactome mapping to %s", path)
    urllib.request.urlretrieve(REACTOME_MAPPING_URL, str(path))


def _ensure_reactome_interactions(path: pathlib.Path) -> None:
    """Download and unzip the Reactome interactions file if missing using known URL.

    Only auto-downloads when the basename is 'FIsInGene_04142025_with_annotations.txt'.
    """
    if path.exists():
        return
    if path.name != REACTOME_INTERACTIONS_TXT:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    zip_path = path.with_suffix(path.suffix + ".zip")
    logging.info("Downloading Reactome interactions zip to %s", zip_path)
    urllib.request.urlretrieve(REACTOME_INTERACTIONS_ZIP_URL, str(zip_path))
    logging.info("Unzipping interactions...")
    with zipfile.ZipFile(str(zip_path), "r") as zf:
        # Extract desired txt; fall back to extracting all and renaming
        members = zf.namelist()
        if REACTOME_INTERACTIONS_TXT in members:
            zf.extract(REACTOME_INTERACTIONS_TXT, path.parent)
        else:
            zf.extractall(path.parent)
    # Ensure the expected target path exists
    extracted = path.parent / REACTOME_INTERACTIONS_TXT
    if extracted.exists() and extracted != path:
        try:
            shutil.move(str(extracted), str(path))
        except Exception:
            pass


def get_disease_genes(pathway2genes, disease_pathways):
    """
    Find all genes that are on disease-specific pathways

    arguments:
    - pathway2genes
    - disease_pathways
    
    returns:
    - list of disease-specific genes
    """
    disease_genes = set()
    for path in disease_pathways:
        genes_on_path = pathway2genes[path]
        logger.info("Found %i genes in %s", len(genes_on_path), path)

        for gene in genes_on_path:
            disease_genes.add(gene)
    
    return(disease_genes)


def get_target_genes(gene2pathways, pathway2genes, target):
    """
    Find all genes that are on the same pathways as target

    arguments:
    - gene2pathways
    - pathway2genes
    - target
    
    returns:
    - list of target-specific genes
    """
    target_genes = set()

    all_paths = gene2pathways[target]
    for path in all_paths:
        genes_on_path = pathway2genes[path]
        for gene in genes_on_path:
            target_genes.add(gene)

    return(target_genes)


def find_overlap(list1, list2):
    """
    Find all genes that are both disease-specific and on the same pathways as target

    arguments:
    - list of disease-specific genes
    - list of target-specific genes
    """
    overlap = list(set(list1) & set(list2))
    return(overlap)


def calculate_scores(genes, gene2pathways, pathway2genes, disease_pathways, target):
    """
    For every gene found in the overlap between disease-specific and target-specific list:
    score = (#disease_paths_with_gene + #target_paths_with_gene) / (#disese_paths + #target_paths)

    arguments:
    - genes: list of genes from find_overlap
    - gene2pathways
    - 
    
    returns:
    - scores: dict, key=gene, value=score
    """
    scores = {}
    for gene in genes:
        # #disease_paths_with_gene
        gene_paths = gene2pathways[gene]
        disease_paths_with_gene = find_overlap(gene_paths, disease_pathways)

        # #target_paths_with_gene
        target_paths = gene2pathways[target]
        target_path_count = 0
        for path in target_paths:
            path_genes = pathway2genes[path]
            if gene in path_genes:
                target_path_count += 1

        scores[gene] = (len(disease_paths_with_gene) + target_path_count) / (len(disease_pathways) + len(target_paths))

    return(scores)


def main(pathway_mapping_file, interactions_file, target):

    logger.info("Parsing gene-to-pathway mapping file")
    (gene2pathways, pathway2genes) = data_parser.parse_pathway_mapping(pathway_mapping_file)

    logger.info("Computing GSEA")
    # Read list of enriched pathways from provided override or legacy default
    ids_path = os.environ.get("GSEA_IDS_FILE_OVERRIDE", "gsea/gsea_output_ids/OT-EFO_0004248_gsea_ids.tsv")
    logger.info("Reading enriched pathway IDs from %s", ids_path)
    disease_pathways = data_parser.parse_disease_pathways(file=ids_path)
    
    # Finding disease-specific and target-specific genes
    disease_genes = get_disease_genes(pathway2genes, disease_pathways)
    logger.info("Found %i disease-specific genes", len(disease_genes))

    target_genes = get_target_genes(gene2pathways, pathway2genes, target)
    logger.info("Found %i genes that are on the same pathways as target", len(target_genes))

    disease_and_target_genes = find_overlap(disease_genes, target_genes)
    logger.info("Found %i genes that are both disease-specific and on the same pathways as target", len(disease_and_target_genes))

    # Scoring
    logger.info("Calculating scores")
    scores = calculate_scores(disease_and_target_genes, gene2pathways, pathway2genes, disease_pathways, target)
    data_parser.scores_to_TSV(scores)

    # Prioritization with network propagation
    logger.info("Parsing interactions file")
    interactions = data_parser.parse_interactions(interactions_file)
    network = networkx.from_edgelist(interactions)
    logger.info("Built network with %i nodes and %i interactions", len(network.nodes()), len(network.edges()))

    # logger.info("Printing interactions")
    # data_parser.interactions_to_TSV(interactions)

    # Visualization


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="TODO"
    )
    
    parser.add_argument('--pathway_mapping_file',
                        help='TODO',
                        type=pathlib.Path,
                        required=True)

    parser.add_argument('--interactions_file',
                        help='TODO',
                        type=pathlib.Path,
                        required=True)

    parser.add_argument('--target',
                        help='TODO',
                        type=str,
                        required=True)

    parser.add_argument('--gsea_ids_file',
                        help='TSV of pathway IDs (column: pathwayId). Defaults to the legacy hardcoded path.',
                        type=pathlib.Path,
                        required=False)

    parser.add_argument('--auto_fetch',
                        help='Automatically download missing Reactome mapping/interactions if using default filenames',
                        action='store_true')

    args = parser.parse_args()

    try:
        # If provided, override the legacy hardcoded GSEA IDs path by temporarily
        # setting the environment variable used in main() to locate the file.
        # For backwards compatibility we keep the default behavior when not set.
        if args.gsea_ids_file is not None:
            # Monkey-patch the module-level default by passing through a kwarg via env
            os.environ["GSEA_IDS_FILE_OVERRIDE"] = str(args.gsea_ids_file)

        # Optionally auto-fetch missing inputs
        if args.auto_fetch:
            _ensure_reactome_mapping(args.pathway_mapping_file)
            _ensure_reactome_interactions(args.interactions_file)

        main(pathway_mapping_file=args.pathway_mapping_file,
             interactions_file=args.interactions_file,
             target=args.target)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)