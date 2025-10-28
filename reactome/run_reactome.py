import os
import sys
import logging
import argparse
import pathlib

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def get_disease_genes(pathway2genes, disease_pathways):
    """
    Find all genes that are on disease-specific pathways

    arguments:
    - pathway2genes: dict, key=pathway, value=list of genes
    - disease_pathways
    
    returns:
    - list of disease-specific genes
    """
    disease_genes = set()

    for path in disease_pathways:
        genes_on_path = pathway2genes[path]

        for gene in genes_on_path:
            disease_genes.add(gene)
    
    return(disease_genes)


def get_target_genes(gene2pathways, pathway2genes, target):
    """
    Find all genes that are on the same pathways as target

    arguments:
    - gene2pathways: dict, key=gene, value=list of pathways with gene
    - pathway2genes: dict, key=pathway, value=list of genes on pathway
    - target: str, gene name of the target of interest
    
    returns:
    - list of target-specific genes
    """
    target_genes = set()

    target_paths = gene2pathways[target]
    for path in target_paths:
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
    For every disease- and target-specific gene:
    score = (#disease_paths_with_gene + #target_paths_with_gene) / (#disese_paths + #target_paths)

    arguments:
    - genes: list of disease- and target-specific genes
    - gene2pathways: dict, key=gene, value=list of pathways with gene
    - pathway2genes: dict, key=pathway, value=list of genes on pathway
    - target: str, gene name of the target of interest
    
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


def main(pathway_mapping_file, disease_pathways_file, target):

    logger.info("Parsing gene-to-pathway mapping file")
    (gene2pathways, pathway2genes) = data_parser.parse_pathway_mapping(pathway_mapping_file)

    logger.info("Computing GSEA")
    disease_pathways = data_parser.parse_disease_pathways(disease_pathways_file)
    logger.info("Found %i disease-specific pathways", len(disease_pathways))
    
    # Finding disease-specific and target-specific genes
    disease_genes = get_disease_genes(pathway2genes, disease_pathways)
    logger.info("Found %i disease-specific genes", len(disease_genes))

    target_genes = get_target_genes(gene2pathways, pathway2genes, target)
    logger.info("Found %i genes that are on the same pathways as target", len(target_genes))

    disease_and_target_genes = find_overlap(disease_genes, target_genes)
    logger.info("Found %i genes that are both disease-specific and on the same pathways as target", len(disease_and_target_genes))

    # Pathway selectivity
    logger.info("Calculating scores")
    scores = calculate_scores(disease_and_target_genes, gene2pathways, pathway2genes, disease_pathways, target)
    data_parser.scores_to_TSV(scores)


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
                        help='Path to pathway mapping file from Reactome, no header',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--disease_pathways_file',
                        help='Path to file with significantly enriched pathways from GSEA, one pathway per line, no header',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--target',
                        help='target gene name',
                        type=str,
                        required=True)

    args = parser.parse_args()

    try:
        main(pathway_mapping_file=args.pathway_mapping_file,
             disease_pathways_file=args.disease_pathways_file,
             target=args.target)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)