import os
import sys
import logging
import argparse
import pathlib

import networkx

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


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


def main(pathway_mapping_file, disease_pathways_file, interactions_file, target):

    logger.info("Parsing gene-to-pathway mapping file")
    (gene2pathways, pathway2genes) = data_parser.parse_pathway_mapping(pathway_mapping_file)

    logger.info("Computing GSEA")
    # Read list of enriched pathways from provided override or legacy default
    disease_pathways = data_parser.parse_disease_pathways(disease_pathways_file)
    
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

    # Network propagation
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
                        help='From Reactome, see GitHub README',
                        type=pathlib.Path,
                        required=True)
    
    parser.add_argument('--disease_pathways_file',
                        help='Significantly enriched pathways from GSEA, one pathway ID per line',
                        type=pathlib.Path,
                        required=True)

    parser.add_argument('--interactions_file',
                        help='From Reactome, see GitHub README',
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
             interactions_file=args.interactions_file,
             target=args.target)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)