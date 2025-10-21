import os
import sys
import logging
import argparse
import pathlib

import networkx

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def get_disease_genes(pathway2genes, enriched_pathways):
    """
    Find all genes that are on disease-specific pathways

    arguments:
    - pathway2genes
    - enriched_pathways
    
    returns:
    - list of disease-specific genes
    """
    disease_genes = set()
    for path in enriched_pathways:
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


def main(pathway_mapping_file, interactions_file):

    logger.info("Parsing gene-to-pathway mapping file")
    (gene2pathways, pathway2genes) = data_parser.parse_pathway_mapping(pathway_mapping_file)

    logger.info("Computing GSEA")
    # For now read list of enriched pathways from a file
    enriched_pathways = data_parser.parse_enriched_pathways(file="/home/kubicaj/open-targets-hackathon/targets-from-pathways/data/gsea_output_ids/OT-EFO_0004248_gsea_ids.tsv")
    
    # Filtering
    disease_genes = get_disease_genes(pathway2genes, enriched_pathways)
    logger.info("Found %i disease-specific genes", len(disease_genes))

    target_genes = get_target_genes(gene2pathways, pathway2genes, target="BTG4")
    logger.info("Found %i genes that are on the same pathways as target", len(target_genes))

    disease_and_target_genes = find_overlap(disease_genes, target_genes)
    logger.info("Found %i genes that are both disease-specific and on the same pathways as target", len(disease_and_target_genes))

    # Prioritization with network propagation
    logger.info("Parsing interactions file")
    interactions = data_parser.parse_interactions(interactions_file)
    network = networkx.from_edgelist(interactions)
    logger.info("Built network with %i nodes and %i interactions", len(network.nodes()), len(network.edges()))

    # logger.info("Printing interactions")
    # data_parser.scores_to_TSV(interactions)

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

    args = parser.parse_args()

    try:
        main(pathway_mapping_file=args.pathway_mapping_file,
             interactions_file=args.interactions_file)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)