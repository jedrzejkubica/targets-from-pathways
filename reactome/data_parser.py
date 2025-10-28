import logging

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_pathway_mapping(pathway_mapping_file):
    """
    Maps genes to pathways and vice versa

    arguments:
    - pathway_mapping_file

    returns:
    - gene2pathways: dict, key=gene, value=list of pathways with gene
    - pathway2genes: dict, key=pathway, value=list of genes on pathway
    """
    gene2pathways = {}
    pathway2genes = {}

    try:
        f = open(pathway_mapping_file, 'r')
    except Exception as e:
        logger.error("Opening provided Reactome mapping file %s: %s", pathway_mapping_file, e)
        raise Exception("Cannot open provided Reactome mapping file")

    for line in f:
        split_line = line.rstrip('\n').split('\t')

        if len(split_line) != 8:
            logger.error("Reactome file %s has bad line (not 8 tab-separated fields): %s",
                         pathway_mapping_file, line)
            raise Exception("Bad line in the Reactome mapping file")

        (ensembleID, geneID, gene, pathwayID, url, pathway_name, evidence, species) = split_line

        if species != "Homo sapiens":
            continue

        gene_name = gene.split(" ")[0]

        if gene_name in gene2pathways:
            gene2pathways[gene_name].append(pathwayID)
        else:
            gene2pathways[gene_name] = [pathwayID]

        if pathwayID in pathway2genes:
            pathway2genes[pathwayID].append(gene_name)
        else:
            pathway2genes[pathwayID] = [gene_name]

    f.close()

    logger.info("Found %i genes on %i pathways", len(gene2pathways), len(pathway2genes))

    return(gene2pathways, pathway2genes)


def parse_disease_pathways(file):
    list = []
    f = open(file, 'r')

    for line in f:
        list.append(line.rstrip('\n'))

    f.close()

    logger.info("Found %i disease-specific pathways", len(list))

    return(list)


def scores_to_TSV(scores):
    '''
    Print scores to stdout in TSV format, 2 columns: gene_name score

    arguments:
    - scores: dict with key=gene, value=score
    '''
    # header
    print("GENE\tSCORE")

    for (gene, score) in sorted(scores.items()):
        print(gene + "\t" + str(score))