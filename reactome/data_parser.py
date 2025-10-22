import logging

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_pathway_mapping(pathway_mapping_file):
    """
    Maps genes to pathways and vice versa

    arguments:
    - pathway_mapping_file, no header

    returns:
    - gene2pathways
    - pathway2genes
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

    logger.info("Found %i genes in %i pathways", len(gene2pathways), len(pathway2genes))

    return(gene2pathways, pathway2genes)


def parse_disease_pathways(file):
    list = []
    f = open(file, 'r')
    # skip header
    f.readline()

    for line in f:
        list.append(line.rstrip('\n'))

    f.close()

    logger.info("Found %i disease-specific pathways", len(list))

    return(list)


def parse_interactions(interactions_file):
    """
    Parses Reactome file 

    arguments:
    - interactions_file: with 5 tab-separated columns

    returns:
    - interactions: list of tuples, directed functional interactions
    """
    interaction_types = set()
    interactions = []
    genes = set()
    
    try:
        f = open(interactions_file, 'r')
    except Exception as e:
        logger.error("Opening provided Reactome file %s: %s", interactions_file, e)
        raise Exception("Cannot open provided Reactome file")
    
    # skip header
    line = f.readline()
    if not line.startswith("Gene1\t"):
        logging.error("Reactome file %s is headerless? expecting headers but got %s",
                      interactions_file, line)
        raise Exception("Reactome file problem")

    for line in f:
        split_line = line.rstrip('\n').split('\t')

        if len(split_line) != 5:
            logger.error("Reactome file %s has bad line (not 5 tab-separated fields): %s",
                         interactions_file, line)
            raise Exception("Bad line in the Reactome file")

        (gene1, gene2, annotation, direction,  score) = split_line

        interaction_types.add(annotation)
        if direction == "->":
            interactions.append((gene1, gene2))
        if direction == "<-":
            interactions.append((gene2, gene1))
        
        genes.add(gene1)
        genes.add(gene2)

    f.close()

    return(interactions)


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


def interactions_to_TSV(interactions):
    '''
    Print interactions to stdout in TSV format, 2 columns: gene1, gene2

    arguments:
    - interactions: list of directed functional interactions
    '''

    for interaction in interactions:
        print(interaction[0] + "\t" + interaction[1])