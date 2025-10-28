import os
import sys
import logging
import argparse
import pathlib

import networkx

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_interactions(interactions_file):
    """
    Parses Reactome TSV file with 5 tab-separated columns: Gene1, Gene2, Annotation, Direction, Score

    creates directed functional interactions

    arguments:
    - interactions_file

    returns:
    - interactions: list of tuples
    """
    interactions = []
    genes = set()
    
    try:
        f = open(interactions_file, 'r')
    except Exception as e:
        logger.error("Opening provided Reactome file %s: %s", interactions_file, e)
        raise Exception("Cannot open provided Reactome file")
    
    # header
    line = f.readline()
    if not line.startswith("Gene1\t"):
        logging.error("Reactome file %s is headerless? expecting headers but got %s",
                      interactions_file, line)
        raise Exception("Reactome file problem")

    for line in f:
        split_line = line.rstrip().split('\t')

        if len(split_line) != 5:
            logger.error("Reactome file %s has bad line (not 5 tab-separated fields): %s",
                         interactions_file, line)
            raise Exception("Bad line in the Reactome file")

        (gene1, gene2, annotation, direction,  score) = split_line

        if direction == "->":
            interactions.append((gene1, gene2))
        if direction == "<-":
            interactions.append((gene2, gene1))
        
        genes.add(gene1)
        genes.add(gene2)

    f.close()

    return(interactions)


def interactions_to_TSV(interactions):
    '''
    Print interactions to stdout in TSV format, 2 columns: gene1, gene2

    arguments:
    - interactions: list of directed functional interactions
    '''
    for interaction in interactions:
        print(interaction[0] + "\t" + interaction[1])


def main(interactions_file):
    logger.info("Parsing interactions file")
    interactions = parse_interactions(interactions_file)

    network = networkx.from_edgelist(interactions)
    logger.info("Built network with %i nodes and %i interactions", len(network.nodes()), len(network.edges()))

    logger.info("Printing interactions")
    interactions_to_TSV(interactions)


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
    
    parser.add_argument('--interactions_file',
                        help='Reactome functional interaction file',
                        type=pathlib.Path,
                        required=True)
    
    args = parser.parse_args()

    try:
        main(interactions_file=args.interactions_file)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)