import os
import logging

import pandas

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_target_parquet(target_parquets_dir):
    """
    Parses target parquet files from Open Targets Platform with 2 columns: id, approvedSymbol

    arguments:
    - target_parquets_dir: directory to Open Targets targets parquets

    return:
    - target2symbol: dict, key=target ID, value=approved symbol
    """
    target2symbol = {}
    parquet_files = []

    for f in os.listdir(target_parquets_dir):
        if f.endswith(".parquet") or f.endswith(".snappy.parquet"):
            parquet_file = os.path.join(target_parquets_dir, f)
            parquet_files.append(parquet_file)

    for file in parquet_files:
        try:
            dfp = pandas.read_parquet(file, columns=["id", "approvedSymbol"])
        except Exception:
            logger.error("Cannot open parquet file %s", file)
        parquet_target2symbol = (dfp[["id", "approvedSymbol"]]
                                .drop_duplicates(subset=["id"])
                                .set_index("id")[ "approvedSymbol"]
                                .to_dict())
        target2symbol.update(parquet_target2symbol)

    logger.info("Found %i target-symbol associations", len(target2symbol))

    return(target2symbol)


def parse_associations_parquet(associations_parquets_dir, disease, datatype):
    """
    Parses associations parquet files from Open Targets Platform with 4 columns:
    disease_id, datatype, score, target_id
    Filter on disease and datatype
    More on the association scores:
    https://platform.opentargets.org/downloads/association_by_datasource_indirect/schema

    arguments:
    - associations_parquets_dir: directory to Open Targets associations parquets

    return:
    - target2score: dict, key=target, value=association score
    """
    target2score = {}
    parquet_files = []

    for f in os.listdir(associations_parquets_dir):
        if f.endswith(".parquet") or f.endswith(".snappy.parquet"):
            parquet_file = os.path.join(associations_parquets_dir, f)
            parquet_files.append(parquet_file)

    for file in parquet_files:
        try:
            dfp = pandas.read_parquet(file, columns=["diseaseId", "targetId", "score", "datatypeId"])
        except Exception:
            logger.error("Cannot open parquet file %s", file)
        
        df_filtered = dfp[(dfp["diseaseId"] == disease) & (dfp["datatypeId"] == datatype)]

        parquet_target2score = (df_filtered[["targetId", "score"]]
                                .drop_duplicates(subset=["targetId"])
                                .set_index("targetId")[ "score"]
                                .to_dict())
        target2score.update(parquet_target2score)

    logger.info("Found %i target-score associations for %s", len(target2score), disease)

    return(target2score)


def build_gsea_input(target2symbol, target2score):
    """
    Create GSEA input as specified here: https://github.com/MaayanLab/blitzgsea

    arguments:
    - target2symbol: dict, key=target, value=approved symbol
    - target2score: dict, key=target, value=association score

    returns:
    - gsea_input: DataFrame with two columns: symbol, score
    """
    symbol2score = {}

    for target in target2score:
        if target in target2symbol:
            symbol = target2symbol[target]
            score = target2score[target]
            symbol2score[symbol] = score

    gsea_input = pandas.DataFrame(list(symbol2score.items()))

    return gsea_input


def parse_gmt_file(gmt_file):
    """
    Parse a GMT file with columns: term_name, description, gene1, gene2, ...

    returns:
    - pathways: dict, key=term_name, key=list of genes
    """
    pathways = {}

    with open(gmt_file, "r") as f:
        for line in f:
            split_line = line.rstrip().split('\t')

            if len(split_line) < 2:
                logger.error("GMT file %s has bad line: %s", gmt_file, split_line)
                raise Exception("Bad line in the GMT file")

            term_name = split_line[0]
            genes = split_line[2:]
            pathways[term_name] = genes

    return pathways


def pathways_to_TSV(pathways):
    '''
    Print significantly enriched pathways to stdout, one pathway ID per line

    arguments:
    - pathways: list of str
    '''

    for pathway in pathways:
        print(pathway)