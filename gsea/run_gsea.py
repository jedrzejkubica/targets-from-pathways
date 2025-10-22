import os
import sys
import argparse
import pathlib
import logging

import pandas

import blitzgsea

import data_parser


# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def run_gsea(gsea_input, library_sets):
    """
    Read gsea_input, run GSEA, and write results.

    arguments:
    - gsea_input
    - gmt_file: path to GMT file (e.g., ReactomePathways_merged.gmt).

    returns:
    - pandas.dataframe with GSEA results
    """
    res_df = blitzgsea.gsea(gsea_input, library_sets).reset_index(names="Term")

    # Include all genes in the pathway for convenience
    res_df["propagated_edge"] = res_df["Term"].apply(
        lambda t: ",".join(library_sets.get(t, [])) if library_sets.get(t) else ""
    )

    # Extract ID from terms like "Term [R-HSA-12345]" and clean the visible name
    term_series = res_df["Term"]
    res_df["ID"] = term_series.str.extract(r"\[([^\]]+)\]", expand=False).fillna("")
    res_df["Term"] = term_series.str.replace(r"\s*\[[^\]]+\]", "", regex=True).str.strip()

    # Normalize leading_edge to a CSV string if provided by blitz.gsea
    if "leading_edge" in res_df.columns:
        res_df["leading_edge"] = res_df["leading_edge"].apply(
            lambda x: ",".join(map(str, x)) if isinstance(x, (list, tuple)) else str(x)
        )

    # Reorder columns for readability
    first_cols = ["Term", "ID"]
    res_df = res_df[first_cols + [c for c in res_df.columns if c not in first_cols]]

    # Ensure FDR (qval) is available: compute from pval if missing
    if "qval" not in res_df.columns and "pval" in res_df.columns:
        p = pandas.to_numeric(res_df["pval"], errors="coerce")
        n = p.shape[0]
        # Sort ascending by p-value
        order = p.sort_values().index
        ranks = pandas.Series(range(1, n + 1), index=order)
        q_raw = (p.loc[order] * n / ranks)
        # Benjamini–Hochberg step-up: cumulative min from bottom
        q_adj = q_raw.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
        res_df["qval"] = q_adj.reindex(res_df.index).fillna(1.0)

    return res_df


def filter_pathways_by_pvalue(gsea_results, pval_threshold, fdr_threshold):
    """
    Filter GSEA results based on p-value or FDR threshold
    
    returns:
     - pathways: list of str, significantly enriched pathways
    """
    # Create p-value and FDR filters
    mask = pandas.Series(True, index=gsea_results.index)
    mask &= gsea_results["pval"] <= pval_threshold
    mask &= gsea_results["qval"] <= fdr_threshold

    gsea_results_filtered = gsea_results.loc[mask, ["ID"]]
    pathways = gsea_results_filtered["ID"].tolist()

    return pathways


def main(target_parquets_dir, associations_parquets_dir, disease, datatype, gmt_file, pval_threshold, fdr_threshold):
    logger.info("Parsing parquet files")
    target2symbol = data_parser.parse_target_parquet(target_parquets_dir)
    target2score = data_parser.parse_associations_parquet(associations_parquets_dir, disease, datatype)

    # Build GSEA input
    logger.info("Building GSEA input")
    gsea_input = data_parser.build_gsea_input(target2symbol, target2score)
    library_sets = data_parser.parse_gmt_file(gmt_file)

    # Run GSEA
    logger.info("Running GSEA")
    gsea_results = run_gsea(gsea_input, library_sets)

    # Filter  pathways by p-value
    pathways = filter_pathways_by_pvalue(gsea_results, pval_threshold, fdr_threshold)

    logger.info("Printing significantly enriched pathways")
    data_parser.pathways_to_TSV(pathways)


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        description=(
            """
            End-to-end pipeline: build GSEA input from associations, run GSEA,"
            "compute Reactome-based disease-target pathway scores.
            """
        )
    )
    parser.add_argument("--target_parquets_dir",
                        type=pathlib.Path,
                        required=True,
                        help="TODO")
    parser.add_argument("--associations_parquets_dir",
                        type=str,
                        default=None,
                        help="TODO")
    parser.add_argument("--disease",
                        type=str,
                        required=True, 
                        help="Disease EFO ID")
    parser.add_argument("--datatype",
                        type=str,
                        help="association datatype filter for GSEA, e.g. genetic_association")
    parser.add_argument("--gmt_file",
                        type=str,
                        required=True,
                        help="Path to Reactome GMT file (e.g., gsea/Reactome_2025/<file>.gmt)")
    parser.add_argument("--pval_threshold",
                        type=float,
                        help="p-value threshold for GSEA")
    parser.add_argument("--fdr_threshold",
                        type=float,
                        help="FDR threshold for GSEA")

    args = parser.parse_args()

    try:
        main(target_parquets_dir=args.target_parquets_dir,
            associations_parquets_dir=args.associations_parquets_dir,
            disease=args.disease,
            datatype=args.datatype,
            gmt_file=args.gmt_file,
            pval_threshold=args.pval_threshold,
            fdr_threshold=args.fdr_threshold)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)


