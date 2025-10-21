import os
import argparse
import pandas as pd
import blitzgsea as blitz


def load_custom_gmt(path: str) -> dict:
    """
    Parse a GMT file into a dict mapping term name to list of genes.
    The format is: term<TAB>description<TAB>gene1<TAB>gene2<...>
    """
    with open(path, "r") as f:
        return {
            parts[0]: parts[2:]
            for line in f
            if (parts := line.strip().split("\t")) and len(parts) > 2
        }


def run_gsea_from_tsv(
    input_tsv: str,
    gmt_file: str,
    output_tsv: str | None = None,
    processes: int = 4,
) -> pd.DataFrame:
    """
    Read TSV with columns 'symbol' and 'globalScore', run GSEA, and write results.

    - input_tsv: path to input TSV with columns 'symbol' and 'globalScore'.
    - gmt_file: path to GMT file (e.g., ReactomePathways_merged.gmt).
    - output_tsv: destination TSV for results; if None, uses <input>_gsea.tsv.
    - processes: number of processes for GSEA.
    """
    library_sets = load_custom_gmt(gmt_file)

    df = pd.read_csv(input_tsv, sep="\t", header=0, index_col=None)

    # blitz.gsea expects a two-column DataFrame with columns [0] (scores) and [1] (genes)
    gsea_df = pd.DataFrame()
    gsea_df[1] = df["symbol"]
    gsea_df[0] = pd.to_numeric(df["globalScore"], errors="coerce")
    gsea_df = gsea_df.dropna(subset=[0])

    res_df = blitz.gsea(gsea_df, library_sets, processes=processes).reset_index(names="Term")

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

    if output_tsv is None:
        output_tsv = f"{os.path.splitext(input_tsv)[0]}_gsea.tsv"

    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)
    res_df.to_csv(output_tsv, sep="\t", index=False)

    return res_df


def filter_pathway_ids_by_pvalue(
    gsea_results: pd.DataFrame,
    pval_threshold: float = 0.05,
    output_ids_tsv: str | None = None,
) -> pd.DataFrame:
    """
    Filter GSEA results by pval <= threshold and return a DataFrame with one
    column 'pathwayId' (renamed from 'ID'). Optionally write it to a TSV.

    - gsea_results: DataFrame produced by run_gsea_from_tsv
    - pval_threshold: inclusive threshold for 'pval' column
    - output_ids_tsv: optional path to save the filtered IDs as TSV
    """
    if "pval" not in gsea_results.columns:
        raise ValueError("Expected 'pval' column in GSEA results")
    if "ID" not in gsea_results.columns:
        raise ValueError("Expected 'ID' column in GSEA results")

    filtered = gsea_results.loc[gsea_results["pval"] <= pval_threshold, ["ID"]].copy()
    filtered = filtered.rename(columns={"ID": "pathwayId"})

    if output_ids_tsv is not None:
        os.makedirs(os.path.dirname(output_ids_tsv) or ".", exist_ok=True)
        filtered.to_csv(output_ids_tsv, sep="\t", index=False)

    return filtered


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run GSEA on TSV of gene scores")
    parser.add_argument("input_tsv", help="TSV with columns 'symbol' and 'globalScore'")
    parser.add_argument("gmt_file", help="Path to GMT file, e.g. ReactomePathways_merged.gmt")
    parser.add_argument(
        "--output",
        dest="output_tsv",
        default=None,
        help="Output TSV path. Defaults to <input>_gsea.tsv",
    )
    parser.add_argument(
        "--pval-threshold",
        dest="pval_threshold",
        type=float,
        default=0.05,
        help="Inclusive p-value threshold to filter pathway IDs (default: 0.05)",
    )
    parser.add_argument(
        "--output-ids",
        dest="output_ids_tsv",
        default=None,
        help="Optional TSV to write filtered pathway IDs (column: pathwayId)",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=4,
        help="Number of processes for GSEA (default: 4)",
    )
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    res = run_gsea_from_tsv(
        input_tsv=args.input_tsv,
        gmt_file=args.gmt_file,
        output_tsv=args.output_tsv,
        processes=args.processes,
    )
    # Apply optional filtering and output of pathway IDs
    filter_pathway_ids_by_pvalue(
        gsea_results=res,
        pval_threshold=args.pval_threshold,
        output_ids_tsv=args.output_ids_tsv,
    )


if __name__ == "__main__":
    main()


