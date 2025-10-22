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

    # Ensure FDR (qval) is available: compute from pval if missing
    if "qval" not in res_df.columns and "pval" in res_df.columns:
        p = pd.to_numeric(res_df["pval"], errors="coerce")
        n = p.shape[0]
        # Sort ascending by p-value
        order = p.sort_values().index
        ranks = pd.Series(range(1, n + 1), index=order)
        q_raw = (p.loc[order] * n / ranks)
        # Benjamini–Hochberg step-up: cumulative min from bottom
        q_adj = q_raw.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
        res_df["qval"] = q_adj.reindex(res_df.index).fillna(1.0)

    if output_tsv is None:
        output_tsv = f"{os.path.splitext(input_tsv)[0]}_gsea.tsv"

    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)
    res_df.to_csv(output_tsv, sep="\t", index=False)

    return res_df


def filter_pathway_ids_by_pvalue(
    gsea_results: pd.DataFrame,
    pval_threshold: float | None = None,
    output_ids_tsv: str | None = None,
    fdr_threshold: float | None = None,
    nes_positive: bool = False,
) -> pd.DataFrame:
    """
    Filter GSEA results and return a DataFrame with one column 'pathwayId'.

    - If fdr_threshold is provided and 'qval' exists, filter by qval <= fdr_threshold.
    - Otherwise filter by pval <= pval_threshold (requires 'pval').
    - If nes_positive is True and 'NES' exists, filter additionally by NES > 0.
    - Writes TSV if output_ids_tsv is provided.
    """
    if "ID" not in gsea_results.columns:
        raise ValueError("Expected 'ID' column in GSEA results")

    df = gsea_results
    mask = pd.Series(True, index=df.index)

    # Apply p-value filter if requested (skip if column missing)
    if pval_threshold is not None:
        if "pval" in df.columns:
            mask &= df["pval"] <= pval_threshold
        else:
            print("[GSEA] Warning: 'pval' column missing; skipping p-value filtering")

    # Apply FDR filter if requested (skip if column missing)
    if fdr_threshold is not None:
        if "qval" in df.columns:
            mask &= df["qval"] <= fdr_threshold
        else:
            print("[GSEA] Warning: 'qval' column missing; skipping FDR filtering")

    if nes_positive and "NES" in df.columns:
        mask &= df["NES"] > 0

    filtered = df.loc[mask, ["ID"]].copy().rename(columns={"ID": "pathwayId"})

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
        default=None,
        help="Optional p-value threshold; if set, keep rows with pval <= threshold",
    )
    parser.add_argument(
        "--fdr-threshold",
        dest="fdr_threshold",
        type=float,
        default=None,
        help="Optional FDR (qval) threshold; if provided and 'qval' exists, overrides p-value",
    )
    parser.add_argument(
        "--nes-positive",
        dest="nes_positive",
        action="store_true",
        help="If set, also require NES > 0 when filtering IDs",
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
        fdr_threshold=args.fdr_threshold,
        nes_positive=args.nes_positive,
    )


if __name__ == "__main__":
    main()


