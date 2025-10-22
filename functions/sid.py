import argparse
import os
from typing import Dict, Tuple, List, Optional

import pandas as pd


# Global caches
GENE_NAME_TO_ID: Optional[Dict[str, str]] = None
GENE_ID_TO_NAME: Optional[Dict[str, str]] = None
DISEASE_NAME_TO_ID: Optional[Dict[str, str]] = None
DISEASE_ASSOC_DATA: Optional[pd.DataFrame] = None


def _load_two_column_tsv_to_map(filepath: str) -> Dict[str, str]:
    """Read a two-column TSV and return mapping col1 -> col2 (as strings)."""
    df = pd.read_csv(filepath, sep="\t", dtype=str)
    if df.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {filepath}")
    # Use the first two columns
    col1, col2 = df.columns[:2]
    # Drop duplicates keeping the first occurrence to ensure one-to-one mapping
    mapped = (
        df[[col1, col2]]
        .dropna()
        .drop_duplicates(subset=[col1])
        .set_index(col1)[col2]
        .to_dict()
    )
    return mapped


def _load_gene_maps(gene_data_file: str = "data/gene_data.txt") -> None:
    global GENE_NAME_TO_ID, GENE_ID_TO_NAME
    if GENE_NAME_TO_ID is not None and GENE_ID_TO_NAME is not None:
        return
    name_to_id = _load_two_column_tsv_to_map(gene_data_file)  # gene_name -> gene_id
    # Invert; if duplicates exist, first occurrence wins
    id_to_name = {gid: gname for gname, gid in name_to_id.items()}
    GENE_NAME_TO_ID = name_to_id
    GENE_ID_TO_NAME = id_to_name


def _load_disease_map(disease_data_file: str = "data/disease_data.txt") -> None:
    global DISEASE_NAME_TO_ID
    if DISEASE_NAME_TO_ID is not None:
        return
    DISEASE_NAME_TO_ID = _load_two_column_tsv_to_map(disease_data_file)  # disease_name -> efo_id


def _load_associations() -> None:
    """Load disease-target association data into memory.

    Prefers parquet files under data/association_by_datatype_indirect/.
    Falls back to TSV at data/input/auto-input/filtered_associations.tsv if present.
    """
    global DISEASE_ASSOC_DATA
    if DISEASE_ASSOC_DATA is not None:
        return

    parquet_dir = "data/association_by_datatype_indirect"
    tsv_fallback = "data/input/auto-input/filtered_associations.tsv"

    if os.path.isdir(parquet_dir):
        parts = [
            os.path.join(parquet_dir, f)
            for f in sorted(os.listdir(parquet_dir))
            if f.endswith(".parquet") or f.endswith(".snappy.parquet")
        ]
        if not parts:
            raise FileNotFoundError(f"No parquet files found in {parquet_dir}")
        # Read only required columns if available
        frames: List[pd.DataFrame] = []
        for p in parts:
            try:
                dfp = pd.read_parquet(p, columns=["diseaseId", "targetId", "score"])  # type: ignore[arg-type]
            except Exception:
                dfp = pd.read_parquet(p)  # type: ignore[arg-type]
            frames.append(dfp)
        df = pd.concat(frames, ignore_index=True)
    elif os.path.isfile(tsv_fallback):
        df = pd.read_csv(tsv_fallback, sep="\t", dtype={"diseaseId": str, "targetId": str})
        missing = {c for c in ["diseaseId", "targetId", "score"] if c not in df.columns}
        if missing:
            raise ValueError(f"Association TSV missing columns: {missing}")
    else:
        raise FileNotFoundError(
            "Could not locate association data. Expected parquet under "
            f"{parquet_dir} or TSV at {tsv_fallback}"
        )

    # Keep only relevant columns
    needed_cols = ["diseaseId", "targetId", "score"]
    present = [c for c in needed_cols if c in df.columns]
    DISEASE_ASSOC_DATA = df[present].copy()


def load_data() -> None:
    """Load all required resources into memory."""
    _load_gene_maps()
    _load_disease_map()
    _load_associations()


def validate_input(gene_name: str, disease_name: str) -> Tuple[str, str]:
    """Validate inputs and return (gene_id, disease_id)."""
    if GENE_NAME_TO_ID is None or DISEASE_NAME_TO_ID is None:
        load_data()
    assert GENE_NAME_TO_ID is not None and DISEASE_NAME_TO_ID is not None
    if gene_name not in GENE_NAME_TO_ID:
        raise ValueError(f"Gene name '{gene_name}' not found in gene map")
    if disease_name not in DISEASE_NAME_TO_ID:
        raise ValueError(f"Disease name '{disease_name}' not found in disease map")
    return GENE_NAME_TO_ID[gene_name], DISEASE_NAME_TO_ID[disease_name]


def get_genes_associated_with_disease(disease_id: str) -> pd.DataFrame:
    """Return targetId and score for genes associated with a disease ID.

    This function preserves the original semantics your colleague described and
    returns the association at the identifier level: columns `targetId`, `score`.
    """
    if DISEASE_ASSOC_DATA is None:
        _load_associations()
    assert DISEASE_ASSOC_DATA is not None
    df = DISEASE_ASSOC_DATA[DISEASE_ASSOC_DATA["diseaseId"] == disease_id]
    df_selected = df[["targetId", "score"]].sort_values(by="score", ascending=False)
    return df_selected.reset_index(drop=True)


def build_gsea_input_for_disease(disease_name: str) -> pd.DataFrame:
    """Create a DataFrame with columns `symbol`, `globalScore` for GSEA input.

    - Maps disease name -> EFO ID
    - Selects associations for that disease
    - Aggregates to the best score per targetId
    - Maps targetId (e.g., ENSG...) -> gene symbol using gene_data.txt
    - Drops entries where symbol mapping is unavailable
    """
    load_data()
    assert DISEASE_ASSOC_DATA is not None and GENE_ID_TO_NAME is not None and DISEASE_NAME_TO_ID is not None

    if disease_name not in DISEASE_NAME_TO_ID:
        raise ValueError(f"Unknown disease name: {disease_name}")

    disease_id = DISEASE_NAME_TO_ID[disease_name]
    subset = DISEASE_ASSOC_DATA.loc[DISEASE_ASSOC_DATA["diseaseId"] == disease_id, ["targetId", "score"]]
    if subset.empty:
        return pd.DataFrame(columns=["symbol", "globalScore"])  # no data for disease

    # Aggregate to one score per target (max score)
    agg = subset.groupby("targetId", as_index=False)["score"].max()

    # Map ENSG -> symbol
    agg["symbol"] = agg["targetId"].map(GENE_ID_TO_NAME)
    agg = agg.dropna(subset=["symbol"])  # require symbol for GSEA

    out = agg[["symbol", "score"]].rename(columns={"score": "globalScore"}).sort_values(
        by="globalScore", ascending=False
    )
    return out.reset_index(drop=True)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare GSEA input (symbol, globalScore) for a disease")
    parser.add_argument("--gene_name", required=True, help="Primary target name (validated, not used in output)")
    parser.add_argument("--disease_name", required=True, help="Disease name as in data/disease_data.txt")
    parser.add_argument(
        "--output",
        dest="output_tsv",
        default=os.path.join("gsea", "input", "auto-input", "final_gsea_input.tsv"),
        help="Output TSV path (default: gsea/input/auto-input/final_gsea_input.tsv)",
    )
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()

    # Validate inputs and load
    load_data()
    validate_input(args.gene_name, args.disease_name)

    gsea_df = build_gsea_input_for_disease(args.disease_name)
    os.makedirs(os.path.dirname(args.output_tsv) or ".", exist_ok=True)
    gsea_df.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
