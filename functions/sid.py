import argparse
import os
from typing import Dict, Tuple, List, Optional

import pandas as pd


# Global caches
GENE_NAME_TO_ID: Optional[Dict[str, str]] = None
GENE_ID_TO_NAME: Optional[Dict[str, str]] = None
DISEASE_NAME_TO_ID: Optional[Dict[str, str]] = None
DISEASE_ASSOC_DATA: Optional[pd.DataFrame] = None

# Configurable data sources (can be set via CLI/env)
ASSOCIATIONS_DIR: Optional[str] = None
ASSOCIATIONS_TSV: Optional[str] = None
TARGET_DIR: Optional[str] = None


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


def _load_gene_maps() -> None:
    """Populate gene symbol/name <-> Ensembl ID maps from Open Targets target parquet.

    Requires TARGET_DIR to point to a directory with target parquet parts containing
    at least the columns 'id' and a symbol column among: approvedSymbol|symbol|gene_symbol.
    """
    global GENE_NAME_TO_ID, GENE_ID_TO_NAME
    if GENE_NAME_TO_ID is not None and GENE_ID_TO_NAME is not None:
        return

    # Resolve target directory from CLI/env config
    target_dir = TARGET_DIR or os.environ.get("SID_TARGET_DIR")
    if not target_dir or not os.path.isdir(target_dir):
        raise FileNotFoundError(
            "Target parquet directory not set or not found. Provide --target_dir or set SID_TARGET_DIR."
        )

    id_to_symbol: Dict[str, str] = {}
    parts = [
        os.path.join(target_dir, f)
        for f in sorted(os.listdir(target_dir))
        if f.endswith(".parquet") or f.endswith(".snappy.parquet")
    ]
    if not parts:
        raise FileNotFoundError(f"No parquet files found in {target_dir}")

    for p in parts:
        try:
            dfp = pd.read_parquet(p)  # type: ignore[arg-type]
        except Exception:
            continue
        cols = {c.lower(): c for c in dfp.columns}
        id_col = cols.get("id")
        sym_col = cols.get("approvedsymbol") or cols.get("symbol") or cols.get("gene_symbol")
        if not id_col or not sym_col:
            continue
        tmp = (
            dfp[[id_col, sym_col]]
            .dropna()
            .drop_duplicates(subset=[id_col])
            .set_index(id_col)[sym_col]
            .to_dict()
        )
        id_to_symbol.update(tmp)

    if not id_to_symbol:
        raise ValueError(
            "Could not build id->symbol map from target parquet. Ensure files have 'id' and a symbol column."
        )

    name_to_id: Dict[str, str] = {}
    for gid, gname in id_to_symbol.items():
        if gname and gid and gname not in name_to_id:
            name_to_id[gname] = gid

    GENE_ID_TO_NAME = id_to_symbol
    GENE_NAME_TO_ID = name_to_id


def _load_disease_map(disease_data_file: str = "data/disease_data.txt") -> None:
    global DISEASE_NAME_TO_ID
    if DISEASE_NAME_TO_ID is not None:
        return
    DISEASE_NAME_TO_ID = _load_two_column_tsv_to_map(disease_data_file)  # disease_name -> efo_id


def _load_associations() -> None:
    """Load disease-target association data into memory using configured sources.

    Priority: ASSOCIATIONS_DIR (parquet directory) → ASSOCIATIONS_TSV (single TSV) → error.
    """
    global DISEASE_ASSOC_DATA
    if DISEASE_ASSOC_DATA is not None:
        return

    parquet_dir = ASSOCIATIONS_DIR or os.environ.get("SID_ASSOCIATIONS_DIR")
    tsv_fallback = ASSOCIATIONS_TSV or os.environ.get("SID_ASSOCIATIONS_TSV")

    if parquet_dir and os.path.isdir(parquet_dir):
        parts = [
            os.path.join(parquet_dir, f)
            for f in sorted(os.listdir(parquet_dir))
            if f.endswith(".parquet") or f.endswith(".snappy.parquet")
        ]
        if not parts:
            raise FileNotFoundError(f"No parquet files found in {parquet_dir}")
        # Read only required columns if available; include optional datatypeId if present
        frames: List[pd.DataFrame] = []
        for p in parts:
            try:
                dfp = pd.read_parquet(
                    p,
                    columns=["diseaseId", "targetId", "score", "datatypeId"],  # type: ignore[arg-type]
                )
            except Exception:
                try:
                    dfp = pd.read_parquet(
                        p,
                        columns=["diseaseId", "targetId", "score"],  # type: ignore[arg-type]
                    )
                except Exception:
                    dfp = pd.read_parquet(p)  # type: ignore[arg-type]
            frames.append(dfp)
        df = pd.concat(frames, ignore_index=True)
    elif tsv_fallback and os.path.isfile(tsv_fallback):
        df = pd.read_csv(tsv_fallback, sep="\t", dtype={"diseaseId": str, "targetId": str})
        missing = {c for c in ["diseaseId", "targetId", "score"] if c not in df.columns}
        if missing:
            raise ValueError(f"Association TSV missing columns: {missing}")
    else:
        raise FileNotFoundError("Association data not found. Provide --associations_dir or --associations_tsv (or set env).")

    # Keep only relevant columns
    needed_cols = ["diseaseId", "targetId", "score", "datatypeId"]
    present = [c for c in needed_cols if c in df.columns]
    DISEASE_ASSOC_DATA = df[present].copy()


def load_data() -> None:
    """Load all required resources into memory."""
    _load_gene_maps()
    _load_disease_map()
    _load_associations()


def configure_data_sources(associations_dir: Optional[str] = None,
                           associations_tsv: Optional[str] = None,
                           target_dir: Optional[str] = None) -> None:
    """Configure data source locations for parquet/TSV inputs."""
    global ASSOCIATIONS_DIR, ASSOCIATIONS_TSV, TARGET_DIR
    ASSOCIATIONS_DIR = associations_dir or ASSOCIATIONS_DIR
    ASSOCIATIONS_TSV = associations_tsv or ASSOCIATIONS_TSV
    TARGET_DIR = target_dir or TARGET_DIR


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


def resolve_disease_id(disease_name: Optional[str], disease_id: Optional[str]) -> str:
    """Resolve a disease ID from either a provided EFO ID or a disease name.

    Prefers disease_id when provided. Falls back to mapping disease_name via
    disease_data.txt. Raises ValueError if neither is available or cannot map.
    """
    load_data()
    if disease_id:
        return disease_id
    if not disease_name:
        raise ValueError("Provide either --disease_id or --disease_name")
    assert DISEASE_NAME_TO_ID is not None
    if disease_name not in DISEASE_NAME_TO_ID:
        raise ValueError(f"Disease name '{disease_name}' not found in disease map")
    return DISEASE_NAME_TO_ID[disease_name]


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


def build_gsea_input_for_disease(disease_name: str, datatype: Optional[str] = None) -> pd.DataFrame:
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
    subset_all = DISEASE_ASSOC_DATA.loc[DISEASE_ASSOC_DATA["diseaseId"] == disease_id]
    # Optional datatype filter (e.g., "genetic_association") if column present
    if datatype is not None and "datatypeId" in subset_all.columns:
        subset_all = subset_all.loc[subset_all["datatypeId"] == datatype]
    subset = subset_all[["targetId", "score"]]
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


def build_gsea_input_for_disease_id(disease_id: str, datatype: Optional[str] = None) -> pd.DataFrame:
    """Create GSEA input DataFrame for a given disease EFO ID.

    Equivalent to build_gsea_input_for_disease but takes an EFO directly.
    """
    load_data()
    assert DISEASE_ASSOC_DATA is not None and GENE_ID_TO_NAME is not None

    subset_all = DISEASE_ASSOC_DATA.loc[DISEASE_ASSOC_DATA["diseaseId"] == disease_id]
    if datatype is not None and "datatypeId" in subset_all.columns:
        subset_all = subset_all.loc[subset_all["datatypeId"] == datatype]
    subset = subset_all[["targetId", "score"]]
    if subset.empty:
        return pd.DataFrame(columns=["symbol", "globalScore"])  # no data for disease

    agg = subset.groupby("targetId", as_index=False)["score"].max()
    agg["symbol"] = agg["targetId"].map(GENE_ID_TO_NAME)
    agg = agg.dropna(subset=["symbol"])  # require symbol for GSEA
    out = agg[["symbol", "score"]].rename(columns={"score": "globalScore"}).sort_values(
        by="globalScore", ascending=False
    )
    return out.reset_index(drop=True)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare GSEA input (symbol, globalScore) for a disease")
    parser.add_argument("--gene_name", required=False, help="Primary target name (optional; validated if provided)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--disease_name", required=False, help="Disease name (requires data/disease_data.txt)")
    group.add_argument("--disease_id", required=False, help="Disease EFO ID (e.g., EFO_0004248)")
    parser.add_argument(
        "--output",
        dest="output_tsv",
        default=os.path.join("gsea", "input", "auto-input", "final_gsea_input.tsv"),
        help="Output TSV path (default: gsea/input/auto-input/final_gsea_input.tsv)",
    )
    parser.add_argument(
        "--datatype",
        dest="datatype",
        default=None,
        help="Optional association datatype filter, e.g. genetic_association",
    )
    # Data source configuration
    parser.add_argument("--associations_dir", default=None, help="Directory of association parquet parts")
    parser.add_argument("--associations_tsv", default=None, help="Single TSV of associations (fallback)")
    parser.add_argument("--target_dir", default=None, help="Directory of target parquet parts")
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()

    # Configure data sources from CLI (env fallbacks are used in loaders)
    configure_data_sources(
        associations_dir=args.associations_dir,
        associations_tsv=args.associations_tsv,
        target_dir=args.target_dir,
    )

    # Validate inputs and load
    load_data()
    # Validate gene if provided
    if args.gene_name:
        if GENE_NAME_TO_ID is None:
            _load_gene_maps()
        assert GENE_NAME_TO_ID is not None
        if args.gene_name not in GENE_NAME_TO_ID:
            raise ValueError(f"Gene name '{args.gene_name}' not found in gene map")

    disease_efo = resolve_disease_id(args.disease_name, args.disease_id)
    gsea_df = build_gsea_input_for_disease_id(disease_efo, datatype=args.datatype)
    os.makedirs(os.path.dirname(args.output_tsv) or ".", exist_ok=True)
    gsea_df.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
