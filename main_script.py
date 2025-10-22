import argparse
import datetime as _dt
import os

from functions import sid
import gsea.run_gsea as run_gsea


def _default_experiment_dir(gene_name: str, disease_name: str) -> str:
    today = _dt.date.today().isoformat()
    safe_disease = disease_name.replace(" ", "_")
    return os.path.join(
        "experiments", f"{safe_disease}_{gene_name}_{today}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "End-to-end pipeline: build GSEA input from associations, run GSEA, "
            "and optionally compute Reactome-based overlap scores."
        )
    )
    parser.add_argument("--gene_name", type=str, required=False, help="Primary target name (optional)")
    parser.add_argument("--disease_name", type=str, required=False, help="Disease name (optional)")
    parser.add_argument("--disease_id", type=str, required=False, help="Disease EFO ID (optional; overrides name)")
    parser.add_argument(
        "--experiment_dir",
        type=str,
        default=None,
        help="Output folder (defaults to experiments/<disease>_<gene>_<date>)",
    )
    parser.add_argument(
        "--gmt_file",
        type=str,
        required=True,
        help="Path to Reactome GMT file (e.g., gsea/Reactome_2025/<file>.gmt)",
    )
    parser.add_argument(
        "--datatype",
        type=str,
        default=None,
        help="Optional association datatype filter, e.g. genetic_association",
    )
    parser.add_argument(
        "--pval_threshold",
        type=float,
        default=None,
        help="Optional p-value threshold; if set, keep rows with pval <= threshold",
    )
    parser.add_argument(
        "--fdr_threshold",
        type=float,
        default=None,
        help="Optional FDR threshold; if set and qval exists, keep qval <= threshold",
    )
    parser.add_argument(
        "--nes_positive",
        action="store_true",
        help="If set, also require NES > 0 when filtering IDs",
    )
    parser.add_argument(
        "--auto_fetch",
        action="store_true",
        help="Auto-download missing Reactome mapping/interactions (default filenames only)",
    )
    # Data source paths forwarded to sid
    parser.add_argument("--associations_dir", type=str, default=None, help="Directory of association parquet parts")
    parser.add_argument("--associations_tsv", type=str, default=None, help="Single TSV of associations (fallback)")
    parser.add_argument("--target_dir", type=str, default=None, help="Directory of target parquet parts")
    parser.add_argument(
        "--run_reactome",
        action="store_true",
        help="If set, run Reactome overlap scoring after GSEA (run.py)",
    )
    parser.add_argument(
        "--pathway_mapping_file",
        type=str,
        default=None,
        help="Reactome Ensembl2Reactome mapping file (required if --run_reactome)",
    )
    parser.add_argument(
        "--interactions_file",
        type=str,
        default=None,
        help="Reactome interactions file (required if --run_reactome)",
    )
    return parser.parse_args()


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def main() -> None:
    args = parse_args()

    exp_dir = args.experiment_dir or _default_experiment_dir(args.gene_name, args.disease_name)
    ensure_dir(exp_dir)

    # 1) Configure data sources and build GSEA input
    sid.configure_data_sources(
        associations_dir=args.associations_dir,
        associations_tsv=args.associations_tsv,
        target_dir=args.target_dir,
    )
    sid.load_data()
    # Optional gene validation
    if args.gene_name:
        try:
            sid.validate_input(args.gene_name, args.disease_name or next(iter({})))
        except Exception:
            # Validate gene name only if provided in map
            if sid.GENE_NAME_TO_ID is None or args.gene_name not in sid.GENE_NAME_TO_ID:
                raise
    # Resolve disease EFO ID and build GSEA input
    disease_efo = sid.resolve_disease_id(args.disease_name, args.disease_id)
    gsea_input_df = sid.build_gsea_input_for_disease_id(disease_efo, datatype=args.datatype)
    gsea_input_tsv = os.path.join(exp_dir, "final_gsea_input.tsv")
    gsea_input_df.to_csv(gsea_input_tsv, sep="\t", index=False)

    # 2) Run GSEA
    gsea_results_tsv = os.path.join(exp_dir, "gsea_results.tsv")
    res_df = run_gsea.run_gsea_from_tsv(
        input_tsv=gsea_input_tsv,
        gmt_file=args.gmt_file,
        output_tsv=gsea_results_tsv,
    )

    # 3) Filter pathway IDs by p-value
    gsea_ids_tsv = os.path.join(exp_dir, "gsea_ids.tsv")
    run_gsea.filter_pathway_ids_by_pvalue(
        gsea_results=res_df,
        pval_threshold=args.pval_threshold,
        output_ids_tsv=gsea_ids_tsv,
        fdr_threshold=args.fdr_threshold,
        nes_positive=args.nes_positive,
    )

    # 4) Optionally run Reactome overlap scoring without symlink (pass IDs directly)
    if args.run_reactome:
        if not args.pathway_mapping_file or not args.interactions_file:
            raise ValueError("--pathway_mapping_file and --interactions_file are required when --run_reactome is set")

        # Defer to run.py for scoring; pass IDs file directly and redirect stdout
        scores_tsv = os.path.join(exp_dir, "scores.tsv")
        cmd = (
            f"python run.py --pathway_mapping_file {args.pathway_mapping_file} "
            f"--interactions_file {args.interactions_file} --target {args.gene_name or ''} "
            f"--gsea_ids_file {gsea_ids_tsv} "
            f"{'--auto_fetch ' if args.auto_fetch else ''}> {scores_tsv}"
        )
        # Use os.system for simplicity to preserve run.py behavior
        rc = os.system(cmd)
        if rc != 0:
            raise RuntimeError("run.py failed; see console output for details")

    print(f"Experiment directory: {exp_dir}")
    print(f" - GSEA input: {gsea_input_tsv}")
    print(f" - GSEA results: {gsea_results_tsv}")
    print(f" - GSEA pathway IDs: {gsea_ids_tsv}")
    if args.run_reactome:
        print(f" - Reactome scores: {os.path.join(exp_dir, 'scores.tsv')}")


if __name__ == "__main__":
    main()

