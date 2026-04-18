import argparse
import sys
import time
from pathlib import Path

import pandas as pd


QTL_COLUMNS = ["tissue", "phenotype_id", "variant_id", "pip", "af", "cs_id"]
FINNGEN_SNP_COLUMNS = [
    "trait",
    "region",
    "v",
    "cs",
    "cs_specific_prob",
    "cs_99",
    "cs_specific_prob_99",
    "chromosome",
    "position",
    "allele1",
    "allele2",
    "maf",
    "beta",
    "p",
    "se",
    "most_severe",
    "gene_most_severe",
]
FINNGEN_CRED_COLUMNS = [
    "trait",
    "region",
    "cs",
    "cs_log10bf",
    "cs_avg_r2",
    "cs_min_r2",
    "low_purity",
    "cs_size",
    "good_cs",
    "cs_id",
    "v",
    "rsid",
    "p",
    "beta",
    "sd",
    "prob",
    "cs_specific_prob",
    "most_severe",
    "gene_most_severe",
]
COLOC_COLUMNS = [
    "tissue",
    "phenotype_id",
    "phenotype_class",
    "qtl_cs_id",
    "qtl_cs_size",
    "qtl_top_variant",
    "qtl_top_pip",
    "finngen_trait",
    "finngen_region",
    "finngen_cs",
    "finngen_cs_size",
    "finngen_top_variant",
    "finngen_top_prob",
    "cs_overlap",
    "clpp",
    "clpa",
    "top_shared_variant",
    "top_shared_qtl_pip",
    "top_shared_finngen_prob",
    "top_shared_clpp",
    "finngen_cs_log10bf",
    "finngen_cs_avg_r2",
    "finngen_cs_min_r2",
    "finngen_low_purity",
    "finngen_good_cs",
    "finngen_lead_variant",
    "finngen_lead_rsid",
    "finngen_lead_p",
    "finngen_lead_beta",
    "finngen_lead_sd",
    "finngen_lead_prob",
    "finngen_lead_cs_specific_prob",
    "finngen_lead_most_severe",
    "finngen_lead_gene_most_severe",
]


def format_duration(seconds: float) -> str:
    seconds = int(max(0, seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours:d}:{minutes:02d}:{seconds:02d}"
    return f"{minutes:d}:{seconds:02d}"


def log_progress(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def parse_bool(value: str) -> bool:
    if value.lower() in {"true", "1", "yes"}:
        return True
    if value.lower() in {"false", "0", "no"}:
        return False
    raise ValueError(f"Expected boolean value, got {value!r}.")


def normalize_qtl_variant(variant_id: pd.Series) -> pd.Series:
    parts = variant_id.str.rsplit("_", n=4, expand=True)
    if parts.shape[1] != 5:
        raise ValueError("QTL variant_id values do not match chr_pos_ref_alt_build format.")
    chrom = parts[0].str.removeprefix("chr")
    return chrom + ":" + parts[1] + ":" + parts[2] + ":" + parts[3]


def normalize_finngen_variant(df: pd.DataFrame) -> pd.Series:
    chrom = df["chromosome"].astype(str).str.removeprefix("chr")
    pos = df["position"].astype(str)
    allele1 = df["allele1"].astype(str)
    allele2 = df["allele2"].astype(str)
    return chrom + ":" + pos + ":" + allele1 + ":" + allele2


def phenotype_class(phenotype_id: pd.Series) -> pd.Series:
    return phenotype_id.str.split(":", n=1).str[0].where(
        phenotype_id.str.contains(":", regex=False),
        "ddp",
    )


def read_qtl_susie(path: str, phenotype_prefixes: list[str]) -> pd.DataFrame:
    qtl = pd.read_csv(path, sep="\t", usecols=QTL_COLUMNS)
    if phenotype_prefixes:
        qtl = qtl.loc[
            qtl["phenotype_id"].apply(
                lambda value: any(str(value).startswith(prefix) for prefix in phenotype_prefixes)
            )
        ].copy()
    if qtl.empty:
        raise ValueError("QTL SuSiE table is empty after phenotype filtering.")

    qtl["variant_key"] = normalize_qtl_variant(qtl["variant_id"].astype(str))
    qtl["phenotype_class"] = phenotype_class(qtl["phenotype_id"].astype(str))
    qtl = qtl.rename(columns={"cs_id": "qtl_cs_id", "pip": "qtl_pip", "variant_id": "qtl_variant_id"})
    dup_cols = ["tissue", "phenotype_id", "qtl_cs_id", "variant_key"]
    if qtl.duplicated(dup_cols).any():
        examples = qtl.loc[qtl.duplicated(dup_cols, keep=False), dup_cols].head(10)
        raise ValueError(f"Duplicate QTL variants within a credible set:\n{examples}")
    return qtl


def read_finngen_snps(paths: list[str], cs_col: str, prob_col: str) -> pd.DataFrame:
    frames = []
    for path in paths:
        df = pd.read_csv(path, sep="\t", usecols=FINNGEN_SNP_COLUMNS)
        df = df.loc[df[cs_col].notna() & (df[cs_col] != -1)].copy()
        if df.empty:
            continue
        df["variant_key"] = normalize_finngen_variant(df)
        df = df.rename(
            columns={
                "trait": "finngen_trait",
                "region": "finngen_region",
                cs_col: "finngen_cs",
                prob_col: "finngen_prob",
                "v": "finngen_variant_id",
                "maf": "finngen_maf",
                "beta": "finngen_beta",
                "p": "finngen_p",
                "se": "finngen_se",
                "most_severe": "finngen_most_severe",
                "gene_most_severe": "finngen_gene_most_severe",
            }
        )
        frames.append(df)
    if not frames:
        return pd.DataFrame(
            columns=[
                "finngen_trait",
                "finngen_region",
                "finngen_cs",
                "finngen_prob",
                "variant_key",
            ]
        )

    snps = pd.concat(frames, ignore_index=True)
    snps["finngen_cs"] = snps["finngen_cs"].astype(int)
    dup_cols = ["finngen_trait", "finngen_region", "finngen_cs", "variant_key"]
    if snps.duplicated(dup_cols).any():
        examples = snps.loc[snps.duplicated(dup_cols, keep=False), dup_cols].head(10)
        raise ValueError(f"Duplicate FinnGen variants within a credible set:\n{examples}")
    return snps


def read_finngen_cred(paths: list[str]) -> pd.DataFrame:
    frames = []
    for path in paths:
        df = pd.read_csv(path, sep="\t", usecols=FINNGEN_CRED_COLUMNS)
        frames.append(df)
    cred = pd.concat(frames, ignore_index=True)
    cred = cred.rename(
        columns={
            "trait": "finngen_trait",
            "region": "finngen_region",
            "cs": "finngen_cs",
            "cs_log10bf": "finngen_cs_log10bf",
            "cs_avg_r2": "finngen_cs_avg_r2",
            "cs_min_r2": "finngen_cs_min_r2",
            "low_purity": "finngen_low_purity",
            "good_cs": "finngen_good_cs",
            "v": "finngen_lead_variant",
            "rsid": "finngen_lead_rsid",
            "p": "finngen_lead_p",
            "beta": "finngen_lead_beta",
            "sd": "finngen_lead_sd",
            "prob": "finngen_lead_prob",
            "cs_specific_prob": "finngen_lead_cs_specific_prob",
            "most_severe": "finngen_lead_most_severe",
            "gene_most_severe": "finngen_lead_gene_most_severe",
        }
    )
    cred["finngen_cs"] = cred["finngen_cs"].astype(int)
    keep = [
        "finngen_trait",
        "finngen_region",
        "finngen_cs",
        "finngen_cs_log10bf",
        "finngen_cs_avg_r2",
        "finngen_cs_min_r2",
        "finngen_low_purity",
        "finngen_good_cs",
        "finngen_lead_variant",
        "finngen_lead_rsid",
        "finngen_lead_p",
        "finngen_lead_beta",
        "finngen_lead_sd",
        "finngen_lead_prob",
        "finngen_lead_cs_specific_prob",
        "finngen_lead_most_severe",
        "finngen_lead_gene_most_severe",
    ]
    cred = cred.loc[:, keep]
    dup_cols = ["finngen_trait", "finngen_region", "finngen_cs"]
    if cred.duplicated(dup_cols).any():
        examples = cred.loc[cred.duplicated(dup_cols, keep=False), dup_cols].head(10)
        raise ValueError(f"Duplicate FinnGen credible-set summaries:\n{examples}")
    return cred


def top_rows(df: pd.DataFrame, group_cols: list[str], value_col: str, prefix: str) -> pd.DataFrame:
    idx = df.groupby(group_cols, observed=True)[value_col].idxmax()
    top = df.loc[idx, group_cols + ["variant_key", value_col]].copy()
    return top.rename(
        columns={
            "variant_key": f"{prefix}_top_variant",
            value_col: f"{prefix}_top_pip" if prefix == "qtl" else f"{prefix}_top_prob",
        }
    )


def compute_qtl_metadata(qtl: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    qtl_group_cols = ["tissue", "phenotype_id", "phenotype_class", "qtl_cs_id"]
    qtl_sizes = (
        qtl.groupby(qtl_group_cols, observed=True)
        .size()
        .rename("qtl_cs_size")
        .reset_index()
    )
    qtl_top = top_rows(qtl, qtl_group_cols, "qtl_pip", "qtl")
    qtl_variant_keys = set(qtl["variant_key"])
    return qtl_sizes, qtl_top, qtl_variant_keys


def compute_coloc(
    qtl: pd.DataFrame,
    finngen: pd.DataFrame,
    cred: pd.DataFrame,
    qtl_sizes: pd.DataFrame,
    qtl_top: pd.DataFrame,
    qtl_variant_keys: set[str],
) -> pd.DataFrame:
    qtl_group_cols = ["tissue", "phenotype_id", "phenotype_class", "qtl_cs_id"]
    finngen_group_cols = ["finngen_trait", "finngen_region", "finngen_cs"]

    if finngen.empty:
        return pd.DataFrame(columns=COLOC_COLUMNS)

    finngen = finngen.loc[finngen["variant_key"].isin(qtl_variant_keys)].copy()
    if finngen.empty:
        return pd.DataFrame(columns=COLOC_COLUMNS)

    finngen_sizes = (
        finngen.groupby(finngen_group_cols, observed=True)
        .size()
        .rename("finngen_cs_size")
        .reset_index()
    )
    finngen_top = top_rows(finngen, finngen_group_cols, "finngen_prob", "finngen")

    shared_keys = set(finngen["variant_key"])
    qtl_overlap = qtl.loc[qtl["variant_key"].isin(shared_keys)].copy()
    if qtl_overlap.empty:
        return pd.DataFrame(columns=COLOC_COLUMNS)

    merged = qtl_overlap.merge(
        finngen,
        on="variant_key",
        how="inner",
        validate="many_to_many",
    )
    if merged.empty:
        return pd.DataFrame(columns=COLOC_COLUMNS)

    merged["variant_clpp"] = merged["qtl_pip"] * merged["finngen_prob"]
    merged["variant_clpa"] = merged[["qtl_pip", "finngen_prob"]].min(axis=1)
    group_cols = qtl_group_cols + finngen_group_cols

    stats = (
        merged.groupby(group_cols, observed=True)
        .agg(
            cs_overlap=("variant_key", "size"),
            clpp=("variant_clpp", "sum"),
            clpa=("variant_clpa", "sum"),
        )
        .reset_index()
    )

    top_shared_idx = merged.groupby(group_cols, observed=True)["variant_clpp"].idxmax()
    top_shared = merged.loc[
        top_shared_idx,
        group_cols + ["variant_key", "qtl_pip", "finngen_prob", "variant_clpp"],
    ].rename(
        columns={
            "variant_key": "top_shared_variant",
            "qtl_pip": "top_shared_qtl_pip",
            "finngen_prob": "top_shared_finngen_prob",
            "variant_clpp": "top_shared_clpp",
        }
    )

    out = stats.merge(qtl_sizes, on=qtl_group_cols, how="left", validate="many_to_one")
    out = out.merge(qtl_top, on=qtl_group_cols, how="left", validate="many_to_one")
    out = out.merge(finngen_sizes, on=finngen_group_cols, how="left", validate="many_to_one")
    out = out.merge(finngen_top, on=finngen_group_cols, how="left", validate="many_to_one")
    out = out.merge(top_shared, on=group_cols, how="left", validate="one_to_one")
    out = out.merge(cred, on=finngen_group_cols, how="left", validate="many_to_one")
    out = out.rename(columns={"qtl_top_pip": "qtl_top_pip"})
    return out.loc[:, COLOC_COLUMNS].sort_values(["clpp", "clpa"], ascending=False)


def write_table(df: pd.DataFrame, path: str, columns: list[str], compression: str | None = None) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if df.empty:
        df = pd.DataFrame(columns=columns)
    else:
        df = df.loc[:, columns]
    df.to_csv(path, sep="\t", index=False, compression=compression)


def summarize_by_trait(df: pd.DataFrame) -> pd.DataFrame:
    columns = ["finngen_trait", "n_coloc_pairs", "n_phenotypes", "n_tissues", "max_clpp", "max_clpa"]
    if df.empty:
        return pd.DataFrame(columns=columns)
    summary = (
        df.groupby("finngen_trait", observed=True)
        .agg(
            n_coloc_pairs=("clpp", "size"),
            n_phenotypes=("phenotype_id", "nunique"),
            n_tissues=("tissue", "nunique"),
            max_clpp=("clpp", "max"),
            max_clpa=("clpa", "max"),
        )
        .reset_index()
    )
    return summary.loc[:, columns]


def summarize_by_phenotype_class(df: pd.DataFrame) -> pd.DataFrame:
    columns = ["phenotype_class", "n_coloc_pairs", "n_phenotypes", "n_traits", "max_clpp", "max_clpa"]
    if df.empty:
        return pd.DataFrame(columns=columns)
    summary = (
        df.groupby("phenotype_class", observed=True)
        .agg(
            n_coloc_pairs=("clpp", "size"),
            n_phenotypes=("phenotype_id", "nunique"),
            n_traits=("finngen_trait", "nunique"),
            max_clpp=("clpp", "max"),
            max_clpa=("clpa", "max"),
        )
        .reset_index()
    )
    return summary.loc[:, columns]


def read_traits(path: str) -> list[str]:
    with open(path, "r") as fh:
        traits = [line.strip() for line in fh if line.strip()]
    if not traits:
        raise ValueError(f"No FinnGen traits listed in {path}.")
    if len(traits) != len(set(traits)):
        raise ValueError(f"Duplicate FinnGen traits listed in {path}.")
    return traits


def filter_coloc(
    df: pd.DataFrame,
    min_cs_overlap: int,
    require_good_cs: bool,
) -> pd.DataFrame:
    df = df.loc[df["cs_overlap"] >= min_cs_overlap].copy()
    if require_good_cs and not df.empty:
        df = df.loc[df["finngen_good_cs"] == True].copy()
    return df


def compute_for_trait(
    qtl: pd.DataFrame,
    qtl_sizes: pd.DataFrame,
    qtl_top: pd.DataFrame,
    qtl_variant_keys: set[str],
    snp_path: str,
    cred_path: str,
    min_cs_overlap: int,
    require_good_cs: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cred = read_finngen_cred([cred_path])

    finngen = read_finngen_snps([snp_path], "cs", "cs_specific_prob")
    primary = compute_coloc(qtl, finngen, cred, qtl_sizes, qtl_top, qtl_variant_keys)
    primary = filter_coloc(primary, min_cs_overlap, require_good_cs)

    finngen_99 = read_finngen_snps([snp_path], "cs_99", "cs_specific_prob_99")
    sensitivity = compute_coloc(qtl, finngen_99, cred, qtl_sizes, qtl_top, qtl_variant_keys)
    sensitivity = filter_coloc(sensitivity, min_cs_overlap, require_good_cs)
    return primary, sensitivity


def finalize_outputs(
    primary: pd.DataFrame,
    sensitivity: pd.DataFrame,
    out_primary: str,
    out_sensitivity: str,
    out_significant: str,
    out_summary_trait: str,
    out_summary_class: str,
    clpp_threshold: float,
    clpa_threshold: float,
) -> None:
    if not primary.empty:
        primary = primary.sort_values(["clpp", "clpa"], ascending=False)
    if not sensitivity.empty:
        sensitivity = sensitivity.sort_values(["clpp", "clpa"], ascending=False)

    significant = primary.loc[
        (primary["clpp"] >= clpp_threshold) | (primary["clpa"] >= clpa_threshold)
    ].copy()

    write_table(primary, out_primary, COLOC_COLUMNS, compression="gzip")
    write_table(sensitivity, out_sensitivity, COLOC_COLUMNS, compression="gzip")
    write_table(significant, out_significant, COLOC_COLUMNS, compression="gzip")
    write_table(
        summarize_by_trait(primary),
        out_summary_trait,
        ["finngen_trait", "n_coloc_pairs", "n_phenotypes", "n_tissues", "max_clpp", "max_clpa"],
    )
    write_table(
        summarize_by_phenotype_class(primary),
        out_summary_class,
        ["phenotype_class", "n_coloc_pairs", "n_phenotypes", "n_traits", "max_clpp", "max_clpa"],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qtl-susie", required=True)
    parser.add_argument("--traits")
    parser.add_argument("--finngen-summary-dir")
    parser.add_argument("--finngen-snps", nargs="+")
    parser.add_argument("--finngen-cred", nargs="+")
    parser.add_argument("--out-primary", required=True)
    parser.add_argument("--out-sensitivity", required=True)
    parser.add_argument("--out-significant")
    parser.add_argument("--out-summary-trait")
    parser.add_argument("--out-summary-class")
    parser.add_argument("--qtl-phenotype-prefixes", default="")
    parser.add_argument("--require-finngen-good-cs", default="true")
    parser.add_argument("--min-cs-overlap", type=int, default=1)
    parser.add_argument("--clpp-threshold", type=float, default=0.01)
    parser.add_argument("--clpa-threshold", type=float, default=0.5)
    parser.add_argument("--progress-every", type=int, default=1)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    phenotype_prefixes = [value for value in args.qtl_phenotype_prefixes.split(",") if value]
    require_good_cs = parse_bool(args.require_finngen_good_cs)
    start_time = time.monotonic()
    progress_every = max(1, args.progress_every)

    if not args.quiet:
        log_progress(f"[0:00] Loading QTL SuSiE table: {args.qtl_susie}")
    qtl = read_qtl_susie(args.qtl_susie, phenotype_prefixes)
    if not args.quiet:
        log_progress(
            f"[{format_duration(time.monotonic() - start_time)}] "
            f"Loaded {len(qtl):,} QTL CS variant rows across "
            f"{qtl['phenotype_id'].nunique():,} phenotypes."
        )
        log_progress(f"[{format_duration(time.monotonic() - start_time)}] Precomputing QTL CS metadata.")
    qtl_sizes, qtl_top, qtl_variant_keys = compute_qtl_metadata(qtl)
    if not args.quiet:
        log_progress(
            f"[{format_duration(time.monotonic() - start_time)}] "
            f"QTL metadata ready for {len(qtl_sizes):,} credible sets."
        )

    if args.traits or args.finngen_summary_dir:
        if not args.traits or not args.finngen_summary_dir:
            raise ValueError("--traits and --finngen-summary-dir must be used together.")
        if not args.out_significant or not args.out_summary_trait or not args.out_summary_class:
            raise ValueError("Summary output paths are required when using --traits.")

        primary_frames = []
        sensitivity_frames = []
        summary_dir = Path(args.finngen_summary_dir)
        traits = read_traits(args.traits)
        if not args.quiet:
            log_progress(
                f"[{format_duration(time.monotonic() - start_time)}] "
                f"Processing {len(traits):,} FinnGen traits from {summary_dir}."
            )
        for i, trait in enumerate(traits, start=1):
            trait_start = time.monotonic()
            snp_path = summary_dir / f"finngen_R12_{trait}.SUSIE_extend.snp.filter.tsv"
            cred_path = summary_dir / f"finngen_R12_{trait}.SUSIE_extend.cred.summary.tsv"
            primary, sensitivity = compute_for_trait(
                qtl,
                qtl_sizes,
                qtl_top,
                qtl_variant_keys,
                str(snp_path),
                str(cred_path),
                args.min_cs_overlap,
                require_good_cs,
            )
            if not primary.empty:
                primary_frames.append(primary)
            if not sensitivity.empty:
                sensitivity_frames.append(sensitivity)
            if not args.quiet and (i == 1 or i == len(traits) or i % progress_every == 0):
                elapsed = time.monotonic() - start_time
                trait_elapsed = time.monotonic() - trait_start
                traits_per_second = i / elapsed if elapsed > 0 else 0
                remaining = (len(traits) - i) / traits_per_second if traits_per_second else 0
                log_progress(
                    f"[{format_duration(elapsed)}] "
                    f"{i:,}/{len(traits):,} traits "
                    f"({100 * i / len(traits):.1f}%) "
                    f"trait={trait} "
                    f"primary_rows={len(primary):,} "
                    f"sensitivity_rows={len(sensitivity):,} "
                    f"trait_time={format_duration(trait_elapsed)} "
                    f"eta={format_duration(remaining)}"
                )

        primary = (
            pd.concat(primary_frames, ignore_index=True)
            if primary_frames
            else pd.DataFrame(columns=COLOC_COLUMNS)
        )
        sensitivity = (
            pd.concat(sensitivity_frames, ignore_index=True)
            if sensitivity_frames
            else pd.DataFrame(columns=COLOC_COLUMNS)
        )
        if not args.quiet:
            log_progress(
                f"[{format_duration(time.monotonic() - start_time)}] "
                f"Writing outputs: primary_rows={len(primary):,}, "
                f"sensitivity_rows={len(sensitivity):,}."
            )
        finalize_outputs(
            primary,
            sensitivity,
            args.out_primary,
            args.out_sensitivity,
            args.out_significant,
            args.out_summary_trait,
            args.out_summary_class,
            args.clpp_threshold,
            args.clpa_threshold,
        )
        if not args.quiet:
            log_progress(f"[{format_duration(time.monotonic() - start_time)}] Done.")
        return

    if not args.finngen_snps or not args.finngen_cred:
        raise ValueError("Either --traits/--finngen-summary-dir or --finngen-snps/--finngen-cred is required.")

    if not args.quiet:
        log_progress(f"[{format_duration(time.monotonic() - start_time)}] Processing single FinnGen trait.")
    primary, sensitivity = compute_for_trait(
        qtl,
        qtl_sizes,
        qtl_top,
        qtl_variant_keys,
        args.finngen_snps[0],
        args.finngen_cred[0],
        args.min_cs_overlap,
        require_good_cs,
    )
    if not args.quiet:
        log_progress(
            f"[{format_duration(time.monotonic() - start_time)}] "
            f"Writing single-trait outputs: primary_rows={len(primary):,}, "
            f"sensitivity_rows={len(sensitivity):,}."
        )
    write_table(primary, args.out_primary, COLOC_COLUMNS, compression="gzip")
    write_table(sensitivity, args.out_sensitivity, COLOC_COLUMNS, compression="gzip")
    if not args.quiet:
        log_progress(f"[{format_duration(time.monotonic() - start_time)}] Done.")


if __name__ == "__main__":
    main()
