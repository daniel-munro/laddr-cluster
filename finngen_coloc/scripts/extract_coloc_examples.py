import argparse
import gzip
import io
import subprocess
from pathlib import Path

import pandas as pd


QTL_COLUMNS = ["tissue", "phenotype_id", "variant_id", "pip", "af", "cs_id"]
COLOC_KEY_COLUMNS = ["tissue", "phenotype_id", "qtl_cs_id", "finngen_trait", "finngen_region", "finngen_cs"]
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


def normalize_qtl_variant(variant_id: pd.Series) -> pd.Series:
    parts = variant_id.astype(str).str.rsplit("_", n=4, expand=True)
    if parts.shape[1] != 5:
        raise ValueError("QTL variant_id values do not match chr_pos_ref_alt_build format.")
    chrom = parts[0].str.removeprefix("chr")
    return chrom + ":" + parts[1] + ":" + parts[2] + ":" + parts[3]


def normalize_finngen_variant(df: pd.DataFrame) -> pd.Series:
    chrom = df["chromosome"].astype(str).str.removeprefix("chr")
    return chrom + ":" + df["position"].astype(str) + ":" + df["allele1"].astype(str) + ":" + df["allele2"].astype(str)


def normalize_sumstat_variant(df: pd.DataFrame) -> pd.Series:
    chrom = df["#chrom"].astype(str).str.removeprefix("chr")
    return chrom + ":" + df["pos"].astype(str) + ":" + df["ref"].astype(str) + ":" + df["alt"].astype(str)


def read_examples(path: str) -> pd.DataFrame:
    examples = pd.read_csv(path, sep="\t", dtype=str)
    if "finngen_trait" not in examples.columns:
        examples = examples.rename(columns={"trait": "finngen_trait"})
    examples["example_id"] = [
        f"example_{i}" for i in range(1, len(examples) + 1)
    ]
    examples["qtl_cs_id"] = examples["qtl_cs_id"].astype(str)
    examples["finngen_cs"] = examples["finngen_cs"].astype(str)
    return examples


def resolve_examples_with_coloc(examples: pd.DataFrame, coloc_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    coloc = pd.read_csv(coloc_path, sep="\t", dtype=str)
    coloc["qtl_cs_id"] = coloc["qtl_cs_id"].astype(str)
    coloc["finngen_cs"] = coloc["finngen_cs"].astype(str)

    resolved = examples.copy()
    for idx, example in resolved.iterrows():
        key_mask = (
            (coloc["tissue"] == example["tissue"])
            & (coloc["qtl_cs_id"] == example["qtl_cs_id"])
            & (coloc["finngen_trait"] == example["finngen_trait"])
            & (coloc["finngen_region"] == example["finngen_region"])
            & (coloc["finngen_cs"] == example["finngen_cs"])
        )
        exact = coloc.loc[key_mask & (coloc["phenotype_id"] == example["phenotype_id"]), "phenotype_id"]
        if len(exact) == 1:
            resolved.loc[idx, "phenotype_id"] = exact.iloc[0]
            continue

        suffix = ":" + example["phenotype_id"]
        suffix_matches = coloc.loc[key_mask & coloc["phenotype_id"].str.endswith(suffix), "phenotype_id"]
        if len(suffix_matches) != 1:
            raise ValueError(
                f"Could not uniquely resolve phenotype_id {example['phenotype_id']!r} "
                f"for {example['example_id']} in coloc table; found {len(suffix_matches)} suffix matches."
            )
        resolved.loc[idx, "phenotype_id"] = suffix_matches.iloc[0]

    annotated = resolved.merge(
        coloc,
        on=COLOC_KEY_COLUMNS,
        how="left",
        validate="one_to_one",
        suffixes=("", "_coloc"),
    )
    missing = annotated["clpp"].isna()
    if missing.any():
        missing_keys = annotated.loc[missing, ["example_id"] + COLOC_KEY_COLUMNS]
        raise ValueError(f"Examples not found in coloc table:\n{missing_keys}")
    return resolved, annotated


def read_qtl_cs(path: str, examples: pd.DataFrame) -> pd.DataFrame:
    qtl = pd.read_csv(path, sep="\t", usecols=QTL_COLUMNS, dtype={"cs_id": str})
    keys = examples.loc[:, ["example_id", "tissue", "phenotype_id", "qtl_cs_id"]].copy()
    qtl = qtl.rename(columns={"cs_id": "qtl_cs_id", "pip": "qtl_pip", "af": "qtl_af"})
    qtl["qtl_cs_id"] = qtl["qtl_cs_id"].astype(str)
    qtl = qtl.merge(keys, on=["tissue", "phenotype_id", "qtl_cs_id"], how="inner", validate="many_to_many")
    if qtl.empty:
        raise ValueError("No QTL credible-set variants matched the selected examples.")
    qtl["variant_key"] = normalize_qtl_variant(qtl["variant_id"])
    return qtl


def read_qtl_nominal(path: Path, example_chrom: str) -> pd.DataFrame:
    nominal = pd.read_parquet(path)
    if "phenotype_id" not in nominal.columns:
        nominal = nominal.reset_index()
    if "phenotype_id" not in nominal.columns:
        raise ValueError(f"No phenotype_id column found in {path}.")
    nominal["variant_key"] = normalize_qtl_variant(nominal["variant_id"])
    nominal["variant_chrom"] = example_chrom
    return nominal


def collect_qtl_region(examples: pd.DataFrame, qtl_cs: pd.DataFrame, nominal_base: str) -> pd.DataFrame:
    frames = []
    for _, example in examples.iterrows():
        chrom = example["finngen_region"].split(":", 1)[0].removeprefix("chr")
        nominal_path = (
            Path(nominal_base)
            / example["tissue"]
            / f"{example['tissue']}_coloc_examples.cis_qtl_pairs.chr{chrom}.parquet"
        )
        nominal = read_qtl_nominal(nominal_path, chrom)
        nominal = nominal.loc[nominal["phenotype_id"] == example["example_id"]].copy()
        if nominal.empty:
            raise ValueError(f"No nominal QTL rows found for {example['example_id']} in {nominal_path}.")

        cs = qtl_cs.loc[qtl_cs["example_id"] == example["example_id"], ["variant_key", "qtl_pip", "qtl_af"]]
        nominal = nominal.merge(cs, on="variant_key", how="left", validate="many_to_one")
        nominal["in_qtl_cs"] = nominal["qtl_pip"].notna()
        nominal = nominal.rename(columns={"phenotype_id": "tensorqtl_phenotype_id"})
        nominal.insert(0, "example_id", example["example_id"])
        nominal.insert(1, "tissue", example["tissue"])
        nominal.insert(2, "qtl_cs_id", example["qtl_cs_id"])
        nominal.insert(3, "phenotype_id", example["phenotype_id"])
        frames.append(nominal)
    return pd.concat(frames, ignore_index=True)


def read_finngen_cs(summary_dir: str, examples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    snp_frames = []
    cred_frames = []
    for trait in examples["finngen_trait"].drop_duplicates():
        snp_path = Path(summary_dir) / f"finngen_R12_{trait}.SUSIE_extend.snp.filter.tsv"
        cred_path = Path(summary_dir) / f"finngen_R12_{trait}.SUSIE_extend.cred.summary.tsv"
        snps = pd.read_csv(snp_path, sep="\t", usecols=FINNGEN_SNP_COLUMNS, dtype={"cs": str})
        cred = pd.read_csv(cred_path, sep="\t", usecols=FINNGEN_CRED_COLUMNS, dtype={"cs": str})
        snp_frames.append(snps)
        cred_frames.append(cred)

    snps = pd.concat(snp_frames, ignore_index=True)
    snps = snps.rename(
        columns={
            "trait": "finngen_trait",
            "region": "finngen_region",
            "cs": "finngen_cs",
            "v": "finngen_variant_id",
            "cs_specific_prob": "finngen_cs_specific_prob",
            "cs_specific_prob_99": "finngen_cs_specific_prob_99",
            "maf": "finngen_maf",
            "beta": "finngen_finemap_beta",
            "p": "finngen_finemap_p",
            "se": "finngen_finemap_se",
            "most_severe": "finngen_most_severe",
            "gene_most_severe": "finngen_gene_most_severe",
        }
    )
    snps["variant_key"] = normalize_finngen_variant(snps)
    snps["finngen_cs"] = snps["finngen_cs"].astype(str)

    cred = pd.concat(cred_frames, ignore_index=True)
    cred = cred.rename(
        columns={
            "trait": "finngen_trait",
            "region": "finngen_region",
            "cs": "finngen_cs",
            "cs_id": "finngen_cs_id",
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
    cred["finngen_cs"] = cred["finngen_cs"].astype(str)
    return snps, cred


def read_tabix_region(path: Path, region: str) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    query_region = region.removeprefix("chr")
    result = subprocess.run(
        ["tabix", str(path), query_region],
        check=True,
        text=True,
        capture_output=True,
    )
    if not result.stdout.strip():
        raise ValueError(f"No FinnGen summary-stat rows found for {path}:{query_region}.")
    return pd.read_csv(io.StringIO(result.stdout), sep="\t", names=header)


def collect_finngen_region(examples: pd.DataFrame, summary_dir: str, sumstats_dir: str) -> pd.DataFrame:
    cs_snps, cred = read_finngen_cs(summary_dir, examples)
    frames = []
    for _, example in examples.iterrows():
        stats_path = Path(sumstats_dir) / f"finngen_R12_{example['finngen_trait']}.gz"
        region_stats = read_tabix_region(stats_path, example["finngen_region"])
        region_stats["variant_key"] = normalize_sumstat_variant(region_stats)

        cs = cs_snps.loc[
            (cs_snps["finngen_trait"] == example["finngen_trait"])
            & (cs_snps["finngen_region"] == example["finngen_region"])
            & (cs_snps["finngen_cs"] == example["finngen_cs"])
        ].copy()
        if cs.empty:
            raise ValueError(f"No FinnGen CS SNPs found for {example['example_id']}.")

        cs_keep = [
            "variant_key",
            "finngen_cs",
            "finngen_variant_id",
            "finngen_cs_specific_prob",
            "cs_99",
            "finngen_cs_specific_prob_99",
            "finngen_maf",
            "finngen_finemap_beta",
            "finngen_finemap_p",
            "finngen_finemap_se",
            "finngen_most_severe",
            "finngen_gene_most_severe",
        ]
        region_stats = region_stats.merge(cs.loc[:, cs_keep], on="variant_key", how="left", validate="many_to_one")
        region_stats["in_finngen_cs"] = region_stats["finngen_cs_specific_prob"].notna()
        region_stats.insert(0, "example_id", example["example_id"])
        region_stats.insert(1, "finngen_trait", example["finngen_trait"])
        region_stats.insert(2, "finngen_region", example["finngen_region"])
        frames.append(region_stats)

    finngen_region = pd.concat(frames, ignore_index=True)
    return finngen_region, cred


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--examples", required=True)
    parser.add_argument("--coloc", required=True)
    parser.add_argument("--qtl-susie", required=True)
    parser.add_argument("--qtl-nominal-base", required=True)
    parser.add_argument("--finngen-summary-dir", required=True)
    parser.add_argument("--finngen-sumstats-dir", required=True)
    parser.add_argument("--out-annotated", required=True)
    parser.add_argument("--out-qtl-region", required=True)
    parser.add_argument("--out-finngen-region", required=True)
    args = parser.parse_args()

    examples = read_examples(args.examples)
    examples, annotated = resolve_examples_with_coloc(examples, args.coloc)
    qtl_cs = read_qtl_cs(args.qtl_susie, examples)
    qtl_region = collect_qtl_region(examples, qtl_cs, args.qtl_nominal_base)
    finngen_region, _ = collect_finngen_region(
        examples,
        args.finngen_summary_dir,
        args.finngen_sumstats_dir,
    )

    Path(args.out_annotated).parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(args.out_annotated, sep="\t", index=False)
    qtl_region.to_csv(args.out_qtl_region, sep="\t", index=False, compression="gzip")
    finngen_region.to_csv(args.out_finngen_region, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()
