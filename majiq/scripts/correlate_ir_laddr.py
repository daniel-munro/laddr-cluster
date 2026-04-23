#!/usr/bin/env python3
"""Correlate MAJIQ retained-intron PSI values with LaDDR phenotypes.

For each gene, this script compares retained-intron PSI vectors against LaDDR
phenotype vectors from the same gene across shared samples. It writes all tested
pairs with Pearson/Spearman correlations and a filtered table of high-confidence
Spearman associations.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


IR_METADATA = {
    "majiq_ir_id",
    "ec_idx",
    "seqid",
    "start",
    "end",
    "strand",
    "gene_id",
    "gene_id_base",
    "gene_name",
    "event_type",
    "is_denovo",
    "event_denovo",
    "ref_exon_start",
    "ref_exon_end",
    "other_exon_start",
    "other_exon_end",
}


def normalize_gene_id(value):
    return str(value).split(".", 1)[0]


def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    valid_idx = np.flatnonzero(valid)
    order = valid_idx[np.argsort(p[valid])]
    ranked = p[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out[order] = np.minimum(adjusted, 1.0)
    return out


def sample_columns(columns, metadata):
    return [col for col in columns if col not in metadata]


def read_groups(path):
    groups = pd.read_csv(path, sep="\t", names=["phenotype_id", "gene_id"])
    groups["gene_id"] = groups["gene_id"].map(normalize_gene_id)
    return dict(zip(groups["phenotype_id"], groups["gene_id"]))


def correlation_pair(x, y, min_shared):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    n = int(valid.sum())
    if n < min_shared:
        return n, np.nan, np.nan, np.nan, np.nan
    xv = x[valid]
    yv = y[valid]
    if np.nanstd(xv) == 0 or np.nanstd(yv) == 0:
        return n, np.nan, np.nan, np.nan, np.nan
    pearson_r, pearson_p = stats.pearsonr(xv, yv)
    spearman_r, spearman_p = stats.spearmanr(xv, yv)
    return n, pearson_r, pearson_p, spearman_r, spearman_p


def correlate(args):
    ir = pd.read_csv(args.ir, sep="\t")
    if "gene_id_base" not in ir.columns:
        ir["gene_id_base"] = ir["gene_id"].map(normalize_gene_id)
    ir_samples = sample_columns(ir.columns, IR_METADATA)

    laddr_header = pd.read_csv(args.laddr_bed, sep="\t", nrows=0)
    laddr_samples = [col for col in laddr_header.columns if col not in {"#chr", "chr", "start", "end", "phenotype_id"}]
    shared_samples = [sample for sample in ir_samples if sample in set(laddr_samples)]
    if len(shared_samples) < args.min_shared_samples:
        raise ValueError(
            f"Only {len(shared_samples)} shared samples between MAJIQ and LaDDR, "
            f"less than --min-shared-samples={args.min_shared_samples}"
        )

    ir_by_gene = {}
    for gene_id, gene_df in ir.groupby("gene_id_base", sort=False):
        ir_by_gene[gene_id] = (
            gene_df.reset_index(drop=True),
            gene_df[shared_samples].to_numpy(dtype=float),
        )
    ir_genes = set(ir_by_gene)

    pheno_to_gene = read_groups(args.laddr_groups)
    usecols = ["#chr", "start", "end", "phenotype_id", *shared_samples]
    records = []
    chunks = pd.read_csv(args.laddr_bed, sep="\t", usecols=usecols, chunksize=args.chunksize)
    for chunk in chunks:
        chunk["gene_id_base"] = chunk["phenotype_id"].map(pheno_to_gene)
        chunk = chunk.loc[chunk["gene_id_base"].isin(ir_genes)]
        if chunk.empty:
            continue
        for gene_id, laddr_gene in chunk.groupby("gene_id_base", sort=False):
            ir_gene, ir_values = ir_by_gene[gene_id]
            laddr_values = laddr_gene[shared_samples].to_numpy(dtype=float)
            for ir_idx, ir_row in ir_gene.iterrows():
                x = ir_values[ir_idx]
                for laddr_idx, (_, laddr_row) in enumerate(laddr_gene.iterrows()):
                    n, pearson_r, pearson_p, spearman_r, spearman_p = correlation_pair(
                        x, laddr_values[laddr_idx], args.min_shared_samples
                    )
                    records.append(
                        {
                            "majiq_ir_id": ir_row["majiq_ir_id"],
                            "ir_gene_id": gene_id,
                            "ir_seqid": ir_row["seqid"],
                            "ir_start": ir_row["start"],
                            "ir_end": ir_row["end"],
                            "ir_strand": ir_row["strand"],
                            "ir_is_denovo": ir_row.get("is_denovo", pd.NA),
                            "ir_event_denovo": ir_row.get("event_denovo", pd.NA),
                            "phenotype_id": laddr_row["phenotype_id"],
                            "laddr_gene_id": gene_id,
                            "n_shared_samples": n,
                            "pearson_r": pearson_r,
                            "pearson_p": pearson_p,
                            "spearman_r": spearman_r,
                            "spearman_p": spearman_p,
                        }
                    )

    out = pd.DataFrame.from_records(records)
    if out.empty:
        out = pd.DataFrame(
            columns=[
                "majiq_ir_id",
                "ir_gene_id",
                "ir_seqid",
                "ir_start",
                "ir_end",
                "ir_strand",
                "ir_is_denovo",
                "ir_event_denovo",
                "phenotype_id",
                "laddr_gene_id",
                "n_shared_samples",
                "pearson_r",
                "pearson_p",
                "spearman_r",
                "spearman_p",
                "pearson_fdr",
                "spearman_fdr",
                "best_laddr_match_for_ir",
            ]
        )
    else:
        out["pearson_fdr"] = bh_fdr(out["pearson_p"].to_numpy())
        out["spearman_fdr"] = bh_fdr(out["spearman_p"].to_numpy())
        out["abs_spearman_r"] = out["spearman_r"].abs()
        best = out.groupby("majiq_ir_id")["abs_spearman_r"].idxmax()
        out["best_laddr_match_for_ir"] = False
        out.loc[best.dropna().astype(int), "best_laddr_match_for_ir"] = True
        out = out.drop(columns=["abs_spearman_r"])

    Path(args.correlations_out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.correlations_out, sep="\t", index=False, compression="infer")
    high = out.loc[
        (out["spearman_r"].abs() >= args.high_abs_spearman)
        & (out["spearman_fdr"] <= args.high_spearman_fdr)
        & (out["n_shared_samples"] >= args.min_shared_samples)
    ].copy()
    high.to_csv(args.high_out, sep="\t", index=False, compression="infer")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Test same-gene correlations between retained-intron PSI values "
            "and LaDDR phenotype values across shared samples."
        )
    )
    parser.add_argument(
        "--ir",
        required=True,
        help="Retained-intron PSI table from extract_ir_psi.py.",
    )
    parser.add_argument(
        "--laddr-bed",
        required=True,
        help="LaDDR phenotype BED table with phenotype_id and sample columns.",
    )
    parser.add_argument(
        "--laddr-groups",
        required=True,
        help="Two-column phenotype-to-gene mapping file used to restrict tests to same-gene pairs.",
    )
    parser.add_argument(
        "--correlations-out",
        required=True,
        help="Output table containing every tested retained-intron/LaDDR phenotype pair.",
    )
    parser.add_argument(
        "--high-out",
        required=True,
        help="Output subset passing the high-correlation Spearman thresholds.",
    )
    parser.add_argument(
        "--min-shared-samples",
        type=int,
        default=100,
        help="Minimum number of samples with non-missing values in both vectors. Default: %(default)s.",
    )
    parser.add_argument(
        "--high-abs-spearman",
        type=float,
        default=0.5,
        help="Minimum absolute Spearman correlation for the high-correlation output. Default: %(default)s.",
    )
    parser.add_argument(
        "--high-spearman-fdr",
        type=float,
        default=0.05,
        help="Maximum Benjamini-Hochberg FDR for Spearman p-values in the high-correlation output. Default: %(default)s.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=5000,
        help="Number of LaDDR phenotype rows to process per pandas chunk. Default: %(default)s.",
    )
    args = parser.parse_args()
    correlate(args)


if __name__ == "__main__":
    main()
