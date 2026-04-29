#!/usr/bin/env python
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RUNS = ["unstranded", "stranded", "stranded_swapped"]
DATASET = "stranded_test"
OUTDIR = ROOT / "summary"


def load_samples(run: str) -> list[str]:
    samples_file = ROOT / "runs" / run / "covg_norm" / DATASET / "samples.txt"
    with open(samples_file) as handle:
        return [line.strip() for line in handle if line.strip()]


def load_phenotypes(run: str) -> pd.DataFrame:
    path = ROOT / "runs" / run / "phenotypes" / f"latent_phenos.{DATASET}.tsv.gz"
    df = pd.read_csv(path, sep="\t")
    df["phenotype_id"] = df["gene_id"].astype(str) + "__" + df["PC"].astype(str)
    return df


def batch_ids(run: str) -> list[int]:
    covg_dir = ROOT / "runs" / run / "covg_norm" / DATASET
    return sorted(int(path.stem.split("_")[1]) for path in covg_dir.glob("batch_*.npy"))


def summarize_coverage(run: str) -> tuple[pd.DataFrame, dict]:
    samples = load_samples(run)
    totals_log2 = np.zeros(len(samples))
    totals_norm = np.zeros(len(samples))
    n_bins = 0
    for batch in batch_ids(run):
        mat = np.load(ROOT / "runs" / run / "covg_norm" / DATASET / f"batch_{batch}.npy")
        totals_log2 += mat.sum(axis=0)
        totals_norm += np.maximum(np.exp2(mat) - 8, 0).sum(axis=0)
        n_bins += mat.shape[0]

    by_sample = pd.DataFrame({
        "analysis": run,
        "sample": samples,
        "n_bins": n_bins,
        "total_log2_covg_norm": totals_log2,
        "mean_log2_covg_norm": totals_log2 / n_bins,
        "approx_total_norm_covg": totals_norm,
        "approx_mean_norm_covg": totals_norm / n_bins,
    })
    aggregate = {
        "n_bins": n_bins,
        "total_log2_covg_norm": totals_log2.sum(),
        "mean_log2_covg_norm": totals_log2.sum() / (n_bins * len(samples)),
        "approx_total_norm_covg": totals_norm.sum(),
        "approx_mean_norm_covg": totals_norm.sum() / (n_bins * len(samples)),
    }
    return by_sample, aggregate


def summarize_run(run: str, phenos: pd.DataFrame, covg_summary: dict) -> dict:
    genes = pd.read_csv(ROOT / "runs" / run / "info" / "genes.tsv", sep="\t")
    median_coverage = float((ROOT / "runs" / run / "info" / "median_coverage.txt").read_text())
    var_per_bin = float((ROOT / "runs" / run / "info" / "var_per_bin.txt").read_text())
    pcs_per_gene = phenos.groupby("gene_id")["PC"].nunique()
    return {
        "analysis": run,
        "n_samples": len(load_samples(run)),
        "n_batches": int((ROOT / "runs" / run / "info" / "n_batches.txt").read_text()),
        "n_genes_in_annotation": genes.shape[0],
        "n_phenotype_rows": phenos.shape[0],
        "n_phenotype_genes": phenos["gene_id"].nunique(),
        "median_pcs_per_gene": pcs_per_gene.median(),
        "max_pcs_per_gene": pcs_per_gene.max(),
        "median_bigwig_sumdata": median_coverage,
        "var_per_bin": var_per_bin,
        **covg_summary,
    }


def summarize_overlaps(phenos_by_run: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    overlap_rows = []
    corr_rows = []
    for left, right in combinations(RUNS, 2):
        left_df = phenos_by_run[left].set_index("phenotype_id")
        right_df = phenos_by_run[right].set_index("phenotype_id")
        common = left_df.index.intersection(right_df.index)
        left_set = set(left_df.index)
        right_set = set(right_df.index)
        overlap_rows.append({
            "analysis_a": left,
            "analysis_b": right,
            "n_phenotypes_a": len(left_set),
            "n_phenotypes_b": len(right_set),
            "n_common_phenotypes": len(common),
            "n_a_only": len(left_set - right_set),
            "n_b_only": len(right_set - left_set),
            "jaccard": len(common) / len(left_set | right_set),
        })

        sample_cols = [col for col in left_df.columns if col.startswith("SRR")]
        per_pheno_corrs = []
        for phenotype_id in common:
            a = left_df.loc[phenotype_id, sample_cols].to_numpy(dtype=float)
            b = right_df.loc[phenotype_id, sample_cols].to_numpy(dtype=float)
            if np.std(a) > 0 and np.std(b) > 0:
                per_pheno_corrs.append(np.corrcoef(a, b)[0, 1])

        all_a = left_df.loc[common, sample_cols].to_numpy(dtype=float).ravel()
        all_b = right_df.loc[common, sample_cols].to_numpy(dtype=float).ravel()
        all_values_corr = np.corrcoef(all_a, all_b)[0, 1] if len(all_a) > 1 else np.nan
        corr_rows.append({
            "analysis_a": left,
            "analysis_b": right,
            "n_common_phenotypes": len(common),
            "all_values_corr": all_values_corr,
            "median_matching_phenotype_corr": np.median(per_pheno_corrs) if per_pheno_corrs else np.nan,
            "median_abs_matching_phenotype_corr": np.median(np.abs(per_pheno_corrs)) if per_pheno_corrs else np.nan,
        })
    return pd.DataFrame(overlap_rows), pd.DataFrame(corr_rows)


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    phenos_by_run = {run: load_phenotypes(run) for run in RUNS}

    coverage_by_sample = []
    run_rows = []
    for run in RUNS:
        by_sample, covg_summary = summarize_coverage(run)
        coverage_by_sample.append(by_sample)
        run_rows.append(summarize_run(run, phenos_by_run[run], covg_summary))

    pd.DataFrame(run_rows).to_csv(OUTDIR / "analysis_summary.tsv", sep="\t", index=False)
    pd.concat(coverage_by_sample, ignore_index=True).to_csv(OUTDIR / "coverage_by_sample.tsv", sep="\t", index=False)

    overlap, corr = summarize_overlaps(phenos_by_run)
    overlap.to_csv(OUTDIR / "phenotype_overlap.tsv", sep="\t", index=False)
    corr.to_csv(OUTDIR / "phenotype_correlations.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
