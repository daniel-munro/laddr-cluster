from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig


def sample_to_individual(sample_id: str) -> str:
    return "-".join(sample_id.split("-")[:2])


def make_fixed_bins(gene: pd.Series, n_bins: int) -> pd.DataFrame:
    edges = np.linspace(int(gene["window_start"]), int(gene["window_end"]), n_bins + 1)
    starts = np.floor(edges[:-1]).astype(int)
    ends = np.floor(edges[1:]).astype(int)
    ends[-1] = int(gene["window_end"])

    return pd.DataFrame(
        {
            "bin_id": [f"bin_{i:04d}" for i in range(n_bins)],
            "seqname": gene["seqname"],
            "start": starts,
            "end": ends,
            "gene_id": gene["gene_id"],
        }
    )


def summarize_bigwig(path: Path, chrom: str, bins: pd.DataFrame) -> list[float]:
    with pyBigWig.open(str(path)) as bw:
        values = [
            bw.stats(chrom, int(row.start), int(row.end), type="mean", exact=True)[0]
            for row in bins.itertuples(index=False)
        ]
    return [0.0 if value is None else value for value in values]


def split_phenotype_id(phenotype_id: str) -> tuple[str, str]:
    phenotype_id = phenotype_id.split(":", 1)[-1]
    gene_id, pc = phenotype_id.split("__", 1)
    return gene_id, pc


def read_laddr_phenotype_row(path: Path, phenotype_id: str) -> pd.Series:
    gene_id, pc = split_phenotype_id(phenotype_id)
    phenotypes = pd.read_csv(path, sep="\t")
    matches = phenotypes.loc[
        (phenotypes["gene_id"] == gene_id)
        & (phenotypes["PC"] == pc)
    ]

    if matches.shape[0] != 1:
        raise ValueError(
            f"Expected exactly one LaDDR phenotype row matching "
            f"gene_id={gene_id!r}, PC={pc!r} in {path}; found {matches.shape[0]}."
        )

    row = matches.iloc[0]
    return pd.Series(
        {
            sample_id: float(value)
            for sample_id, value in row.items()
            if sample_id not in {"gene_id", "PC"}
        },
        name=phenotype_id,
    )


def make_sample_manifest(
    samples_path: Path,
    coverage_manifest_path: Path,
    coverage_dir: Path,
    tissue: str,
) -> pd.DataFrame:
    samples = pd.read_csv(samples_path, sep="\t", header=None, names=["sample_id"])
    samples["individual_id"] = samples["sample_id"].map(sample_to_individual)

    manifest = pd.read_csv(
        coverage_manifest_path,
        sep="\t",
        header=None,
        names=["tissue", "sample_id", "path"],
    )
    sample_manifest = samples.merge(
        manifest.loc[manifest["tissue"] == tissue, ["sample_id", "path"]],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
    sample_manifest["path"] = sample_manifest["path"].map(lambda path: coverage_dir / path)
    return sample_manifest


def assign_deciles(phenotype: pd.Series) -> pd.Series:
    ranked = phenotype.rank(method="first")
    deciles = pd.qcut(ranked, 10, labels=False) + 1
    return deciles.astype(int)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--phenotype-id", required=True)
    parser.add_argument("--genes", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--coverage-manifest", required=True, type=Path)
    parser.add_argument("--coverage-dir", required=True, type=Path)
    parser.add_argument("--phenotypes", required=True, type=Path)
    parser.add_argument("--n-bins", required=True, type=int)
    parser.add_argument("--out-bins", required=True, type=Path)
    parser.add_argument("--out-decile-median-coverage", required=True, type=Path)
    args = parser.parse_args()

    genes = pd.read_csv(args.genes, sep="\t")
    gene = genes.loc[genes["gene_id"] == args.gene_id].iloc[0]
    bins = make_fixed_bins(gene, args.n_bins)
    bins.insert(1, "tissue", args.tissue)

    sample_manifest = make_sample_manifest(
        samples_path=args.samples,
        coverage_manifest_path=args.coverage_manifest,
        coverage_dir=args.coverage_dir,
        tissue=args.tissue,
    )
    phenotype = read_laddr_phenotype_row(args.phenotypes, args.phenotype_id)

    sample_manifest["phenotype_value"] = sample_manifest["sample_id"].map(phenotype)
    sample_manifest = sample_manifest.dropna(subset=["phenotype_value"]).reset_index(drop=True)
    if sample_manifest.empty:
        raise ValueError(f"No samples have phenotype values for {phenotype.name}.")
    sample_manifest["decile"] = assign_deciles(sample_manifest["phenotype_value"])

    raw_coverage = pd.DataFrame(
        [
            summarize_bigwig(path, gene["seqname"], bins)
            for path in sample_manifest["path"]
        ],
        columns=bins["bin_id"],
    )
    raw_coverage["decile"] = sample_manifest["decile"].to_numpy()

    decile_medians = raw_coverage.groupby("decile", sort=True).median(numeric_only=True)
    decile_counts = sample_manifest.groupby("decile", sort=True).size().rename("n_samples")
    decile_medians = decile_medians.join(decile_counts)
    summary = decile_medians.reset_index().melt(
        id_vars=["decile", "n_samples"],
        var_name="bin_id",
        value_name="median_coverage",
    )
    summary = summary.merge(bins, on="bin_id", how="left", validate="many_to_one")
    summary.insert(2, "phenotype_id", phenotype.name)
    summary = summary[
        [
            "tissue",
            "gene_id",
            "phenotype_id",
            "decile",
            "n_samples",
            "bin_id",
            "seqname",
            "start",
            "end",
            "median_coverage",
        ]
    ]

    args.out_bins.parent.mkdir(parents=True, exist_ok=True)
    args.out_decile_median_coverage.parent.mkdir(parents=True, exist_ok=True)
    bins.to_csv(args.out_bins, sep="\t", index=False, compression="gzip")
    summary.to_csv(
        args.out_decile_median_coverage,
        sep="\t",
        index=False,
        compression="gzip",
    )


if __name__ == "__main__":
    main()
