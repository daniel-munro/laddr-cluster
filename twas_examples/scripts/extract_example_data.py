from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig


def sample_to_individual(sample_id):
    return "-".join(sample_id.split("-")[:2])


def parse_variant_id(variant_id):
    chrom, pos, ref, alt, build = variant_id.split("_")
    return {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "build": build,
    }


def format_genotype_counts(genotypes, ref, alt):
    genotype = pd.Series(pd.NA, index=genotypes.index, dtype="string")
    genotype_label = pd.Series(pd.NA, index=genotypes.index, dtype="string")

    for alt_count, gt, label in [
        (0, "0/0", f"{ref}/{ref}"),
        (1, "0/1", f"{ref}/{alt}"),
        (2, "1/1", f"{alt}/{alt}"),
    ]:
        mask = genotypes["alt_allele_count"] == alt_count
        genotype.loc[mask] = gt
        genotype_label.loc[mask] = label

    genotypes["genotype"] = genotype
    genotypes["genotype_label"] = genotype_label
    return genotypes


def make_fixed_bins(gene, n_bins):
    edges = np.linspace(gene["window_start"], gene["window_end"], n_bins + 1)
    starts = np.floor(edges[:-1]).astype(int)
    ends = np.floor(edges[1:]).astype(int)
    ends[-1] = int(gene["window_end"])

    return pd.DataFrame(
        {
            "raw_bin_id": [f"raw_bin_{i:03d}" for i in range(n_bins)],
            "seqname": gene["seqname"],
            "start": starts,
            "end": ends,
            "gene_id": gene["gene_id"],
        }
    )


def summarize_bigwig(path, chrom, bins):
    with pyBigWig.open(str(path)) as bw:
        values = [
            bw.stats(chrom, int(row.start), int(row.end), type="mean", exact=True)[0]
            for row in bins.itertuples(index=False)
        ]
    return [0.0 if value is None else value for value in values]


out_dir = Path(snakemake.output.coverage).parent
out_dir.mkdir(parents=True, exist_ok=True)

tissue = snakemake.params.tissue
gene_id = snakemake.params.gene_id
variant_id = snakemake.params.variant_id
variant = parse_variant_id(variant_id)
raw_coverage_n_bins = int(snakemake.params.raw_coverage_bins)

samples = pd.read_csv(snakemake.input.samples, sep="\t", header=None, names=["sample_id"])
samples["individual_id"] = samples["sample_id"].map(sample_to_individual)

bins = pd.read_csv(snakemake.input.bins, sep="\t")
gene_mask = bins["gene_id"] == gene_id
gene_bins = bins.loc[gene_mask].reset_index(drop=True).copy()
gene_bins.insert(0, "bin_id", [f"bin_{i:03d}" for i in range(gene_bins.shape[0])])

covg = np.load(snakemake.input.covg, mmap_mode="r")
gene_covg = np.power(2.0, covg[gene_mask.to_numpy(), :]) - 8.0

if gene_covg.shape[1] != samples.shape[0]:
    raise ValueError(
        f"Coverage sample count for {tissue} does not match samples.txt: "
        f"{gene_covg.shape[1]} vs {samples.shape[0]}"
    )

coverage = pd.concat(
    [
        samples.reset_index(drop=True),
        pd.DataFrame(gene_covg.T, columns=gene_bins["bin_id"]),
    ],
    axis=1,
)
coverage.insert(2, "tissue", tissue)
coverage.insert(3, "gene_id", gene_id)

genes = pd.read_csv(snakemake.input.genes, sep="\t")
gene = genes.loc[genes["gene_id"] == gene_id].iloc[0]
raw_coverage_bins = make_fixed_bins(gene, raw_coverage_n_bins)
raw_coverage_bins.insert(1, "tissue", tissue)

manifest = pd.read_csv(
    snakemake.input.coverage_manifest,
    sep="\t",
    header=None,
    names=["tissue", "sample_id", "path"],
)
sample_manifest = samples[["sample_id", "individual_id"]].merge(
    manifest.loc[manifest["tissue"] == tissue, ["sample_id", "path"]],
    on="sample_id",
    how="left",
    validate="one_to_one",
)
sample_manifest["path"] = sample_manifest["path"].map(
    lambda path: Path(snakemake.params.coverage_dir) / path
)

raw_coverage_values = [
    summarize_bigwig(path, gene["seqname"], raw_coverage_bins)
    for path in sample_manifest["path"]
]
raw_coverage = pd.concat(
    [
        sample_manifest[["sample_id", "individual_id"]].reset_index(drop=True),
        pd.DataFrame(raw_coverage_values, columns=raw_coverage_bins["raw_bin_id"]),
    ],
    axis=1,
)
raw_coverage.insert(2, "tissue", tissue)
raw_coverage.insert(3, "gene_id", gene_id)

raw = pd.read_csv(snakemake.input.genotypes, sep="\t")
ref_count_col = f'{variant_id}_{variant["ref"]}'

individual_genotypes = raw.loc[:, ["IID", ref_count_col]].rename(
    columns={"IID": "individual_id", ref_count_col: "ref_allele_count"}
)
individual_genotypes["ref_allele_count"] = pd.to_numeric(
    individual_genotypes["ref_allele_count"], errors="coerce"
).astype("Int64")
individual_genotypes["alt_allele_count"] = (
    2 - individual_genotypes["ref_allele_count"]
).astype("Int64")
individual_genotypes = format_genotype_counts(
    individual_genotypes,
    ref=variant["ref"],
    alt=variant["alt"],
)
individual_genotypes.insert(1, "variant_id", variant_id)
individual_genotypes.insert(2, "chrom", variant["chrom"])
individual_genotypes.insert(3, "pos", variant["pos"])
individual_genotypes.insert(4, "ref_allele", variant["ref"])
individual_genotypes.insert(5, "alt_allele", variant["alt"])

gene_bins.insert(1, "tissue", tissue)

gene_bins.to_csv(snakemake.output.bins, sep="\t", index=False, compression="gzip")
coverage.to_csv(snakemake.output.coverage, sep="\t", index=False, compression="gzip")
raw_coverage_bins.to_csv(
    snakemake.output.raw_coverage_bins,
    sep="\t",
    index=False,
    compression="gzip",
)
raw_coverage.to_csv(
    snakemake.output.raw_coverage,
    sep="\t",
    index=False,
    compression="gzip",
)
individual_genotypes.to_csv(
    snakemake.output.individual_genotypes,
    sep="\t",
    index=False,
    compression="gzip",
)
