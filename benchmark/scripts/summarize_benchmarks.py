#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, required=True)
    parser.add_argument("--pantry-threads", type=int, required=True)
    parser.add_argument("--bigwig-threads", type=int, required=True)
    parser.add_argument("--output-tsv", required=True)
    return parser.parse_args()


def read_benchmark(path):
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one row in benchmark file: {path}")
    return rows[0]


def to_float(value):
    if value in (None, "", "NA", "nan"):
        return None
    return float(value)


def round_or_blank(value):
    if value is None:
        return ""
    return f"{value:.3f}"


def benchmark_seconds(row):
    for key in ("s", "seconds"):
        if key in row:
            return to_float(row[key])
    raise KeyError(f"No runtime column found in benchmark row: {row}")


def benchmark_max_rss_mb(row):
    for key in ("max_rss", "max_rss_mb"):
        if key in row:
            value = to_float(row[key])
            if value is None:
                return None
            if key == "max_rss":
                return value / 1024
            return value
    return None


def sum_component(paths):
    total_seconds = 0.0
    max_rss_mb = None
    for path in paths:
        row = read_benchmark(path)
        total_seconds += benchmark_seconds(row)
        rss = benchmark_max_rss_mb(row)
        if rss is not None:
            max_rss_mb = rss if max_rss_mb is None else max(max_rss_mb, rss)
    return total_seconds, max_rss_mb


def one_component(path):
    row = read_benchmark(path)
    return benchmark_seconds(row), benchmark_max_rss_mb(row)


def glob_component(base_dir, *patterns):
    paths = []
    for pattern in patterns:
        paths.extend(sorted(base_dir.glob(pattern)))
    return paths


def write_table(path, rows, delimiter):
    fieldnames = [
        "method",
        "component",
        "samples",
        "threads",
        "runtime_seconds",
        "max_rss_mb",
        "notes",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()
    benchmark_dir = Path("benchmarks")

    alignment_paths = glob_component(
        benchmark_dir,
        "alignment/star_align/*.tsv",
        "alignment/index_bam/*.tsv",
    )
    bigwig_paths = glob_component(benchmark_dir, "laddr/bigwig/*.tsv")

    pantry_components = [
        ("alt_TSS/alt_polyA", glob_component(benchmark_dir, "pantry/alt_TSS_polyA/**/*.tsv")),
        ("expression/isoforms", glob_component(benchmark_dir, "pantry/expression/**/*.tsv")),
        ("splicing", glob_component(benchmark_dir, "pantry/splicing/**/*.tsv")),
        ("stability", glob_component(benchmark_dir, "pantry/stability/**/*.tsv")),
        ("bam shrinking", glob_component(benchmark_dir, "alignment/shrink_bam/*.tsv")),
    ]
    laddr_components = [
        ("laddr setup", [benchmark_dir / "laddr" / "setup.tsv"]),
        ("laddr binning", [benchmark_dir / "laddr" / "binning.tsv"]),
        ("laddr coverage", [benchmark_dir / "laddr" / "coverage.tsv"]),
        ("laddr fit", [benchmark_dir / "laddr" / "fit.tsv"]),
        ("laddr transform", glob_component(benchmark_dir, "laddr/transform/*.tsv")),
    ]

    rows = []

    alignment_seconds, alignment_rss = sum_component(alignment_paths)
    rows.append(
        {
            "method": "shared",
            "component": "alignment",
            "samples": args.samples,
            "threads": "",
            "runtime_seconds": round_or_blank(alignment_seconds),
            "max_rss_mb": round_or_blank(alignment_rss),
            "notes": "Sum of STAR alignment and BAM indexing across samples.",
        }
    )

    bigwig_seconds, bigwig_rss = sum_component(bigwig_paths)
    rows.append(
        {
            "method": "LaDDR",
            "component": "bamCoverage",
            "samples": args.samples,
            "threads": args.bigwig_threads,
            "runtime_seconds": round_or_blank(bigwig_seconds),
            "max_rss_mb": round_or_blank(bigwig_rss),
            "notes": "Sum of per-sample bigWig generation from Pantry BAMs.",
        }
    )

    pantry_total = 0.0
    pantry_total_rss = None
    for component, paths in pantry_components:
        seconds, rss = sum_component(paths)
        pantry_total += seconds
        if rss is not None:
            pantry_total_rss = rss if pantry_total_rss is None else max(pantry_total_rss, rss)
        rows.append(
            {
                "method": "Pantry",
                "component": component,
                "samples": args.samples,
                "threads": args.pantry_threads,
                "runtime_seconds": round_or_blank(seconds),
                "max_rss_mb": round_or_blank(rss),
                "notes": "",
            }
        )

    rows.append(
        {
            "method": "Pantry",
            "component": "Pantry total",
            "samples": args.samples,
            "threads": args.pantry_threads,
            "runtime_seconds": round_or_blank(pantry_total),
            "max_rss_mb": round_or_blank(pantry_total_rss),
            "notes": "Sum of Pantry alignment-adjacent shrinking plus modality-group benchmark runs.",
        }
    )

    laddr_total = 0.0
    laddr_total_rss = None
    for component, paths in laddr_components:
        seconds, rss = sum_component(paths)
        laddr_total += seconds
        if rss is not None:
            laddr_total_rss = rss if laddr_total_rss is None else max(laddr_total_rss, rss)
        rows.append(
            {
                "method": "LaDDR",
                "component": component,
                "samples": args.samples,
                "threads": "",
                "runtime_seconds": round_or_blank(seconds),
                "max_rss_mb": round_or_blank(rss),
                "notes": "",
            }
        )

    rows.append(
        {
            "method": "LaDDR",
            "component": "LaDDR total",
            "samples": args.samples,
            "threads": "",
            "runtime_seconds": round_or_blank(laddr_total),
            "max_rss_mb": round_or_blank(laddr_total_rss),
            "notes": "Sum of LaDDR bamCoverage, setup, binning, coverage, fit, and transform.",
        }
    )

    write_table(args.output_tsv, rows, "\t")


if __name__ == "__main__":
    main()
