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
            return value
    return None


def benchmark_cpu_time_seconds(row):
    for key in ("cpu_time",):
        if key in row:
            return to_float(row[key])
    return None


def sum_component(paths):
    total_seconds = 0.0
    total_cpu_seconds = 0.0
    max_rss_mb = None
    for path in paths:
        row = read_benchmark(path)
        total_seconds += benchmark_seconds(row)
        cpu_time = benchmark_cpu_time_seconds(row)
        if cpu_time is not None:
            total_cpu_seconds += cpu_time
        rss = benchmark_max_rss_mb(row)
        if rss is not None:
            max_rss_mb = rss if max_rss_mb is None else max(max_rss_mb, rss)
    return total_seconds, total_cpu_seconds, max_rss_mb


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
        "max_threads",
        "runtime_seconds",
        "cpu_time_seconds",
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
        ("alt_TSS/alt_polyA", 16, glob_component(benchmark_dir, "pantry/alt_TSS_polyA/**/*.tsv")),
        ("expression/isoforms", 16, glob_component(benchmark_dir, "pantry/expression/**/*.tsv")),
        ("splicing", 1, glob_component(benchmark_dir, "pantry/splicing/**/*.tsv")),
        ("stability", 8, glob_component(benchmark_dir, "pantry/stability/**/*.tsv")),
        ("bam shrinking", 1, glob_component(benchmark_dir, "pantry/shrink_bam/*.tsv")),
    ]
    laddr_components = [
        ("bamCoverage", 16, glob_component(benchmark_dir, "laddr/bigwig/*.tsv")),
        ("laddr setup", 1, [benchmark_dir / "laddr" / "setup.tsv"]),
        ("laddr binning", 1, glob_component(benchmark_dir, "laddr/binning/*.tsv")),
        ("laddr coverage", 1, glob_component(benchmark_dir, "laddr/coverage/*.tsv")),
        ("laddr fit", 1, glob_component(benchmark_dir, "laddr/fit/*.tsv")),
        ("laddr transform", 1, glob_component(benchmark_dir, "laddr/transform/*.tsv")),
    ]

    rows = []

    alignment_seconds, alignment_cpu_seconds, alignment_rss = sum_component(alignment_paths)
    rows.append(
        {
            "method": "shared",
            "component": "alignment",
            "samples": args.samples,
            "max_threads": 16,
            "runtime_seconds": round_or_blank(alignment_seconds),
            "cpu_time_seconds": round_or_blank(alignment_cpu_seconds),
            "max_rss_mb": round_or_blank(alignment_rss),
            "notes": "Sum of STAR alignment and BAM indexing across samples.",
        }
    )

    pantry_total = 0.0
    pantry_total_cpu = 0.0
    pantry_total_rss = None
    for component, max_threads, paths in pantry_components:
        seconds, cpu_seconds, rss = sum_component(paths)
        pantry_total += seconds
        pantry_total_cpu += cpu_seconds
        if rss is not None:
            pantry_total_rss = rss if pantry_total_rss is None else max(pantry_total_rss, rss)
        rows.append(
            {
                "method": "Pantry",
                "component": component,
                "samples": args.samples,
                "max_threads": max_threads,
                "runtime_seconds": round_or_blank(seconds),
                "cpu_time_seconds": round_or_blank(cpu_seconds),
                "max_rss_mb": round_or_blank(rss),
                "notes": "",
            }
        )

    rows.append(
        {
            "method": "Pantry",
            "component": "Pantry total",
            "samples": args.samples,
            "max_threads": 16,
            "runtime_seconds": round_or_blank(pantry_total),
            "cpu_time_seconds": round_or_blank(pantry_total_cpu),
            "max_rss_mb": round_or_blank(pantry_total_rss),
            "notes": "Sum of Pantry alignment-adjacent shrinking plus modality-group benchmark runs.",
        }
    )

    laddr_total = 0.0
    laddr_total_cpu = 0.0
    laddr_total_rss = None
    for component, max_threads, paths in laddr_components:
        seconds, cpu_seconds, rss = sum_component(paths)
        laddr_total += seconds
        laddr_total_cpu += cpu_seconds
        if rss is not None:
            laddr_total_rss = rss if laddr_total_rss is None else max(laddr_total_rss, rss)
        rows.append(
            {
                "method": "LaDDR",
                "component": component,
                "samples": args.samples,
                "max_threads": max_threads,
                "runtime_seconds": round_or_blank(seconds),
                "cpu_time_seconds": round_or_blank(cpu_seconds),
                "max_rss_mb": round_or_blank(rss),
                "notes": "",
            }
        )

    rows.append(
        {
            "method": "LaDDR",
            "component": "LaDDR total",
            "samples": args.samples,
            "max_threads": 16,
            "runtime_seconds": round_or_blank(laddr_total),
            "cpu_time_seconds": round_or_blank(laddr_total_cpu),
            "max_rss_mb": round_or_blank(laddr_total_rss),
            "notes": "Sum of LaDDR bamCoverage, setup, binning, coverage, fit, and transform.",
        }
    )

    write_table(args.output_tsv, rows, "\t")


if __name__ == "__main__":
    main()
