#!/usr/bin/env python3
"""Summarize the MAJIQ retained-intron and LaDDR comparison outputs.

The TSV output is a machine-readable metric/value table. The Markdown output is
a short human-readable report describing what each count represents in the
retained-intron, correlation, and evidence-annotation stages.
"""

import argparse
from pathlib import Path

import pandas as pd


METRIC_DESCRIPTIONS = {
    "retained_introns": (
        "Distinct retained-intron events extracted from the MAJIQ PSI table."
    ),
    "retained_intron_genes": (
        "Distinct Ensembl genes represented by those retained-intron events."
    ),
    "ir_laddr_correlation_pairs": (
        "Same-gene retained-intron/LaDDR phenotype pairs tested for correlation."
    ),
    "retained_introns_with_laddr_match": (
        "Retained introns with at least one same-gene LaDDR phenotype tested."
    ),
    "high_correlation_pairs": (
        "Pairs passing the configured absolute Spearman and Spearman FDR thresholds."
    ),
    "high_correlation_retained_introns": (
        "Distinct retained introns represented among high-correlation pairs."
    ),
    "high_correlation_laddr_phenotypes": (
        "Distinct LaDDR phenotypes represented among high-correlation pairs."
    ),
    "high_pairs_with_independent_xqtl": (
        "High-correlation pairs whose LaDDR phenotype has independent cis-xQTL evidence."
    ),
    "high_pairs_with_twas_weight": (
        "High-correlation pairs whose LaDDR phenotype has an available TWAS weight."
    ),
    "high_pairs_with_twas_hit": (
        "High-correlation pairs whose LaDDR phenotype is also a significant TWAS hit."
    ),
}


def bool_sum(df, col):
    if col not in df.columns:
        return 0
    return int(df[col].fillna(False).astype(bool).sum())


def metric_lookup(metrics):
    return dict(metrics)


def pct(numerator, denominator):
    if denominator == 0:
        return "NA"
    return f"{100 * numerator / denominator:.1f}%"


def metric_table(metrics):
    lines = ["| Metric | Value | Description |", "| --- | ---: | --- |"]
    for metric, value in metrics:
        lines.append(f"| `{metric}` | {value} | {METRIC_DESCRIPTIONS[metric]} |")
    return lines


def markdown_report(metrics):
    values = metric_lookup(metrics)
    retained_introns = values["retained_introns"]
    high_pairs = values["high_correlation_pairs"]
    tested_pairs = values["ir_laddr_correlation_pairs"]

    lines = [
        "# MAJIQ Retained-Intron Comparison",
        "",
        "This report summarizes retained introns extracted from MAJIQ PSI output, "
        "their same-gene correlations with LaDDR phenotypes, and overlap between "
        "the high-correlation pairs and external xQTL/TWAS evidence.",
        "",
        "## Overview",
        "",
        f"- {values['retained_introns']} retained introns from "
        f"{values['retained_intron_genes']} genes were available for comparison.",
        f"- {values['retained_introns_with_laddr_match']} retained introns "
        f"({pct(values['retained_introns_with_laddr_match'], retained_introns)}) "
        "had at least one same-gene LaDDR phenotype tested.",
        f"- {values['high_correlation_pairs']} of {tested_pairs} tested pairs "
        f"({pct(values['high_correlation_pairs'], tested_pairs)}) passed the "
        "high-correlation thresholds.",
        f"- High-correlation pairs involved "
        f"{values['high_correlation_retained_introns']} retained introns and "
        f"{values['high_correlation_laddr_phenotypes']} LaDDR phenotypes.",
        "",
        "## Evidence Overlap",
        "",
        f"- Independent xQTL evidence: {values['high_pairs_with_independent_xqtl']} "
        f"high-correlation pairs ({pct(values['high_pairs_with_independent_xqtl'], high_pairs)}).",
        f"- Available TWAS weights: {values['high_pairs_with_twas_weight']} "
        f"high-correlation pairs ({pct(values['high_pairs_with_twas_weight'], high_pairs)}).",
        f"- Significant TWAS hits: {values['high_pairs_with_twas_hit']} "
        f"high-correlation pairs ({pct(values['high_pairs_with_twas_hit'], high_pairs)}).",
        "",
        "## Metric Definitions",
        "",
        *metric_table(metrics),
    ]
    return "\n".join(lines) + "\n"


def summarize(args):
    ir = pd.read_csv(args.ir, sep="\t")
    corr = pd.read_csv(args.correlations, sep="\t")
    high = pd.read_csv(args.high, sep="\t")

    metrics = [
        ("retained_introns", ir["majiq_ir_id"].nunique()),
        ("retained_intron_genes", ir["gene_id_base"].nunique()),
        ("ir_laddr_correlation_pairs", len(corr)),
        ("retained_introns_with_laddr_match", corr["majiq_ir_id"].nunique() if len(corr) else 0),
        ("high_correlation_pairs", len(high)),
        ("high_correlation_retained_introns", high["majiq_ir_id"].nunique() if len(high) else 0),
        ("high_correlation_laddr_phenotypes", high["phenotype_id"].nunique() if len(high) else 0),
        ("high_pairs_with_independent_xqtl", bool_sum(high, "has_independent_xqtl")),
        ("high_pairs_with_twas_weight", bool_sum(high, "has_twas_weight")),
        ("high_pairs_with_twas_hit", bool_sum(high, "has_twas_hit")),
    ]
    summary = pd.DataFrame(metrics, columns=["metric", "value"])
    Path(args.summary_tsv).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_tsv, sep="\t", index=False)

    Path(args.summary_md).write_text(markdown_report(metrics))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Write machine-readable and Markdown summaries for the MAJIQ "
            "retained-intron/LaDDR comparison."
        )
    )
    parser.add_argument(
        "--ir",
        required=True,
        help="Retained-intron PSI table from extract_ir_psi.py.",
    )
    parser.add_argument(
        "--correlations",
        required=True,
        help="All retained-intron/LaDDR correlation pairs from correlate_ir_laddr.py.",
    )
    parser.add_argument(
        "--high",
        required=True,
        help="Annotated high-correlation table from annotate_hits.py.",
    )
    parser.add_argument(
        "--summary-tsv",
        required=True,
        help="Output metric/value TSV summary.",
    )
    parser.add_argument(
        "--summary-md",
        required=True,
        help="Output Markdown report with metric interpretation.",
    )
    args = parser.parse_args()
    summarize(args)


if __name__ == "__main__":
    main()
