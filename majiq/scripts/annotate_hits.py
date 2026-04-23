#!/usr/bin/env python3
"""Annotate high retained-intron/LaDDR correlations with QTL and TWAS evidence.

The input is the filtered high-correlation table from correlate_ir_laddr.py. This
script adds boolean columns indicating whether each LaDDR phenotype has an
independent cis xQTL, a TWAS weight, and a significant TWAS hit.
"""

import argparse
import tarfile
from pathlib import Path

import pandas as pd


def read_phenotype_ids(path):
    df = pd.read_csv(path, sep="\t")
    return set(df["phenotype_id"].astype(str))


def twas_weight_ids(path):
    ids = set()
    with tarfile.open(path, "r:bz2") as tf:
        for member in tf.getmembers():
            name = Path(member.name).name
            if name.endswith(".wgt.RDat"):
                ids.add(name[: -len(".wgt.RDat")])
    return ids


def twas_hit_ids(path):
    df = pd.read_csv(path, sep="\t")
    id_col = "phenotype_id" if "phenotype_id" in df.columns else "ID" if "ID" in df.columns else None
    if id_col is None:
        raise ValueError(f"TWAS hits file {path} must contain phenotype_id or ID")
    return set(df[id_col].astype(str))


def annotate(args):
    high = pd.read_csv(args.input, sep="\t")
    phenotype_ids = high["phenotype_id"].astype(str)
    high["has_independent_xqtl"] = phenotype_ids.isin(read_phenotype_ids(args.xqtl_independent))
    high["has_twas_weight"] = phenotype_ids.isin(twas_weight_ids(args.twas_weights_archive))
    twas_hits = twas_hit_ids(args.twas_hits)
    high["has_twas_hit"] = phenotype_ids.isin(twas_hits)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    high.to_csv(args.output, sep="\t", index=False, compression="infer")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Add xQTL and TWAS evidence flags to high retained-intron/LaDDR "
            "correlation pairs."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="High-correlation retained-intron/LaDDR table from correlate_ir_laddr.py.",
    )
    parser.add_argument(
        "--xqtl-independent",
        required=True,
        help="Significant independent cis-QTL table with a phenotype_id column.",
    )
    parser.add_argument(
        "--twas-weights-archive",
        required=True,
        help="tar.bz2 archive containing TWAS weight files named <phenotype_id>.wgt.RDat.",
    )
    parser.add_argument(
        "--twas-hits",
        required=True,
        help="Significant TWAS hit table with phenotype_id or ID.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Annotated high-correlation output table.",
    )
    args = parser.parse_args()
    annotate(args)


if __name__ == "__main__":
    main()
