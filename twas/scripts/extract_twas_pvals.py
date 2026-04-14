import argparse
from pathlib import Path

import pandas as pd


TRAIT_DF = pd.read_csv("input/gwas/gwas_metadata.txt", sep="\t")
TRAITS = TRAIT_DF["Tag"].tolist()
KDP_MODALITIES = ["alt_polyA", "alt_TSS", "expression", "isoforms", "splicing", "stability"]


def ddp_file(trait: str) -> str:
    return f"output/geuvadis-full-Geuvadis/fusion.geuvadis-full-Geuvadis.{trait}.tsv"


def kdp_file(modality: str, trait: str) -> str:
    return f"../../pantry/twas/output/GEUVADIS/{modality}/fusion.GEUVADIS.{modality}.{trait}.tsv"


def extract_ddp() -> pd.DataFrame:
    dfs = []
    for trait in TRAITS:
        df = pd.read_csv(ddp_file(trait), sep="\t", usecols=["ID", "TWAS.P"])
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def extract_kdp() -> pd.DataFrame:
    dfs = []
    for modality in KDP_MODALITIES:
        for trait in TRAITS:
            df = pd.read_csv(kdp_file(modality, trait), sep="\t", usecols=["ID", "TWAS.P"])
            dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["ddp", "kdp"], required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if args.dataset == "ddp":
        df = extract_ddp()
    else:
        df = extract_kdp()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()
