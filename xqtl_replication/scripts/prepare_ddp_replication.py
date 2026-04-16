import argparse
from pathlib import Path

import pandas as pd


OUTPUT_COLUMNS = [
    "phenotype_id",
    "variant_id",
    "gtex_pval_nominal",
    "gtex_slope",
    "gtex_slope_se",
    "gtex_rank",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--pairs-out", required=True)
    parser.add_argument("--phenotypes-out", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    df = df.drop_duplicates(subset=["phenotype_id", "variant_id"]).copy()
    df = df.rename(
        columns={
            "pval_nominal": "gtex_pval_nominal",
            "slope": "gtex_slope",
            "slope_se": "gtex_slope_se",
            "rank": "gtex_rank",
        }
    )
    df = df.loc[:, OUTPUT_COLUMNS]

    pairs_out = Path(args.pairs_out)
    pairs_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(pairs_out, sep="\t", index=False, compression="gzip")

    phenotype_ids = pd.Series(df["phenotype_id"].drop_duplicates().sort_values())
    phenotype_ids.to_csv(args.phenotypes_out, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
