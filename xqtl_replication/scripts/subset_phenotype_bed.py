import argparse

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--output-bed", required=True)
    args = parser.parse_args()

    phenotype_ids = pd.read_csv(args.phenotypes, sep="\t", header=None, dtype=str)[0]
    phenotype_set = set(phenotype_ids)

    bed = pd.read_csv(args.input_bed, sep="\t", dtype={"#chr": str, "phenotype_id": str})
    keep = bed["phenotype_id"].isin(phenotype_set)
    subset = bed.loc[keep].copy()

    if subset.empty:
        raise ValueError("None of the requested phenotypes were found in the BED.")

    subset.to_csv(args.output_bed, sep="\t", index=False)


if __name__ == "__main__":
    main()
