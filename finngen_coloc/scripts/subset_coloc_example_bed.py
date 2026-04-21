import argparse
from pathlib import Path

import pandas as pd


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, span = region.split(":", 1)
    start, end = span.split("-", 1)
    return chrom, int(start), int(end)


def resolve_phenotype_id(phenotype_id: str, bed_ids: pd.Series) -> str:
    if phenotype_id in set(bed_ids):
        return phenotype_id

    suffix = ":" + phenotype_id
    matches = bed_ids.loc[bed_ids.str.endswith(suffix)].drop_duplicates()
    if len(matches) != 1:
        raise ValueError(
            f"Could not uniquely resolve phenotype_id {phenotype_id!r} in phenotype BED; "
            f"found {len(matches)} suffix matches."
        )
    return matches.iloc[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--examples", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--output-bed", required=True)
    args = parser.parse_args()

    examples = pd.read_csv(args.examples, sep="\t", dtype=str)
    examples = examples.loc[examples["tissue"] == args.tissue].copy()
    if examples.empty:
        raise ValueError(f"No examples found for tissue {args.tissue}.")
    examples["example_id"] = [f"example_{i}" for i in examples.index + 1]

    bed = pd.read_csv(args.input_bed, sep="\t", dtype={"#chr": str, "phenotype_id": str})
    bed_ids = bed["phenotype_id"].astype(str)
    subset_rows = []
    for _, example in examples.iterrows():
        phenotype_id = resolve_phenotype_id(example["phenotype_id"], bed_ids)
        bed_row = bed.loc[bed["phenotype_id"] == phenotype_id].copy()
        if len(bed_row) != 1:
            raise ValueError(f"Expected one BED row for {phenotype_id}, found {len(bed_row)}.")

        chrom, start, end = parse_region(example["finngen_region"])
        bed_row["#chr"] = chrom
        bed_row["start"] = start
        bed_row["end"] = end
        bed_row["phenotype_id"] = example["example_id"]
        subset_rows.append(bed_row)

    subset = pd.concat(subset_rows, ignore_index=True)

    Path(args.output_bed).parent.mkdir(parents=True, exist_ok=True)
    subset.to_csv(args.output_bed, sep="\t", index=False)


if __name__ == "__main__":
    main()
