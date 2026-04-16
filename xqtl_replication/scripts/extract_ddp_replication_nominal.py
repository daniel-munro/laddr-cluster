import argparse
from pathlib import Path

import pandas as pd


GTEX_COLUMNS = [
    "phenotype_id",
    "variant_id",
    "gtex_pval_nominal",
    "gtex_slope",
    "gtex_slope_se",
    "gtex_rank",
]
GEUVADIS_COLUMNS = [
    "phenotype_id",
    "variant_id",
    "geuvadis_pval_nominal",
    "geuvadis_pval_nominal_threshold",
    "geuvadis_slope",
    "geuvadis_slope_se",
]
OUTPUT_COLUMNS = [
    "phenotype_id",
    "variant_id",
    "gtex_pval_nominal",
    "gtex_slope",
    "gtex_slope_se",
    "gtex_rank",
    "geuvadis_pval_nominal",
    "geuvadis_pval_nominal_threshold",
    "geuvadis_slope",
    "geuvadis_slope_se",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--nominal-dir", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    pairs = pairs.rename(
        columns={
            "ddp_pval_nominal": "gtex_pval_nominal",
            "ddp_slope": "gtex_slope",
            "ddp_slope_se": "gtex_slope_se",
            "ddp_rank": "gtex_rank",
        }
    )
    pairs = pairs.loc[:, GTEX_COLUMNS].copy()
    pairs["pair_id"] = pairs["phenotype_id"] + "\t" + pairs["variant_id"]
    pair_ids = set(pairs["pair_id"])

    nominal_dir = Path(args.nominal_dir)
    parquet_paths = sorted(nominal_dir.glob(f"{args.prefix}.cis_qtl_pairs.*.parquet"))
    if not parquet_paths:
        raise FileNotFoundError(
            f"No cis_nominal parquet files found in {nominal_dir} for prefix {args.prefix}."
        )

    nominal_dfs = []
    for path in parquet_paths:
        df = pd.read_parquet(path)
        df["pair_id"] = df["phenotype_id"] + "\t" + df["variant_id"]
        df = df.loc[df["pair_id"].isin(pair_ids)].copy()
        if not df.empty:
            nominal_dfs.append(df)

    if not nominal_dfs:
        raise ValueError("None of the requested DDP variant-phenotype pairs were found in cis_nominal output.")

    nominal = pd.concat(nominal_dfs, ignore_index=True)
    if nominal["pair_id"].duplicated().any():
        raise ValueError("Duplicate phenotype_id/variant_id pairs found in cis_nominal output.")

    nominal = nominal.rename(
        columns={
            "pval_nominal": "geuvadis_pval_nominal",
            "slope": "geuvadis_slope",
            "slope_se": "geuvadis_slope_se",
        }
    )
    nominal["group_id"] = nominal["phenotype_id"].str.split("__", n=1).str[0]

    thresholds = pd.read_csv(
        args.thresholds,
        sep="\t",
        usecols=["group_id", "pval_nominal_threshold"],
    ).rename(columns={"pval_nominal_threshold": "geuvadis_pval_nominal_threshold"})
    if thresholds["group_id"].duplicated().any():
        raise ValueError("Duplicate group_id values found in Geuvadis cis-QTL threshold file.")

    nominal = nominal.merge(thresholds, on="group_id", how="left", validate="many_to_one")
    missing_threshold = nominal["geuvadis_pval_nominal_threshold"].isna()
    if missing_threshold.any():
        missing_groups = sorted(nominal.loc[missing_threshold, "group_id"].unique())
        raise ValueError(f"Missing Geuvadis nominal p-value thresholds for groups: {missing_groups[:10]}")

    nominal = nominal.loc[:, GEUVADIS_COLUMNS + ["pair_id"]]

    out = pairs.merge(
        nominal.drop(columns="pair_id"),
        on=["phenotype_id", "variant_id"],
        how="inner",
        validate="one_to_one",
    )
    if out.empty:
        raise ValueError("None of the requested pairs were found in cis_nominal output.")

    out = out.loc[:, OUTPUT_COLUMNS]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()
