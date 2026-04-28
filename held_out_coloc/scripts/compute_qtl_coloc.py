import argparse
from pathlib import Path

import pandas as pd


SUSIE_COLUMNS = [
    "analysis",
    "held_out_modality",
    "phenotype_id",
    "variant_id",
    "pip",
    "af",
    "cs_id",
    "gene_id",
]

COLOC_COLUMNS = [
    "held_out_modality",
    "comparison",
    "gene_id",
    "left_phenotype_id",
    "left_phenotype_class",
    "left_cs_id",
    "left_cs_size",
    "left_top_variant",
    "left_top_pip",
    "right_phenotype_id",
    "right_phenotype_class",
    "right_cs_id",
    "right_cs_size",
    "right_top_variant",
    "right_top_pip",
    "cs_overlap",
    "clpp",
    "clpa",
    "top_shared_variant",
    "top_shared_left_pip",
    "top_shared_right_pip",
    "top_shared_clpp",
]

SUMMARY_COLUMNS = [
    "held_out_modality",
    "n_original_cs",
    "n_baseline_coloc_cs",
    "n_held_out_coloc_cs",
    "n_held_out_only_cs",
    "baseline_coloc_fraction",
    "held_out_coloc_fraction",
    "held_out_only_fraction",
    "estimated_recaptured_excess",
]


def parse_modalities(value: str) -> list[str]:
    modalities = [item for item in value.split(",") if item]
    if not modalities:
        raise ValueError("At least one modality is required.")
    if len(modalities) != len(set(modalities)):
        raise ValueError(f"Duplicate modalities: {modalities}")
    return modalities


def phenotype_class(phenotype_id: pd.Series) -> pd.Series:
    return phenotype_id.str.split(":", n=1).str[0]


def gene_id_from_phenotype(phenotype_id: pd.Series) -> pd.Series:
    payload = phenotype_id.str.split(":", n=1).str[1]
    parsed = payload.str.split("__", n=1).str[0]
    if parsed.isna().any():
        examples = phenotype_id.loc[parsed.isna()].head(10).tolist()
        raise ValueError(f"Could not parse gene IDs from phenotype_id values: {examples}")
    return parsed


def normalize_variant(variant_id: pd.Series) -> pd.Series:
    parts = variant_id.str.rsplit("_", n=4, expand=True)
    if parts.shape[1] != 5:
        raise ValueError("variant_id values do not match chr_pos_ref_alt_build format.")
    chrom = parts[0].str.removeprefix("chr")
    return chrom + ":" + parts[1] + ":" + parts[2] + ":" + parts[3]


def read_susie(path: str) -> pd.DataFrame:
    susie = pd.read_csv(path, sep="\t", usecols=SUSIE_COLUMNS, dtype={"cs_id": str})
    if susie.empty:
        raise ValueError(f"SuSiE table is empty: {path}")
    expected_gene_id = gene_id_from_phenotype(susie["phenotype_id"].astype(str))
    mismatch = susie["gene_id"].astype(str) != expected_gene_id.astype(str)
    if mismatch.any():
        examples = susie.loc[mismatch, ["phenotype_id", "gene_id"]].head(10)
        raise ValueError(f"gene_id does not match phenotype_id in {path}:\n{examples}")

    susie["phenotype_class"] = phenotype_class(susie["phenotype_id"].astype(str))
    susie["variant_key"] = normalize_variant(susie["variant_id"].astype(str))
    susie["pip"] = susie["pip"].astype(float)
    dup_cols = ["analysis", "held_out_modality", "phenotype_id", "cs_id", "variant_key"]
    if susie.duplicated(dup_cols).any():
        examples = susie.loc[susie.duplicated(dup_cols, keep=False), dup_cols].head(10)
        raise ValueError(f"Duplicate variants within a credible set in {path}:\n{examples}")
    return susie


def prepare_side(df: pd.DataFrame, side: str) -> pd.DataFrame:
    rename = {
        "phenotype_id": f"{side}_phenotype_id",
        "phenotype_class": f"{side}_phenotype_class",
        "cs_id": f"{side}_cs_id",
        "pip": f"{side}_pip",
        "variant_id": f"{side}_variant_id",
        "af": f"{side}_af",
    }
    keep = ["gene_id", "phenotype_id", "phenotype_class", "cs_id", "variant_key", "variant_id", "pip", "af"]
    return df.loc[:, keep].rename(columns=rename)


def top_rows(df: pd.DataFrame, group_cols: list[str], pip_col: str, prefix: str) -> pd.DataFrame:
    idx = df.groupby(group_cols, observed=True)[pip_col].idxmax()
    return df.loc[idx, group_cols + ["variant_key", pip_col]].rename(
        columns={"variant_key": f"{prefix}_top_variant", pip_col: f"{prefix}_top_pip"}
    )


def compute_coloc(left: pd.DataFrame, right: pd.DataFrame, modality: str, comparison: str) -> pd.DataFrame:
    left = prepare_side(left, "left")
    right = prepare_side(right, "right")
    if left.empty:
        raise ValueError(f"No baseline credible sets for modality {modality}.")
    if right.empty:
        return pd.DataFrame(columns=COLOC_COLUMNS)

    left_group_cols = ["gene_id", "left_phenotype_id", "left_phenotype_class", "left_cs_id"]
    right_group_cols = ["gene_id", "right_phenotype_id", "right_phenotype_class", "right_cs_id"]
    group_cols = left_group_cols + ["right_phenotype_id", "right_phenotype_class", "right_cs_id"]

    shared_right = right.loc[right["variant_key"].isin(set(left["variant_key"]))].copy()
    if shared_right.empty:
        return pd.DataFrame(columns=COLOC_COLUMNS)
    shared_left = left.loc[left["variant_key"].isin(set(shared_right["variant_key"]))].copy()
    merged = shared_left.merge(shared_right, on=["gene_id", "variant_key"], how="inner", validate="many_to_many")
    if merged.empty:
        return pd.DataFrame(columns=COLOC_COLUMNS)

    merged["variant_clpp"] = merged["left_pip"] * merged["right_pip"]
    merged["variant_clpa"] = merged[["left_pip", "right_pip"]].min(axis=1)
    stats = (
        merged.groupby(group_cols, observed=True)
        .agg(
            cs_overlap=("variant_key", "size"),
            clpp=("variant_clpp", "sum"),
            clpa=("variant_clpa", "sum"),
        )
        .reset_index()
    )

    left_sizes = left.groupby(left_group_cols, observed=True).size().rename("left_cs_size").reset_index()
    right_sizes = right.groupby(right_group_cols, observed=True).size().rename("right_cs_size").reset_index()
    left_top = top_rows(left, left_group_cols, "left_pip", "left")
    right_top = top_rows(right, right_group_cols, "right_pip", "right")
    top_shared_idx = merged.groupby(group_cols, observed=True)["variant_clpp"].idxmax()
    top_shared = merged.loc[
        top_shared_idx,
        group_cols + ["variant_key", "left_pip", "right_pip", "variant_clpp"],
    ].rename(
        columns={
            "variant_key": "top_shared_variant",
            "left_pip": "top_shared_left_pip",
            "right_pip": "top_shared_right_pip",
            "variant_clpp": "top_shared_clpp",
        }
    )

    out = stats.merge(left_sizes, on=left_group_cols, how="left", validate="many_to_one")
    out = out.merge(right_sizes, on=right_group_cols, how="left", validate="many_to_one")
    out = out.merge(left_top, on=left_group_cols, how="left", validate="many_to_one")
    out = out.merge(right_top, on=right_group_cols, how="left", validate="many_to_one")
    out = out.merge(top_shared, on=group_cols, how="left", validate="one_to_one")
    out.insert(0, "comparison", comparison)
    out.insert(0, "held_out_modality", modality)
    return out.loc[:, COLOC_COLUMNS].sort_values(["clpp", "clpa"], ascending=False)


def original_cs_keys(df: pd.DataFrame) -> pd.Series:
    return (
        df["left_phenotype_id"].astype(str)
        + "\t"
        + df["left_cs_id"].astype(str)
    )


def summarize(significant: pd.DataFrame, original: pd.DataFrame, modalities: list[str]) -> pd.DataFrame:
    rows = []
    for modality in modalities:
        original_keys = set(
            original.loc[original["phenotype_class"] == modality, "phenotype_id"].astype(str)
            + "\t"
            + original.loc[original["phenotype_class"] == modality, "cs_id"].astype(str)
        )
        modality_hits = significant.loc[significant["held_out_modality"] == modality].copy()
        baseline = modality_hits.loc[modality_hits["comparison"] == "baseline_control"]
        held_out = modality_hits.loc[modality_hits["comparison"] == "held_out"]
        baseline_keys = set(original_cs_keys(baseline)) if not baseline.empty else set()
        held_out_keys = set(original_cs_keys(held_out)) if not held_out.empty else set()
        held_out_only = held_out_keys - baseline_keys
        n_original = len(original_keys)
        rows.append(
            {
                "held_out_modality": modality,
                "n_original_cs": n_original,
                "n_baseline_coloc_cs": len(baseline_keys),
                "n_held_out_coloc_cs": len(held_out_keys),
                "n_held_out_only_cs": len(held_out_only),
                "baseline_coloc_fraction": len(baseline_keys) / n_original if n_original else 0,
                "held_out_coloc_fraction": len(held_out_keys) / n_original if n_original else 0,
                "held_out_only_fraction": len(held_out_only) / n_original if n_original else 0,
                "estimated_recaptured_excess": len(held_out_keys) - len(baseline_keys),
            }
        )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def write_table(df: pd.DataFrame, path: str, columns: list[str], compression: str | None = None) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if df.empty:
        df = pd.DataFrame(columns=columns)
    else:
        df = df.loc[:, columns]
    df.to_csv(path, sep="\t", index=False, compression=compression)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-susie", required=True)
    parser.add_argument("--held-out-susie", required=True)
    parser.add_argument("--modalities", required=True)
    parser.add_argument("--out-all", required=True)
    parser.add_argument("--out-significant", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--clpp-threshold", type=float, default=0.01)
    parser.add_argument("--clpa-threshold", type=float, default=0.5)
    args = parser.parse_args()

    modalities = parse_modalities(args.modalities)
    baseline = read_susie(args.baseline_susie)
    held_out = read_susie(args.held_out_susie)

    baseline_latent = baseline.loc[baseline["phenotype_class"] == "latent_residual"].copy()
    if baseline_latent.empty:
        raise ValueError("No baseline latent_residual credible sets found.")

    frames = []
    for modality in modalities:
        original = baseline.loc[baseline["phenotype_class"] == modality].copy()
        baseline_control = compute_coloc(original, baseline_latent, modality, "baseline_control")
        if not baseline_control.empty:
            frames.append(baseline_control)

        held_out_latent = held_out.loc[
            (held_out["held_out_modality"] == modality)
            & (held_out["phenotype_class"] == "latent_residual")
        ].copy()
        held_out_coloc = compute_coloc(original, held_out_latent, modality, "held_out")
        if not held_out_coloc.empty:
            frames.append(held_out_coloc)

    coloc = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=COLOC_COLUMNS)
    significant = coloc.loc[
        (coloc["clpp"] >= args.clpp_threshold) | (coloc["clpa"] >= args.clpa_threshold)
    ].copy()
    summary = summarize(significant, baseline, modalities)

    write_table(coloc, args.out_all, COLOC_COLUMNS, compression="gzip")
    write_table(significant, args.out_significant, COLOC_COLUMNS, compression="gzip")
    write_table(summary, args.out_summary, SUMMARY_COLUMNS)


if __name__ == "__main__":
    main()
