import tarfile
from argparse import Namespace
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.annotate_hits import annotate
from scripts.correlate_ir_laddr import correlate
from scripts.extract_ir_psi import extract_ir


def test_extract_ir_psi_explodes_semicolon_rows(tmp_path):
    inp = tmp_path / "majiq.psi.tsv"
    out = tmp_path / "ir.tsv.gz"
    inp.write_text(
        "\t".join(
            [
                "ec_idx",
                "seqid",
                "strand",
                "gene_id",
                "gene_name",
                "event_type",
                "start",
                "end",
                "is_intron",
                "S1 raw_psi_mean",
                "S2 raw_psi_mean",
            ]
        )
        + "\n"
        + "0;1\tchr1\t+\tENSG1\tG1\tsource\t100;200\t150;250\tTrue;False\t0.2\t0.4\n"
    )
    extract_ir(inp, out)
    df = pd.read_csv(out, sep="\t")
    assert df["majiq_ir_id"].tolist() == ["IR:ENSG1:chr1:100:150:+:0"]
    assert df["S1"].tolist() == [0.2]


def test_correlate_same_gene_only(tmp_path):
    ir = tmp_path / "ir.tsv.gz"
    laddr = tmp_path / "laddr.bed.gz"
    groups = tmp_path / "groups.txt"
    corr = tmp_path / "corr.tsv.gz"
    high = tmp_path / "high.tsv.gz"

    pd.DataFrame(
        {
            "majiq_ir_id": ["IR:ENSG1:chr1:100:150:+:0"],
            "seqid": ["chr1"],
            "start": [100],
            "end": [150],
            "strand": ["+"],
            "gene_id": ["ENSG1"],
            "gene_id_base": ["ENSG1"],
            "S1": [0.0],
            "S2": [1.0],
            "S3": [2.0],
        }
    ).to_csv(ir, sep="\t", index=False, compression="infer")
    pd.DataFrame(
        {
            "#chr": ["chr1", "chr1"],
            "start": [100, 100],
            "end": [101, 101],
            "phenotype_id": ["ENSG1__PC1", "ENSG2__PC1"],
            "S1": [0.0, 2.0],
            "S2": [1.0, 1.0],
            "S3": [2.0, 0.0],
        }
    ).to_csv(laddr, sep="\t", index=False, compression="infer")
    groups.write_text("ENSG1__PC1\tENSG1\nENSG2__PC1\tENSG2\n")

    correlate(
        Namespace(
            ir=str(ir),
            phenotype_bed=str(laddr),
            phenotype_groups=str(groups),
            phenotype_name="LaDDR",
            correlations_out=str(corr),
            high_out=str(high),
            min_shared_samples=3,
            high_abs_spearman=0.5,
            high_spearman_fdr=0.05,
            chunksize=1,
        )
    )
    df = pd.read_csv(corr, sep="\t")
    assert df["phenotype_id"].tolist() == ["ENSG1__PC1"]
    assert df["spearman_r"].iloc[0] == 1.0


def test_annotate_hits(tmp_path):
    high = tmp_path / "high.tsv.gz"
    indep = tmp_path / "indep.tsv.gz"
    twas_hits = tmp_path / "twas_hits.tsv"
    archive = tmp_path / "weights.tar.bz2"
    out = tmp_path / "annotated.tsv.gz"

    pd.DataFrame({"phenotype_id": ["ENSG1__PC1"]}).to_csv(
        high, sep="\t", index=False, compression="infer"
    )
    pd.DataFrame({"phenotype_id": ["ENSG1__PC1"]}).to_csv(
        indep, sep="\t", index=False, compression="infer"
    )
    pd.DataFrame({"phenotype_id": ["ENSG1__PC1"]}).to_csv(
        twas_hits, sep="\t", index=False
    )
    weight_file = tmp_path / "ENSG1__PC1.wgt.RDat"
    weight_file.write_text("dummy\n")
    with tarfile.open(archive, "w:bz2") as tf:
        tf.add(weight_file, arcname="geuvadis-residual-Geuvadis-latent/ENSG1__PC1.wgt.RDat")

    annotate(
        Namespace(
            input=str(high),
            xqtl_independent=str(indep),
            twas_weights_archive=str(archive),
            twas_hits=str(twas_hits),
            output=str(out),
        )
    )
    df = pd.read_csv(out, sep="\t")
    assert bool(df["has_independent_xqtl"].iloc[0])
    assert bool(df["has_twas_weight"].iloc[0])
    assert bool(df["has_twas_hit"].iloc[0])
