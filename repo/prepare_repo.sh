set -e

mkdir -p repo/{info,covariates,geuvadis,laddr_models,rna_phenotypes,qtls,twas_model_summaries}

rsync -a ../data/gtex/tissues.gtex.txt repo/info/tissues-gtex.txt
rsync -av ../../pantry/repo/info/LDREF_b38ids_chr.tar.bz2 repo/info/

echo "=== GEUVADIS ==="

rsync -a ../pheast/geuvadis/input/phenotypes/geuvadis-full-Geuvadis-latent.bed.gz repo/geuvadis/ddp-GEUVADIS.bed.gz
rsync -a ../pheast/geuvadis/input/phenotypes/geuvadis-residual-Geuvadis-latent.bed.gz repo/geuvadis/rddp-GEUVADIS.bed.gz

rsync -a ../pheast/geuvadis/intermediate/covar/geuvadis-full-Geuvadis-latent.covar.tsv repo/geuvadis/ddp-GEUVADIS.covar.tsv
rsync -a ../pheast/geuvadis/intermediate/covar/geuvadis-residual-Geuvadis-cross_latent.covar.tsv repo/geuvadis/hybrid-GEUVADIS.covar.tsv
rsync -a ../pheast/geuvadis/intermediate/covar/geuvadis-full-Geuvadis-latent.covar.plink.tsv repo/geuvadis/ddp-GEUVADIS.covar.plink.tsv
rsync -a ../pheast/geuvadis/intermediate/covar/geuvadis-residual-Geuvadis-latent.covar.plink.tsv repo/geuvadis/rddp-GEUVADIS.covar.plink.tsv

rsync -a ../pheast/geuvadis/output/qtl/geuvadis-full-Geuvadis-latent.cis_qtl.txt.gz repo/geuvadis/ddp-GEUVADIS.cis_qtl.txt.gz
rsync -a ../pheast/geuvadis/output/qtl/geuvadis-residual-Geuvadis-cross_latent.cis_qtl.txt.gz repo/geuvadis/hybrid-GEUVADIS.cis_qtl.txt.gz
rsync -a ../pheast/geuvadis/output/qtl/geuvadis-full-Geuvadis-latent.cis_independent_qtl.txt.gz repo/geuvadis/ddp-GEUVADIS.cis_independent_qtl.txt.gz
rsync -a ../pheast/geuvadis/output/qtl/geuvadis-residual-Geuvadis-cross_latent.cis_independent_qtl.txt.gz repo/geuvadis/hybrid-GEUVADIS.cis_independent_qtl.txt.gz

rsync -a ../pheast/geuvadis/intermediate/twas/geuvadis-full-Geuvadis-latent.profile repo/geuvadis/ddp-GEUVADIS.profile
rsync -a ../pheast/geuvadis/intermediate/twas/geuvadis-residual-Geuvadis-latent.profile repo/geuvadis/rddp-GEUVADIS.profile

gzip -c ../twas/output/twas_hits.geuvadis-full-Geuvadis.tsv > repo/geuvadis/ddp-GEUVADIS.twas_hits.tsv.gz
gzip -c ../twas/output/twas_hits.geuvadis-residual-Geuvadis.tsv > repo/geuvadis/rddp-GEUVADIS.twas_hits.tsv.gz

echo "=== GTEx ==="
cat ../data/gtex/tissues.gtex.txt | while read tissue; do
    echo $tissue
    rsync -a ../pheast/gtextcga-full/input/phenotypes/$tissue-latent.bed.gz repo/rna_phenotypes/ddp-$tissue.bed.gz
    rsync -a ../pheast/gtex-residual/input/phenotypes/$tissue-latent.bed.gz repo/rna_phenotypes/rddp-$tissue.bed.gz

    rsync -a ../pheast/gtextcga-full/intermediate/covar/$tissue-latent.covar.tsv repo/covariates/ddp-$tissue.covar.tsv
    rsync -a ../pheast/gtex-residual/intermediate/covar/$tissue-cross_latent.covar.tsv repo/covariates/hybrid-$tissue.covar.tsv
    rsync -a ../pheast/gtextcga-full/intermediate/covar/$tissue-latent.covar.plink.tsv repo/covariates/ddp-$tissue.covar.plink.tsv
    rsync -a ../pheast/gtex-residual/intermediate/covar/$tissue-latent.covar.plink.tsv repo/covariates/rddp-$tissue.covar.plink.tsv

    rsync -a ../pheast/gtextcga-full/intermediate/twas/$tissue-latent.profile repo/twas_model_summaries/ddp-$tissue.profile
    rsync -a ../pheast/gtex-residual/intermediate/twas/$tissue-latent.profile repo/twas_model_summaries/rddp-$tissue.profile

    rsync -a ../pheast/gtextcga-full/output/qtl/$tissue-latent.cis_qtl.txt.gz repo/qtls/ddp-$tissue.cis_qtl.txt.gz
    rsync -a ../pheast/gtex-residual/output/qtl/$tissue-cross_latent.cis_qtl.txt.gz repo/qtls/hybrid-$tissue.cis_qtl.txt.gz
    rsync -a ../pheast/gtextcga-full/output/qtl/$tissue-latent.cis_independent_qtl.txt.gz repo/qtls/ddp-$tissue.cis_independent_qtl.txt.gz
    rsync -a ../pheast/gtex-residual/output/qtl/$tissue-cross_latent.cis_independent_qtl.txt.gz repo/qtls/hybrid-$tissue.cis_independent_qtl.txt.gz
done

echo "=== GTEx TWAS ==="
head -n1 ../twas/output/twas_hits.gtextcga-full-ADPSBQ.tsv | sed 's/^/TISSUE\t/' > repo/twas_hits.gtex-ddp.tsv
head -n1 ../twas/output/twas_hits.gtex-residual-ADPSBQ.tsv | sed 's/^/TISSUE\t/' > repo/twas_hits.gtex-rddp.tsv
cat ../data/gtex/tissues.gtex.txt | while read tissue; do
    tail -n+2 ../twas/output/twas_hits.gtextcga-full-$tissue.tsv | sed "s/^/$tissue\t/" >> repo/twas_hits.gtex-ddp.tsv
    tail -n+2 ../twas/output/twas_hits.gtex-residual-$tissue.tsv | sed "s/^/$tissue\t/" >> repo/twas_hits.gtex-rddp.tsv
done
gzip -f repo/twas_hits.gtex-ddp.tsv
gzip -f repo/twas_hits.gtex-rddp.tsv

echo "=== LaDDR models ==="
tar -cjf repo/laddr_models/info_and_bins.tar.bz2 --owner=0 --group=0 -C ../phenos/gtextcga-full info gene_bins
tar -cjf repo/laddr_models/models-ddp.tar.bz2 --owner=0 --group=0 -C ../phenos/gtextcga-full models
tar -cjf repo/laddr_models/models-rddp.tar.bz2 --owner=0 --group=0 -C ../phenos/gtex-residual models

tree repo > file_tree.txt

echo "=== Done ==="
