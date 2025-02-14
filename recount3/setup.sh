set -e

Rscript get_bigwig_urls.R

mkdir -p bigwig
cd bigwig
for url in $(cat ../samples.tsv | tail -n +2 | cut -f 3); do
    wget $url
done
cd ..

latent-rna init brain-blood --template snakemake

# Create coverage manifest file from samples.tsv
cat samples.tsv | tail -n +2 | awk -F'\t' '{split($3,a,"/"); print "brain-blood\t"$1"\t../bigwig/"a[length(a)]}' > brain-blood/coverage_manifest.tsv
