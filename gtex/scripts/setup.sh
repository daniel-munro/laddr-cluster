set -e

mkdir -p data
Rscript scripts/get_bigwig_urls.R

mkdir -p data/bigwig
for tissue in $(zcat data/bigwigs.tsv.gz | tail -n +2 | cut -f 2 | sort | uniq); do
    echo "Downloading bigwig files for $tissue"
    mkdir -p "data/bigwig/$tissue"
    for url in $(zcat data/bigwigs.tsv.gz | tail -n +2 | awk -v t="$tissue" '$2 == t {print $3}'); do
        wget -P "data/bigwig/$tissue" $url
    done
done
cd ..

# Initialize project to train models on 5 GTEx tissues and on all tissues
latent-rna init gtex5 --template snakemake
latent-rna init gtex --template snakemake

# Create coverage manifest files
# zcat data/bigwigs.gtex5.tsv.gz | tail -n +2 | awk -F'\t' '{split($3,a,"/"); print $2"\t"$1"\t../data/bigwig/"$2"/"a[length(a)]}' > gtex5/coverage_manifest.tsv
# zcat data/bigwigs.tsv.gz | tail -n +2 | awk -F'\t' '{split($3,a,"/"); print $2"\t"$1"\t../data/bigwig/"$2"/"a[length(a)]}' > gtex/coverage_manifest.tsv
python scripts/create_coverage_manifest.py \
    -i data/bigwigs.tsv.gz \
    -g ../../pantry/GTEx/geno/ids.txt \
    -o gtex/coverage_manifest.tsv
# Filter gtex/coverage_manifest.tsv to get gtex5/coverage_manifest.tsv
awk 'NR==FNR {tissues[$1]=1; next} $1 in tissues' data/tissues.gtex5.txt gtex/coverage_manifest.tsv > gtex5/coverage_manifest.tsv

# (Edit configs)

cd gtex5
latent-rna setup
cd ../gtex
latent-rna setup
