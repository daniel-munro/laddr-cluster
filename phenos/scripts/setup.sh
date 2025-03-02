set -e

mkdir -p data
Rscript scripts/get_bigwig_urls.R

mkdir -p ../data/gtex/bigwig
for tissue in $(zcat ../data/gtex/bigwigs.tsv.gz | tail -n +2 | cut -f 2 | sort | uniq); do
    echo "Downloading bigwig files for $tissue"
    mkdir -p "../data/gtex/bigwig/$tissue"
    for url in $(zcat ../data/gtex/bigwigs.tsv.gz | tail -n +2 | awk -v t="$tissue" '$2 == t {print $3}'); do
        wget -P "../data/gtex/bigwig/$tissue" $url
    done
done
cd ..

# Initialize project to train models on 5 GTEx tissues and on all tissues
latent-rna init gtex5-full --template snakemake
latent-rna init gtex-full --template snakemake

# Create coverage manifest files
python scripts/create_coverage_manifest.py \
    -i ../data/gtex/bigwigs.tsv.gz \
    -s <(cut -f1 ../../pantry/GTEx/samples_tissues.txt) \
    -o gtex5-full/coverage_manifest.tsv
# Filter gtex/coverage_manifest.tsv to get gtex5/coverage_manifest.tsv
awk 'NR==FNR {tissues[$1]=1; next} $1 in tissues' ../data/gtex/tissues.gtex5.txt gtex-full/coverage_manifest.tsv > gtex5-full/coverage_manifest.tsv

# (Edit configs)

cd gtex5-full
latent-rna setup
cd ../gtex-full
latent-rna setup

# (Run phenotyping using Snakemake)
