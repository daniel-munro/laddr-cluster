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

mkdir -p ../data/tcga/bigwig
# for study in $(zcat ../data/tcga/bigwigs.tsv.gz | tail -n +2 | cut -f 2 | sort | uniq); do
for study in $(cat ../data/tcga/studies.tcga5.txt); do
    echo "Downloading bigwig files for $study"
    mkdir -p "../data/tcga/bigwig/$study"
    for url in $(zcat ../data/tcga/bigwigs.tsv.gz | tail -n +2 | awk -v s="$study" '$2 == s {print $3}'); do
        wget -P "../data/tcga/bigwig/$study" $url
    done
done

# Initialize project to train models on 5 GTEx tissues and on all tissues, likewise for TCGA
latent-rna init gtex5-full --template snakemake
latent-rna init gtex-full --template snakemake
latent-rna init tcga5-full --template snakemake
latent-rna init tcga-full --template snakemake

# Create coverage manifest files
python scripts/create_coverage_manifest.py \
    -i ../data/gtex/bigwigs.tsv.gz \
    -s <(cut -f1 ../../pantry/GTEx/samples_tissues.txt) \
    --filter-gtex \
    -d ../../data/gtex/bigwig \
    -o gtex-full/coverage_manifest.tsv
# Filter gtex-full/coverage_manifest.tsv to get gtex5-full/coverage_manifest.tsv
awk 'NR==FNR {tissues[$1]=1; next} $1 in tissues' ../data/gtex/tissues.gtex5.txt gtex-full/coverage_manifest.tsv > gtex5-full/coverage_manifest.tsv

python scripts/create_coverage_manifest.py \
    -i ../data/tcga/bigwigs.tsv.gz \
    -d ../../data/tcga/bigwig \
    -o tcga-full/coverage_manifest.tsv
# Filter tcga-full/coverage_manifest.tsv to get tcga5-full/coverage_manifest.tsv
awk 'NR==FNR {studies[$1]=1; next} $1 in studies' ../data/tcga/studies.tcga5.txt tcga-full/coverage_manifest.tsv > tcga5-full/coverage_manifest.tsv

# (Edit configs)

cd gtex5-full && latent-rna setup && cd ..
cd gtex-full && latent-rna setup && cd ..
cd tcga5-full && latent-rna setup && cd ..
cd tcga-full && latent-rna setup && cd ..

# (Run phenotyping using Snakemake)
