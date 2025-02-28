set -e

GTF=../../ref/human/gencode.v47.annotation.gtf.gz

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
latent-rna init phenos-gtex5 --template snakemake
latent-rna init phenos-gtex --template snakemake

# Create coverage manifest files
python scripts/create_coverage_manifest.py \
    -i data/bigwigs.tsv.gz \
    -s <(cut -f1 ../../pantry/GTEx/samples_tissues.txt) \
    -o phenos-gtex/coverage_manifest.tsv
# Filter gtex/coverage_manifest.tsv to get gtex5/coverage_manifest.tsv
awk 'NR==FNR {tissues[$1]=1; next} $1 in tissues' data/tissues.gtex5.txt phenos-gtex/coverage_manifest.tsv > phenos-gtex5/coverage_manifest.tsv

# (Edit configs)

cd phenos-gtex5
latent-rna setup
cd ../phenos-gtex
latent-rna setup

# (Run phenotyping using Snakemake)

############
## Pheast ##
############

mkdir -p pheast-gtex5
# (Copy Pheast template dir into pheast-gtex5 for each tissue, remove example data, edit configs)

for tissue in $(cat data/tissues.gtex5.txt); do
    echo "Processing $tissue"
    cp scripts/qtl.smk pheast-gtex5/$tissue/steps/
    ## Prepare latent phenotypes
    mkdir -p pheast-gtex5/$tissue/input/phenotypes
    python3 ~/tools/Pantry/Project/scripts/assemble_bed.py \
        --type latent \
        --input phenos-gtex5/phenotypes/latent_phenos.$tissue.tsv.gz \
        --ref_anno $GTF \
        --output phenos-gtex5/phenotypes/latent_phenos.$tissue.bed

    python3 scripts/prepare_gtex_bed.py \
        --bed phenos-gtex5/phenotypes/latent_phenos.$tissue.bed \
        --individuals ../../pantry/GTEx/geno/ids.txt \
        --out pheast-gtex5/$tissue/input/phenotypes/latent_full.bed
    bgzip pheast-gtex5/$tissue/input/phenotypes/latent_full.bed
    tabix -p bed pheast-gtex5/$tissue/input/phenotypes/latent_full.bed.gz

    ## Get latent phenotype groups
    zcat pheast-gtex5/$tissue/input/phenotypes/latent_full.bed.gz \
        | tail -n +2 \
        | cut -f4 \
        | awk '{{ g=$1; sub(/:.*$/, "", g); print $1 "\t" g }}' \
        > pheast-gtex5/$tissue/input/phenotypes/latent_full.phenotype_groups.txt

    ## Get samples with genotypes
    zcat pheast-gtex5/$tissue/input/phenotypes/latent_full.bed.gz \
        | head -n1 \
        | cut -f5- \
        | sed 's/\t/\n/g' \
        > pheast-gtex5/$tissue/input/samples.txt
done

# (Run Pheast using Snakemake)
