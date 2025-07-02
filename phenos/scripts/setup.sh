set -e

mkdir -p data
Rscript scripts/get_bigwig_urls.R

mkdir -p ../data/gtex/bigwig
for tissue in $(zcat ../data/gtex/bigwig/bigwigs.tsv.gz | tail -n +2 | cut -f 2 | sort | uniq); do
    echo "Downloading bigwig files for $tissue"
    mkdir -p "../data/gtex/bigwig/$tissue"
    for url in $(zcat ../data/gtex/bigwig/bigwigs.tsv.gz | tail -n +2 | awk -v t="$tissue" '$2 == t {print $3}'); do
        wget --no-verbose -P "../data/gtex/bigwig/$tissue" $url
    done
done

mkdir -p ../data/tcga/bigwig
for study in $(zcat ../data/tcga/bigwig/bigwigs.tsv.gz | tail -n +2 | cut -f 2 | sort | uniq); do
    echo "Downloading bigwig files for $study"
    mkdir -p "../data/tcga/bigwig/$study"
    for url in $(zcat ../data/tcga/bigwig/bigwigs.tsv.gz | tail -n +2 | awk -v s="$study" '$2 == s {print $3}'); do
        wget --no-verbose -P "../data/tcga/bigwig/$study" $url
    done
done

# Initialize projects to train models on 5 GTEx tissues, all tissues, and all tissues plus TCGA studies
latent-rna init gtex5-full --template snakemake
latent-rna init gtex-full --template snakemake
latent-rna init gtextcga-full --template snakemake

# Create coverage manifest files
python scripts/create_coverage_manifest.py \
    -i ../data/gtex/bigwigs.tsv.gz \
    -s <(cut -f1 ../../pantry/GTEx/samples_tissues.txt) \
    --collection gtex \
    -o gtex-full/coverage_manifest.tsv
# Filter gtex-full/coverage_manifest.tsv to get gtex5-full/coverage_manifest.tsv
awk 'NR==FNR {tissues[$1]=1; next} $1 in tissues' ../data/gtex/tissues.gtex5.txt gtex-full/coverage_manifest.tsv > gtex5-full/coverage_manifest.tsv

python scripts/create_coverage_manifest.py \
    -i ../data/tcga/bigwigs.tsv.gz \
    --collection tcga \
    -o gtextcga-full/coverage_manifest.tcga.tsv

# Combine GTEx and TCGA manifests with appropriate prefixes
awk 'BEGIN {OFS="\t"} {print $1, $2, "gtex/bigwig/" $3}' gtex-full/coverage_manifest.tsv > gtextcga-full/coverage_manifest.tsv
awk 'BEGIN {OFS="\t"} {print $1, $2, "tcga/bigwig/" $3}' gtextcga-full/coverage_manifest.tcga.tsv >> gtextcga-full/coverage_manifest.tsv

# (Edit configs)

cd gtex5-full && latent-rna setup && cd ..
cd gtex-full && latent-rna setup && cd ..
cd gtextcga-full && latent-rna setup && cd ..

# (Run phenotyping using Snakemake)

#############
# Geuvadis ##
#############

# Generate bigWig files from BAM files
threads=16
mkdir -p ../data/geuvadis/bigwig
for sample in $(cat todo.txt); do
    echo "Generating bigWig file for $sample"
    bamCoverage -b ../data/geuvadis/align/intermediate/star_out/$sample.Aligned.sortedByCoord.out.bam -o ../data/geuvadis/bigwig/$sample.bw -of bigwig --binSize 1 -p $threads
done

latent-rna init geuvadis-full --template snakemake

cat ../data/geuvadis/samples.txt | awk '{ print "Geuvadis\t" $1 "\t" $1 ".bw" }' > geuvadis-full/coverage_manifest.tsv

