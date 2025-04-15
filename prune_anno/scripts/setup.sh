set -e

############
## Latent ##
############

# Gencode annotation for latent
python scripts/prune_annotations.py \
    --gtf ../../ref/human/gencode.v47.annotation.gtf.gz \
    --output_prefix data/gencode.v47.annotation \
    --remove_percent 20 \
    --steps 3

cat data/gencode.v47.annotation.pruned_40.gtf \
    | awk '$3 == "gene"' \
    | awk -F'\t' '{ match($9, /gene_id "([^"]+)"/, a); print a[1] }' \
    > data/gencode.v47.annotation.pruned_40.genes.txt

for percent in 40 60 80; do
    latent-rna init latent-${percent} --template snakemake
    grep '^BRNCTXB' ../phenos/gtex5-full/coverage_manifest.tsv > latent-${percent}/coverage_manifest.tsv
done

# Edit each config file to use the new gtf file

# Run `latent-rna setup` in each directory

############
## Pantry ##
############

## Ensembl annotation for Pantry
python scripts/prune_annotations.py \
    --gtf ../../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.gtf.gz \
    --output_prefix data/Homo_sapiens.GRCh38.113.chr \
    --remove_percent 20 \
    --steps 3

cat data/Homo_sapiens.GRCh38.113.chr.pruned_40.gtf \
    | awk '$3 == "gene"' \
    | awk -F'\t' '{ match($9, /gene_id "([^"]+)"/, a); print a[1] }' \
    > data/Homo_sapiens.GRCh38.113.chr.pruned_40.genes.txt

cat ../../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.gtf.gz \
    | awk '$3 == "gene"' \
    | grep 'gene_biotype "protein_coding"' \
    | awk -F'\t' '{ match($9, /gene_id "([^"]+)"/, a); print a[1] }' \
    > data/Homo_sapiens.GRCh38.113.chr.pcg.txt

## Extract fastqs (temporary directory enables parallelization)
extract_fastq() {
    local sample=$1
    local bamdir=/data/hps/assoc/private/gdml/from_rss/scripps_data/gtex/v8/bams
    local bam=${bamdir}/${sample}.Aligned.sortedByCoord.out.patched.md.bam
    local tmpdir=data/fastq/tmp_${sample}
    echo $sample
    python3 scripts/run_SamToFastq.py $bam \
        --prefix $sample \
        --output_dir ${tmpdir} \
        --jar ~/tools/picard.jar
    mv ${tmpdir}/* data/fastq/
    rmdir ${tmpdir}
}
export -f extract_fastq
parallel --halt now,fail=1 -j 8 extract_fastq :::: data/samples.txt

## Create fastq map file
awk '{print $1 "_1.fastq.gz\t" $1 "_2.fastq.gz\t" $1}' data/samples.txt > data/fastq_map.txt

# Edit each config file to use the new gtf file
