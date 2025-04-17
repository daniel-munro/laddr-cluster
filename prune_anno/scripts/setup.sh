set -e

## Convert chrom names in Ensembl GTF to match those in bigWig files (using Ensembl for compatibility with txrevise)
zcat ../../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.gtf.gz \
    | sed 's/^MT/M/' \
    | sed '/^[^#]/s/^/chr/' \
    > data/Homo_sapiens.GRCh38.113.chr.pruned_100.gtf

python scripts/prune_annotations.py \
    --gtf data/Homo_sapiens.GRCh38.113.chr.pruned_100.gtf \
    --canonical ../../ref/human/knownCanonical.txt.gz \
    --output_prefix data/Homo_sapiens.GRCh38.113.chr \
    --remove_percent 20 \
    --steps 5

grep '^BRNCTXB' ../phenos/gtex5-full/coverage_manifest.tsv > data/coverage_manifest.tsv

############
## Latent ##
############

for prune in 0 20 40 60 80 100; do
    latent-rna init latent-${prune} --template snakemake
done

## Edit each config file to use the new gtf file

## Run `latent-rna setup` in each directory

############
## Pantry ##
############

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

for prune in 0 20 40 60 80 100; do
    rsync -av ~/tools/Pantry/phenotyping/ pantry-${prune} --exclude input --exclude intermediate --exclude output --exclude .snakemake
done

## Edit each config file to use the new gtf file
