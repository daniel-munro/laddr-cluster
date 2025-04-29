set -e

read_length=75
threads=16
gtf=../../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.chrom.gtf

mkdir -p data/fastq align data/bigwig latent

shuf -n 100 ../../data/geuvadis/samples.txt | sort > data/samples.txt

#####################
## FASTQs and BAMs ##
#####################

for n in 50 75; do
    mkdir -p data/fastq/${n}bp
    while read sample; do
        seqtk trimfq -b 0 -e $((read_length-n)) ../../data/geuvadis/fastq/${sample}_1.fastq.gz | gzip > data/fastq/${n}bp/${sample}_R1.fastq.gz
        seqtk trimfq -b 0 -e $((read_length-n)) ../../data/geuvadis/fastq/${sample}_2.fastq.gz | gzip > data/fastq/${n}bp/${sample}_R2.fastq.gz
    done < data/samples.txt
done

awk '{print $1 "_R1.fastq.gz\t" $1 "_R2.fastq.gz\t" $1}' data/samples.txt > data/fastq/fastq_map.pe.txt
awk '{print $1 "_R1.fastq.gz\t" $1}' data/samples.txt > data/fastq/fastq_map.se1.txt
awk '{print $1 "_R2.fastq.gz\t" $1}' data/samples.txt > data/fastq/fastq_map.se2.txt

for sim in pe-50 pe-75 se1-50 se1-75 se2-50 se2-75; do
    rsync -av ~/tools/Pantry/phenotyping/ align/${sim} --exclude input --exclude intermediate --exclude output --exclude .snakemake
done

## (Edit each config to point to correct inputs and edit Snakefile to run alignment only)

## (Run Pantry)

#######################
## Latent phenotypes ##
#######################

for sim in pe-50 pe-75 se1-50 se1-75 se2-50 se2-75; do
    echo $sim
    while read sample; do
        bamCoverage -b align/${sim}/intermediate/star_out/${sample}.Aligned.sortedByCoord.out.bam -o data/bigwig/${sim}/${sample}.bw -of bigwig --binSize 1 -p $threads
    done < data/samples.txt
    latent-rna init latent/${sim} --template snakemake
    awk -v sim=$sim '{print sim "\t" $1 "\t" $1 ".bw" }' data/samples.txt > latent/${sim}/coverage_manifest.tsv
done

## (For each, edit config.yaml, run latent-rna setup, and run latent phenotyping)

##########
## QTLs ##
##########

rsync -av ~/tools/Pantry/phenotyping/ qtl --exclude input --exclude intermediate --exclude output --exclude .snakemake
mkdir -p qtl/input/phenotypes

## (Edit config)

for sim in pe-50 pe-75 se1-50 se1-75 se2-50 se2-75; do
    python3 ~/tools/Pantry/phenotyping/scripts/assemble_bed.py \
        --type latent \
        --input latent/${sim}/phenotypes/latent_phenos.${sim}.tsv.gz \
        --ref_anno $gtf \
        --output qtl/input/phenotypes/${sim}.unnorm.bed

    python3 ~/tools/Pantry/phenotyping/scripts/normalize_phenotypes.py \
        --input qtl/input/phenotypes/${sim}.unnorm.bed \
        --output qtl/input/phenotypes/${sim}.bed
    bgzip qtl/input/phenotypes/${sim}.bed

    tabix -p bed qtl/input/phenotypes/${sim}.bed.gz

    ## Get latent phenotype groups
    zcat qtl/input/phenotypes/${sim}.bed.gz \
        | tail -n +2 \
        | cut -f4 \
        | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
        > qtl/input/phenotypes/${sim}.phenotype_groups.txt
done

## (Run Pheast)
