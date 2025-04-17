set -e

tissue=BRNCTXB
projdir=qtl

rsync -av ~/tools/Pantry/pheast/ $projdir --exclude input --exclude intermediate --exclude output --exclude .snakemake
cp ../pheast/scripts/qtl.smk $projdir/steps/
cp scripts/config_pheast.yml $projdir/config.yml
mkdir -p $projdir/input/phenotypes_original/unnorm

############
## Latent ##
############

gtf=../../ref/human/gencode.v47.annotation.gtf.gz

for prune in 0 40 100; do
    echo "Processing latent-${prune}"

    python3 ~/tools/Pantry/phenotyping/scripts/assemble_bed.py \
        --type latent \
        --input latent-${prune}/phenotypes/latent_phenos.${tissue}.tsv.gz \
        --ref_anno $gtf \
        --output $projdir/input/phenotypes_original/unnorm/latent-${prune}.unnorm.bed

    python3 ~/tools/Pantry/phenotyping/scripts/normalize_phenotypes.py \
        --input $projdir/input/phenotypes_original/unnorm/latent-${prune}.unnorm.bed \
        --output $projdir/input/phenotypes_original/latent-${prune}.bed
    bgzip $projdir/input/phenotypes_original/unnorm/latent-${prune}.unnorm.bed
    bgzip $projdir/input/phenotypes_original/latent-${prune}.bed

    ## Get latent phenotype groups
    zcat $projdir/input/phenotypes_original/latent-${prune}.bed.gz \
        | tail -n +2 \
        | cut -f4 \
        | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
        > $projdir/input/phenotypes_original/latent-${prune}.phenotype_groups.txt
done

############
## Pantry ##
############

for prune in 0 20 40 60 80 100; do
    echo "Processing pantry-${prune}"

    cd pantry-${prune}
    bash scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability
    cd ..

    (zcat pantry-${prune}/output/cross_modality.bed.gz | head -n1; \
        zcat pantry-${prune}/output/cross_modality.bed.gz | tail -n+2 | sed 's/^/chr/') \
        | bgzip > $projdir/input/phenotypes_original/pantry-${prune}.bed.gz
    cp pantry-${prune}/output/cross_modality.phenotype_groups.txt $projdir/input/phenotypes_original/pantry-${prune}.phenotype_groups.txt
done

# Convert to individual IDs and filter samples
python3 ~/tools/Pantry/pheast/scripts/prepare_phenotypes.py \
    --indir $projdir/input/phenotypes_original \
    --out $projdir/input/phenotypes \
    --map ../data/gtex/sample_individual_map.tsv \
    --individuals ../../pantry/GTEx/geno/ids.txt

for type in latent; do # pantry; do
    for prune in 0 40 100; do
        tabix -p bed $projdir/input/phenotypes/${type}-${prune}.bed.gz
    done
done

## Get samples with genotypes
zcat $projdir/input/phenotypes/latent-0.bed.gz \
    | head -n1 \
    | cut -f5- \
    | sed 's/\t/\n/g' \
    > $projdir/input/samples.txt

