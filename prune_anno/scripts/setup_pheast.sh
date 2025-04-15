set -e

# prune is e.g. 40, 60, 80
prune=$1
gtf=../../ref/human/gencode.v47.annotation.gtf.gz
tissue=BRNCTXB

mkdir -p pheast

############
## Latent ##
############

projdir="pheast/latent-${prune}"

# Confirm projdir does not already exist
if [ -d "$projdir" ]; then
    echo "Error: $projdir already exists"
    exit 1
fi
rsync -av ~/tools/Pantry/pheast/ $projdir --exclude input --exclude intermediate --exclude output
cp ../pheast/scripts/qtl.smk $projdir/steps/
cp scripts/config_latent.yml $projdir/config.yml

# (Remove example data, edit config)
rm -rf $projdir/.snakemake

echo "Processing latent-${prune}"

mkdir -p $projdir/input/phenotypes_original

python3 ~/tools/Pantry/phenotyping/scripts/assemble_bed.py \
    --type latent \
    --input latent-${prune}/phenotypes/latent_phenos.${tissue}.tsv.gz \
    --ref_anno $gtf \
    --output $projdir/input/phenotypes_original/latent_full.unnorm.bed

python3 ~/tools/Pantry/phenotyping/scripts/normalize_phenotypes.py \
    --input $projdir/input/phenotypes_original/latent_full.unnorm.bed \
    --output $projdir/input/phenotypes_original/latent_full.bed
bgzip $projdir/input/phenotypes_original/latent_full.bed

# Convert to individual IDs and filter samples
python3 ~/tools/Pantry/pheast/scripts/prepare_phenotypes.py \
    --indir $projdir/input/phenotypes_original \
    --out $projdir/input/phenotypes \
    --map ../data/gtex/sample_individual_map.tsv \
    --individuals ../../pantry/GTEx/geno/ids.txt

tabix -p bed $projdir/input/phenotypes/latent_full.bed.gz

## Get latent phenotype groups
zcat $projdir/input/phenotypes/latent_full.bed.gz \
    | tail -n +2 \
    | cut -f4 \
    | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
    > $projdir/input/phenotypes/latent_full.phenotype_groups.txt

## Get samples with genotypes
zcat $projdir/input/phenotypes/latent_full.bed.gz \
    | head -n1 \
    | cut -f5- \
    | sed 's/\t/\n/g' \
    > $projdir/input/samples.txt

# (Run Pheast using Snakemake)

############
## Pantry ##
############

projdir="pheast/pantry-${prune}"

# Confirm projdir does not already exist
if [ -d "$projdir" ]; then
    echo "Error: $projdir already exists"
    exit 1
fi
rsync -av ~/tools/Pantry/pheast/ $projdir --exclude input --exclude intermediate --exclude output
cp ../pheast/scripts/qtl.smk $projdir/steps/
cp scripts/config_pantry.yml $projdir/config.yml

# (Remove example data, edit config)
rm -rf $projdir/.snakemake

echo "Processing pantry-${prune}"

cd pantry-${prune}
bash scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability
cd ..

mkdir -p $projdir/input/phenotypes_original
(zcat pantry-${prune}/output/cross_modality.bed.gz | head -n1; \
    zcat pantry-${prune}/output/cross_modality.bed.gz | tail -n+2 | sed 's/^/chr/') \
    | bgzip > $projdir/input/phenotypes_original/cross_modality.bed.gz
cp pantry-${prune}/output/cross_modality.phenotype_groups.txt $projdir/input/phenotypes_original/

# Convert to individual IDs and filter samples
python3 ~/tools/Pantry/pheast/scripts/prepare_phenotypes.py \
    --indir $projdir/input/phenotypes_original \
    --out $projdir/input/phenotypes \
    --map ../data/gtex/sample_individual_map.tsv \
    --individuals ../../pantry/GTEx/geno/ids.txt

rm -r $projdir/input/phenotypes_original

tabix -p bed $projdir/input/phenotypes/cross_modality.bed.gz

## Get samples with genotypes
zcat $projdir/input/phenotypes/cross_modality.bed.gz \
    | head -n1 \
    | cut -f5- \
    | sed 's/\t/\n/g' \
    > $projdir/input/samples.txt

# (Run Pheast using Snakemake)
