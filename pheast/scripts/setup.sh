set -e

group=$1
mode=$2
tissue=$3
#gtf=../../ref/human/gencode.v47.annotation.gtf.gz
# Use GTF without gene ID versions
gtf=../../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.chrom.gtf

projdir="${group}-${mode}-${tissue}"

# Confirm projdir does not already exist
if [ -d "$projdir" ]; then
    echo "Error: $projdir already exists"
    exit 1
fi
rsync -av ~/tools/Pantry/pheast/ $projdir --exclude input --exclude intermediate --exclude output --exclude .snakemake
cp scripts/covariates.smk ${projdir}/steps/
cp scripts/qtl.smk ${projdir}/steps/

# (Edit config)

echo "Processing $tissue"

mkdir -p ${projdir}/input/phenotypes_original

python3 ~/tools/Pantry/phenotyping/scripts/assemble_bed.py \
    --type latent \
    --input ../phenos/${group}-${mode}/phenotypes/latent_phenos.${tissue}.tsv.gz \
    --ref_anno $gtf \
    --output ${projdir}/input/phenotypes_original/latent_${mode}.unnorm.bed

python3 ~/tools/Pantry/phenotyping/scripts/normalize_phenotypes.py \
    --input ${projdir}/input/phenotypes_original/latent_${mode}.unnorm.bed \
    --output ${projdir}/input/phenotypes_original/latent_${mode}.bed
bgzip ${projdir}/input/phenotypes_original/latent_${mode}.bed

# Convert to individual IDs and filter samples
python3 ~/tools/Pantry/pheast/scripts/prepare_phenotypes.py \
    --indir ${projdir}/input/phenotypes_original \
    --out ${projdir}/input/phenotypes \
    --map ../data/gtex/sample_individual_map.tsv \
    --individuals ../../pantry/GTEx/geno/ids.txt

tabix -p bed ${projdir}/input/phenotypes/latent_${mode}.bed.gz

## Get latent phenotype groups
zcat ${projdir}/input/phenotypes/latent_${mode}.bed.gz \
    | tail -n +2 \
    | cut -f4 \
    | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
    > ${projdir}/input/phenotypes/latent_${mode}.phenotype_groups.txt

## Get samples with genotypes
zcat ${projdir}/input/phenotypes/latent_${mode}.bed.gz \
    | head -n1 \
    | cut -f5- \
    | sed 's/\t/\n/g' \
    > ${projdir}/input/samples.txt

# (Run Pheast using Snakemake)

##############
## Residual ##
##############

## If mode is 'residual', then also generate cross-modality phenotypes
if [ "$mode" = "residual" ]; then
    echo "Generating cross-modality phenotypes"
    cd ${projdir}/input
    mkdir output
    for modality in alt_polyA alt_TSS expression isoforms splicing stability; do
        ln -s ../../../../data/gtex/pantry_phenos/${tissue}/${modality}.bed.gz output/${modality}.bed.gz
    done
    ln -s ../phenotypes_original/latent_${mode}.bed.gz output/latent_residual.bed.gz
    bash ~/tools/Pantry/phenotyping/scripts/combine_modalities.sh \
        alt_polyA alt_TSS expression isoforms splicing stability latent_residual
    cd ../..
    for modality in alt_polyA alt_TSS expression isoforms splicing stability latent_residual; do
        rm ${projdir}/input/output/${modality}.bed.gz
    done

    python3 ~/tools/Pantry/pheast/scripts/prepare_phenotypes.py \
        --indir ${projdir}/input/output \
        --out ${projdir}/input/phenotypes2 \
        --map ../data/gtex/sample_individual_map.tsv \
        --individuals ../../pantry/GTEx/geno/ids.txt
    mv -i ${projdir}/input/phenotypes2/* ${projdir}/input/phenotypes/
    rmdir ${projdir}/input/phenotypes2
    rm -r ${projdir}/input/output
fi
