set -e

rsync -av ~/tools/Pantry/pheast/ qtl --exclude input --exclude intermediate --exclude output --exclude .snakemake
cp ../pheast/scripts/covariates.smk qtl/steps/
cp ../pheast/scripts/qtl.smk qtl/steps/

# (Run Snakemake to set up input files)

# (Edit config and run Pheast)

