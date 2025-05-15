set -e

for projdir in geuvadis-full geuvadis-residual gtex5-full gtex5-residual gtextcga-full gtex-residual; do
    # Confirm projdir does not already exist
    if [ -d "$projdir" ]; then
        echo "Error: $projdir already exists"
        exit 1
    fi
    rsync -av ~/tools/Pantry/pheast/ $projdir --exclude input --exclude intermediate --exclude output --exclude .snakemake
    cp scripts/covariates.smk ${projdir}/steps/
    cp scripts/qtl.smk ${projdir}/steps/
    cp scripts/twas.smk ${projdir}/steps/
done

# (Run Snakemake to set up input files)

# (Edit configs and run Pheast)

