# To run another trial of all parameters:
# Copy config and directories that don't need to be recomputed
trial=trial1
mkdir $trial
rsync -av \
    --exclude="*/covg_norm" \
    --exclude="*/gene_bins" \
    --exclude="*/models" \
    --exclude="*/phenotypes" \
    trial0/binning-* $trial/
rsync -av \
    --exclude="*/models" \
    --exclude="*/phenotypes" \
    trial0/model-* $trial/

## To avoid irrelevant variability in model performance, use the same info/bins/covg files from trial0 for all trials.
## And to properly include variability in the adaptive binning, regenerate info for all a2 and a3 binning runs each trial.
