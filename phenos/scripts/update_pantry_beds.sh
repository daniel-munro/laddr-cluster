set -e

# mkdir data/gtex/pantry_phenos

for tissue in $(cat data/gtex/tissues.all54.txt); do
    echo $tissue
	python phenos/scripts/update_pantry_beds.py \
	    ../pantry/GTEx/phenos/$tissue/output/ \
	    ../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.chrom.gtf \
	    data/gtex/pantry_phenos/$tissue
	bgzip data/gtex/pantry_phenos/$tissue/*.bed
done
