set -e

# mkdir data/gtex/pantry_phenos

for tissue in $(cat todo.txt); do
    echo $tissue
	python pheast/scripts/update_pantry_beds.py \
	    ../pantry/GTEx/phenos/$tissue/output/ \
	    ../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.chrom.gtf \
	    data/gtex/pantry_phenos/$tissue
	bgzip data/gtex/pantry_phenos/$tissue/*.bed
done
