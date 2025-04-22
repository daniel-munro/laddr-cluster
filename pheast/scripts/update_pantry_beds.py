#!/usr/bin/env python3

import sys
import pandas as pd
import gzip
from pathlib import Path

def parse_gene_id(phenotype_id):
    """Extract gene ID from phenotype ID format '{gene_id}__{etc}'."""
    return phenotype_id.split('__')[0]

def update_bed_file(input_file, output_file, tss_info):
    """Update phenotype IDs and TSS coordinates in a BED file."""
    with gzip.open(input_file, 'rt') as f:
        lines = f.readlines()
    
    updated_lines = []
    for line in lines:
        fields = line.strip().split('\t')
        if len(fields) < 4:
            continue
            
        # Update phenotype ID format
        phenotype_id = fields[3]
        if '.' in phenotype_id:
            phenotype_id = phenotype_id.replace('.', '__', 1)
        elif ':' in phenotype_id:
            phenotype_id = phenotype_id.replace(':', '__', 1)
        # Update TSS coordinates if gene is in tss_info
        gene_id = parse_gene_id(phenotype_id)
        if gene_id in tss_info.index:
            tss = tss_info.loc[gene_id, 'tss']
            strand = tss_info.loc[gene_id, 'strand']
            if strand == '+':
                fields[1] = str(tss)  # tss gives coordinate before TSS base
                fields[2] = str(tss + 1)
            else:  # strand == '-'
                fields[1] = str(tss - 1)  # tss gives coordinate after TSS base
                fields[2] = str(tss)
        
        fields[3] = phenotype_id
        updated_lines.append('\t'.join(fields) + '\n')
    
    with open(output_file, 'w') as f:
        f.writelines(updated_lines)

def main():
    if len(sys.argv) != 4:
        print("Usage: python update_pantry_beds.py <input_dir> <tss_info.tsv> <output_dir>")
        sys.exit(1)
    
    input_dir = Path(sys.argv[1])
    tss_info_file = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])
    
    output_dir.mkdir(exist_ok=False)
    
    tss_info = pd.read_csv(tss_info_file, sep='\t', index_col='gene_id')
    
    # Process each BED file
    modalities = ['alt_polyA', 'alt_TSS', 'expression', 'isoforms', 'splicing', 'stability']
    for modality in modalities:
        bed_file = input_dir / f"{modality}.bed.gz"
        output_file = output_dir / f"{modality}.bed"
        print(f"Processing {bed_file}...")
        update_bed_file(bed_file, output_file, tss_info)
        print(f"Created {output_file}")

if __name__ == '__main__':
    main() 
