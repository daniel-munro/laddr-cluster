#!/usr/bin/env python3

import sys
import pandas as pd
import gzip
from pathlib import Path
from gtfparse import read_gtf

def parse_gene_id(phenotype_id):
    """Extract gene ID from phenotype ID format '{gene_id}__{etc}'."""
    return phenotype_id.split('__')[0]

def tss_info_from_gtf(gtf_file):
    """Parse GTF file to extract gene information and calculate TSS coordinates."""
    gtf_df = read_gtf(gtf_file, result_type='pandas')
    gene_df = gtf_df[gtf_df['feature'] == 'gene'].copy()
    tss_info = {}
    for _, row in gene_df.iterrows():
        gene_id = row['gene_id']
        strand = row['strand']
        tss = int(row['start']) if strand == '+' else int(row['end'])
        tss_info[gene_id] = tss
    return tss_info

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
        if gene_id in tss_info:
            tss = tss_info[gene_id]
            fields[1] = str(tss - 1)  # convert to 0-based coordinates
            fields[2] = str(tss)
        
        fields[3] = phenotype_id
        updated_lines.append('\t'.join(fields) + '\n')
    
    with open(output_file, 'w') as f:
        f.writelines(updated_lines)

def main():
    if len(sys.argv) != 4:
        print("Usage: python update_pantry_beds.py <input_dir> <gtf_file> <output_dir>")
        sys.exit(1)
    
    input_dir = Path(sys.argv[1])
    gtf_file = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])
    
    output_dir.mkdir(exist_ok=False)
    
    print(f"Parsing GTF file {gtf_file}...")
    tss_info = tss_info_from_gtf(gtf_file)
    print(f"Found {len(tss_info)} genes in GTF file")
    
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
