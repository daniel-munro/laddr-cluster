"""Successively prune a proportion of the isoforms in a gtf file.

For example, remove a random 20% of the isoforms (and remove any genes with no
remaining isoforms). Then remove another random 20% (of the original amount),
and so on. These will be used to run Pantry and latent-rna, and the results will
demonstrate the effect of isoform sparsity on the results of each method.
"""

import argparse
from collections import defaultdict
import gtfparse
import gzip
import numpy as np
import random
import re
import polars as pl

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', type=str, required=True, help='Path to the gtf file')
    parser.add_argument('--output_prefix', type=str, required=True, help='Prefix for the output gtf files')
    parser.add_argument('--remove_percent', type=int, required=True, help='Percentage (out of 100) of original isoforms to remove at each step')
    parser.add_argument('--steps', type=int, default=3, help='Number of pruning steps to perform')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    args = parser.parse_args()

    # Set random seed for reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)

    gtf = gtfparse.read_gtf(args.gtf)
    
    transcripts = gtf.filter(pl.col('feature') == 'transcript')
    
    gene_to_transcripts = defaultdict(list)
    for row in transcripts.iter_rows(named=True):
        gene_to_transcripts[row['gene_id']].append(row['transcript_id'])
    
    all_transcript_ids = transcripts['transcript_id'].to_list()
    
    total_transcripts = len(all_transcript_ids)
    n_remove_per_step = int(total_transcripts * args.remove_percent / 100)
    
    remaining_transcript_ids = set(all_transcript_ids)
    
    for step in range(args.steps):
        transcripts_to_remove = random.sample(list(remaining_transcript_ids), n_remove_per_step)
        remaining_transcript_ids = remaining_transcript_ids - set(transcripts_to_remove)
        
        # Find genes that still have at least one transcript
        remaining_genes = set()
        for gene, gene_transcripts in gene_to_transcripts.items():
            if any(tid in remaining_transcript_ids for tid in gene_transcripts):
                remaining_genes.add(gene)
        
        # Save the pruned GTF file (filter original lines to avoid having to convert back to GTF format)
        percent_remaining = 100 - (step + 1) * args.remove_percent
        output_file = f"{args.output_prefix}.pruned_{percent_remaining}.gtf"
        # gtf_formatted.write_csv(output_file, separator='\t', include_header=False)
        open_fn = gzip.open if args.gtf.endswith('.gz') else open
        mode = 'rt' if args.gtf.endswith('.gz') else 'r'
        with open_fn(args.gtf, mode) as f:
            with open(output_file, 'w') as out:
                for line in f:
                    gene_id = None
                    match = re.search(r'gene_id "([^"]+)"', line)
                    if match:
                        gene_id = match.group(1)
                    transcript_id = None
                    match = re.search(r'transcript_id "([^"]+)"', line)
                    if match:
                        transcript_id = match.group(1)
                    if gene_id is None or gene_id in remaining_genes:
                        if transcript_id is None or transcript_id in remaining_transcript_ids:
                            out.write(line)
        
        print(f"Step {step+1}: Removed {n_remove_per_step} transcripts ({len(remaining_transcript_ids)} remaining, {percent_remaining}% of original)")
        print(f"  - {len(remaining_genes)} genes with at least one transcript")
        print(f"  - Output saved to {output_file}")

if __name__ == "__main__":
    main()
