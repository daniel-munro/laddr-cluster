"""Successively prune a proportion of the non-canonical isoforms in a gtf file.

For example, remove a random 20% of the non-canonical isoforms. Then remove
another random 20% (of the original amount), and so on. These will be used to
run Pantry and latent-rna, and the results will demonstrate the effect of
isoform sparsity on the results of each method.
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
    parser.add_argument('--canonical', type=str, required=True, help='Path to the canonical transcript file')
    parser.add_argument('--output_prefix', type=str, required=True, help='Prefix for the output gtf files')
    parser.add_argument('--remove_percent', type=int, required=True, help='Percentage (out of 100) of original non-canonical isoforms to remove at each step')
    parser.add_argument('--steps', type=int, default=3, help='Number of pruning steps to perform. 100 - remove_percent * steps should be 0 or greater.')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    gtf = gtfparse.read_gtf(args.gtf)
    canonical_df = pl.read_csv(args.canonical, separator='\t', has_header=False)
    canonical_ids = set(tx.split('.')[0] for tx in canonical_df.get_column('column_5'))

    transcripts = gtf.filter(pl.col('feature') == 'transcript')
    transcripts = transcripts.with_columns(
        pl.col('transcript_id').str.split('.').list.first().is_in(canonical_ids).alias('canonical')
    )
    
    n_canonical_per_gene = transcripts.group_by('gene_id').agg(
        pl.col('canonical').sum().alias('n_canonical')
    )
    assert (n_canonical_per_gene.get_column('n_canonical') == 1).all(), \
        "Some genes do not have exactly one canonical transcript"
    
    noncanonical = transcripts.filter(pl.col('canonical').not_())
    noncanonical_ids = noncanonical.get_column('transcript_id').to_list()
    random.shuffle(noncanonical_ids)
    
    n_remove_per_step = int(len(noncanonical_ids) * args.remove_percent / 100)
    to_remove = set()
    
    for step in range(args.steps):
        percent_remaining = 100 - (step + 1) * args.remove_percent
        # If 0% should remain, remove all rather than leaving a few due to rounding
        if args.remove_percent == 0:
            to_remove.update(noncanonical_ids)
            noncanonical_ids = []
        else:
            to_remove.update(noncanonical_ids[:n_remove_per_step])
            noncanonical_ids = noncanonical_ids[n_remove_per_step:]
        
        # Save the pruned GTF file (filter original lines to avoid having to convert back to GTF format)
        output_file = f"{args.output_prefix}.pruned_{percent_remaining}.gtf"
        open_fn = gzip.open if args.gtf.endswith('.gz') else open
        mode = 'rt' if args.gtf.endswith('.gz') else 'r'
        with open_fn(args.gtf, mode) as f:
            with open(output_file, 'w') as out:
                for line in f:
                    transcript_id = None
                    match = re.search(r'transcript_id "([^"]+)"', line)
                    if match:
                        transcript_id = match.group(1)
                    if transcript_id is None or transcript_id not in to_remove:
                        out.write(line)
        
        print(f"Step {step+1}: Removed {n_remove_per_step} transcripts ({len(noncanonical_ids)} remaining, {percent_remaining}% of original)", flush=True)
        print(f"  - Output saved to {output_file}", flush=True)

if __name__ == "__main__":
    main()
