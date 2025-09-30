#!/usr/bin/env python3
"""Compute pairwise Pearson correlations of TWAS Z-scores across traits per tissue.

Requirements (relative to script dir):
    ../data/gtex/tissues.gtex.txt  (one tissue per line)
    ../input/gwas/gwas_metadata.txt (has column Tag with trait IDs)
    ../output/gtextcga-full-{tissue}/fusion.gtextcga-full-{tissue}.{trait}.tsv (has TWAS.Z column)

Assumes: all fusion files for a tissue share identical phenotype ordering so correlation can be done on aligned rows without checking IDs.
Output: TSV columns: tissue, trait1, trait2, r
"""

from pathlib import Path
import argparse
import itertools
import sys
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Pairwise TWAS Z correlations per tissue")
    base = Path(__file__).parent
    p.add_argument('--tissues', type=Path, default=Path('../data/gtex/tissues.gtex.txt'))
    p.add_argument('--gwas-metadata', type=Path, default=Path('input/gwas/gwas_metadata.txt'))
    p.add_argument('--output-root', type=Path, default=Path('output'))
    p.add_argument('--out', type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    tissues = [t for t in args.tissues.read_text().splitlines() if t.strip()]
    traits = pd.read_csv(args.gwas_metadata, sep='\t', dtype=str)['Tag'].dropna().unique().tolist()

    rows = []
    for tissue in tissues:
        print(tissue, flush=True)
        cols = {}
        prefix = args.output_root / f"gtextcga-full-{tissue}"
        for trait in traits:
            fp = prefix / f"fusion.gtextcga-full-{tissue}.{trait}.tsv"
            df = pd.read_csv(fp, sep='\t', usecols=['TWAS.Z'])
            cols[trait] = df['TWAS.Z']
        if len(cols) < 2:
            continue
        z = pd.DataFrame(cols)
        corr = z.corr(method='pearson')
        for t1, t2 in itertools.combinations(corr.columns, 2):
            rows.append({'tissue': tissue, 'trait1': t1, 'trait2': t2, 'r': corr.at[t1, t2]})

    out_df = pd.DataFrame(rows, columns=['tissue', 'trait1', 'trait2', 'r'])
    out_df.to_csv(args.out, sep='\t', index=False)

if __name__ == '__main__':
    main()

