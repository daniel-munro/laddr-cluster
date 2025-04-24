import argparse
import pandas as pd
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, required=True, help="Input TSV with sample, dataset, url columns")
    parser.add_argument("-s", type=Path, help="File with sample IDs to filter to (optional)")
    parser.add_argument("--collection", choices=['gtex', 'tcga'], help="Either 'gtex' or 'tcga', used to determine how to filter bigWigs")
    parser.add_argument("-o", type=Path, required=True, help="Output manifest file")
    return parser.parse_args()

def main():
    args = parse_args()
    
    df = pd.read_csv(args.i, sep='\t')
    if args.s:
        sample_ids = pd.read_csv(args.s, header=None, names=['id']).id.tolist()
        df = df[df['sample'].isin(sample_ids)]
    df['file_name'] = df.url.apply(lambda x: Path(x).name)
    
    if args.collection == 'gtex':
        # Store total sample count before filtering
        sample_count = df['sample'].nunique()
        # Filter to only include .1.ALL.bw files when multiple files exist per sample
        # (The additional bigWigs are often identical)
        df = df[df.file_name.str.endswith('.1.ALL.bw')]
        # Assert that we haven't lost any samples
        filtered_sample_count = df['sample'].nunique()
        assert sample_count == filtered_sample_count, f"Missing .1.ALL.bw files for {sample_count - filtered_sample_count} samples"
    else:
        # For non-GTEx samples, keep first file when multiple exist per sample
        df = df.groupby('sample').first().reset_index()
        df = df.sort_values(['dataset', 'sample'])
    # Construct the bigwig paths
    df['bigwig_path'] = df.apply(lambda x: Path(x.dataset) / x.file_name, axis=1)
    
    manifest = df[['dataset', 'sample', 'bigwig_path']]
    manifest.to_csv(args.o, sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()
