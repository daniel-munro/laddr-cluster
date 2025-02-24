import argparse
import pandas as pd
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input TSV with sample, tissue, url columns")
    parser.add_argument("-g", help="File with GTEx individual IDs for individuals with genotypes to include")
    parser.add_argument("-o", help="Output manifest file")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Read individual IDs included in the genotypes
    individual_ids = pd.read_csv(args.g, header=None, names=['individual']).individual.tolist()
    
    df = pd.read_csv(args.i, sep='\t')
    
    df['individual'] = df['sample'].str.split('-').str[:2].str.join('-')
    # Only include individuals with genotypes
    df = df[df.individual.isin(individual_ids)]
    
    df['file_name'] = df.url.apply(lambda x: Path(x).name)
    
    # Store total sample count before filtering
    sample_count = df['sample'].nunique()
    # Filter to only include .1.ALL.bw files when multiple files exist per sample
    # (The additional bigWigs are often identical)
    df = df[df.file_name.str.endswith('.1.ALL.bw')]
    # Assert that we haven't lost any samples
    filtered_sample_count = df['sample'].nunique()
    assert sample_count == filtered_sample_count, f"Missing .1.ALL.bw files for {sample_count - filtered_sample_count} samples"
    
    # Construct the bigwig paths
    df['bigwig_path'] = df.apply(
        lambda x: f"../data/bigwig/{x.tissue}/{x.file_name}",
        axis=1
    )
    
    manifest = df[['tissue', 'sample', 'bigwig_path']]
    manifest.to_csv(args.o, sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()
