import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Convert GTEx sample IDs to individual IDs and filter to genotyped IDs.')
parser.add_argument('--bed', type=str, help='BED file with sample IDs as columnn headers 5-end')
parser.add_argument('--individuals', type=str, help='File with list of individual IDs to keep (can include individuals not in BED file that will be ignored)')
parser.add_argument('--out', type=str, help='Output BED file')
args = parser.parse_args()

df = pd.read_csv(args.bed, sep='\t', dtype=str)

## Convert GTEx sample IDs to individual IDs
ids = df.columns[4:]
# Strip off everything from 2nd dash onwards in sample IDs:
ids = [id.split('-')[0] + '-' + id.split('-')[1] for id in ids]
# Check for duplicates
unique_ids = set(ids)
if len(ids) != len(unique_ids):
    duplicates = [id for id in unique_ids if ids.count(id) > 1]
    print("Found duplicate individual IDs:", duplicates)
    raise AssertionError("Duplicate individual IDs found")
df.columns = list(df.columns[:4]) + ids

## Subset to genotyped individuals
with open(args.individuals, 'r') as f:
    individuals = set(f.read().splitlines())
keep = [col for i, col in enumerate(df.columns) if i < 4 or col in individuals]
df = df[keep]

df.to_csv(args.out, sep='\t', index=False)
