"""List genes with models, along with their batch and number of PCs"""

import sys
import os
import pickle
import pandas as pd
from pathlib import Path
from tqdm import tqdm

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <project_dir> <output_file>")
        sys.exit(1)
    project_dir = Path(sys.argv[1])
    output_file = sys.argv[2]

    # Read number of batches
    n_batches_file = project_dir / 'info' / 'n_batches.txt'
    with open(n_batches_file) as f:
        n_batches = int(f.read().strip())
 
    rows = []
    for batch in tqdm(range(n_batches), desc="Processing batches"):
        models_file = project_dir / 'models' / f'models_batch_{batch}.pickle'
        with open(models_file, 'rb') as f:
            models = pickle.load(f)
        for gene_id, model in models['models'].items():
            # Number of PCs: for PCA, n_components; for FPCA, len(explained_variance_ratio_)
            if hasattr(model.model, 'n_components_'):
                n_pcs = model.model.n_components_
            elif hasattr(model.model, 'explained_variance_ratio_'):
                n_pcs = len(model.model.explained_variance_ratio_)
            else:
                n_pcs = None
            rows.append({'gene_id': gene_id, 'batch': batch, 'n_phenotypes': n_pcs})

    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()
