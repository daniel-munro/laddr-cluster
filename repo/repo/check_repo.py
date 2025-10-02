from pathlib import Path

tissues_gtex = [
    'ADPSBQ', 'ADPVSC', 'ADRNLG', 'ARTAORT', 'ARTCRN', 'ARTTBL', 'BREAST',
    'BRNACC', 'BRNAMY', 'BRNCDT', 'BRNCHA', 'BRNCHB', 'BRNCTXA', 'BRNCTXB',
    'BRNHPP', 'BRNHPT', 'BRNNCC', 'BRNPTM', 'BRNSNG', 'BRNSPC', 'CLNSGM',
    'CLNTRN', 'ESPGEJ', 'ESPMCS', 'ESPMSL', 'FIBRBLS', 'HRTAA', 'HRTLV',
    'KDNCTX', 'LCL', 'LIVER', 'LUNG', 'MSCLSK', 'NERVET', 'OVARY',
    'PNCREAS', 'PRSTTE', 'PTTARY', 'SKINNS', 'SKINS', 'SLVRYG', 'SNTTRM',
    'SPLEEN', 'STMACH', 'TESTIS', 'THYROID', 'UTERUS', 'VAGINA', 'WHLBLD',
]

files = []

## Metadata and processed results
files_info = ['gwas_metadata.txt', 'LDREF_b38ids_chr.tar.bz2', 'pcg_and_lncrna.tsv', 'tissues-gtex.txt']
files += [f'info/{file}' for file in files_info]
files_processed = [
]
files += [f'processed/{file}' for file in files_processed]

# Geuvadis files
files_geuvadis = [
    'ddp-GEUVADIS.bed.gz',
    'ddp-GEUVADIS.cis_independent_qtl.txt.gz',
    'ddp-GEUVADIS.cis_qtl.txt.gz',
    'ddp-GEUVADIS.covar.plink.tsv',
    'ddp-GEUVADIS.covar.tsv',
    'ddp-GEUVADIS.profile',
    'ddp-GEUVADIS.twas_hits.tsv.gz',
    'hybrid-GEUVADIS.cis_independent_qtl.txt.gz',
    'hybrid-GEUVADIS.cis_qtl.txt.gz',
    'hybrid-GEUVADIS.covar.tsv',
    'qtls.geuvadis-ddp_hybrid_kdp.tsv.gz',
    'rddp-GEUVADIS.bed.gz',
    'rddp-GEUVADIS.covar.plink.tsv',
    'rddp-GEUVADIS.profile',
    'rddp-GEUVADIS.twas_hits.tsv.gz',
]
files += [f'geuvadis/{file}' for file in files_geuvadis]

# LaDDR models
files_models = ['info_and_bins.tar.bz2', 'models-ddp.tar.bz2', 'models-rddp.tar.bz2']
files += [f'laddr_models/{file}' for file in files_models]

files += [
    'qtls.gtex-ddp.tsv.gz',
    'qtls.gtex-hybrid.tsv.gz',
    'qtls.gtex-kdp.tsv.gz',
    'twas_hits.gtex-ddp.tsv.gz',
    'twas_hits.gtex-rddp.tsv.gz',
]

for tissue in tissues_gtex:
    files.append(f'covariates/ddp-{tissue}.covar.tsv')
    files.append(f'covariates/hybrid-{tissue}.covar.tsv')
    files.append(f'covariates/ddp-{tissue}.covar.plink.tsv')
    files.append(f'covariates/rddp-{tissue}.covar.plink.tsv')
    files.append(f'rna_phenotypes/ddp-{tissue}.bed.gz')
    files.append(f'rna_phenotypes/rddp-{tissue}.bed.gz')
    files.append(f'qtls/ddp-{tissue}.cis_qtl.txt.gz')
    files.append(f'qtls/hybrid-{tissue}.cis_qtl.txt.gz')
    files.append(f'qtls/ddp-{tissue}.cis_independent_qtl.txt.gz')
    files.append(f'qtls/hybrid-{tissue}.cis_independent_qtl.txt.gz')
    files.append(f'twas_model_summaries/ddp-{tissue}.profile')
    files.append(f'twas_model_summaries/rddp-{tissue}.profile')

files = [Path(file) for file in files]
for file in files:
    assert file.exists(), file
print(f'All {len(files)} files present.')

files += [Path(file) for file in ['README.md', 'check_repo.py']]
for file in Path('.').glob('**/*'):
    if file.is_file() and file not in files:
        print(f'Extra file: {file}')
