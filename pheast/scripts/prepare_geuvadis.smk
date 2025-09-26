"""Process latent phenotyping output and set up for Pheast"""

# Use GTF without gene ID versions
ref_anno = '../../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.chrom.gtf'
modalities = ['alt_polyA', 'alt_TSS', 'expression', 'isoforms', 'splicing', 'stability']
file_types = ['bed.gz', 'bed.gz.tbi', 'phenotype_groups.txt']

localrules:
    latent_pheno_groups,
    prepare_phenotypes_geuvadis,
    index_bed,

rule all:
    input:
        expand('geuvadis/input/phenotypes/geuvadis-full-Geuvadis-latent.{file_type}', file_type=file_types),
        expand('geuvadis/input/phenotypes/geuvadis-residual-Geuvadis-latent.{file_type}', file_type=file_types),
        expand('geuvadis/input/phenotypes/geuvadis-residual-Geuvadis-cross_latent.{file_type}', file_type=file_types),
        expand('geuvadis/input/phenotypes/pantry-Geuvadis-cross_pantry.{file_type}', file_type=file_types),
        expand('geuvadis/input/phenotypes/gtextcga-full-Geuvadis-latent.{file_type}', file_type=file_types),
        expand('geuvadis/input/phenotypes/gtex-residual-Geuvadis-latent.{file_type}', file_type=file_types),
        expand('geuvadis/input/phenotypes/gtex-residual-Geuvadis-cross_latent.{file_type}', file_type=file_types),

rule assemble_latent_bed:
    """Add gene info to latent phenotypes to generate BED file"""
    input:
        phenos = '../phenos/{group}-{mode}/phenotypes/latent_phenos.Geuvadis.tsv.gz',
        ref_anno = ref_anno,
    output:
        bed = 'geuvadis/input/phenotypes/unnorm/{group}-{mode}-Geuvadis-latent.unnorm.bed',
    resources:
        mem_mb = 32000,
    shell:
        """
        python3 ~/tools/Pantry/phenotyping/scripts/assemble_bed.py \
            --type latent \
            --input {input.phenos} \
            --ref_anno {input.ref_anno} \
            --output {output.bed}
        """

rule normalize_latent:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed = 'geuvadis/input/phenotypes/unnorm/{group}-{mode}-Geuvadis-latent.unnorm.bed',
    output:
        bed = 'geuvadis/input/phenotypes/{group}-{mode}-Geuvadis-latent.bed.gz',
    params:
        bed = 'geuvadis/input/phenotypes/{group}-{mode}-Geuvadis-latent.bed',
    resources:
        mem_mb = 48000,
    shell:
        """
        python3 ~/tools/Pantry/phenotyping/scripts/normalize_phenotypes.py \
            --input {input.bed} \
            --output {params.bed}
        bgzip {params.bed}
        """

rule latent_pheno_groups:
    """Group latent phenotypes by gene"""
    input:
        bed = 'geuvadis/input/phenotypes/{group}-{mode}-Geuvadis-latent.bed.gz',
    output:
        groups = 'geuvadis/input/phenotypes/{group}-{mode}-Geuvadis-latent.phenotype_groups.txt',
    shell:
        """
        zcat {input.bed} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
            > {output.groups}
        """

rule combine_modalities_hybrid:
    """Combine latent and Pantry phenotypes"""
    input:
        latent_bed = 'geuvadis/input/phenotypes/{group}-residual-Geuvadis-latent.bed.gz',
        pantry_beds = expand('../../pantry/phenos/Geuvadis/output/{modality}.bed.gz', modality = modalities),
    output:
        bed = 'geuvadis/input/phenotypes/{group}-residual-Geuvadis-cross_latent.bed.gz',
        groups = 'geuvadis/input/phenotypes/{group}-residual-Geuvadis-cross_latent.phenotype_groups.txt',
    params:
        tmp_dir = 'geuvadis/input/phenotypes/tmp_{group}',
        modalities = modalities,
    shell:
        """
        set -e
        mkdir {params.tmp_dir}
        mkdir {params.tmp_dir}/output
        ln -s ../../{wildcards.group}-residual-Geuvadis-latent.bed.gz {params.tmp_dir}/output/latent_residual.bed.gz
        for modality in {params.modalities}; do
            ln -s ../../../../../../../pantry/phenos/Geuvadis/output/$modality.bed.gz {params.tmp_dir}/output/$modality.bed.gz
        done
        cd {params.tmp_dir}
        bash ~/tools/Pantry/phenotyping/scripts/combine_modalities.sh {params.modalities} latent_residual
        cd -
        mv {params.tmp_dir}/output/cross_modality.bed.gz {output.bed}
        mv {params.tmp_dir}/output/cross_modality.phenotype_groups.txt {output.groups}
        rm -r {params.tmp_dir}
        """

rule combine_modalities_pantry:
    """Combine Pantry phenotypes for residual analysis comparison"""
    input:
        pantry_beds = expand('../../pantry/phenos/Geuvadis/output/{modality}.bed.gz', modality = modalities),
    output:
        bed = 'geuvadis/input/phenotypes/pantry-Geuvadis-cross_pantry.bed.gz',
        groups = 'geuvadis/input/phenotypes/pantry-Geuvadis-cross_pantry.phenotype_groups.txt',
    params:
        tmp_dir = 'geuvadis/input/phenotypes/tmp_pantry',
        modalities = modalities,
    shell:
        """
        set -e
        mkdir {params.tmp_dir}
        mkdir {params.tmp_dir}/output
        for modality in {params.modalities}; do
            ln -s ../../../../../../../pantry/phenos/Geuvadis/output/$modality.bed.gz {params.tmp_dir}/output/$modality.bed.gz
        done
        cd {params.tmp_dir}
        bash ~/tools/Pantry/phenotyping/scripts/combine_modalities.sh {params.modalities}
        cd -
        mv {params.tmp_dir}/output/cross_modality.bed.gz {output.bed}
        mv {params.tmp_dir}/output/cross_modality.phenotype_groups.txt {output.groups}
        rm -r {params.tmp_dir}
        """

rule index_bed:
    """Index a BED file."""
    input:
        bed = 'geuvadis/input/phenotypes/{file}.bed.gz',
    output:
        tbi = 'geuvadis/input/phenotypes/{file}.bed.gz.tbi'
    shell:
        'tabix {input.bed}'
