"""Process latent phenotyping output and set up for Pheast"""

# Use GTF without gene ID versions
ref_anno = '../../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.chrom.gtf'
modalities = ['alt_polyA', 'alt_TSS', 'expression', 'isoforms', 'splicing', 'stability']
file_types = ['bed.gz', 'bed.gz.tbi', 'phenotype_groups.txt']

with open('../data/gtex/tissues.gtex5.txt') as f:
    tissues_gtex5 = f.read().splitlines()
with open('../data/gtex/tissues.gtex.txt') as f:
    tissues_gtex = f.read().splitlines()

localrules:
    latent_pheno_groups,
    index_bed,

rule all:
    input:
        expand('gtex5-full/input/phenotypes/{tissue}-latent.{file_type}', tissue=tissues_gtex, file_type=file_types),
        expand('gtex5-residual/input/phenotypes/{tissue}-latent.{file_type}', tissue=tissues_gtex5, file_type=file_types),
        expand('gtex5-residual/input/phenotypes/{tissue}-cross_latent.{file_type}', tissue=tissues_gtex5, file_type=file_types),
        expand('gtex-full/input/phenotypes/{tissue}-latent.{file_type}', tissue=tissues_gtex, file_type=file_types),
        expand('gtex-residual/input/phenotypes/{tissue}-latent.{file_type}', tissue=tissues_gtex, file_type=file_types),
        expand('gtex-residual/input/phenotypes/{tissue}-cross_latent.{file_type}', tissue=tissues_gtex, file_type=file_types),
        expand('gtextcga-full/input/phenotypes/{tissue}-latent.{file_type}', tissue=tissues_gtex, file_type=file_types),

rule assemble_latent_bed:
    """Add gene info to latent phenotypes to generate BED file"""
    input:
        phenos = '../phenos/{group}-{mode}/phenotypes/latent_phenos.{tissue}.tsv.gz',
        ref_anno = ref_anno,
    output:
        bed = '{group}-{mode}/input/phenotypes_original/unnorm/{tissue}-latent.unnorm.bed',
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
        bed = '{group}-{mode}/input/phenotypes_original/unnorm/{tissue}-latent.unnorm.bed',
    output:
        bed = '{group}-{mode}/input/phenotypes_original/{tissue}-latent.bed.gz',
    params:
        bed = '{group}-{mode}/input/phenotypes_original/{tissue}-latent.bed',
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
        bed = '{group}-{mode}/input/phenotypes_original/{tissue}-latent.bed.gz',
    output:
        groups = '{group}-{mode}/input/phenotypes_original/{tissue}-latent.phenotype_groups.txt',
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
        latent_bed = '{group}-residual/input/phenotypes_original/{tissue}-latent.bed.gz',
        pantry_beds = expand('../../pantry/phenos/{{tissue}}/output/{modality}.bed.gz', modality = modalities),
    output:
        bed = '{group}-residual/input/phenotypes_original/{tissue}-cross_latent.bed.gz',
        groups = '{group}-residual/input/phenotypes_original/{tissue}-cross_latent.phenotype_groups.txt',
    params:
        tmp_dir = '{group}-residual/input/phenotypes_original/tmp_{tissue}',
        modalities = modalities,
    shell:
        """
        set -e
        mkdir {params.tmp_dir}
        mkdir {params.tmp_dir}/output
        ln -s ../../{wildcards.tissue}-latent.bed.gz {params.tmp_dir}/output/latent_residual.bed.gz
        for modality in {params.modalities}; do
            ln -s ../../../../../../../pantry/phenos/{wildcards.tissue}/output/$modality.bed.gz {params.tmp_dir}/output/$modality.bed.gz
        done
        cd {params.tmp_dir}
        bash ~/tools/Pantry/phenotyping/scripts/combine_modalities.sh \
            {params.modalities} latent_residual
        cd -
        mv {params.tmp_dir}/output/cross_modality.bed.gz {output.bed}
        mv {params.tmp_dir}/output/cross_modality.phenotype_groups.txt {output.groups}
        rm -r {params.tmp_dir}
        """

rule prepare_phenotypes_full:
    """Convert to individual IDs and filter samples"""
    input:
        latent_bed = '{group}-full/input/phenotypes_original/{tissue}-latent.bed.gz',
        latent_groups = '{group}-full/input/phenotypes_original/{tissue}-latent.phenotype_groups.txt',
        sample_map = '../data/gtex/sample_individual_map.tsv',
        individuals = '../../pantry/GTEx/geno/ids.txt',
    output:
        latent_bed = '{group,gtex.*}-full/input/phenotypes/{tissue}-latent.bed.gz',
        latent_groups = '{group,gtex.*}-full/input/phenotypes/{tissue}-latent.phenotype_groups.txt',
    params:
        tmp_in_dir = '{group}-full/input/phenotypes_original/tmp_in_{tissue}',
        tmp_out_dir = '{group}-full/input/phenotypes_original/tmp_out_{tissue}',
        pheno_dir = '{group}-full/input/phenotypes',
    resources:
        mem_mb = 32000,
    shell:
        """
        mkdir {params.tmp_in_dir}
        ln -s ../{wildcards.tissue}-latent.bed.gz {params.tmp_in_dir}/{wildcards.tissue}-latent.bed.gz
        ln -s ../{wildcards.tissue}-latent.phenotype_groups.txt {params.tmp_in_dir}/{wildcards.tissue}-latent.phenotype_groups.txt
        python3 ~/tools/Pantry/pheast/scripts/prepare_phenotypes.py \
            --indir {params.tmp_in_dir} \
            --out {params.tmp_out_dir} \
            --map {input.sample_map} \
            --individuals {input.individuals}
        mv -i {params.tmp_out_dir}/* {params.pheno_dir}
        rm -r {params.tmp_in_dir}
        rmdir {params.tmp_out_dir}
        """

rule prepare_phenotypes_residual:
    """Convert to individual IDs and filter samples"""
    input:
        latent_bed = '{group}-residual/input/phenotypes_original/{tissue}-latent.bed.gz',
        latent_groups = '{group}-residual/input/phenotypes_original/{tissue}-latent.phenotype_groups.txt',
        cross_bed = '{group}-residual/input/phenotypes_original/{tissue}-cross_latent.bed.gz',
        cross_groups = '{group}-residual/input/phenotypes_original/{tissue}-cross_latent.phenotype_groups.txt',
        sample_map = '../data/gtex/sample_individual_map.tsv',
        individuals = '../../pantry/GTEx/geno/ids.txt',
    output:
        latent_bed = '{group,gtex.*}-residual/input/phenotypes/{tissue}-latent.bed.gz',
        latent_groups = '{group,gtex.*}-residual/input/phenotypes/{tissue}-latent.phenotype_groups.txt',
        cross_bed = '{group,gtex.*}-residual/input/phenotypes/{tissue}-cross_latent.bed.gz',
        cross_groups = '{group,gtex.*}-residual/input/phenotypes/{tissue}-cross_latent.phenotype_groups.txt',
    params:
        tmp_in_dir = '{group}-residual/input/phenotypes_original/tmp_in_{tissue}',
        tmp_out_dir = '{group}-residual/input/phenotypes_original/tmp_out_{tissue}',
        pheno_dir = '{group}-residual/input/phenotypes',
    resources:
        mem_mb = 64000,
    shell:
        """
        mkdir {params.tmp_in_dir}
        ln -s ../{wildcards.tissue}-latent.bed.gz {params.tmp_in_dir}/{wildcards.tissue}-latent.bed.gz
        ln -s ../{wildcards.tissue}-latent.phenotype_groups.txt {params.tmp_in_dir}/{wildcards.tissue}-latent.phenotype_groups.txt
        ln -s ../{wildcards.tissue}-cross_latent.bed.gz {params.tmp_in_dir}/{wildcards.tissue}-cross_latent.bed.gz
        ln -s ../{wildcards.tissue}-cross_latent.phenotype_groups.txt {params.tmp_in_dir}/{wildcards.tissue}-cross_latent.phenotype_groups.txt
        python3 ~/tools/Pantry/pheast/scripts/prepare_phenotypes.py \
            --indir {params.tmp_in_dir} \
            --out {params.tmp_out_dir} \
            --map {input.sample_map} \
            --individuals {input.individuals}
        mv -i {params.tmp_out_dir}/* {params.pheno_dir}
        rm -r {params.tmp_in_dir}
        rmdir {params.tmp_out_dir}
        """

rule index_bed:
    """Index a BED file."""
    input:
        bed = '{group}-{mode}/input/phenotypes/{file}.bed.gz',
    output:
        tbi = '{group}-{mode}/input/phenotypes/{file}.bed.gz.tbi'
    shell:
        'tabix {input.bed}'
