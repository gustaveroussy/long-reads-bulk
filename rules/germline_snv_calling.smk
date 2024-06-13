"""
##########################################################################
These rules make the SNV Calling for germline variants
##########################################################################
"""
wildcard_constraints:
    bam_name = '|'.join([x for x in BAMQC_BAM_NAME]),
    claire3_model = '|'.join([x for x in NAME_CLAIR3_MODEL])

"""
This rule makes the SNV Calling by clair3 with various models
"""

def clair3_input_bam(wildcards):
    index = BAMQC_BAM_NAME.index(wildcards.bam_name)
    return BAMQ_SYMLINK_FILES[index]

def clair3_input_model_path(wildcards):
    index = NAME_CLAIR3_MODEL.index(wildcards.claire3_model)
    return config["clair3"]["model"][index]


rule clair3:
    input:
        bam_file = clair3_input_bam,
        fa_ref = config["references"]["genome"],
        clair3_path = clair3_input_model_path
    output:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/{bam_name}_merge_output.vcf.gz"),
        vcf_tbi = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/{bam_name}_merge_output.vcf.gz.tbi"),
        full_al_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/{bam_name}_full_alignment.vcf.gz"),
        full_al_vcf_tbi = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/{bam_name}_full_alignment.vcf.gz.tbi"),
        pil_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/{bam_name}_pileup.vcf.gz"),
        pil_vcf_tbi = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/{bam_name}_pileup.vcf.gz.tbi")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_CLAIR3
    shell:
        """
        run_clair3.sh --bam_fn={input.bam_file} --ref_fn={input.fa_ref} --threads={threads} --platform="ont" --model_path={input.clair3_path} --output={OUTPUT_DIR}/SNV_Calling/clair3/{wildcards.claire3_model}/{wildcards.bam_name}/{wildcards.bam_name}_ --include_all_ctgs

        """

    
"""
This rule makes the SNV Calling by pepper_margin_deepvariant
"""

def pepper_margin_deepvariant_input_bam(wildcards):
    index = BAMQC_BAM_NAME.index(wildcards.bam_name)
    return BAMQ_SYMLINK_FILES[index]

rule pepper_margin_deepvariant:
    input:
        bam_file = pepper_margin_deepvariant_input_bam,
        fa_ref = config["references"]["genome"]
    output:
        os.path.normpath(OUTPUT_DIR + "/SNV_Calling/pepper_margin_deepvariant/{bam_name}/{bam_name}_XXX.vcf.gz") # to do
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        path_fa_ref = os.path.dirname(config["references"]["genome"])
    shell:
        """
        singularity exec --no-home --bind {OUTPUT_DIR},{params.path_fa_ref} {SING_ENV_PEPPER_DEEPVARIANT} run_pepper_margin_deepvariant call_variant \
        --bam {input.bam_file} \
        --fasta {input.fa_ref} \
        --output_dir {OUTPUT_DIR}/SNV_Calling/pepper_margin_deepvariant/{wildcards.bam_name}/ \
        --output_prefix {wildcards.bam_name}_ \
        --threads {threads} \
        --sample_name {wildcards.bam_name} \
        --ont_r10_q20

        """
