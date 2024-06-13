"""
##########################################################################
These rules make the SNV Calling for somatic variants 
##########################################################################
"""
wildcard_constraints:
    sample_id = '|'.join([x for x in BAMQC_SAMPLE_NAME])

"""
This rule makes the SNV Calling by clair3 with various models
"""

def clairs_input_normal_bam(wildcards):
    index = BAMQC_BAM_NAME.index(wildcards.sample_id + "_normal")
    return BAMQ_SYMLINK_FILES[index]

def clairs_input_tumor_bam(wildcards):
    index = BAMQC_BAM_NAME.index(wildcards.sample_id + "_tumor")
    return BAMQ_SYMLINK_FILES[index]

rule clairs:
    input:
        normal_bam_file = clairs_input_normal_bam,
        tumor_bam_file = clairs_input_tumor_bam,
        fa_ref = config["references"]["genome"],
    output:
        snv_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_id}/{sample_id}_snv.vcf.gz"),
        indel_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_id}/{sample_id}_indel.vcf.gz")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        path_fa_ref = os.path.dirname(config["references"]["genome"]),
        model = config["clairs"]["model"]
    shell:
        """
        singularity exec \
          -B {OUTPUT_DIR},{params.path_fa_ref} \
          /mnt/beegfs/userdata/n_rabearivelo/containers/clairs_latest.sif \
          /opt/bin/run_clairs \
          --tumor_bam_fn {input.tumor_bam_file} \
          --normal_bam_fn {input.normal_bam_file} \
          --ref_fn {input.fa_ref} \
          --threads {threads} \
          --platform {params.model} \
          --output_dir {OUTPUT_DIR} \
          --output_prefix {wildcards.sample_id}_snv \
          --indel_output_prefix {wildcards.sample_id}_indel \
          --sample_name {wildcards.sample_id} \
          --include_all_ctgs \
          --remove_intermediate_dir \
          --conda_prefix /opt/conda/envs/clairs \
          --enable_indel_calling
  
        """

 