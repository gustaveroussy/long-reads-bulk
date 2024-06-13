"""
##########################################################################
These rules make the SNV Annotation
##########################################################################
"""
wildcard_constraints:
    bam_name = '|'.join([x for x in BAMQC_BAM_NAME]),
    claire3_model = '|'.join([x for x in NAME_CLAIR3_MODEL]),
    path_calling_tool_params = "clair3.+|pepper_margin_deepvariant|clairs",
    filter = '|'.join([x for x in SNPSIFT_FILTERS_NAMES]),
    compl = "snv|indel",

"""
This rule makes the annotation of SNV by snpEff
"""

def snpeff_annotation_input_vcf(wildcards):
    print(str(wildcards.path_calling_tool_params))
    if config["variant_calling_mode"] == "germline" :
        if wildcards.path_calling_tool_params.startswith("clair3") :
            return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.bam_name) + "/" + str(wildcards.bam_name) + "_merge_output.vcf.gz")
        elif wildcards.path_calling_tool_params.startswith("pepper_margin_deepvariant") :
            return os.path.normpath(OUTPUT_DIR + "/SNV_Calling" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.bam_name) + "/" + str(wildcards.bam_name) + "_XXX.vcf.gz")
#    elif config["variant_calling_mode"] == "somatic" :
        #os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{wildcards.sample_id}/{wildcards.sample_id}_snv.vcf.gz"),
        #os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{wildcards.sample_id}/{wildcards.sample_id}_indel.vcf.gz")
        #return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{wildcards.sample_id}/{wildcards.sample_id}_{wildcards.compl}.vcf.gz")

rule snpeff_annotation:
    input:
        vcf_file = snpeff_annotation_input_vcf
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/snpEff/{bam_name}_annotated.vcf"),
        snpEff_stat_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/snpEff/{bam_name}_annotated_stats.csv"),
        snpEff_summary_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/snpEff/{bam_name}_annotated_summary.html")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        genome_name = config["references"]["genome_name"]
    shell:
        """
        singularity exec --no-home --bind {OUTPUT_DIR} {SING_ENV_SNPEFF} \
        java -jar /home/snpEff/snpEff.jar \ # to do
        {params.genome_name} \
        {input.vcf_file} \
        -csvStats {output.snpEff_stat_file} \
        -stats {output.snpEff_summary_file} > {output.annotated_vcf_file}
   
        """


"""
This rule makes the annotation of SNV by snpSift with dbsnp
"""

rule snpsift_annotation_dbsnp:
    input:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/snpEff/{bam_name}_annotated.vcf")
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/snpEff/{bam_name}_annotated_dbnsfp.vcf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        database = config["references"]["dbsnp"]
    shell:
        """
        singularity exec --no-home --bind {OUTPUT_DIR} {SING_ENV_SNPEFF} \
        java -jar /home/snpEff/SnpSift.jar \ # to do
        dbnsfp -v -db {params.database} \
        {input.annotated_vcf_file} > {output.annotated_vcf_file}

        """

"""
This rule makes the annotation of SNV by snpSift with clinvar
"""

rule snpsift_annotation_clinvar:
    input:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/snpEff/{bam_name}_annotated_dbnsfp.vcf")
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/snpEff/{bam_name}_annotated_dbnsfp_clinvar.vcf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        database = config["references"]["clinvar"]
    shell:
        """
        singularity exec --no-home --bind {OUTPUT_DIR} {SING_ENV_SNPEFF} \
        java -jar /home/snpEff/SnpSift.jar \ # to do
        annotate -v {params.database} \
        {input.annotated_vcf_file} > {output.annotated_vcf_file}

        """


"""
This rule makes the filtering of SNV by snpSift
"""
def snpsift_filter_params(wildcards):
    index = SNPSIFT_FILTERS_NAMES.index(wildcards.filter)
    return SNPSIFT_FILTERS[index]

rule snpsift_filtering:
    input:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/snpEff/{bam_name}_annotated.vcf")
    output:
        filtered_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/snpEff/{bam_name}_annotated_{filter}.vcf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        filter = snpsift_filter_params
    shell:
        """
        echo 'cat {input.annotated_vcf_file} | \
        java -jar /home/snpEff/SnpSift.jar \ # to do
        filter {params.filter} > {output.filtered_vcf_file}' | singularity exec --no-home --bind {OUTPUT_DIR} {SING_ENV_SNPEFF} bash

        """
