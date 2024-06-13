"""
##########################################################################
These rules make the SNV Phasing
##########################################################################
"""
wildcard_constraints:
    bam_name = '|'.join([x for x in BAMQC_BAM_NAME]),
    claire3_model = '|'.join([x for x in NAME_CLAIR3_MODEL]),
    path_calling_tool_params = "clair3.+|pepper_margin_deepvariant",
    filter = '|'.join([x for x in SNPSIFT_FILTERS_NAMES])

"""
This rule makes the phasing of SNV by whatshap
"""

def phasing_input_bam(wildcards):
    index = BAMQC_BAM_NAME.index(wildcards.bam_name)
    return BAMQ_SYMLINK_FILES[index]
    
def phasing_input_vcf(wildcards):
    if config["variant_calling_mode"] == "germline" :
        if wildcards.path_calling_tool_params.startswith("clair3") :
            return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.bam_name) + "/" + str(wildcards.bam_name) + "_merge_output.vcf.gz")
        elif wildcards.path_calling_tool_params.startswith("pepper_margin_deepvariant") :
            return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.bam_name) + "/" + str(wildcards.bam_name) + "_XXX.vcf.gz")
#    elif config["variant_calling_mode"] == "somatic" :
#        return rules.clairs.output.vcf_file
    

rule phasing:
    input:
        bam_file = phasing_input_bam,
        vcf_file = phasing_input_vcf
    output:
        phased_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/whatshap/{bam_name}_phased.vcf.gz")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_WHATSHAP
    params:
        genome_fasta = config["references"]["genome"]
    shell:
        """
        whatshap phase --output {output.phased_vcf_file} --reference {params.genome_fasta} --mapping-quality 20 --ignore-read-groups {input.vcf_file} {input.bam_file}
   
        """


"""
This rule makes the stat of phasing of SNV by whatshap
"""
def input_vcf_gz(wildcards):
    if config["variant_calling_mode"] == "germline" :
        return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.bam_name) + "/whatshap/" + str(wildcards.bam_name) + "_phased.vcf.gz")
#    elif config["variant_calling_mode"] == "somatic" :
#        return rules.clairs.output.vcf_file
   
rule phasing_stat:
    input:
        vcf_file = input_vcf_gz
    output:
        phased_stat_txt = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/whatshap/{bam_name}_phasing_stats.txt"),
        phased_stat_tsv = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/whatshap/{bam_name}_phasing_stats.tsv"),
        phased_block_tsv = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/whatshap/{bam_name}_phasing_haplotype_blocks.tsv"),
        phased_block_gtf = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/whatshap/{bam_name}_phasing_haplotype_blocks.gtf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_WHATSHAP
    shell:
        """
        whatshap stats --tsv{output.phased_stat_tsv} --block-list={output.phased_block_tsv} --gtf={output.phased_block_gtf} {input.vcf_file} > {output.phased_stat_txt}
   
        """

"""
This rule makes the vcf index file
"""
rule tabix_vcf:
    input:
        vcf_gz_file = input_vcf_gz
    output:
        vcf_index = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/whatshap/{bam_name}_phased.vcf.gz.tbi")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_WHATSHAP
    shell:
        """
        tabix {input.vcf_gz_file}
        """

"""
This rule makes the haplotagging of phasing of SNV by whatshap
"""

def haplotagging_input_index(wildcards):
    if config["variant_calling_mode"] == "germline" :
        return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.bam_name) + "/whatshap/" + str(wildcards.bam_name) + "_phased.vcf.gz.tbi")
#    elif config["variant_calling_mode"] == "somatic" :
#        return 

rule phasing_haplotagging:
    input:
        bam_file = phasing_input_bam,
        vcf_gz_file = input_vcf_gz,
        vcf_index_file = haplotagging_input_index
    output:
        haplotag_list_tsv = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/whatshap/{bam_name}_phasing_haplotag_list.tsv"),
        haplotag_bam = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{bam_name}/whatshap/{bam_name}_haplotagged.bam")
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_WHATSHAP
    params:
        genome_fasta = config["references"]["genome"]
    shell:
        """
        whatshap haplotag --output-threads={threads} --ignore-read-groups --output-haplotag-list {output.haplotag_list_tsv} --output {output.haplotag_bam} --reference {params.genome_fasta} {input.vcf_gz_file} {input.bam_file}

        """

