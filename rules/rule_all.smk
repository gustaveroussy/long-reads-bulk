"""
##########################################################################
This function make the input of rule all
##########################################################################
"""

def get_targets():
  targets = {}
  if config["steps"]["bam_qc"]:
      targets["bam_qc"]=[
        #qualimap
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/qualimapReport.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/raw_data_qualimapReport/genome_fraction_coverage.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/genome_results.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/raw_data_qualimapReport/coverage_histogram.txt"), bam_name=BAMQC_BAM_NAME),
        #mosdepth
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q15.mosdepth.global.dist.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q15.mosdepth.region.dist.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q15.mosdepth.summary.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q15.per-base.bed.gz"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q15.per-base.bed.gz.csi"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q15.regions.bed.gz"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q15.regions.bed.gz.csi"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q20.mosdepth.global.dist.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q20.mosdepth.region.dist.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q20.mosdepth.summary.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q20.per-base.bed.gz"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q20.per-base.bed.gz.csi"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q20.regions.bed.gz"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_for_CNV_Q20.regions.bed.gz.csi"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}.mosdepth.global.dist.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}.mosdepth.region.dist.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}.mosdepth.summary.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}.per-base.bed.gz"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}.per-base.bed.gz.csi"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}.regions.bed.gz"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}.regions.bed.gz.csi"), bam_name=BAMQC_BAM_NAME),
        #nanoplot
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_dot.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_dot.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_kde.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_kde.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_dot.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_dot.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_kde.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_kde.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_dot.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_dot.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_kde.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_kde.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_dot.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_dot.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_kde.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_kde.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_NanoPlot_20240404_2323.log"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_NanoPlot-report.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_NanoStats.txt"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedHistogramReadlength.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedHistogramReadlength.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedLogTransformed_HistogramReadlength.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedLogTransformed_HistogramReadlength.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_dot.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_dot.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_kde.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_kde.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_dot.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_dot.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_kde.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_kde.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedHistogramReadlength.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedHistogramReadlength.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedLogTransformed_HistogramReadlength.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedLogTransformed_HistogramReadlength.png"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Yield_By_Length.html"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Yield_By_Length.png"), bam_name=BAMQC_BAM_NAME)
        ]
  if config["variant_calling_mode"] == "germline" and config["steps"]["snv_calling"]:
      targets["germline_snv_calling"]=[
        #clair3
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/{bam_name}_merge_output.vcf.gz"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/{bam_name}_merge_output.vcf.gz.tbi"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL)
        #pepper_margin_deepvariant
        ]
  if config["variant_calling_mode"] == "somatic" and config["steps"]["snv_calling"]:
      targets["somatic_snv_calling"]=[
        expand( os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_id}/{sample_id}_snv.vcf.gz"), sample_id=BAMQC_SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_id}/{sample_id}_indel.vcf.gz"), sample_id=BAMQC_SAMPLE_NAME)
        ]
  if config["steps"]["snv_calling"] and config["variant_calling_mode"] == "germline": #to change for only: if config["steps"]["snv_calling"]:
      targets["snv_annotation"]=[
        #clair3 & snpEff & snpSift
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/snpEff/{bam_name}_annotated.vcf"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/snpEff/{bam_name}_annotated_stats.csv"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/snpEff/{bam_name}_annotated_summary.html"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/snpEff/{bam_name}_annotated_dbnsfp.vcf"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/snpEff/{bam_name}_annotated_dbnsfp_clinvar.vcf"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/snpEff/{bam_name}_annotated_{filter}.vcf"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL, filter=SNPSIFT_FILTERS_NAMES)
        #pepper_margin_deepvariant & snpEff & snpSift
        #clairs & snpEff & snpSift
        ]
  if config["steps"]["phasing"] and config["variant_calling_mode"] == "germline": #to change for only: if config["steps"]["phasing"]:
      targets["phasing"]=[
        #clair3
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/whatshap/{bam_name}_phased.vcf.gz"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/whatshap/{bam_name}_phasing_stats.txt"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/whatshap/{bam_name}_phasing_stats.tsv"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/whatshap/{bam_name}_phasing_haplotype_blocks.tsv"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/whatshap/{bam_name}_phasing_haplotype_blocks.gtf"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/whatshap/{bam_name}_phased.vcf.gz.tbi"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/whatshap/{bam_name}_phasing_haplotag_list.tsv"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL),
        expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{bam_name}/whatshap/{bam_name}_haplotagged.bam"), bam_name=BAMQC_BAM_NAME, claire3_model=NAME_CLAIR3_MODEL)    
        #pepper_margin_deepvariant
        ]
  if config["variant_calling_mode"] == "germline" and config["steps"]["sv_calling"]:
      targets["germline_sv_calling"]=[
        #sniffles
        expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/{bam_name}/{bam_name}_SV.vcf"), bam_name=BAMQC_BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/{bam_name}/{bam_name}_SV.snf"), bam_name=BAMQC_BAM_NAME),
        #sniffles2-plot
        ]
      if len(BAMQC_BAM_NAME) > 1 :
        targets["germline_sv_calling"].append(os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/all_samples_SV.vcf"))
        
        
  return targets