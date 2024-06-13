import pandas as pd
import os

sys.stderr.write("\n############################################################ \n")
sys.stderr.write("\n\n\t Bulk long reads pipeline \n\n")
sys.stderr.write("\n############################################################ \n\n")

### parameters ###################################################################################################################################

sys.stderr.write("\n#################### Setting Parameters ####################\n\n")

#pipeline directory
PIPELINE_DIR = workflow.snakefile
PIPELINE_DIR = PIPELINE_DIR.replace("/Snakefile", "")

#output directory (working directory)
OUTPUT_DIR = os.getcwd()

#set enviroements and parameters of steps
if config["steps"]["bam_qc"]:
  CONDA_ENV_QUALITYMAP = PIPELINE_DIR + "/envs/conda/qualimap.yaml"
  CONDA_ENV_MOSDEPTH = PIPELINE_DIR + "/envs/conda/mosdepth.yaml"
  SING_ENV_NANOPLOT = PIPELINE_DIR + "/envs/singularity/nanoplot_1.42.0.simg"
  
#if config["steps"]["fastq_qc"]:

if config["steps"]["snv_calling"]:
  CONDA_ENV_CLAIR3 = PIPELINE_DIR + "/envs/conda/clair3.yaml"
  SING_ENV_PEPPER_DEEPVARIANT = PIPELINE_DIR + "/envs/singularity/pepper_deepvariant_r0.8.sif"
  NAME_CLAIR3_MODEL = [os.path.basename(model) for model in config["clair3"]["model"]]
  SING_ENV_SNPEFF = PIPELINE_DIR + "/envs/singularity/snpeff_5.2c.simg"
  SNPSIFT_FILTERS = config["snpsift"]["filters"]
  SNPSIFT_FILTERS_NAMES = [filter.replace("'", "").replace("(", "").replace(")", "").replace(" ", "").replace(">", "sup").replace("<", "inf").replace("=", "eq") for filter in SNPSIFT_FILTERS] #to do: find a way to keep "()" if therer are more thant 2 parentheses. Ex: "((QUAL >= 10) && "(QUAL <= 30)) || (FILTER = 'PASS')""

#if config["steps"]["sv_calling"]:
  CONDA_ENV_SNIFFLES = PIPELINE_DIR + "/envs/conda/sniffles.yaml"

if config["steps"]["phasing"]:
  CONDA_ENV_WHATSHAP = PIPELINE_DIR + "/envs/conda/whatshap.yaml"



sys.stderr.write("Parameters validated.\n")

sys.stderr.write("\n################### Checking Design File ###################\n\n")

#read design file
design=pd.read_table(config["design"],sep=",")
#design=pd.read_table("/mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/script/design_bam_germline.tsv",sep=",")

#check if input_format and variant_calling_mode are agree with design format
if config["input_format"] == "bam" and config["variant_calling_mode"] == "germline" : format_design = ['sample_id', 'bam_file']
if config["input_format"] == "bam" and config["variant_calling_mode"] == "somatic" : format_design = ['sample_id', 'bam_file_tumor', 'bam_file_normal']
if config["input_format"] == "fastq" and config["variant_calling_mode"] == "germline" : format_design = ['sample_id', 'upstream_fastq_file', 'downstream_fastq_file']
if config["input_format"] == "fastq" and config["variant_calling_mode"] == "somatic" : format_design = ['sample_id', 'upstream_fastq_file_tumor', 'downstream_fastq_file_tumor', 'upstream_fastq_file_normal', 'downstream_fastq_file_normal']
if set(format_design).issubset(design.columns):
  sys.stderr.write("Design file well formated.\n")
  design=design[format_design]
else: sys.exit("Error in the format of the design file: missing column(s).")

#check if all sample_id are differents
if not len(set(design["sample_id"])) == len(design["sample_id"]): sys.exit("Error: All sample_id have to be different! Check your design file.")

#check if all files are differents on each line
for i in range(0, len(design["sample_id"]), 1):
  if not len(set(design.iloc[i])) == len(design.iloc[i]): sys.exit("Error: All files have to be different! Check your design file.")

#make the data structure
BAMQC_SAMPLE_NAME =  []
BAMQC_BAM_NAME =  []
BAMQ_ORIG_FILE = []
BAMQ_SYMLINK_FILES = []
for col in design.drop(columns="sample_id").columns.tolist():
  for line in range(0,len(design[col].tolist()),1):
    BAMQC_SAMPLE_NAME.append(design["sample_id"].iloc[line])
    BAMQ_ORIG_FILE.append(design[col].iloc[line])
    if col == "bam_file" : SUPPL_NAME=""
    if col == "bam_file_tumor": SUPPL_NAME="_tumor"
    if col == "bam_file_normal": SUPPL_NAME="_normal"
    BAMQC_BAM_NAME.append(design["sample_id"].iloc[line] + SUPPL_NAME)
    BAMQ_SYMLINK_FILES.append(OUTPUT_DIR + "/symlink_input/" + design["sample_id"].iloc[line] + SUPPL_NAME + ".bam")
    


all_files=design.drop(columns="sample_id").stack().tolist()



sys.stderr.write("\n########################### Run ############################\n\n")

### rule all ###################################################################################################################################

include: "rules/rule_all.smk"
rule all:
    input:
        **get_targets()
    message:
        "Pipeline finished!"

### real rules ###################################################################################################################################


# Include rules files
if config["steps"]["bam_qc"]:
    include: "rules/bam_qc.smk"

if config["steps"]["snv_calling"] and config["variant_calling_mode"] == "germline":
    include: "rules/germline_snv_calling.smk"
    include: "rules/snv_annotation.smk"

if config["steps"]["sv_calling"] and config["variant_calling_mode"] == "germline":
    include: "rules/germline_sv_calling.smk"

if config["steps"]["snv_calling"] and config["variant_calling_mode"] == "somatic":
    include: "rules/somatic_snv_calling.smk"
#    include: "rules/snv_annotation.smk"

if config["steps"]["phasing"] and config["variant_calling_mode"] == "germline":
    include: "rules/phasing.smk"


