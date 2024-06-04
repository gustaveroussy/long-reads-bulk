#!/bin/bash
################################################################
## Variant Calling Tools Testing
## Using : sbatch run_variant_calling.sh
################################################################
 
#SBATCH --job-name=3582_1000NG_vs_HG38_HAC_clair3_default_10threads_1CPU_16Gmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=10
#SBATCH --mem=16G
#SBATCH --partition=longq

# Input information
#WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"

#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3700_R10"
#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#SAMPLE_NAME="3582_R10_250NG"
#SAMPLE_NAME="3582_R10_1000NG"

#INPUT="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.sorted.bam"
#INPUT="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_T2T-CHM13.sorted.bam"

# R9 model - Guppy 6.5.7 - SUP & HAC : https://github.com/HKU-BAL/Clair3/blob/main/docs/guppy5_20220113.md "use the sup model for hac data. But use sup data for the best performance." + You can use the Guppy5 model r941_prom_sup_g5014 for Guppy6 data.
# Clair3 included model
#CLAIR3_MODEL="r941_prom_sup_g5014"

# R10 model - Guppy 6.5.7 - HAC
# Rerio model (considered deprecated since it's using Guppy), use the latest when analyzing Dorado basecalled data
#CLAIR3_MODEL="/mnt/beegfs/userdata/n_rabearivelo/rerio/clair3_models/r1041_e82_400bps_hac_g632"
# Latest R10 model from Rerio for SUP / downloaded on 2024-02-25
CLAIR3_MODEL="/mnt/beegfs/userdata/n_rabearivelo/rerio/clair3_models/r1041_e82_400bps_sup_v430"
CLAIR3_THREADS="10"

# Reference genomes (FASTA)
#HG38="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta"
#T2T="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/T2T-CHM13v2.0/Homo_sapiens-GCA_009914755.4-unmasked.fa"

# Clair3 - SNP & Long indels Calling
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate clair3

# Default required parameters
#Included Clair3 models
#run_clair3.sh --bam_fn=${INPUT} --ref_fn=${HG38} --threads=${CLAIR3_THREADS} --platform="ont" --model_path="${CONDA_PREFIX}/bin/models/${CLAIR3_MODEL}" --output="${WDIR}/data_output/${SAMPLE_NAME}/clair3/" --include_all_ctgs
#Rerio models
run_clair3.sh --bam_fn=${BAM} --ref_fn=${ref_fasta} --threads=${CLAIR3_THREADS} --platform="ont" --model_path="${CLAIR3_MODEL}" --output="${WDIR}/data_output/dorado/${SAMPLE_NAME}/clair3_HAC_${REF}/" --include_all_ctgs

#--threads: Max threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The chunks being processed simultaneously is ceil($threads/4)*3. 3 is the overloading factor.
#If 1 cpu/task & 16G memory + CLAIR3_THREADS="20" : [WARNING] Threads setting exceeds maximum available threads 2, set threads=2
#"We find that Clair3 is sensitive to high read depth: variant calling performance can suffer when read depth is excessive." https://labs.epi2me.io/giab-2023.05/

# Add long indels detection (EXPERIMENTAL)
#run_clair3.sh --bam_fn=${INPUT} --ref_fn=${HG38} --threads=${CLAIR3_THREADS} --platform="ont" --model_path="${CONDA_PREFIX}/bin/models/${CLAIR3_MODEL}" --output="${WDIR}/data_output/clair3/" --enable_long_indel

conda deactivate
