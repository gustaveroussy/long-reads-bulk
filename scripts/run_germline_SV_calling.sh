#!/bin/bash
################################################################
## Using : sbatch run_structural_variant_calling.sh
################################################################

#SBATCH --job-name=SV_calling_250NG_hg38
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --mem=20G
#SBATCH --partition=shortq
#SBATCH --output=SV_calling_250NG_hg38-%x.%j.out
#SBATCH --error=SV_calling_250NG_hg38-%x.%j.err
 
# Input information
#WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
#SAMPLE_NAME="3582_R10_250NG"
#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.sorted.bam"
#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_T2T-CHM13.sorted.bam"
#REF="hg38"

#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_T2T-CHM13.chr20.sorted.bam"
#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.chr20.sorted.bam"

mkdir -p ${WDIR}/data_output/dorado/${SAMPLE_NAME}/sniffles/

source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate sniffles

#sniffles -i ${BAM} -v ${WDIR}/data_output/${SAMPLE_NAME}/sniffles/${SAMPLE_NAME}_hg38_SV.vcf 
sniffles --allow-overwrite -i ${BAM} -v ${WDIR}/data_output/dorado/${SAMPLE_NAME}/sniffles/${SAMPLE_NAME}_${REF}_SV.vcf

#sniffles -i ${BAM} -v ${WDIR}/data_output/${SAMPLE_NAME}/sniffles/${SAMPLE_NAME}_hg38_SV.chr20.vcf
#sniffles -i ${BAM} -v ${WDIR}/data_output/${SAMPLE_NAME}/sniffles/${SAMPLE_NAME}_T2T-CHM13_SV.chr20.vcf

conda deactivate