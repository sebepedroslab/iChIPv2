#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name star_align
#SBATCH --qos=long
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1-12
#SBATCH --output=/no_backup/asebe/smontgomery/logs/%x.%a.out
#SBATCH --error=/no_backup/asebe/smontgomery/logs/%x.%a.err
##################################
# make bash behave more robustly #
##################################
# set -e
# set -u
# set -o pipefail
# === begin ENVIRONMENT SETUP ===
#1. Species; chromsize calculated from scaffold/chromosome sizes in fasta.fai
##Loop through species
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]]; then
    species=Acas #Acanthamoeba castellani; 14 marks (H3,input)
    chromsize=42019824
elif [[ ${SLURM_ARRAY_TASK_ID} == 2 ]]; then
    species=Atha #Arabidopsis thaliana; 11 marks (H3,input)
    chromsize=119667750
elif [[ ${SLURM_ARRAY_TASK_ID} == 3 ]]; then
    species=Ddis #Dictyostelium discoideum; 14 marks (H3, input)
    chromsize=34134454
elif [[ ${SLURM_ARRAY_TASK_ID} == 4 ]]; then
    species=Ngru #Naegleria gruberii; 10 marks (H3, input)
    chromsize=40964085
elif [[ ${SLURM_ARRAY_TASK_ID} == 5 ]]; then
    species=Nvec #Nematostella vectensis; 14 marks (H3, input)
    chromsize=269418438
elif [[ ${SLURM_ARRAY_TASK_ID} == 6 ]]; then
    species=Ppat #Physcomitrium patens; 11 marks (H3, input)
    chromsize=471852792
elif [[ ${SLURM_ARRAY_TASK_ID} == 7 ]]; then
    species=Spun #Spizellomyces punctatus; 14 marks (H3,input)
    chromsize=24131112
elif [[ ${SLURM_ARRAY_TASK_ID} == 8 ]]; then
    species=Tthe #Tetrahymena thermophila; 13 marks (H3, input)
    chromsize=103014375
elif [[ ${SLURM_ARRAY_TASK_ID} == 9 ]]; then
    species=Bnat #Bigelowiella natans; 14 marks (H3, input)
    chromsize=94701163
elif [[ ${SLURM_ARRAY_TASK_ID} == 10 ]]; then
    species=Gthe #Guillardia theta; 11 marks (H3, input)
    chromsize=87266873
elif [[ ${SLURM_ARRAY_TASK_ID} == 11 ]]; then
    species=Cfra #Creolimax fragrantissima; 10 marks (H3, input)
    chromsize=44821703
elif [[ ${SLURM_ARRAY_TASK_ID} == 12 ]]; then
    species=Scer #Saccharomyces cerevisiae; 11 marks (H3, input)
    chromsize=12157105
fi
#2. list of sample names on separate lines (eg copied from excel)
sample_list=/users/asebe/smontgomery/input_files/${species}_rnaseq.txt
#3. set directory containing sample folders
work_folder=/no_backup/asebe/smontgomery/chip/${species}
#4. STAR genome directory
star_dir=/users/asebe/smontgomery/genomes/${species}/STAR_index

## Load the required modules
eval "$(conda shell.bash hook)"
conda activate /users/asebe/cnavarrete/miniconda3/envs/ichip/

##Make relevant folders and move to working folder
mkdir -p ${work_folder}
cd $work_folder
mkdir /no_backup/asebe/smontgomery/chip/${species}/rnaseq/
#set tmp file for deeptools
export TMPDIR=/nfs/scratch01/asebe/smontgomery/deeptool_tmp
# make sure the directory exists
mkdir -p $TMPDIR

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===
##MODIFY WHEN DATA ARRIVES
cat /no_backup/asebe/smontgomery/chip/${species}/fastq/SRR11235958.fastq \
        /no_backup/asebe/smontgomery/chip/${species}/fastq/SRR11235959.fastq \
        /no_backup/asebe/smontgomery/chip/${species}/fastq/SRR11236317.fastq \
        /no_backup/asebe/smontgomery/chip/${species}/fastq/SRR11236318.fastq > /no_backup/asebe/smontgomery/chip/${species}/fastq/rnaseq.fastq
STAR --runThreadN 6 \
    --genomeDir ${star_dir} \
    --readFilesIn /no_backup/asebe/smontgomery/chip/${species}/fastq/rnaseq.fastq \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 2 \
    --limitBAMsortRAM 300000000000 \
    --outFilterMismatchNmax 3 \
    --alignIntronMax 5000 \
    --genomeLoad LoadAndRemove \
    --outSAMunmapped None \
    --readNameSeparator ' ' \
    --quantMode GeneCounts \
    --outFileNamePrefix /no_backup/asebe/smontgomery/chip/${species}/rnaseq/${species}_STAR_

##Generate bigwig file
#5. set bam file of RNA-seq reads
if [[ ${species} == "Nvec" ]]; then
    inputbam=/no_backup/asebe/aelek/RNAseq_nvec/Nvec_adult_polyp_Aligned.sortedByCoord.out.bam
else 
    inputbam=/no_backup/asebe/smontgomery/chip/${species}/rnaseq/${species}_STAR_Aligned.sortedByCoord.out.bam
fi
bamCoverage \
    --bam ${inputbam} \
    --outFileName bw/${species}_geneExp.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --ignoreDuplicates \
    --normalizeUsing RPKM \
    --numberOfProcessors 8

##Clean up files
gzip /no_backup/asebe/smontgomery/chip/${species}/fastq/rnaseq.fastq
rm /no_backup/asebe/smontgomery/chip/${species}/fastq/SRR11235958.fastq \
    /no_backup/asebe/smontgomery/chip/${species}/fastq/SRR11235959.fastq \
    /no_backup/asebe/smontgomery/chip/${species}/fastq/SRR11236317.fastq \
    /no_backup/asebe/smontgomery/chip/${species}/fastq/SRR11236318.fastq 