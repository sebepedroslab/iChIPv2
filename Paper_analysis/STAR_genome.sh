#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name star_genome
#SBATCH --qos=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1-12
#SBATCH --output=/users/asebe/smontgomery/tmp/logs/%x.%a.out
#SBATCH --error=/users/asebe/smontgomery/tmp/logs/%x.%a.err
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
#2. set directory containing sample folders
work_folder=/users/asebe/smontgomery/genomes/${species}
#3. GTF files
gtffile=/users/asebe/xgraubove/genomes/data/${species}_long.annot.gtf
#4. Genome file
genomefile=/users/asebe/xgraubove/genomes/data/${species}_gDNA.fasta

## Load the required modules
eval "$(conda shell.bash hook)"
conda activate /users/asebe/cnavarrete/miniconda3/envs/ichip/

##Make relevant folders and move to working folder
mkdir -p ${work_folder}
cd $work_folder
mkdir profiles
#set tmp file for deeptools
export TMPDIR=$work_folder/deeptool_tmp
# make sure the directory exists
mkdir -p $TMPDIR

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===

STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir /users/asebe/smontgomery/genomes/${species}/STAR_index \
    --genomeFastaFiles ${genomefile} \
    --sjdbGTFfile ${gtffile} \
    --genomeSAindexNbases 12
