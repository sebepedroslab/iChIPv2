#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name chromhmm_binarize
#SBATCH --qos=vshort
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=8000
#SBATCH --array=1-12
#SBATCH --output=/no_backup/asebe/smontgomery/logs/%x.%a.out
#SBATCH --error=/no_backup/asebe/smontgomery/logs/%x.%a.err
##################################
# make bash behave more robustly #
##################################
set -e
set -u
set -o pipefail
# === begin ENVIRONMENT SETUP ===
#1. Species; chromsize calculated from scaffold/chromosome sizes in fasta.fai
##Loop through species
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]]; then
    species=Acas #Acanthamoeba castellani; 14 marks (H3,input)
    chromsize=43956874
elif [[ ${SLURM_ARRAY_TASK_ID} == 2 ]]; then
    species=Atha #Arabidopsis thaliana; 11 marks (H3,input)
    chromsize=131559676
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
    chromsize=481355555
elif [[ ${SLURM_ARRAY_TASK_ID} == 7 ]]; then
    species=Spun #Spizellomyces punctatus; 14 marks (H3,input)
    chromsize=24131112
elif [[ ${SLURM_ARRAY_TASK_ID} == 8 ]]; then
    species=Tthe #Tetrahymena thermophila; 13 marks (H3, input)
    chromsize=103349470
elif [[ ${SLURM_ARRAY_TASK_ID} == 9 ]]; then
    species=Bnat #Bigelowiella natans; 14 marks (H3, input)
    chromsize=94701163
elif [[ ${SLURM_ARRAY_TASK_ID} == 10 ]]; then
    species=Gthe #Guillardia theta; 11 marks (H3, input)
    chromsize=87145349
elif [[ ${SLURM_ARRAY_TASK_ID} == 11 ]]; then
    species=Cfra #Creolimax fragrantissima; 10 marks (H3, input)
    chromsize=44821703
elif [[ ${SLURM_ARRAY_TASK_ID} == 12 ]]; then
    species=Scer #Saccharomyces cerevisiae; 11 marks (H3, input)
    chromsize=12071326
fi
#2. set directory containing sample folders
work_folder=/no_backup/asebe/smontgomery/chip/${species}

##Load modules
eval "$(conda shell.bash hook)"
conda activate /users/asebe/cnavarrete/miniconda3/envs/ichip/

##Make relevant folders and move to working folder
mkdir -p $work_folder/
cd $work_folder
export TMPDIR=/nfs/scratch01/asebe/smontgomery/tmp
mkdir -p $TMPDIR
mkdir bambinary

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===

##Binarize bam files
java -mx4000M -jar /users/asebe/smontgomery/ChromHMM/ChromHMM.jar BinarizeBam -b 200 -paired -gzip /users/asebe/smontgomery/ChromHMM/CHROMSIZES/${species}.txt bam/ /users/asebe/smontgomery/input_files/${species}_chromHMM_files.txt bambinary/

