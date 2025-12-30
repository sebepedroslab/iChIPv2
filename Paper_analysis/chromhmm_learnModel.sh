#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name chromhmm_learn
#SBATCH --qos=short
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=8000
#SBATCH --array=2-50
#SBATCH --output=/no_backup/asebe/smontgomery/logs/%x.%a.out
#SBATCH --error=/no_backup/asebe/smontgomery/logs/%x.%a.err
##################################
# make bash behave more robustly #
##################################
set -e
set -u
set -o pipefail
# === begin ENVIRONMENT SETUP ===
#1. Species; chromsize calculated from scaffold/chromosome sizes in fasta.fai; comment/uncomment lines to run script per species
species=Acas #Acanthamoeba castellani; 14 marks (H3,input)
chromsize=43956874
# species=Atha #Arabidopsis thaliana; 11 marks (H3,input)
# chromsize=131559676
# species=Bnat #Bigelowiela natans; 14 marks (H3, input)
# chromsize=94701163
# species=Cfra #Creolimax fra; 11 marks (H3, input)
# chromsize=44821703
# species=Ddis #Dictyostelium dis; 14 marks (H3, input)
# chromsize=34134454
# species=Gthe #Guiardia theta; 11 marks (H3, input)
# chromsize=87145349
# species=Ngru #Naegleria gruberi; 11 marks (H3, input)
# chromsize=40964085
# species=Nvec #Nematostella vec; 14 marks (H3, input)
# chromsize=269418438
# species=Ppat #Physcomitrium patens; 11 marks (H3, input)
# chromsize=481355555
# species=Scer #Saccharomyces cerevisiae; 11 marks (H3, input)
# chromsize=12071326
# species=Spun #Spizellomyces punctatus; 14 marks (H3,input)
# chromsize=24131112
# species=Tthe #Tetrahymena therm; 11 marks (H3, input)
# chromsize=103349470
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

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===

##Run ChromHMM for X states
java -mx1600M -jar /users/asebe/smontgomery/ChromHMM/ChromHMM.jar LearnModel bambinary/ state_${SLURM_ARRAY_TASK_ID}/ ${SLURM_ARRAY_TASK_ID} ${species}
