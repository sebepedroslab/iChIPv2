#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name log2r
#SBATCH --qos=short
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=8000
#SBATCH --array=2-13
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
species=Acas #Acanthamoeba castellani; 14 marks (H3,input)
chromsize=43956874
# species=Atha #Arabidopsis thaliana; 11 marks (H3,input)
# chromsize=131559676
# species=Bnat #Bigelowiela natans; 14 marks (H3, input)
# chromsize=94701163
# species=Cfra #Creolimax fra; 10 marks (H3, input)
# chromsize=44821703
# species=Ddis #Dictyostelium dis; 14 marks (H3, input)
# chromsize=34134454
# species=Gthe #Guiardia theta; 11 marks (H3, input)
# chromsize=87145349
# species=Ngru #Naegleria gruberii; 10 marks (H3, input)
# chromsize=40964085
# species=Nvec #Nematostella vec; 14 marks (H3, input)
# chromsize=269418438
# species=Ppat #Physcomitrium patens; 11 marks (H3, input)
# chromsize=481355555
# species=Scer #Saccharomyces cerevisiae; 11 marks (H3, input)
# chromsize=12071326
# species=Spun #Spizellomyces punctatus; 14 marks (H3,input)
# chromsize=24131112
# species=Tthe #Tetrahymena therm; 13 marks (H3, input)
# chromsize=103349470
#2. list of sample names on separate lines (eg copied from excel)
sample_list=/users/asebe/smontgomery/input_files/${species}_files.txt
#3. set directory containing sample folders
work_folder=/no_backup/asebe/smontgomery/chip/${species}
#4. line in sample file input sample name is on
inputline=1

##File names
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`
INPUT=`sed -n "${inputline} p" $sample_list | awk '{print $1}'`

## Load the required modules
eval "$(conda shell.bash hook)"
conda activate /users/asebe/cnavarrete/miniconda3/envs/ichip/

##Make relevant folders and move to working folder
mkdir -p ${work_folder}
cd $work_folder
mkdir log2r
#set tmp file for deeptools
export TMPDIR=/nfs/scratch01/asebe/smontgomery/deeptool_tmp
# make sure the directory exists
mkdir -p $TMPDIR

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===

# #generate a log2 ratio coverage file with aligned BAM files
bigwigCompare -b1 bw/${NAME}.bw -b2 bw/${species}_input.bw --operation log2 --skipZeroOverZero -of bigwig --binSize=10 -o log2r/${NAME}.log2r.input.bw

