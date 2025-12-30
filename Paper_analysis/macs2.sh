#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name macs2
#SBATCH --qos=shorter
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=8000
#SBATCH --array=1-14
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
#2. list of sample names on separate lines (eg copied from excel)
sample_list=/users/asebe/smontgomery/input_files/${species}_files.txt
#3. set directory containing sample folders
work_folder=/no_backup/asebe/smontgomery/chip/${species}
#4. Data type "SE" or "PE"
datatype="PE"

##File names
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`

## Load the required modules
eval "$(conda shell.bash hook)"
conda activate /users/asebe/cnavarrete/miniconda3/envs/ichip/
module load MACS2/2.1.0.20150731-foss-2018a-Python-2.7.14

##Make relevant folders and move to working folder
mkdir -p ${work_folder}/macs2
cd $work_folder
#set tmp file for deeptools
export TMPDIR=/nfs/scratch01/asebe/smontgomery/deeptool_tmp
# make sure the directory exists
mkdir -p $TMPDIR

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===
#################################################################################################################################################
##Call peaks 
tmpvar=`echo "${NAME}" | cut -d '_' -f2`
if [ ${tmpvar} == H3 ] || [ ${tmpvar} == input ]; then
    ##Call H3 and input peaks to create blacklist 
    macs2 callpeak -t bam/${NAME}.bam --broad --max-gap 300 -f BAMPE -g ${chromsize} -n ${NAME} --outdir macs2/
else
    ##Narrow
    macs2 callpeak -t bam/${NAME}.bam -c bam/${species}_input.bam --max-gap 300 -f BAMPE -g ${chromsize} -n ${NAME} --outdir macs2/
    ##Broad
    macs2 callpeak -t bam/${NAME}.bam -c bam/${species}_input.bam --broad --max-gap 300 -f BAMPE -g ${chromsize} -n ${NAME} --outdir macs2/
fi

