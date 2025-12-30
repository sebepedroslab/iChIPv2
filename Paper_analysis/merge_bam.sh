#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name mergebam
#SBATCH --qos=short
#SBATCH --time=06:00:00
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
# set -e
# set -u
# set -o pipefail
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
FILE1=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $2}'`
FILE2=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $3}'`
FILE3=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $4}'`
FILE4=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $5}'`
FILE5=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $6}'`
FILE6=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $7}'`

## Load the required modules
eval "$(conda shell.bash hook)"
conda activate /users/asebe/cnavarrete/miniconda3/envs/ichip/
module load SAMtools/1.20-GCC-13.2.0

##Make relevant folders and move to working folder
mkdir -p $work_folder/
cd $work_folder
export TMPDIR=/nfs/scratch01/asebe/smontgomery/tmp
mkdir -p $TMPDIR
mkdir inputbam
mkdir bam
mkdir bw

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===

##Merge bam files for replicates

if [[ -f bam_files/${species}/${FILE6}.q30.rmdup.bam ]]; then
	samtools merge -f bam/${NAME}.bam bam_files/${species}/${FILE1}.q30.rmdup.bam bam_files/${species}/${FILE2}.q30.rmdup.bam bam_files/${species}/${FILE3}.q30.rmdup.bam bam_files/${species}/${FILE4}.q30.rmdup.bam bam_files/${species}/${FILE5}.q30.rmdup.bam bam_files/${species}/${FILE6}.q30.rmdup.bam
elif [[ -f bam_files/${species}/${FILE5}.q30.rmdup.bam ]]; then
	samtools merge -f bam/${NAME}.bam bam_files/${species}/${FILE1}.q30.rmdup.bam bam_files/${species}/${FILE2}.q30.rmdup.bam bam_files/${species}/${FILE3}.q30.rmdup.bam bam_files/${species}/${FILE4}.q30.rmdup.bam bam_files/${species}/${FILE5}.q30.rmdup.bam
elif [[ -f bam_files/${species}/${FILE4}.q30.rmdup.bam ]]; then
	samtools merge -f bam/${NAME}.bam bam_files/${species}/${FILE1}.q30.rmdup.bam bam_files/${species}/${FILE2}.q30.rmdup.bam bam_files/${species}/${FILE3}.q30.rmdup.bam bam_files/${species}/${FILE4}.q30.rmdup.bam
elif [[ -f bam_files/${species}/${FILE3}.q30.rmdup.bam ]]; then
	samtools merge -f bam/${NAME}.bam bam_files/${species}/${FILE1}.q30.rmdup.bam bam_files/${species}/${FILE2}.q30.rmdup.bam bam_files/${species}/${FILE3}.q30.rmdup.bam
elif [[ -f bam_files/${species}/${FILE2}.q30.rmdup.bam ]]; then
	samtools merge -f bam/${NAME}.bam bam_files/${species}/${FILE1}.q30.rmdup.bam bam_files/${species}/${FILE2}.q30.rmdup.bam
elif [[ -f bam_files/${species}/${FILE1}.q30.rmdup.bam ]]; then
	cp bam_files/${species}/${FILE1}.q30.rmdup.bam bam/${NAME}.bam
else
	echo "You're working on thin ice, bud!"
fi
samtools index bam/${NAME}.bam
if [ ${datatype} == PE ]; then
	bamCoverage -b bam/${NAME}.bam -o bw/${NAME}.bw --effectiveGenomeSize ${chromsize} --binSize=10 --extendReads --normalizeUsing RPGC
elif [ ${datatype} == SE ]; then
	bamCoverage -b bam/${NAME}.bam -o bw/${NAME}.bw --effectiveGenomeSize ${chromsize} --binSize=10 --normalizeUsing RPGC
fi