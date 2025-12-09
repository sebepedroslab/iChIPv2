#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name EDTA
#SBATCH --qos=marathon
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100000
#SBATCH --array=1-12
#SBATCH --output=/users/asebe/smontgomery/tmp/logs/%x.%a.out
#SBATCH --error=/users/asebe/smontgomery/tmp/logs/%x.%a.err
##################################
# make bash behave more robustly #
##################################
set -e
set -u
set -o pipefail
##Adapted from Shuangyang and Elin to make gtf files for TE annotations
##Made complicated by Sean
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
#3. set genome and annotation files
fastafile=/users/asebe/xgraubove/genomes/data/${species}_gDNA.fasta
cdsfile=/users/asebe/xgraubove/genomes/data/${species}_long.cds.fasta

eval "$(conda shell.bash hook)"
# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===

##Create renamed fasta file
mkdir -p ${work_folder}/EDTA2
cd ${work_folder}/EDTA2
/users/asebe/xgraubove/miniconda3/bin/bioawk -c fastx '{ print ">prov_"NR"\n"$2 }' ${fastafile} > renamed_genome.${species}.fasta
/users/asebe/xgraubove/miniconda3/bin/bioawk -c fastx '{ print "prov_"NR"\t"$1  }' ${fastafile} > renamed_genome.${species}.dict.txt

##Run EDTA2
conda activate /users/asebe/smontgomery/.conda/envs/EDTA2_env/envs/EDTA2.2/
perl /users/asebe/smontgomery/EDTA2/EDTA.pl --genome renamed_genome.${species}.fasta --overwrite 1 --anno 1 --force 1 --cds ${cdsfile} --threads 32 #--evaluate 1

##Rename scaffolds
for gfi in *.gff3 ; do
	awk 'BEGIN{OFS="\t"} NR==FNR { l[$1]=$2;next} { for(i = 1; i <= 1; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' renamed_genome.${species}.dict.txt ${gfi} > ${gfi%%.gff3}.original_scaffolds.gff
done

cp ${work_folder}/EDTA2/renamed_genome.${species}.fasta.mod.EDTA.TEanno.original_scaffolds.gff ${work_folder}/${species}.TE.all.gff3
