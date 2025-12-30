#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name chromhmm_makeRef
#SBATCH --qos=test
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
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
    species=Ngru #Naegleria gruberi; 11 marks (H3, input)
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
    species=Tthe #Tetrahymena thermophila; 11 marks (H3, input)
    chromsize=103349470
elif [[ ${SLURM_ARRAY_TASK_ID} == 9 ]]; then
    species=Bnat #Bigelowiella natans; 14 marks (H3, input)
    chromsize=94701163
elif [[ ${SLURM_ARRAY_TASK_ID} == 10 ]]; then
    species=Gthe #Guillardia theta; 11 marks (H3, input)
    chromsize=87145349
elif [[ ${SLURM_ARRAY_TASK_ID} == 11 ]]; then
    species=Cfra #Creolimax fragrantissima; 11 marks (H3, input)
    chromsize=44821703
elif [[ ${SLURM_ARRAY_TASK_ID} == 12 ]]; then
    species=Scer #Saccharomyces cerevisiae; 11 marks (H3, input)
    chromsize=12071326
fi
#2. set directory containing sample folders
work_folder=/users/asebe/smontgomery/

##Load modules
eval "$(conda shell.bash hook)"
conda activate /users/asebe/cnavarrete/miniconda3/envs/ichip/

##Make relevant folders and move to working folder
mkdir -p $work_folder/
cd $work_folder
export TMPDIR=/nfs/scratch01/asebe/smontgomery/deeptool_tmp
mkdir -p $TMPDIR

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===

##Make CHROMSIZES file
if [[ -f genomes/${species}/${species}.chromsizes.txt ]]; then
	cp genomes/${species}/${species}.chromsizes.txt ChromHMM/CHROMSIZES/${species}.txt
else
	cut -f1,2 /users/asebe/cnavarrete/proj/iChIP/genome_indexes/${species}_genome_data/${species}_gDNA.fasta.fai > ChromHMM/CHROMSIZES/${species}.txt
	cp ChromHMM/CHROMSIZES/${species}.txt genomes/${species}/${species}.chromsizes.txt 
fi

##Make COORDS files for Gene, Exon, TES, TSS, TEall, DNA TEs, DNA TE fragments, intact DNA TEs, RT TEs, LINEs, RT TE fragments, intact RT TEs, SINEs
mkdir ChromHMM/COORDS/${species}
rm ChromHMM/COORDS/${species}/*.bed.gz
##Gene
cut -f1,2,3 genomes/${species}/${species}.gene.PCG.bed > ChromHMM/COORDS/${species}/Gene.${species}.bed
gzip --force ChromHMM/COORDS/${species}/Gene.${species}.bed
##Exon
cut -f1,4,5 genomes/${species}/${species}.gene.exon.gtf > ChromHMM/COORDS/${species}/Exon.${species}.bed
gzip --force ChromHMM/COORDS/${species}/Exon.${species}.bed
##TES
awk -v OFS='\t' '{if ( $6 == "-" ) print $1,$2,$2+1 ; else print $1,$3-1,$3 }' genomes/${species}/${species}.gene.PCG.bed > ChromHMM/COORDS/${species}/TES.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TES.${species}.bed
##TSS
awk -v OFS='\t' '{if ( $6 == "-" ) print $1,$3-1,$3 ; else print $1,$2,$2+1 }' genomes/${species}/${species}.gene.PCG.bed > ChromHMM/COORDS/${species}/TSS.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TSS.${species}.bed
##TEall
cut -f1,2,3 genomes/${species}/${species}.TE.all.bed > ChromHMM/COORDS/${species}/TEall.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TEall.${species}.bed
##TE DNA
cut -f1,2,3 genomes/${species}/${species}.TE.DNA.bed > ChromHMM/COORDS/${species}/TEDNA.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TEDNA.${species}.bed
##TE DNA fragment
cut -f1,2,3 genomes/${species}/${species}.TE.DNA.fragment.bed > ChromHMM/COORDS/${species}/TEDNAfragment.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TEDNAfragment.${species}.bed
##TE DNA intact
cut -f1,2,3 genomes/${species}/${species}.TE.DNA.intact.bed > ChromHMM/COORDS/${species}/TEDNAintact.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TEDNAintact.${species}.bed
##TE intact
cut -f1,2,3 genomes/${species}/${species}.TE.intact.bed > ChromHMM/COORDS/${species}/TEintact.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TEintact.${species}.bed
##TE RT
cut -f1,2,3 genomes/${species}/${species}.TE.RT.bed > ChromHMM/COORDS/${species}/TERT.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TERT.${species}.bed
##TE RT LINE
cut -f1,2,3 genomes/${species}/${species}.TE.RT.LINE.bed > ChromHMM/COORDS/${species}/TERTLINE.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TERTLINE.${species}.bed
##TE RT LTR fragment
cut -f1,2,3 genomes/${species}/${species}.TE.RT.LTR.fragment.bed > ChromHMM/COORDS/${species}/TERTLTRfragment.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TERTLTRfragment.${species}.bed
##TE RT LTR intact
cut -f1,2,3 genomes/${species}/${species}.TE.RT.LTR.intact.bed > ChromHMM/COORDS/${species}/TERTLTRintact.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TERTLTRintact.${species}.bed
##TE RT SINE
cut -f1,2,3 genomes/${species}/${species}.TE.RT.SINE.bed > ChromHMM/COORDS/${species}/TERTSINE.${species}.bed
gzip --force ChromHMM/COORDS/${species}/TERTSINE.${species}.bed


##Make ANCHORFILES for TES, TSS
mkdir ChromHMM/ANCHORFILES/${species}
##TES
awk -v OFS='\t' '{if ( $6 == "-" ) print $1,$2,$6 ; else print $1,$3-1,$6 }' genomes/${species}/${species}.gene.PCG.bed > ChromHMM/ANCHORFILES/${species}/TES.${species}.txt
gzip --force ChromHMM/ANCHORFILES/${species}/TES.${species}.txt
##TSS
awk -v OFS='\t' '{if ( $6 == "-" ) print $1,$3-1,$6 ; else print $1,$2,$6 }' genomes/${species}/${species}.gene.PCG.bed > ChromHMM/ANCHORFILES/${species}/TSS.${species}.txt
gzip --force ChromHMM/ANCHORFILES/${species}/TSS.${species}.txt
