#!/bin/bash
##################
# slurm settings #
##################
#SBATCH --job-name plot_genome
#SBATCH --qos=short
#SBATCH --time=02:00:00
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
work_folder=/no_backup/asebe/smontgomery/chip/${species}
#3. Gene bed file
genebed=/users/asebe/smontgomery/genomes/${species}/${species}.gene.PCG.bed
#4. Transposon bed files
TEnogene=/users/asebe/smontgomery/genomes/${species}/${species}.TE.nogene.bed

## Load the required modules
eval "$(conda shell.bash hook)"
conda activate /users/asebe/cnavarrete/miniconda3/envs/ichip/

##Make relevant folders and move to working folder
mkdir -p ${work_folder}
cd $work_folder
mkdir profiles
#set tmp file for deeptools
export TMPDIR=/nfs/scratch01/asebe/smontgomery/deeptool_tmp
# make sure the directory exists
mkdir -p $TMPDIR

# === end ENVIRONMENT SETUP ===

# === begin RUN THE PROGRAM ===
####################################################################################################################################################################################################################################
##Spun
if  [[ ${species} == "Spun" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Spun_H3K9ac.log2r.input.bw \
        log2r/Spun_H3K27ac.log2r.input.bw \
        log2r/Spun_H4K16ac.log2r.input.bw \
        log2r/Spun_H3K4me2.log2r.input.bw \
        log2r/Spun_H3K4me3.log2r.input.bw \
        log2r/Spun_H3K36me3.log2r.input.bw \
        log2r/Spun_H3K79me1.log2r.input.bw \
        log2r/Spun_H3K79me2.log2r.input.bw \
        log2r/Spun_H3K79me3.log2r.input.bw \
        log2r/Spun_H3K27me3.log2r.input.bw \
        log2r/Spun_H3K9me1.log2r.input.bw \
        log2r/Spun_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_3/cluster_x/g' | sed 's/cluster_5/cluster_3/g' | sed 's/cluster_6/cluster_5/g' | sed 's/cluster_4/cluster_6/g' | sed 's/cluster_x/cluster_4/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Spun_H3K9ac.log2r.input.bw \
        log2r/Spun_H3K27ac.log2r.input.bw \
        log2r/Spun_H4K16ac.log2r.input.bw \
        log2r/Spun_H3K4me2.log2r.input.bw \
        log2r/Spun_H3K4me3.log2r.input.bw \
        log2r/Spun_H3K36me3.log2r.input.bw \
        log2r/Spun_H3K79me1.log2r.input.bw \
        log2r/Spun_H3K79me2.log2r.input.bw \
        log2r/Spun_H3K79me3.log2r.input.bw \
        log2r/Spun_H3K27me3.log2r.input.bw \
        log2r/Spun_H3K9me1.log2r.input.bw \
        log2r/Spun_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Spun_H3K9ac.log2r.input.bw \
        log2r/Spun_H3K27ac.log2r.input.bw \
        log2r/Spun_H4K16ac.log2r.input.bw \
        log2r/Spun_H3K4me2.log2r.input.bw \
        log2r/Spun_H3K4me3.log2r.input.bw \
        log2r/Spun_H3K36me3.log2r.input.bw \
        log2r/Spun_H3K79me1.log2r.input.bw \
        log2r/Spun_H3K79me2.log2r.input.bw \
        log2r/Spun_H3K79me3.log2r.input.bw \
        log2r/Spun_H3K27me3.log2r.input.bw \
        log2r/Spun_H3K9me1.log2r.input.bw \
        log2r/Spun_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_6' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Spun_H3K9ac.log2r.input.bw \
        log2r/Spun_H3K27ac.log2r.input.bw \
        log2r/Spun_H4K16ac.log2r.input.bw \
        log2r/Spun_H3K4me2.log2r.input.bw \
        log2r/Spun_H3K4me3.log2r.input.bw \
        log2r/Spun_H3K36me3.log2r.input.bw \
        log2r/Spun_H3K79me1.log2r.input.bw \
        log2r/Spun_H3K79me2.log2r.input.bw \
        log2r/Spun_H3K79me3.log2r.input.bw \
        log2r/Spun_H3K27me3.log2r.input.bw \
        log2r/Spun_H3K9me1.log2r.input.bw \
        log2r/Spun_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Spun_H3K9ac.log2r.input.bw \
        log2r/Spun_H3K27ac.log2r.input.bw \
        log2r/Spun_H4K16ac.log2r.input.bw \
        log2r/Spun_H3K4me2.log2r.input.bw \
        log2r/Spun_H3K4me3.log2r.input.bw \
        log2r/Spun_H3K36me3.log2r.input.bw \
        log2r/Spun_H3K79me1.log2r.input.bw \
        log2r/Spun_H3K79me2.log2r.input.bw \
        log2r/Spun_H3K79me3.log2r.input.bw \
        log2r/Spun_H3K27me3.log2r.input.bw \
        log2r/Spun_H3K9me1.log2r.input.bw \
        log2r/Spun_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Spun_H3K9ac.log2r.input.bw \
        log2r/Spun_H3K27ac.log2r.input.bw \
        log2r/Spun_H4K16ac.log2r.input.bw \
        log2r/Spun_H3K4me2.log2r.input.bw \
        log2r/Spun_H3K4me3.log2r.input.bw \
        log2r/Spun_H3K36me3.log2r.input.bw \
        log2r/Spun_H3K79me1.log2r.input.bw \
        log2r/Spun_H3K79me2.log2r.input.bw \
        log2r/Spun_H3K79me3.log2r.input.bw \
        log2r/Spun_H3K27me3.log2r.input.bw \
        log2r/Spun_H3K9me1.log2r.input.bw \
        log2r/Spun_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt
    computeMatrix reference-point \
            --referencePoint TSS \
            -b 5000 -a 5000 \
            -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt  \
            -S \
            log2r/Spun_H3K9ac.log2r.input.bw \
            log2r/Spun_H3K27ac.log2r.input.bw \
            log2r/Spun_H4K16ac.log2r.input.bw \
            log2r/Spun_H3K4me2.log2r.input.bw \
            log2r/Spun_H3K4me3.log2r.input.bw \
            log2r/Spun_H3K36me3.log2r.input.bw \
            log2r/Spun_H3K79me1.log2r.input.bw \
            log2r/Spun_H3K79me2.log2r.input.bw \
            log2r/Spun_H3K79me3.log2r.input.bw \
            log2r/Spun_H3K27me3.log2r.input.bw \
            log2r/Spun_H3K9me1.log2r.input.bw \
            log2r/Spun_H3K9me3.log2r.input.bw \
            --skipZeros --missingDataAsZero \
            -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.high.gz 
        plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.high.gz    \
            --perGroup \
            --plotHeight 18 --plotWidth 13 \
            --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
            --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
            --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.high.pdf

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Spun_H3K9ac.log2r.input.bw \
        log2r/Spun_H3K27ac.log2r.input.bw \
        log2r/Spun_H4K16ac.log2r.input.bw \
        log2r/Spun_H3K4me2.log2r.input.bw \
        log2r/Spun_H3K4me3.log2r.input.bw \
        log2r/Spun_H3K36me3.log2r.input.bw \
        log2r/Spun_H3K79me1.log2r.input.bw \
        log2r/Spun_H3K79me2.log2r.input.bw \
        log2r/Spun_H3K79me3.log2r.input.bw \
        log2r/Spun_H3K27me3.log2r.input.bw \
        log2r/Spun_H3K9me1.log2r.input.bw \
        log2r/Spun_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --kmeans 3 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k3.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k3.pdf 
    head -1 profiles/log2r.${species}.TEnogene.region.heatmap.k3.txt > profiles/log2r.${species}.TEnogene.region.heatmap.k3.reorder.txt
    grep -v '#' profiles/log2r.${species}.TEnogene.region.heatmap.k3.txt | sed 's/cluster_2/cluster_x/g' | sed 's/cluster_3/cluster_2/g' | sed 's/cluster_x/cluster_3/g' | sort -s -k13 >> profiles/log2r.${species}.TEnogene.region.heatmap.k3.reorder.txt
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k3.reorder.txt \
        -S \
        log2r/Spun_H3K9ac.log2r.input.bw \
        log2r/Spun_H3K27ac.log2r.input.bw \
        log2r/Spun_H4K16ac.log2r.input.bw \
        log2r/Spun_H3K4me2.log2r.input.bw \
        log2r/Spun_H3K4me3.log2r.input.bw \
        log2r/Spun_H3K36me3.log2r.input.bw \
        log2r/Spun_H3K79me1.log2r.input.bw \
        log2r/Spun_H3K79me2.log2r.input.bw \
        log2r/Spun_H3K79me3.log2r.input.bw \
        log2r/Spun_H3K27me3.log2r.input.bw \
        log2r/Spun_H3K9me1.log2r.input.bw \
        log2r/Spun_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k3.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k3.pdf 


#############################################################################################################################################################################################################################################################################
##Cfra
elif  [[ ${species} == "Cfra" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Cfra_H3K9ac.log2r.input.bw \
        log2r/Cfra_H3K27ac.log2r.input.bw \
        log2r/Cfra_H4K16ac.log2r.input.bw \
        log2r/Cfra_H3K4me2.log2r.input.bw \
        log2r/Cfra_H3K4me3.log2r.input.bw \
        log2r/Cfra_H3K36me3.log2r.input.bw \
        log2r/Cfra_H3K79me1.log2r.input.bw \
        log2r/Cfra_H3K79me2.log2r.input.bw \
        log2r/Cfra_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_3/cluster_1/g' | sed 's/cluster_2/cluster_3/g' | sed 's/cluster_x/cluster_2/g' | sed 's/cluster_4/cluster_x/g' | sed 's/cluster_5/cluster_4/g' | sed 's/cluster_x/cluster_5/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Cfra_H3K9ac.log2r.input.bw \
        log2r/Cfra_H3K27ac.log2r.input.bw \
        log2r/Cfra_H4K16ac.log2r.input.bw \
        log2r/Cfra_H3K4me2.log2r.input.bw \
        log2r/Cfra_H3K4me3.log2r.input.bw \
        log2r/Cfra_H3K36me3.log2r.input.bw \
        log2r/Cfra_H3K79me1.log2r.input.bw \
        log2r/Cfra_H3K79me2.log2r.input.bw \
        log2r/Cfra_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Cfra_H3K9ac.log2r.input.bw \
        log2r/Cfra_H3K27ac.log2r.input.bw \
        log2r/Cfra_H4K16ac.log2r.input.bw \
        log2r/Cfra_H3K4me2.log2r.input.bw \
        log2r/Cfra_H3K4me3.log2r.input.bw \
        log2r/Cfra_H3K36me3.log2r.input.bw \
        log2r/Cfra_H3K79me1.log2r.input.bw \
        log2r/Cfra_H3K79me2.log2r.input.bw \
        log2r/Cfra_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_6' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Cfra_H3K9ac.log2r.input.bw \
        log2r/Cfra_H3K27ac.log2r.input.bw \
        log2r/Cfra_H4K16ac.log2r.input.bw \
        log2r/Cfra_H3K4me2.log2r.input.bw \
        log2r/Cfra_H3K4me3.log2r.input.bw \
        log2r/Cfra_H3K36me3.log2r.input.bw \
        log2r/Cfra_H3K79me1.log2r.input.bw \
        log2r/Cfra_H3K79me2.log2r.input.bw \
        log2r/Cfra_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Cfra_H3K9ac.log2r.input.bw \
        log2r/Cfra_H3K27ac.log2r.input.bw \
        log2r/Cfra_H4K16ac.log2r.input.bw \
        log2r/Cfra_H3K4me2.log2r.input.bw \
        log2r/Cfra_H3K4me3.log2r.input.bw \
        log2r/Cfra_H3K36me3.log2r.input.bw \
        log2r/Cfra_H3K79me1.log2r.input.bw \
        log2r/Cfra_H3K79me2.log2r.input.bw \
        log2r/Cfra_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Cfra_H3K9ac.log2r.input.bw \
        log2r/Cfra_H3K27ac.log2r.input.bw \
        log2r/Cfra_H4K16ac.log2r.input.bw \
        log2r/Cfra_H3K4me2.log2r.input.bw \
        log2r/Cfra_H3K4me3.log2r.input.bw \
        log2r/Cfra_H3K36me3.log2r.input.bw \
        log2r/Cfra_H3K79me1.log2r.input.bw \
        log2r/Cfra_H3K79me2.log2r.input.bw \
        log2r/Cfra_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt  \
        -S \
        log2r/Cfra_H3K9ac.log2r.input.bw \
        log2r/Cfra_H3K27ac.log2r.input.bw \
        log2r/Cfra_H4K16ac.log2r.input.bw \
        log2r/Cfra_H3K4me2.log2r.input.bw \
        log2r/Cfra_H3K4me3.log2r.input.bw \
        log2r/Cfra_H3K36me3.log2r.input.bw \
        log2r/Cfra_H3K79me1.log2r.input.bw \
        log2r/Cfra_H3K79me2.log2r.input.bw \
        log2r/Cfra_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.high.pdf

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Cfra_H3K9ac.log2r.input.bw \
        log2r/Cfra_H3K27ac.log2r.input.bw \
        log2r/Cfra_H4K16ac.log2r.input.bw \
        log2r/Cfra_H3K4me2.log2r.input.bw \
        log2r/Cfra_H3K4me3.log2r.input.bw \
        log2r/Cfra_H3K36me3.log2r.input.bw \
        log2r/Cfra_H3K79me1.log2r.input.bw \
        log2r/Cfra_H3K79me2.log2r.input.bw \
        log2r/Cfra_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 \
        --kmeans 5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k5.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k5.pdf 
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k5.txt\
        -S \
        log2r/Cfra_H3K9ac.log2r.input.bw \
        log2r/Cfra_H3K27ac.log2r.input.bw \
        log2r/Cfra_H4K16ac.log2r.input.bw \
        log2r/Cfra_H3K4me2.log2r.input.bw \
        log2r/Cfra_H3K4me3.log2r.input.bw \
        log2r/Cfra_H3K36me3.log2r.input.bw \
        log2r/Cfra_H3K79me1.log2r.input.bw \
        log2r/Cfra_H3K79me2.log2r.input.bw \
        log2r/Cfra_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k5.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k5.pdf 


#################################################################################################################################################
##Acas
elif  [[ ${species} == "Acas" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Acas_H3K9ac.log2r.input.bw \
        log2r/Acas_H3K27ac.log2r.input.bw \
        log2r/Acas_H4K16ac.log2r.input.bw \
        log2r/Acas_H3K4me2.log2r.input.bw \
        log2r/Acas_H3K4me3.log2r.input.bw \
        log2r/Acas_H3K36me3.log2r.input.bw \
        log2r/Acas_H3K79me1.log2r.input.bw \
        log2r/Acas_H3K79me2.log2r.input.bw \
        log2r/Acas_H3K79me3.log2r.input.bw \
        log2r/Acas_H3K27me3.log2r.input.bw \
        log2r/Acas_H3K9me1.log2r.input.bw \
        log2r/Acas_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_4/cluster_x/g' | sed 's/cluster_6/cluster_4/g' | sed 's/cluster_x/cluster_6/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Acas_H3K9ac.log2r.input.bw \
        log2r/Acas_H3K27ac.log2r.input.bw \
        log2r/Acas_H4K16ac.log2r.input.bw \
        log2r/Acas_H3K4me2.log2r.input.bw \
        log2r/Acas_H3K4me3.log2r.input.bw \
        log2r/Acas_H3K36me3.log2r.input.bw \
        log2r/Acas_H3K79me1.log2r.input.bw \
        log2r/Acas_H3K79me2.log2r.input.bw \
        log2r/Acas_H3K79me3.log2r.input.bw \
        log2r/Acas_H3K27me3.log2r.input.bw \
        log2r/Acas_H3K9me1.log2r.input.bw \
        log2r/Acas_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_3' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Acas_H3K9ac.log2r.input.bw \
        log2r/Acas_H3K27ac.log2r.input.bw \
        log2r/Acas_H4K16ac.log2r.input.bw \
        log2r/Acas_H3K4me2.log2r.input.bw \
        log2r/Acas_H3K4me3.log2r.input.bw \
        log2r/Acas_H3K36me3.log2r.input.bw \
        log2r/Acas_H3K79me1.log2r.input.bw \
        log2r/Acas_H3K79me2.log2r.input.bw \
        log2r/Acas_H3K79me3.log2r.input.bw \
        log2r/Acas_H3K27me3.log2r.input.bw \
        log2r/Acas_H3K9me1.log2r.input.bw \
        log2r/Acas_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_6' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Acas_H3K9ac.log2r.input.bw \
        log2r/Acas_H3K27ac.log2r.input.bw \
        log2r/Acas_H4K16ac.log2r.input.bw \
        log2r/Acas_H3K4me2.log2r.input.bw \
        log2r/Acas_H3K4me3.log2r.input.bw \
        log2r/Acas_H3K36me3.log2r.input.bw \
        log2r/Acas_H3K79me1.log2r.input.bw \
        log2r/Acas_H3K79me2.log2r.input.bw \
        log2r/Acas_H3K79me3.log2r.input.bw \
        log2r/Acas_H3K27me3.log2r.input.bw \
        log2r/Acas_H3K9me1.log2r.input.bw \
        log2r/Acas_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Acas_H3K9ac.log2r.input.bw \
        log2r/Acas_H3K27ac.log2r.input.bw \
        log2r/Acas_H4K16ac.log2r.input.bw \
        log2r/Acas_H3K4me2.log2r.input.bw \
        log2r/Acas_H3K4me3.log2r.input.bw \
        log2r/Acas_H3K36me3.log2r.input.bw \
        log2r/Acas_H3K79me1.log2r.input.bw \
        log2r/Acas_H3K79me2.log2r.input.bw \
        log2r/Acas_H3K79me3.log2r.input.bw \
        log2r/Acas_H3K27me3.log2r.input.bw \
        log2r/Acas_H3K9me1.log2r.input.bw \
        log2r/Acas_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Acas_H3K9ac.log2r.input.bw \
        log2r/Acas_H3K27ac.log2r.input.bw \
        log2r/Acas_H4K16ac.log2r.input.bw \
        log2r/Acas_H3K4me2.log2r.input.bw \
        log2r/Acas_H3K4me3.log2r.input.bw \
        log2r/Acas_H3K36me3.log2r.input.bw \
        log2r/Acas_H3K79me1.log2r.input.bw \
        log2r/Acas_H3K79me2.log2r.input.bw \
        log2r/Acas_H3K79me3.log2r.input.bw \
        log2r/Acas_H3K27me3.log2r.input.bw \
        log2r/Acas_H3K9me1.log2r.input.bw \
        log2r/Acas_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_3' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Acas_H3K9ac.log2r.input.bw \
        log2r/Acas_H3K27ac.log2r.input.bw \
        log2r/Acas_H4K16ac.log2r.input.bw \
        log2r/Acas_H3K4me2.log2r.input.bw \
        log2r/Acas_H3K4me3.log2r.input.bw \
        log2r/Acas_H3K36me3.log2r.input.bw \
        log2r/Acas_H3K79me1.log2r.input.bw \
        log2r/Acas_H3K79me2.log2r.input.bw \
        log2r/Acas_H3K79me3.log2r.input.bw \
        log2r/Acas_H3K27me3.log2r.input.bw \
        log2r/Acas_H3K9me1.log2r.input.bw \
        log2r/Acas_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt > profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt
    grep -v '#' profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_6/cluster_1/g' | sed 's/cluster_5/cluster_6/g' | sed 's/cluster_3/cluster_5/g' | sed 's/cluster_x/cluster_3/g' | sed 's/cluster_2/cluster_x/g' | sed 's/cluster_4/cluster_2/g' | sed 's/cluster_x/cluster_4/g' | sort -s -k13 >> profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt\
        -S \
        log2r/Acas_H3K9ac.log2r.input.bw \
        log2r/Acas_H3K27ac.log2r.input.bw \
        log2r/Acas_H4K16ac.log2r.input.bw \
        log2r/Acas_H3K4me2.log2r.input.bw \
        log2r/Acas_H3K4me3.log2r.input.bw \
        log2r/Acas_H3K36me3.log2r.input.bw \
        log2r/Acas_H3K79me1.log2r.input.bw \
        log2r/Acas_H3K79me2.log2r.input.bw \
        log2r/Acas_H3K79me3.log2r.input.bw \
        log2r/Acas_H3K27me3.log2r.input.bw \
        log2r/Acas_H3K9me1.log2r.input.bw \
        log2r/Acas_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k6.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k6.pdf 


#################################################################################################################################################
##Atha
elif  [[ ${species} == "Atha" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Atha_H3K9ac.log2r.input.bw \
        log2r/Atha_H3K27ac.log2r.input.bw \
        log2r/Atha_H4K16ac.log2r.input.bw \
        log2r/Atha_H3K4me2.log2r.input.bw \
        log2r/Atha_H3K4me3.log2r.input.bw \
        log2r/Atha_H3K36me3.log2r.input.bw \
        log2r/Atha_H3K27me3.log2r.input.bw \
        log2r/Atha_H3K9me1.log2r.input.bw \
        log2r/Atha_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_4/cluster_x/g' | sed 's/cluster_6/cluster_4/g' | sed 's/cluster_5/cluster_6/g' | sed 's/cluster_x/cluster_5/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Atha_H3K9ac.log2r.input.bw \
        log2r/Atha_H3K27ac.log2r.input.bw \
        log2r/Atha_H4K16ac.log2r.input.bw \
        log2r/Atha_H3K4me2.log2r.input.bw \
        log2r/Atha_H3K4me3.log2r.input.bw \
        log2r/Atha_H3K36me3.log2r.input.bw \
        log2r/Atha_H3K27me3.log2r.input.bw \
        log2r/Atha_H3K9me1.log2r.input.bw \
        log2r/Atha_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_2' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Atha_H3K9ac.log2r.input.bw \
        log2r/Atha_H3K27ac.log2r.input.bw \
        log2r/Atha_H4K16ac.log2r.input.bw \
        log2r/Atha_H3K4me2.log2r.input.bw \
        log2r/Atha_H3K4me3.log2r.input.bw \
        log2r/Atha_H3K36me3.log2r.input.bw \
        log2r/Atha_H3K27me3.log2r.input.bw \
        log2r/Atha_H3K9me1.log2r.input.bw \
        log2r/Atha_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_5' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Atha_H3K9ac.log2r.input.bw \
        log2r/Atha_H3K27ac.log2r.input.bw \
        log2r/Atha_H4K16ac.log2r.input.bw \
        log2r/Atha_H3K4me2.log2r.input.bw \
        log2r/Atha_H3K4me3.log2r.input.bw \
        log2r/Atha_H3K36me3.log2r.input.bw \
        log2r/Atha_H3K27me3.log2r.input.bw \
        log2r/Atha_H3K9me1.log2r.input.bw \
        log2r/Atha_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Atha_H3K9ac.log2r.input.bw \
        log2r/Atha_H3K27ac.log2r.input.bw \
        log2r/Atha_H4K16ac.log2r.input.bw \
        log2r/Atha_H3K4me2.log2r.input.bw \
        log2r/Atha_H3K4me3.log2r.input.bw \
        log2r/Atha_H3K36me3.log2r.input.bw \
        log2r/Atha_H3K27me3.log2r.input.bw \
        log2r/Atha_H3K9me1.log2r.input.bw \
        log2r/Atha_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Atha_H3K9ac.log2r.input.bw \
        log2r/Atha_H3K27ac.log2r.input.bw \
        log2r/Atha_H4K16ac.log2r.input.bw \
        log2r/Atha_H3K4me2.log2r.input.bw \
        log2r/Atha_H3K4me3.log2r.input.bw \
        log2r/Atha_H3K36me3.log2r.input.bw \
        log2r/Atha_H3K27me3.log2r.input.bw \
        log2r/Atha_H3K9me1.log2r.input.bw \
        log2r/Atha_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_2' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Atha_H3K9ac.log2r.input.bw \
        log2r/Atha_H3K27ac.log2r.input.bw \
        log2r/Atha_H4K16ac.log2r.input.bw \
        log2r/Atha_H3K4me2.log2r.input.bw \
        log2r/Atha_H3K4me3.log2r.input.bw \
        log2r/Atha_H3K36me3.log2r.input.bw \
        log2r/Atha_H3K27me3.log2r.input.bw \
        log2r/Atha_H3K9me1.log2r.input.bw \
        log2r/Atha_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 7 8 9 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt > profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt  | sed 's/cluster_3/cluster_x/g' | sed 's/cluster_4/cluster_3/g' | sed 's/cluster_5/cluster_4/g' | sed 's/cluster_6/cluster_5/g' | sed 's/cluster_x/cluster_6/g' | sort -s -k13 >> profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt 
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt \
        -S \
        log2r/Atha_H3K9ac.log2r.input.bw \
        log2r/Atha_H3K27ac.log2r.input.bw \
        log2r/Atha_H4K16ac.log2r.input.bw \
        log2r/Atha_H3K4me2.log2r.input.bw \
        log2r/Atha_H3K4me3.log2r.input.bw \
        log2r/Atha_H3K36me3.log2r.input.bw \
        log2r/Atha_H3K27me3.log2r.input.bw \
        log2r/Atha_H3K9me1.log2r.input.bw \
        log2r/Atha_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k6.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 7 8 9 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k6.pdf 



################################################################################################################################################
##Scer
elif  [[ ${species} == "Scer" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Scer_H3K9ac.log2r.input.bw \
        log2r/Scer_H3K27ac.log2r.input.bw \
        log2r/Scer_H4K16ac.log2r.input.bw \
        log2r/Scer_H3K4me2.log2r.input.bw \
        log2r/Scer_H3K4me3.log2r.input.bw \
        log2r/Scer_H3K36me3.log2r.input.bw \
        log2r/Scer_H3K79me1.log2r.input.bw \
        log2r/Scer_H3K79me2.log2r.input.bw \
        log2r/Scer_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_3/cluster_1/g' | sed 's/cluster_2/cluster_3/g' | sed 's/cluster_4/cluster_2/g' | sed 's/cluster_x/cluster_4/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Scer_H3K9ac.log2r.input.bw \
        log2r/Scer_H3K27ac.log2r.input.bw \
        log2r/Scer_H4K16ac.log2r.input.bw \
        log2r/Scer_H3K4me2.log2r.input.bw \
        log2r/Scer_H3K4me3.log2r.input.bw \
        log2r/Scer_H3K36me3.log2r.input.bw \
        log2r/Scer_H3K79me1.log2r.input.bw \
        log2r/Scer_H3K79me2.log2r.input.bw \
        log2r/Scer_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_2' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Scer_H3K9ac.log2r.input.bw \
        log2r/Scer_H3K27ac.log2r.input.bw \
        log2r/Scer_H4K16ac.log2r.input.bw \
        log2r/Scer_H3K4me2.log2r.input.bw \
        log2r/Scer_H3K4me3.log2r.input.bw \
        log2r/Scer_H3K36me3.log2r.input.bw \
        log2r/Scer_H3K79me1.log2r.input.bw \
        log2r/Scer_H3K79me2.log2r.input.bw \
        log2r/Scer_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_6' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Scer_H3K9ac.log2r.input.bw \
        log2r/Scer_H3K27ac.log2r.input.bw \
        log2r/Scer_H4K16ac.log2r.input.bw \
        log2r/Scer_H3K4me2.log2r.input.bw \
        log2r/Scer_H3K4me3.log2r.input.bw \
        log2r/Scer_H3K36me3.log2r.input.bw \
        log2r/Scer_H3K79me1.log2r.input.bw \
        log2r/Scer_H3K79me2.log2r.input.bw \
        log2r/Scer_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Scer_H3K9ac.log2r.input.bw \
        log2r/Scer_H3K27ac.log2r.input.bw \
        log2r/Scer_H4K16ac.log2r.input.bw \
        log2r/Scer_H3K4me2.log2r.input.bw \
        log2r/Scer_H3K4me3.log2r.input.bw \
        log2r/Scer_H3K36me3.log2r.input.bw \
        log2r/Scer_H3K79me1.log2r.input.bw \
        log2r/Scer_H3K79me2.log2r.input.bw \
        log2r/Scer_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Scer_H3K9ac.log2r.input.bw \
        log2r/Scer_H3K27ac.log2r.input.bw \
        log2r/Scer_H4K16ac.log2r.input.bw \
        log2r/Scer_H3K4me2.log2r.input.bw \
        log2r/Scer_H3K4me3.log2r.input.bw \
        log2r/Scer_H3K36me3.log2r.input.bw \
        log2r/Scer_H3K79me1.log2r.input.bw \
        log2r/Scer_H3K79me2.log2r.input.bw \
        log2r/Scer_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_2' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt

 
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Scer_H3K9ac.log2r.input.bw \
        log2r/Scer_H3K27ac.log2r.input.bw \
        log2r/Scer_H4K16ac.log2r.input.bw \
        log2r/Scer_H3K4me2.log2r.input.bw \
        log2r/Scer_H3K4me3.log2r.input.bw \
        log2r/Scer_H3K36me3.log2r.input.bw \
        log2r/Scer_H3K79me1.log2r.input.bw \
        log2r/Scer_H3K79me2.log2r.input.bw \
        log2r/Scer_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 7 8 9 \
        --kmeans 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k2.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k2.pdf 
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k2.txt\
        -S \
        log2r/Scer_H3K9ac.log2r.input.bw \
        log2r/Scer_H3K27ac.log2r.input.bw \
        log2r/Scer_H4K16ac.log2r.input.bw \
        log2r/Scer_H3K4me2.log2r.input.bw \
        log2r/Scer_H3K4me3.log2r.input.bw \
        log2r/Scer_H3K36me3.log2r.input.bw \
        log2r/Scer_H3K79me1.log2r.input.bw \
        log2r/Scer_H3K79me2.log2r.input.bw \
        log2r/Scer_H3K79me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k2.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 7 8 9 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k2.pdf 


#################################################################################################################################################
##Nvec
elif  [[ ${species} == "Nvec" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Nvec_H3K9ac.log2r.input.bw \
        log2r/Nvec_H3K27ac.log2r.input.bw \
        log2r/Nvec_H4K16ac.log2r.input.bw \
        log2r/Nvec_H3K4me2.log2r.input.bw \
        log2r/Nvec_H3K4me3.log2r.input.bw \
        log2r/Nvec_H3K36me3.log2r.input.bw \
        log2r/Nvec_H3K79me1.log2r.input.bw \
        log2r/Nvec_H3K79me2.log2r.input.bw \
        log2r/Nvec_H3K79me3.log2r.input.bw \
        log2r/Nvec_H3K27me3.log2r.input.bw \
        log2r/Nvec_H3K9me1.log2r.input.bw \
        log2r/Nvec_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_5/cluster_1/g' | sed 's/cluster_6/cluster_5/g' | sed 's/cluster_4/cluster_6/g' | sed 's/cluster_x/cluster_4/g' | sed 's/cluster_2/cluster_x/g' | sed 's/cluster_3/cluster_2/g' | sed 's/cluster_x/cluster_3/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Nvec_H3K9ac.log2r.input.bw \
        log2r/Nvec_H3K27ac.log2r.input.bw \
        log2r/Nvec_H4K16ac.log2r.input.bw \
        log2r/Nvec_H3K4me2.log2r.input.bw \
        log2r/Nvec_H3K4me3.log2r.input.bw \
        log2r/Nvec_H3K36me3.log2r.input.bw \
        log2r/Nvec_H3K79me1.log2r.input.bw \
        log2r/Nvec_H3K79me2.log2r.input.bw \
        log2r/Nvec_H3K79me3.log2r.input.bw \
        log2r/Nvec_H3K27me3.log2r.input.bw \
        log2r/Nvec_H3K9me1.log2r.input.bw \
        log2r/Nvec_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Nvec_H3K9ac.log2r.input.bw \
        log2r/Nvec_H3K27ac.log2r.input.bw \
        log2r/Nvec_H4K16ac.log2r.input.bw \
        log2r/Nvec_H3K4me2.log2r.input.bw \
        log2r/Nvec_H3K4me3.log2r.input.bw \
        log2r/Nvec_H3K36me3.log2r.input.bw \
        log2r/Nvec_H3K79me1.log2r.input.bw \
        log2r/Nvec_H3K79me2.log2r.input.bw \
        log2r/Nvec_H3K79me3.log2r.input.bw \
        log2r/Nvec_H3K27me3.log2r.input.bw \
        log2r/Nvec_H3K9me1.log2r.input.bw \
        log2r/Nvec_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_5' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Nvec_H3K9ac.log2r.input.bw \
        log2r/Nvec_H3K27ac.log2r.input.bw \
        log2r/Nvec_H4K16ac.log2r.input.bw \
        log2r/Nvec_H3K4me2.log2r.input.bw \
        log2r/Nvec_H3K4me3.log2r.input.bw \
        log2r/Nvec_H3K36me3.log2r.input.bw \
        log2r/Nvec_H3K79me1.log2r.input.bw \
        log2r/Nvec_H3K79me2.log2r.input.bw \
        log2r/Nvec_H3K79me3.log2r.input.bw \
        log2r/Nvec_H3K27me3.log2r.input.bw \
        log2r/Nvec_H3K9me1.log2r.input.bw \
        log2r/Nvec_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Nvec_H3K9ac.log2r.input.bw \
        log2r/Nvec_H3K27ac.log2r.input.bw \
        log2r/Nvec_H4K16ac.log2r.input.bw \
        log2r/Nvec_H3K4me2.log2r.input.bw \
        log2r/Nvec_H3K4me3.log2r.input.bw \
        log2r/Nvec_H3K36me3.log2r.input.bw \
        log2r/Nvec_H3K79me1.log2r.input.bw \
        log2r/Nvec_H3K79me2.log2r.input.bw \
        log2r/Nvec_H3K79me3.log2r.input.bw \
        log2r/Nvec_H3K27me3.log2r.input.bw \
        log2r/Nvec_H3K9me1.log2r.input.bw \
        log2r/Nvec_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Nvec_H3K9ac.log2r.input.bw \
        log2r/Nvec_H3K27ac.log2r.input.bw \
        log2r/Nvec_H4K16ac.log2r.input.bw \
        log2r/Nvec_H3K4me2.log2r.input.bw \
        log2r/Nvec_H3K4me3.log2r.input.bw \
        log2r/Nvec_H3K36me3.log2r.input.bw \
        log2r/Nvec_H3K79me1.log2r.input.bw \
        log2r/Nvec_H3K79me2.log2r.input.bw \
        log2r/Nvec_H3K79me3.log2r.input.bw \
        log2r/Nvec_H3K27me3.log2r.input.bw \
        log2r/Nvec_H3K9me1.log2r.input.bw \
        log2r/Nvec_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt  \
        -S \
        log2r/Nvec_H3K9ac.log2r.input.bw \
        log2r/Nvec_H3K27ac.log2r.input.bw \
        log2r/Nvec_H4K16ac.log2r.input.bw \
        log2r/Nvec_H3K4me2.log2r.input.bw \
        log2r/Nvec_H3K4me3.log2r.input.bw \
        log2r/Nvec_H3K36me3.log2r.input.bw \
        log2r/Nvec_H3K79me1.log2r.input.bw \
        log2r/Nvec_H3K79me2.log2r.input.bw \
        log2r/Nvec_H3K79me3.log2r.input.bw \
        log2r/Nvec_H3K27me3.log2r.input.bw \
        log2r/Nvec_H3K9me1.log2r.input.bw \
        log2r/Nvec_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.high.pdf

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Nvec_H3K9ac.log2r.input.bw \
        log2r/Nvec_H3K27ac.log2r.input.bw \
        log2r/Nvec_H4K16ac.log2r.input.bw \
        log2r/Nvec_H3K4me2.log2r.input.bw \
        log2r/Nvec_H3K4me3.log2r.input.bw \
        log2r/Nvec_H3K36me3.log2r.input.bw \
        log2r/Nvec_H3K79me1.log2r.input.bw \
        log2r/Nvec_H3K79me2.log2r.input.bw \
        log2r/Nvec_H3K79me3.log2r.input.bw \
        log2r/Nvec_H3K27me3.log2r.input.bw \
        log2r/Nvec_H3K9me1.log2r.input.bw \
        log2r/Nvec_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt  > profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt  | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_5/cluster_1/g' | sed 's/cluster_6/cluster_5/g' | sed 's/cluster_x/cluster_6/g' | sort -s -k13 >> profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt 
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt \
        -S \
        log2r/Nvec_H3K9ac.log2r.input.bw \
        log2r/Nvec_H3K27ac.log2r.input.bw \
        log2r/Nvec_H4K16ac.log2r.input.bw \
        log2r/Nvec_H3K4me2.log2r.input.bw \
        log2r/Nvec_H3K4me3.log2r.input.bw \
        log2r/Nvec_H3K36me3.log2r.input.bw \
        log2r/Nvec_H3K79me1.log2r.input.bw \
        log2r/Nvec_H3K79me2.log2r.input.bw \
        log2r/Nvec_H3K79me3.log2r.input.bw \
        log2r/Nvec_H3K27me3.log2r.input.bw \
        log2r/Nvec_H3K9me1.log2r.input.bw \
        log2r/Nvec_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k6.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k6.pdf 

#################################################################################################################################################
##Ppat
elif  [[ ${species} == "Ppat" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Ppat_H3K9ac.log2r.input.bw \
        log2r/Ppat_H3K27ac.log2r.input.bw \
        log2r/Ppat_H4K16ac.log2r.input.bw \
        log2r/Ppat_H3K4me2.log2r.input.bw \
        log2r/Ppat_H3K4me3.log2r.input.bw \
        log2r/Ppat_H3K36me3.log2r.input.bw \
        log2r/Ppat_H3K27me3.log2r.input.bw \
        log2r/Ppat_H3K9me1.log2r.input.bw \
        log2r/Ppat_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_2/cluster_1/g' | sed 's/cluster_4/cluster_2/g' | sed 's/cluster_3/cluster_4/g' | sed 's/cluster_x/cluster_3/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Ppat_H3K9ac.log2r.input.bw \
        log2r/Ppat_H3K27ac.log2r.input.bw \
        log2r/Ppat_H4K16ac.log2r.input.bw \
        log2r/Ppat_H3K4me2.log2r.input.bw \
        log2r/Ppat_H3K4me3.log2r.input.bw \
        log2r/Ppat_H3K36me3.log2r.input.bw \
        log2r/Ppat_H3K27me3.log2r.input.bw \
        log2r/Ppat_H3K9me1.log2r.input.bw \
        log2r/Ppat_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Ppat_H3K9ac.log2r.input.bw \
        log2r/Ppat_H3K27ac.log2r.input.bw \
        log2r/Ppat_H4K16ac.log2r.input.bw \
        log2r/Ppat_H3K4me2.log2r.input.bw \
        log2r/Ppat_H3K4me3.log2r.input.bw \
        log2r/Ppat_H3K36me3.log2r.input.bw \
        log2r/Ppat_H3K27me3.log2r.input.bw \
        log2r/Ppat_H3K9me1.log2r.input.bw \
        log2r/Ppat_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_6' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Ppat_H3K9ac.log2r.input.bw \
        log2r/Ppat_H3K27ac.log2r.input.bw \
        log2r/Ppat_H4K16ac.log2r.input.bw \
        log2r/Ppat_H3K4me2.log2r.input.bw \
        log2r/Ppat_H3K4me3.log2r.input.bw \
        log2r/Ppat_H3K36me3.log2r.input.bw \
        log2r/Ppat_H3K27me3.log2r.input.bw \
        log2r/Ppat_H3K9me1.log2r.input.bw \
        log2r/Ppat_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Ppat_H3K9ac.log2r.input.bw \
        log2r/Ppat_H3K27ac.log2r.input.bw \
        log2r/Ppat_H4K16ac.log2r.input.bw \
        log2r/Ppat_H3K4me2.log2r.input.bw \
        log2r/Ppat_H3K4me3.log2r.input.bw \
        log2r/Ppat_H3K36me3.log2r.input.bw \
        log2r/Ppat_H3K27me3.log2r.input.bw \
        log2r/Ppat_H3K9me1.log2r.input.bw \
        log2r/Ppat_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Ppat_H3K9ac.log2r.input.bw \
        log2r/Ppat_H3K27ac.log2r.input.bw \
        log2r/Ppat_H4K16ac.log2r.input.bw \
        log2r/Ppat_H3K4me2.log2r.input.bw \
        log2r/Ppat_H3K4me3.log2r.input.bw \
        log2r/Ppat_H3K36me3.log2r.input.bw \
        log2r/Ppat_H3K27me3.log2r.input.bw \
        log2r/Ppat_H3K9me1.log2r.input.bw \
        log2r/Ppat_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Ppat_H3K9ac.log2r.input.bw \
        log2r/Ppat_H3K27ac.log2r.input.bw \
        log2r/Ppat_H4K16ac.log2r.input.bw \
        log2r/Ppat_H3K4me2.log2r.input.bw \
        log2r/Ppat_H3K4me3.log2r.input.bw \
        log2r/Ppat_H3K36me3.log2r.input.bw \
        log2r/Ppat_H3K27me3.log2r.input.bw \
        log2r/Ppat_H3K9me1.log2r.input.bw \
        log2r/Ppat_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 7 8 9 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt > profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt
    grep -v '#' profiles/log2r.${species}.TEnogene.region.heatmap.k6.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_3/cluster_1/g' | sed 's/cluster_2/cluster_3/g' | sed 's/cluster_x/cluster_2/g' | sort -s -k13 >> profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k6.reorder.txt \
        -S \
        log2r/Ppat_H3K9ac.log2r.input.bw \
        log2r/Ppat_H3K27ac.log2r.input.bw \
        log2r/Ppat_H4K16ac.log2r.input.bw \
        log2r/Ppat_H3K4me2.log2r.input.bw \
        log2r/Ppat_H3K4me3.log2r.input.bw \
        log2r/Ppat_H3K36me3.log2r.input.bw \
        log2r/Ppat_H3K27me3.log2r.input.bw \
        log2r/Ppat_H3K9me1.log2r.input.bw \
        log2r/Ppat_H3K9me2.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k6.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 7 8 9 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me2 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k6.pdf 

################################################################################################################################################
##Tthe
elif  [[ ${species} == "Tthe" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Tthe_H3K9ac.log2r.input.bw \
        log2r/Tthe_H3K27ac.log2r.input.bw \
        log2r/Tthe_H4K16ac.log2r.input.bw \
        log2r/Tthe_H3K4me2.log2r.input.bw \
        log2r/Tthe_H3K4me3.log2r.input.bw \
        log2r/Tthe_H3K36me3.log2r.input.bw \
        log2r/Tthe_H3K79me1.log2r.input.bw \
        log2r/Tthe_H3K27me3.log2r.input.bw \
        log2r/Tthe_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Tthe_H3K9ac.log2r.input.bw \
        log2r/Tthe_H3K27ac.log2r.input.bw \
        log2r/Tthe_H4K16ac.log2r.input.bw \
        log2r/Tthe_H3K4me2.log2r.input.bw \
        log2r/Tthe_H3K4me3.log2r.input.bw \
        log2r/Tthe_H3K36me3.log2r.input.bw \
        log2r/Tthe_H3K79me1.log2r.input.bw \
        log2r/Tthe_H3K27me3.log2r.input.bw \
        log2r/Tthe_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_3' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Tthe_H3K9ac.log2r.input.bw \
        log2r/Tthe_H3K27ac.log2r.input.bw \
        log2r/Tthe_H4K16ac.log2r.input.bw \
        log2r/Tthe_H3K4me2.log2r.input.bw \
        log2r/Tthe_H3K4me3.log2r.input.bw \
        log2r/Tthe_H3K36me3.log2r.input.bw \
        log2r/Tthe_H3K79me1.log2r.input.bw \
        log2r/Tthe_H3K27me3.log2r.input.bw \
        log2r/Tthe_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_5' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Tthe_H3K9ac.log2r.input.bw \
        log2r/Tthe_H3K27ac.log2r.input.bw \
        log2r/Tthe_H4K16ac.log2r.input.bw \
        log2r/Tthe_H3K4me2.log2r.input.bw \
        log2r/Tthe_H3K4me3.log2r.input.bw \
        log2r/Tthe_H3K36me3.log2r.input.bw \
        log2r/Tthe_H3K79me1.log2r.input.bw \
        log2r/Tthe_H3K27me3.log2r.input.bw \
        log2r/Tthe_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Tthe_H3K9ac.log2r.input.bw \
        log2r/Tthe_H3K27ac.log2r.input.bw \
        log2r/Tthe_H4K16ac.log2r.input.bw \
        log2r/Tthe_H3K4me2.log2r.input.bw \
        log2r/Tthe_H3K4me3.log2r.input.bw \
        log2r/Tthe_H3K36me3.log2r.input.bw \
        log2r/Tthe_H3K79me1.log2r.input.bw \
        log2r/Tthe_H3K27me3.log2r.input.bw \
        log2r/Tthe_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Tthe_H3K9ac.log2r.input.bw \
        log2r/Tthe_H3K27ac.log2r.input.bw \
        log2r/Tthe_H4K16ac.log2r.input.bw \
        log2r/Tthe_H3K4me2.log2r.input.bw \
        log2r/Tthe_H3K4me3.log2r.input.bw \
        log2r/Tthe_H3K36me3.log2r.input.bw \
        log2r/Tthe_H3K79me1.log2r.input.bw \
        log2r/Tthe_H3K27me3.log2r.input.bw \
        log2r/Tthe_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_3' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Tthe_H3K9ac.log2r.input.bw \
        log2r/Tthe_H3K27ac.log2r.input.bw \
        log2r/Tthe_H4K16ac.log2r.input.bw \
        log2r/Tthe_H3K4me2.log2r.input.bw \
        log2r/Tthe_H3K4me3.log2r.input.bw \
        log2r/Tthe_H3K36me3.log2r.input.bw \
        log2r/Tthe_H3K79me1.log2r.input.bw \
        log2r/Tthe_H3K27me3.log2r.input.bw \
        log2r/Tthe_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 8 9 \
        --kmeans 1 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k1.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k1.pdf 
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k1.txt\
        -S \
        log2r/Tthe_H3K9ac.log2r.input.bw \
        log2r/Tthe_H3K27ac.log2r.input.bw \
        log2r/Tthe_H4K16ac.log2r.input.bw \
        log2r/Tthe_H3K4me2.log2r.input.bw \
        log2r/Tthe_H3K4me3.log2r.input.bw \
        log2r/Tthe_H3K36me3.log2r.input.bw \
        log2r/Tthe_H3K79me1.log2r.input.bw \
        log2r/Tthe_H3K27me3.log2r.input.bw \
        log2r/Tthe_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k1.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 8 9 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me3 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k1.pdf 

#################################################################################################################################################
##Ddis
elif  [[ ${species} == "Ddis" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Ddis_H3K9ac.log2r.input.bw \
        log2r/Ddis_H3K27ac.log2r.input.bw \
        log2r/Ddis_H4K16ac.log2r.input.bw \
        log2r/Ddis_H3K4me2.log2r.input.bw \
        log2r/Ddis_H3K4me3.log2r.input.bw \
        log2r/Ddis_H3K36me3.log2r.input.bw \
        log2r/Ddis_H3K79me1.log2r.input.bw \
        log2r/Ddis_H3K79me2.log2r.input.bw \
        log2r/Ddis_H3K79me3.log2r.input.bw \
        log2r/Ddis_H3K27me3.log2r.input.bw \
        log2r/Ddis_H3K9me1.log2r.input.bw \
        log2r/Ddis_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_5/cluster_1/g' | sed 's/cluster_4/cluster_5/g' | sed 's/cluster_6/cluster_4/g' | sed 's/cluster_2/cluster_6/g' | sed 's/cluster_3/cluster_2/g' | sed 's/cluster_x/cluster_3/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Ddis_H3K9ac.log2r.input.bw \
        log2r/Ddis_H3K27ac.log2r.input.bw \
        log2r/Ddis_H4K16ac.log2r.input.bw \
        log2r/Ddis_H3K4me2.log2r.input.bw \
        log2r/Ddis_H3K4me3.log2r.input.bw \
        log2r/Ddis_H3K36me3.log2r.input.bw \
        log2r/Ddis_H3K79me1.log2r.input.bw \
        log2r/Ddis_H3K79me2.log2r.input.bw \
        log2r/Ddis_H3K79me3.log2r.input.bw \
        log2r/Ddis_H3K27me3.log2r.input.bw \
        log2r/Ddis_H3K9me1.log2r.input.bw \
        log2r/Ddis_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Ddis_H3K9ac.log2r.input.bw \
        log2r/Ddis_H3K27ac.log2r.input.bw \
        log2r/Ddis_H4K16ac.log2r.input.bw \
        log2r/Ddis_H3K4me2.log2r.input.bw \
        log2r/Ddis_H3K4me3.log2r.input.bw \
        log2r/Ddis_H3K36me3.log2r.input.bw \
        log2r/Ddis_H3K79me1.log2r.input.bw \
        log2r/Ddis_H3K79me2.log2r.input.bw \
        log2r/Ddis_H3K79me3.log2r.input.bw \
        log2r/Ddis_H3K27me3.log2r.input.bw \
        log2r/Ddis_H3K9me1.log2r.input.bw \
        log2r/Ddis_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_6' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Ddis_H3K9ac.log2r.input.bw \
        log2r/Ddis_H3K27ac.log2r.input.bw \
        log2r/Ddis_H4K16ac.log2r.input.bw \
        log2r/Ddis_H3K4me2.log2r.input.bw \
        log2r/Ddis_H3K4me3.log2r.input.bw \
        log2r/Ddis_H3K36me3.log2r.input.bw \
        log2r/Ddis_H3K79me1.log2r.input.bw \
        log2r/Ddis_H3K79me2.log2r.input.bw \
        log2r/Ddis_H3K79me3.log2r.input.bw \
        log2r/Ddis_H3K27me3.log2r.input.bw \
        log2r/Ddis_H3K9me1.log2r.input.bw \
        log2r/Ddis_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Ddis_H3K9ac.log2r.input.bw \
        log2r/Ddis_H3K27ac.log2r.input.bw \
        log2r/Ddis_H4K16ac.log2r.input.bw \
        log2r/Ddis_H3K4me2.log2r.input.bw \
        log2r/Ddis_H3K4me3.log2r.input.bw \
        log2r/Ddis_H3K36me3.log2r.input.bw \
        log2r/Ddis_H3K79me1.log2r.input.bw \
        log2r/Ddis_H3K79me2.log2r.input.bw \
        log2r/Ddis_H3K79me3.log2r.input.bw \
        log2r/Ddis_H3K27me3.log2r.input.bw \
        log2r/Ddis_H3K9me1.log2r.input.bw \
        log2r/Ddis_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Ddis_H3K9ac.log2r.input.bw \
        log2r/Ddis_H3K27ac.log2r.input.bw \
        log2r/Ddis_H4K16ac.log2r.input.bw \
        log2r/Ddis_H3K4me2.log2r.input.bw \
        log2r/Ddis_H3K4me3.log2r.input.bw \
        log2r/Ddis_H3K36me3.log2r.input.bw \
        log2r/Ddis_H3K79me1.log2r.input.bw \
        log2r/Ddis_H3K79me2.log2r.input.bw \
        log2r/Ddis_H3K79me3.log2r.input.bw \
        log2r/Ddis_H3K27me3.log2r.input.bw \
        log2r/Ddis_H3K9me1.log2r.input.bw \
        log2r/Ddis_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_1' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Ddis_H3K9ac.log2r.input.bw \
        log2r/Ddis_H3K27ac.log2r.input.bw \
        log2r/Ddis_H4K16ac.log2r.input.bw \
        log2r/Ddis_H3K4me2.log2r.input.bw \
        log2r/Ddis_H3K4me3.log2r.input.bw \
        log2r/Ddis_H3K36me3.log2r.input.bw \
        log2r/Ddis_H3K79me1.log2r.input.bw \
        log2r/Ddis_H3K79me2.log2r.input.bw \
        log2r/Ddis_H3K79me3.log2r.input.bw \
        log2r/Ddis_H3K27me3.log2r.input.bw \
        log2r/Ddis_H3K9me1.log2r.input.bw \
        log2r/Ddis_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --kmeans 5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k5.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k5.pdf 
    head -1 profiles/log2r.${species}.TEnogene.region.heatmap.k5.txt  > profiles/log2r.${species}.TEnogene.region.heatmap.k5.reorder.txt 
    grep -v '#' profiles/log2r.${species}.TEnogene.region.heatmap.k5.txt  | sed 's/cluster_2/cluster_x/g' | sed 's/cluster_3/cluster_2/g' | sed 's/cluster_x/cluster_3/g' | sort -s -k13 >> profiles/log2r.${species}.TEnogene.region.heatmap.k5.reorder.txt 
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k5.reorder.txt \
        -S \
        log2r/Ddis_H3K9ac.log2r.input.bw \
        log2r/Ddis_H3K27ac.log2r.input.bw \
        log2r/Ddis_H4K16ac.log2r.input.bw \
        log2r/Ddis_H3K4me2.log2r.input.bw \
        log2r/Ddis_H3K4me3.log2r.input.bw \
        log2r/Ddis_H3K36me3.log2r.input.bw \
        log2r/Ddis_H3K79me1.log2r.input.bw \
        log2r/Ddis_H3K79me2.log2r.input.bw \
        log2r/Ddis_H3K79me3.log2r.input.bw \
        log2r/Ddis_H3K27me3.log2r.input.bw \
        log2r/Ddis_H3K9me1.log2r.input.bw \
        log2r/Ddis_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k5.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k5.pdf 

#################################################################################################################################################
##Ngru
elif  [[ ${species} == "Ngru" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Ngru_H3K9ac.log2r.input.bw \
        log2r/Ngru_H3K27ac.log2r.input.bw \
        log2r/Ngru_H4K16ac.log2r.input.bw \
        log2r/Ngru_H3K4me2.log2r.input.bw \
        log2r/Ngru_H3K4me3.log2r.input.bw \
        log2r/Ngru_H3K36me3.log2r.input.bw \
        log2r/Ngru_H3K27me3.log2r.input.bw \
        log2r/Ngru_H3K9me1.log2r.input.bw \
        log2r/Ngru_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_3/cluster_x/g' | sed 's/cluster_5/cluster_3/g' | sed 's/cluster_6/cluster_5/g' | sed 's/cluster_x/cluster_6/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Ngru_H3K9ac.log2r.input.bw \
        log2r/Ngru_H3K27ac.log2r.input.bw \
        log2r/Ngru_H4K16ac.log2r.input.bw \
        log2r/Ngru_H3K4me2.log2r.input.bw \
        log2r/Ngru_H3K4me3.log2r.input.bw \
        log2r/Ngru_H3K36me3.log2r.input.bw \
        log2r/Ngru_H3K27me3.log2r.input.bw \
        log2r/Ngru_H3K9me1.log2r.input.bw \
        log2r/Ngru_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_2' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Ngru_H3K9ac.log2r.input.bw \
        log2r/Ngru_H3K27ac.log2r.input.bw \
        log2r/Ngru_H4K16ac.log2r.input.bw \
        log2r/Ngru_H3K4me2.log2r.input.bw \
        log2r/Ngru_H3K4me3.log2r.input.bw \
        log2r/Ngru_H3K36me3.log2r.input.bw \
        log2r/Ngru_H3K27me3.log2r.input.bw \
        log2r/Ngru_H3K9me1.log2r.input.bw \
        log2r/Ngru_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_6' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Ngru_H3K9ac.log2r.input.bw \
        log2r/Ngru_H3K27ac.log2r.input.bw \
        log2r/Ngru_H4K16ac.log2r.input.bw \
        log2r/Ngru_H3K4me2.log2r.input.bw \
        log2r/Ngru_H3K4me3.log2r.input.bw \
        log2r/Ngru_H3K36me3.log2r.input.bw \
        log2r/Ngru_H3K27me3.log2r.input.bw \
        log2r/Ngru_H3K9me1.log2r.input.bw \
        log2r/Ngru_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Ngru_H3K9ac.log2r.input.bw \
        log2r/Ngru_H3K27ac.log2r.input.bw \
        log2r/Ngru_H4K16ac.log2r.input.bw \
        log2r/Ngru_H3K4me2.log2r.input.bw \
        log2r/Ngru_H3K4me3.log2r.input.bw \
        log2r/Ngru_H3K36me3.log2r.input.bw \
        log2r/Ngru_H3K27me3.log2r.input.bw \
        log2r/Ngru_H3K9me1.log2r.input.bw \
        log2r/Ngru_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Ngru_H3K9ac.log2r.input.bw \
        log2r/Ngru_H3K27ac.log2r.input.bw \
        log2r/Ngru_H4K16ac.log2r.input.bw \
        log2r/Ngru_H3K4me2.log2r.input.bw \
        log2r/Ngru_H3K4me3.log2r.input.bw \
        log2r/Ngru_H3K36me3.log2r.input.bw \
        log2r/Ngru_H3K27me3.log2r.input.bw \
        log2r/Ngru_H3K9me1.log2r.input.bw \
        log2r/Ngru_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_2' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Ngru_H3K9ac.log2r.input.bw \
        log2r/Ngru_H3K27ac.log2r.input.bw \
        log2r/Ngru_H4K16ac.log2r.input.bw \
        log2r/Ngru_H3K4me2.log2r.input.bw \
        log2r/Ngru_H3K4me3.log2r.input.bw \
        log2r/Ngru_H3K36me3.log2r.input.bw \
        log2r/Ngru_H3K27me3.log2r.input.bw \
        log2r/Ngru_H3K9me1.log2r.input.bw \
        log2r/Ngru_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 7 8 \
        --kmeans 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k2.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k2.pdf 
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k2.txt\
        -S \
        log2r/Ngru_H3K9ac.log2r.input.bw \
        log2r/Ngru_H3K27ac.log2r.input.bw \
        log2r/Ngru_H4K16ac.log2r.input.bw \
        log2r/Ngru_H3K4me2.log2r.input.bw \
        log2r/Ngru_H3K4me3.log2r.input.bw \
        log2r/Ngru_H3K36me3.log2r.input.bw \
        log2r/Ngru_H3K27me3.log2r.input.bw \
        log2r/Ngru_H3K9me1.log2r.input.bw \
        log2r/Ngru_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k2.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 7 8 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k2.pdf 

#################################################################################################################################################
##Bnat
elif  [[ ${species} == "Bnat" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Bnat_H3K9ac.log2r.input.bw \
        log2r/Bnat_H3K27ac.log2r.input.bw \
        log2r/Bnat_H4K16ac.log2r.input.bw \
        log2r/Bnat_H3K4me2.log2r.input.bw \
        log2r/Bnat_H3K4me3.log2r.input.bw \
        log2r/Bnat_H3K36me3.log2r.input.bw \
        log2r/Bnat_H3K79me1.log2r.input.bw \
        log2r/Bnat_H3K79me2.log2r.input.bw \
        log2r/Bnat_H3K79me3.log2r.input.bw \
        log2r/Bnat_H3K27me3.log2r.input.bw \
        log2r/Bnat_H3K9me1.log2r.input.bw \
        log2r/Bnat_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 6 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k6.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k6.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k6.txt > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_2/cluster_1/g' | sed 's/cluster_4/cluster_2/g' | sed 's/cluster_6/cluster_4/g' | sed 's/cluster_x/cluster_6/g' | sed 's/cluster_3/cluster_x/g' | sed 's/cluster_5/cluster_3/g' | sed 's/cluster_x/cluster_5/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt  \
        -S \
        log2r/Bnat_H3K9ac.log2r.input.bw \
        log2r/Bnat_H3K27ac.log2r.input.bw \
        log2r/Bnat_H4K16ac.log2r.input.bw \
        log2r/Bnat_H3K4me2.log2r.input.bw \
        log2r/Bnat_H3K4me3.log2r.input.bw \
        log2r/Bnat_H3K36me3.log2r.input.bw \
        log2r/Bnat_H3K79me1.log2r.input.bw \
        log2r/Bnat_H3K79me2.log2r.input.bw \
        log2r/Bnat_H3K79me3.log2r.input.bw \
        log2r/Bnat_H3K27me3.log2r.input.bw \
        log2r/Bnat_H3K9me1.log2r.input.bw \
        log2r/Bnat_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_2' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt  \
        -S \
        log2r/Bnat_H3K9ac.log2r.input.bw \
        log2r/Bnat_H3K27ac.log2r.input.bw \
        log2r/Bnat_H4K16ac.log2r.input.bw \
        log2r/Bnat_H3K4me2.log2r.input.bw \
        log2r/Bnat_H3K4me3.log2r.input.bw \
        log2r/Bnat_H3K36me3.log2r.input.bw \
        log2r/Bnat_H3K79me1.log2r.input.bw \
        log2r/Bnat_H3K79me2.log2r.input.bw \
        log2r/Bnat_H3K79me3.log2r.input.bw \
        log2r/Bnat_H3K27me3.log2r.input.bw \
        log2r/Bnat_H3K9me1.log2r.input.bw \
        log2r/Bnat_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.txt | egrep '#|cluster_6' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt  \
        -S \
        log2r/Bnat_H3K9ac.log2r.input.bw \
        log2r/Bnat_H3K27ac.log2r.input.bw \
        log2r/Bnat_H4K16ac.log2r.input.bw \
        log2r/Bnat_H3K4me2.log2r.input.bw \
        log2r/Bnat_H3K4me3.log2r.input.bw \
        log2r/Bnat_H3K36me3.log2r.input.bw \
        log2r/Bnat_H3K79me1.log2r.input.bw \
        log2r/Bnat_H3K79me2.log2r.input.bw \
        log2r/Bnat_H3K79me3.log2r.input.bw \
        log2r/Bnat_H3K27me3.log2r.input.bw \
        log2r/Bnat_H3K9me1.log2r.input.bw \
        log2r/Bnat_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.highlow.txt  \
        -S \
        log2r/Bnat_H3K9ac.log2r.input.bw \
        log2r/Bnat_H3K27ac.log2r.input.bw \
        log2r/Bnat_H4K16ac.log2r.input.bw \
        log2r/Bnat_H3K4me2.log2r.input.bw \
        log2r/Bnat_H3K4me3.log2r.input.bw \
        log2r/Bnat_H3K36me3.log2r.input.bw \
        log2r/Bnat_H3K79me1.log2r.input.bw \
        log2r/Bnat_H3K79me2.log2r.input.bw \
        log2r/Bnat_H3K79me3.log2r.input.bw \
        log2r/Bnat_H3K27me3.log2r.input.bw \
        log2r/Bnat_H3K9me1.log2r.input.bw \
        log2r/Bnat_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt  \
        -S \
        log2r/Bnat_H3K9ac.log2r.input.bw \
        log2r/Bnat_H3K27ac.log2r.input.bw \
        log2r/Bnat_H4K16ac.log2r.input.bw \
        log2r/Bnat_H3K4me2.log2r.input.bw \
        log2r/Bnat_H3K4me3.log2r.input.bw \
        log2r/Bnat_H3K36me3.log2r.input.bw \
        log2r/Bnat_H3K79me1.log2r.input.bw \
        log2r/Bnat_H3K79me2.log2r.input.bw \
        log2r/Bnat_H3K79me3.log2r.input.bw \
        log2r/Bnat_H3K27me3.log2r.input.bw \
        log2r/Bnat_H3K9me1.log2r.input.bw \
        log2r/Bnat_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k6.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.txt | egrep '#|cluster_2' > profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k6.reorder.no5kb.high.txt  \
        -S \
        log2r/Bnat_H3K9ac.log2r.input.bw \
        log2r/Bnat_H3K27ac.log2r.input.bw \
        log2r/Bnat_H4K16ac.log2r.input.bw \
        log2r/Bnat_H3K4me2.log2r.input.bw \
        log2r/Bnat_H3K4me3.log2r.input.bw \
        log2r/Bnat_H3K36me3.log2r.input.bw \
        log2r/Bnat_H3K79me1.log2r.input.bw \
        log2r/Bnat_H3K79me2.log2r.input.bw \
        log2r/Bnat_H3K79me3.log2r.input.bw \
        log2r/Bnat_H3K27me3.log2r.input.bw \
        log2r/Bnat_H3K9me1.log2r.input.bw \
        log2r/Bnat_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k6.reorder.no5kb.high.pdf

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Bnat_H3K9ac.log2r.input.bw \
        log2r/Bnat_H3K27ac.log2r.input.bw \
        log2r/Bnat_H4K16ac.log2r.input.bw \
        log2r/Bnat_H3K4me2.log2r.input.bw \
        log2r/Bnat_H3K4me3.log2r.input.bw \
        log2r/Bnat_H3K36me3.log2r.input.bw \
        log2r/Bnat_H3K79me1.log2r.input.bw \
        log2r/Bnat_H3K79me2.log2r.input.bw \
        log2r/Bnat_H3K79me3.log2r.input.bw \
        log2r/Bnat_H3K27me3.log2r.input.bw \
        log2r/Bnat_H3K9me1.log2r.input.bw \
        log2r/Bnat_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --kmeans 5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k5.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k5.pdf 
    head -1 profiles/log2r.${species}.TEnogene.region.heatmap.k5.txt > profiles/log2r.${species}.TEnogene.region.heatmap.k5.reorder.txt
    grep -v '#' profiles/log2r.${species}.TEnogene.region.heatmap.k5.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_3/cluster_1/g' | sed 's/cluster_4/cluster_3/g' | sed 's/cluster_2/cluster_4/g' | sed 's/cluster_x/cluster_2/g' | sort -s -k13 >> profiles/log2r.${species}.TEnogene.region.heatmap.k5.reorder.txt
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k5.reorder.txt \
        -S \
        log2r/Bnat_H3K9ac.log2r.input.bw \
        log2r/Bnat_H3K27ac.log2r.input.bw \
        log2r/Bnat_H4K16ac.log2r.input.bw \
        log2r/Bnat_H3K4me2.log2r.input.bw \
        log2r/Bnat_H3K4me3.log2r.input.bw \
        log2r/Bnat_H3K36me3.log2r.input.bw \
        log2r/Bnat_H3K79me1.log2r.input.bw \
        log2r/Bnat_H3K79me2.log2r.input.bw \
        log2r/Bnat_H3K79me3.log2r.input.bw \
        log2r/Bnat_H3K27me3.log2r.input.bw \
        log2r/Bnat_H3K9me1.log2r.input.bw \
        log2r/Bnat_H3K9me3.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#D448BA" "#70015D" "#4120A9" "#41B6C4" "#3262AB" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k5.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 10 12 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K79me2 H3K79me3 H3K27me3 H3K9me1 H3K9me3 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k5.pdf 


#################################################################################################################################################
##Gthe
elif  [[ ${species} == "Gthe" ]]; then
    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${genebed} \
        -S \
        log2r/Gthe_H3K9ac.log2r.input.bw \
        log2r/Gthe_H3K27ac.log2r.input.bw \
        log2r/Gthe_H4K16ac.log2r.input.bw \
        log2r/Gthe_H3K4me2.log2r.input.bw \
        log2r/Gthe_H3K4me3.log2r.input.bw \
        log2r/Gthe_H3K36me3.log2r.input.bw \
        log2r/Gthe_H3K79me1.log2r.input.bw \
        log2r/Gthe_H3K27me3.log2r.input.bw \
        log2r/Gthe_H3K9me1.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.gz 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 1 2 3 4 5 6 \
        --kmeans 8 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        --outFileSortedRegions profiles/log2r.${species}.gene.noregion.heatmap.k8.txt \
        -out profiles/log2r.${species}.gene.noregion.heatmap.k8.pdf 
    head -1 profiles/log2r.${species}.gene.noregion.heatmap.k8.txt > profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.txt 
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k8.txt | sed 's/cluster_3/cluster_x/g' | sed 's/cluster_7/cluster_3/g' | sed 's/cluster_8/cluster_7/g' | sed 's/cluster_6/cluster_8/g' | sed 's/cluster_x/cluster_6/g' | sed 's/cluster_4/cluster_x/g' | sed 's/cluster_5/cluster_4/g' | sed 's/cluster_x/cluster_5/g' | sort -s -k13 >> profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.txt 
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.txt  \
        -S \
        log2r/Gthe_H3K9ac.log2r.input.bw \
        log2r/Gthe_H3K27ac.log2r.input.bw \
        log2r/Gthe_H4K16ac.log2r.input.bw \
        log2r/Gthe_H3K4me2.log2r.input.bw \
        log2r/Gthe_H3K4me3.log2r.input.bw \
        log2r/Gthe_H3K36me3.log2r.input.bw \
        log2r/Gthe_H3K79me1.log2r.input.bw \
        log2r/Gthe_H3K27me3.log2r.input.bw \
        log2r/Gthe_H3K9me1.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#41B6C4" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k8.reorder.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -2 --zMax 2 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k8.reorder.sortgene.pdf 
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k8.reorder.sortgene.1.5.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.txt | egrep '#|cluster_4' > profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.high.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.high.txt  \
        -S \
        log2r/Gthe_H3K9ac.log2r.input.bw \
        log2r/Gthe_H3K27ac.log2r.input.bw \
        log2r/Gthe_H4K16ac.log2r.input.bw \
        log2r/Gthe_H3K4me2.log2r.input.bw \
        log2r/Gthe_H3K4me3.log2r.input.bw \
        log2r/Gthe_H3K36me3.log2r.input.bw \
        log2r/Gthe_H3K79me1.log2r.input.bw \
        log2r/Gthe_H3K27me3.log2r.input.bw \
        log2r/Gthe_H3K9me1.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.high.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#41B6C4" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k8.reorder.high.pdf
    cat profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.txt | egrep '#|cluster_8' > profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.low.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.low.txt  \
        -S \
        log2r/Gthe_H3K9ac.log2r.input.bw \
        log2r/Gthe_H3K27ac.log2r.input.bw \
        log2r/Gthe_H4K16ac.log2r.input.bw \
        log2r/Gthe_H3K4me2.log2r.input.bw \
        log2r/Gthe_H3K4me3.log2r.input.bw \
        log2r/Gthe_H3K36me3.log2r.input.bw \
        log2r/Gthe_H3K79me1.log2r.input.bw \
        log2r/Gthe_H3K27me3.log2r.input.bw \
        log2r/Gthe_H3K9me1.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.low.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#41B6C4" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k8.reorder.low.pdf
    grep -v '#' profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.low.txt | cat profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.high.txt - > profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.highlow.txt
    computeMatrix reference-point \
        --referencePoint TSS \
        -b 500 -a 3000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.highlow.txt  \
        -S \
        log2r/Gthe_H3K9ac.log2r.input.bw \
        log2r/Gthe_H3K27ac.log2r.input.bw \
        log2r/Gthe_H4K16ac.log2r.input.bw \
        log2r/Gthe_H3K4me2.log2r.input.bw \
        log2r/Gthe_H3K4me3.log2r.input.bw \
        log2r/Gthe_H3K36me3.log2r.input.bw \
        log2r/Gthe_H3K79me1.log2r.input.bw \
        log2r/Gthe_H3K27me3.log2r.input.bw \
        log2r/Gthe_H3K9me1.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.highlow.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#41B6C4" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k8.reorder.highlow.pdf

    computeMatrix reference-point \
        --referencePoint TSS \
        -b 5000 -a 5000 \
        -R profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.no5kb.txt  \
        -S \
        log2r/Gthe_H3K9ac.log2r.input.bw \
        log2r/Gthe_H3K27ac.log2r.input.bw \
        log2r/Gthe_H4K16ac.log2r.input.bw \
        log2r/Gthe_H3K4me2.log2r.input.bw \
        log2r/Gthe_H3K4me3.log2r.input.bw \
        log2r/Gthe_H3K36me3.log2r.input.bw \
        log2r/Gthe_H3K79me1.log2r.input.bw \
        log2r/Gthe_H3K27me3.log2r.input.bw \
        log2r/Gthe_H3K9me1.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz 
    plotProfile -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#41B6C4" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.plotProfile.k8.reorder.no5kb.pdf
    plotHeatmap -m profiles/log2r.${species}.gene.noregion.TSS.reorder.no5kb.gz   \
        --sortUsing region_length \
        --linesAtTickMarks \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        -out profiles/log2r.${species}.gene.noregion.TSS.heatmap.k8.reorder.sortgene.no5kb.pdf 
    cat profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.no5kb.txt | egrep '#|cluster_4' > profiles/log2r.${species}.gene.noregion.heatmap.k8.reorder.no5kb.high.txt

    computeMatrix scale-regions \
        -b 0 -a 0 \
        -R ${TEnogene} \
        -S \
        log2r/Gthe_H3K9ac.log2r.input.bw \
        log2r/Gthe_H3K27ac.log2r.input.bw \
        log2r/Gthe_H4K16ac.log2r.input.bw \
        log2r/Gthe_H3K4me2.log2r.input.bw \
        log2r/Gthe_H3K4me3.log2r.input.bw \
        log2r/Gthe_H3K36me3.log2r.input.bw \
        log2r/Gthe_H3K79me1.log2r.input.bw \
        log2r/Gthe_H3K27me3.log2r.input.bw \
        log2r/Gthe_H3K9me1.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 8 9 \
        --kmeans 8 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        --outFileSortedRegions profiles/log2r.${species}.TEnogene.region.heatmap.k8.txt \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k8.pdf 
    head -1 profiles/log2r.${species}.TEnogene.region.heatmap.k8.txt > profiles/log2r.${species}.TEnogene.region.heatmap.k8.reorder.txt
    grep -v '#' profiles/log2r.${species}.TEnogene.region.heatmap.k8.txt | sed 's/cluster_1/cluster_x/g' | sed 's/cluster_4/cluster_1/g' | sed 's/cluster_6/cluster_4/g' | sed 's/cluster_2/cluster_6/g' | sed 's/cluster_5/cluster_2/g' | sed 's/cluster_3/cluster_5/g' | sed 's/cluster_7/cluster_3/g' | sed 's/cluster_x/cluster_7/g' | sort -s -k13 >> profiles/log2r.${species}.TEnogene.region.heatmap.k8.reorder.txt
    computeMatrix scale-regions \
        -b 500 -a 500 \
        -R profiles/log2r.${species}.TEnogene.region.heatmap.k8.reorder.txt \
        -S \
        log2r/Gthe_H3K9ac.log2r.input.bw \
        log2r/Gthe_H3K27ac.log2r.input.bw \
        log2r/Gthe_H4K16ac.log2r.input.bw \
        log2r/Gthe_H3K4me2.log2r.input.bw \
        log2r/Gthe_H3K4me3.log2r.input.bw \
        log2r/Gthe_H3K36me3.log2r.input.bw \
        log2r/Gthe_H3K79me1.log2r.input.bw \
        log2r/Gthe_H3K27me3.log2r.input.bw \
        log2r/Gthe_H3K9me1.log2r.input.bw \
        --skipZeros --missingDataAsZero \
        -o profiles/log2r.${species}.TEnogene.region.gz 
    plotProfile -m profiles/log2r.${species}.TEnogene.region.gz    \
        --perGroup \
        --plotHeight 18 --plotWidth 13 \
        --colors "#94FA7E" "#257604" "#112B0A" "#FFE04D" "#F46D43" "#B82300" "#F7BEEE" "#4120A9" "#41B6C4" \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        --numPlotsPerRow 1 \
        -out profiles/log2r.${species}.TEnogene.region.plotProfile.k8.pdf
    plotHeatmap -m profiles/log2r.${species}.TEnogene.region.gz   \
        --colorMap coolwarm \
        --whatToShow 'heatmap and colorbar' \
        --sortUsingSamples 8 9 \
        --zMin -1.5 --zMax 1.5 \
        --samplesLabel H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3 H3K36me3 H3K79me1 H3K27me3 H3K9me1 \
        -out profiles/log2r.${species}.TEnogene.region.heatmap.k8.pdf 
fi