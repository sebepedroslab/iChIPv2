# === begin RUN THE PROGRAM ===
# #################################################################################################################################################

for species in Bnat Acas Atha Cfra Ddis Gthe Ngru Nvec Ppat Scer Spun Tthe; do \
    cat /home/smontgomery/cluster/smontgomery/genomes/${species}/${species}.gene.all.bed \
        | grep 'protein_coding' \
        | awk -v OFS='\t' '{if ($6 == "+") print $1,$2,$3,$4,1,$6; else print $1,$2,$3,$4,"-1",$6}' \
        | awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' \
        | sort -k1,1 -k2,2n \
        | awk -v OFS='\t' '{print $0,"gene"}' \
        | bedtools merge -i stdin -c 2,4,5 -o count_distinct,mean,collapse \
        | awk '{OFS="\t"}{if ($4 ==1) print $1,$2,$3,$5}' > /home/smontgomery/cluster/smontgomery/genomes/${species}/${species}.gene.all.score.unique.bed;
    /home/smontgomery/Documents/ucscTools/bedGraphToBigWig /home/smontgomery/cluster/smontgomery/genomes/${species}/${species}.gene.all.score.unique.bed /home/smontgomery/cluster/smontgomery/genomes/${species}/${species}.chromsizes.txt /home/smontgomery/cluster/smontgomery/genomes/${species}/${species}.gene.all.score.bw;

    for mark in H3K9ac H3K27ac H4K16ac H3K4me2 H3K4me3; do \
        cat /home/smontgomery/no_backup/smontgomery/chip/${species}_all/macs2/${species}_${mark}_peaks.narrowPeak \
            | awk '{OFS="\t"}{if ($5 > 100) print $0}' > /home/smontgomery/no_backup/smontgomery/chip/${species}_all/macs2/${species}_${mark}_peaks.narrowPeak.filt100;

        computeMatrix scale-regions \
            -b 1000 -a 1000 \
            -R /home/smontgomery/no_backup/smontgomery/chip/${species}_all/macs2/${species}_${mark}_peaks.narrowPeak.filt100 \
            -S /home/smontgomery/cluster/smontgomery/genomes/${species}/${species}.gene.all.score.bw \
            --skipZeros --missingDataAsZero \
            -o /home/smontgomery/no_backup/smontgomery/chip/${species}_all/macs2/${species}_${mark}_peaks.narrowPeak.filt100.matrix.gz; 
        plotHeatmap -m /home/smontgomery/no_backup/smontgomery/chip/${species}_all/macs2/${species}_${mark}_peaks.narrowPeak.filt100.matrix.gz \
            --colorMap coolwarm \
            --whatToShow 'heatmap and colorbar' \
            --kmeans 10 \
            --samplesLabel Genes \
            --outFileSortedRegions /home/smontgomery/no_backup/smontgomery/chip/${species}_all/macs2/${species}_${mark}_peaks.narrowPeak.filt100.regions.k10.txt \
            -out /home/smontgomery/no_backup/smontgomery/chip/${species}_all/macs2/${species}_${mark}_peaks.narrowPeak.filt100.heatmap.k10.pdf;
    done;
done 

sort -k1,1 -k2,2n no_backup/smontgomery/chip/Bnat_all/macs2/Bnat_H3K4me2_peaks.narrowPeak.filt100.regions.k10.txt \
    | bedtools closest -a stdin -b <(sort -k1,1 -k3,3n cluster/smontgomery/genomes/Bnat/Bnat_expr_sorted.tsv | cut -f1,3,4,6,11,5,7) -k 2 \
    | grep 'cluster_2' \
    | grep -v '-' \
    | sort -k1,1 -k2,2n -k15,15n \
    | awk -v OFS='\t' '{print $1,$2,$3,$4,$13,$20}' \
    | bedtools merge -i stdin -c 6 -o collapse \
    | grep ',' > no_backup/smontgomery/chip/Bnat_all/macs2/Bnat_H3K4me2_peaks.narrowPeak.filt100.regions.k10.RNAcounts.cluster_2.txt
sort -k1,1 -k2,2n no_backup/smontgomery/chip/Bnat_all/macs2/Bnat_H3K4me2_peaks.narrowPeak.filt100.regions.k10.txt \
    | bedtools closest -a stdin -b <(sort -k1,1 -k3,3n cluster/smontgomery/genomes/Bnat/Bnat_expr_sorted.tsv | cut -f1,3,4,6,11,5,7) -k 2 \
    | grep 'cluster_9' \
    | grep -v '+' \
    | sort -k1,1 -k2,2n -k15,15n \
    | awk -v OFS='\t' '{print $1,$2,$3,$4,$13,$20}' \
    | bedtools merge -i stdin -c 6 -o collapse \
    | grep ',' > no_backup/smontgomery/chip/Bnat_all/macs2/Bnat_H3K4me2_peaks.narrowPeak.filt100.regions.k10.RNAcounts.cluster_9.txt

