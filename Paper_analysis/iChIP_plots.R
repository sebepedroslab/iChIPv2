# install.packages("umap")
# install.packages("dbscan")
# install.packages("cowplot")
# install.packages("gridGraphics")
# remotes::install_github('VPetukhov/ggrastr')
library("umap")
library("dbscan")
library("cowplot")
library("gridGraphics")
library("GenomicRanges")
library("plyranges")
library("ggforce")
library("scico")
library("vcd")
library("rtracklayer")
library("ggplot2")
library("Gviz")
library("clustree")
library("cluster")
library("ggrepel")
library("karyoploteR")
label_position <- function(species) {
    if (length(unique(species)) == 1) {
        position <- as.character(unique(species))
    } else {
        position <- NA
    }
    return(position)
}

##Source functions
source("/home/smontgomery/OneDrive/Scripts/iChIP2_paper/functions.R")

##Species names
specieslist=data.frame(species=c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru") ,statenumber=c(16,9,11,16,14,9,12,10,9,11,14,12))
nullList=c("Atha_State12_8","Atha_State12_12","Acas_State14_6","Bnat_State14_13","Cfra_State9_4","Ddis_State9_7","Gthe_State9_4","Ngru_State12_11","Nvec_State16_2","Nvec_State16_6","Ppat_State10_1","Scer_State11_2","Spun_State16_9","Spun_State16_10","Tthe_State11_9")

##Set colours
cols <- c("H3K9ac"="#94FA7E","H3K27ac"="#257604","H4K16ac"="#112B0A","H3K4me2"="#FFE04D","H3K4me3"="#F46D43","H3K36me3"="#B82300","H3K79me1"="#F7BEEE","H3K79me2"="#D448BA","H3K79me3"="#70015D","H3K27me3"="#4120A9","H3K9me1"="#41B6C4","H3K9me2"="#3262AB","H3K9me3"="#3262AB","H3"="grey","input"="black")
##Probably need to add more colours and classes
TEcols=c("DNA"="steelblue","TIR"="steelblue1","Helitron"="steelblue4","TE_DNA_fragment"="cyan","TE_DNA_intact"="darkblue","DNAauto"="aquamarine","DNAnona"="lightgreen",
       "RT"="firebrick3","LTR"="firebrick","TE_RT_LTR_fragment"="pink","TE_RT_LTR_intact"="darkred","TE_RT_LINE"="orange","TE_RT_SINE"="orange3","LINE"="orange","SINE"="orange3","PLE"="yellow4",
       "MITE"="purple","TE_RT_DIRS"="grey","DIRS"="violet","NA"="grey")

##Set directories
raw_folder="/home/smontgomery/cluster/xgraubove/genomes/data/"
genome_folder="/home/smontgomery/cluster/smontgomery/genomes/"
chip_folder="/home/smontgomery/no_backup/smontgomery/chip/"
plot_folder="/home/smontgomery/OneDrive/iChIPv2_SAM_plots/00Revision/"


###############################################################################################################################################################################################################################################################
##Plot genome sizes, TE percentage, gene densities
##Calculate genome coverage
##DONE
tegenomeperc <- data.frame()
for (species in specieslist$species){
    grange = get_our_grange(species=species)
    tes=grange[mcols(grange)$source == "EDTA"]
    if (species == "Acas"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Atha"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Bnat"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k5.reorder.txt"))
    } else if (species == "Cfra"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k5.txt"))
    } else if (species == "Ddis"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k5.reorder.txt"))
    } else if (species == "Gthe"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k8.reorder.txt"))
    } else if (species == "Ngru"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k2.txt"))
    } else if (species == "Nvec"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Ppat"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Scer"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k2.txt"))
    } else if (species == "Spun"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k3.reorder.txt"))
    } else if (species == "Tthe"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k1.txt"))
    }
    
    if (species == "Gthe"){
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.reorder.txt"))
    } else {
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.reorder.txt"))
    }
    tegenomeperc[species,1] <- species
    tegenomeperc[species,2] <- sum(seqlengths(seqinfo(grange)),na.rm = TRUE)
    tegenomeperc[species,3] <- sum(teregions$V11)/sum(seqlengths(seqinfo(grange)),na.rm = TRUE)*100
    tegenomeperc[species,4] <- sum(generegions$V11)/sum(seqlengths(seqinfo(grange)),na.rm = TRUE)*100
}
colnames(tegenomeperc) <- c("species","genomesize","teperc","geneperc")
##Calculate gene density
binsize <- 10000
gene_density_summary <- data.frame()
for (species in specieslist$species){
    gene_density <- data.frame()
    grange = get_our_grange(species=species)
    TSSgrange <- add_TSS(grange)
    for (chr in levels(seqnames(GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1))))) {
        chr_end <- seqlengths(GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))[chr]
        tmp <- data.frame(table(cut(as.data.frame(TSSgrange[seqnames(TSSgrange)==chr])[,2],breaks=c(seq(0,chr_end,binsize),chr_end))))
        tmp <-
            tidyr::separate(tmp,
                            1,
                            into = c("Start", "End"),
                            sep = ",")
        tmp$Start <- as.numeric(gsub("\\(", "", tmp$Start)) + 1
        tmp$End <- as.numeric(gsub("\\]", "", tmp$End))
        
        # add chr info
        tmp$Chr <- chr
        
        # reorder columns
        tmp <- tmp[c(4, 1:3)]
        colnames(tmp) <- c("Chr", "Start", "End", "Count")
        
        # modify last end_pos to chr_end
        tmp[nrow(tmp), "End"] <- chr_end
        gene_density <- (rbind(gene_density,tmp))
    }
    gene_density_summary[species,1] <- species
    gene_density_summary[species,2] <- mean(gene_density$Count)
    gene_density_summary[species,3] <- median(gene_density$Count)
}
colnames(gene_density_summary) <- c("species","mean","median")
##Plots
pdf(paste0(plot_folder,"Fig2_draft/genomesize.pdf"),width=10,height=10)
ggplot(melt(tegenomeperc,id.vars = "species")%>%filter(variable=="genomesize"),aes(y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))),x=variable,fill=value/1000000)) +
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", high = "purple3",name="Mb",limits=c(0,500)) +
    theme_classic()+ # minimal theme
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=18),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          strip.text.x = element_text(size=18))
dev.off()
pdf(paste0(plot_folder,"Fig2_draft/TEperc.pdf"),width=10,height=10)
ggplot(melt(tegenomeperc,id.vars = "species")%>%filter(variable=="teperc"),aes(y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))),x=variable,fill=value)) +
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", mid = "#000000", high="#000000",midpoint = 33,limits=c(0,68),name=NULL) +
    theme_classic()+ # minimal theme
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=18),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          strip.text.x = element_text(size=18))
dev.off()
pdf(paste0(plot_folder,"Fig2_draft/geneperc.pdf"),width=10,height=10)
ggplot(melt(tegenomeperc,id.vars = "species")%>%filter(variable=="geneperc"),aes(y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))),x=variable,fill=value)) +
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", mid = "#000000", high="#000000",midpoint = 33,limits=c(0,68),name=NULL) +
    theme_classic()+ # minimal theme
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=18),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          strip.text.x = element_text(size=18))
dev.off()
pdf(paste0(plot_folder,"Fig2_draft/genedensity.pdf"),width=10,height=10)
ggplot(gene_density_summary,aes(y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))),x=1,fill=mean)) +
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", high = "black",name="Genes per 10kb") +
    theme_classic()+ # minimal theme
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=18),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          strip.text.x = element_text(size=18))
dev.off()

#########################################################################################################################################################################################################################################################################
##Percentage of genome without marks - peak-based - and peak size distributions
##DONE
null_genome_perc <- data.frame()
totPeakSizes <- data.frame()
totGeneRegions <- data.frame()
crm=c("Pt","Mt","Mito","ChrM","ChrmC","chrmt")
for (species in specieslist$species){
    grange = get_our_grange(species=species)
    genome_size <- sum(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2)
    if (species == "Gthe"){
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.reorder.txt"))
    } else {
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.reorder.txt"))
    }
    generegions$species <- species
    totGeneRegions <- rbind(totGeneRegions,generegions)
    H3K4me3peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K4me3_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(H3K4me3peaksbroad)[1:3] <- c("chrom","start","end")
    H3K4me3peaksbroad <- c(GRanges(H3K4me3peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K4me3peaksbroad=H3K4me3peaksbroad[!(seqnames(H3K4me3peaksbroad) %in% c(crm))]
    H3K4me3peaksbroad <- H3K4me3peaksbroad%>%mutate(size=width(H3K4me3peaksbroad))
    H3K4me3peaksbroad$numGene <- countOverlaps(H3K4me3peaksbroad,grange[mcols(grange)$type=="protein_coding"])

    H3K36me3peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K36me3_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(H3K36me3peaksbroad)[1:3] <- c("chrom","start","end")
    H3K36me3peaksbroad <- c(GRanges(H3K36me3peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K36me3peaksbroad=H3K36me3peaksbroad[!(seqnames(H3K36me3peaksbroad) %in% c(crm))]
    H3K36me3peaksbroad <- H3K36me3peaksbroad%>%mutate(size=width(H3K36me3peaksbroad))
    H3K36me3peaksbroad$numGene <- countOverlaps(H3K36me3peaksbroad,grange[mcols(grange)$type=="protein_coding"])

    H3K4me2peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K4me2_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(H3K4me2peaksbroad)[1:3] <- c("chrom","start","end")
    H3K4me2peaksbroad <- c(GRanges(H3K4me2peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K4me2peaksbroad=H3K4me2peaksbroad[!(seqnames(H3K4me2peaksbroad) %in% c(crm))]
    H3K4me2peaksbroad <- H3K4me2peaksbroad%>%mutate(size=width(H3K4me2peaksbroad))
    H3K4me2peaksbroad$numGene <- countOverlaps(H3K4me2peaksbroad,grange[mcols(grange)$type=="protein_coding"])

    H3K9acpeaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K9ac_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(H3K9acpeaksbroad)[1:3] <- c("chrom","start","end")
    H3K9acpeaksbroad <- c(GRanges(H3K9acpeaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K9acpeaksbroad=H3K9acpeaksbroad[!(seqnames(H3K9acpeaksbroad) %in% c(crm))]
    H3K9acpeaksbroad <- H3K9acpeaksbroad%>%mutate(size=width(H3K9acpeaksbroad))
    H3K9acpeaksbroad$numGene <- countOverlaps(H3K9acpeaksbroad,grange[mcols(grange)$type=="protein_coding"])

    H3K27acpeaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K27ac_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(H3K27acpeaksbroad)[1:3] <- c("chrom","start","end")
    H3K27acpeaksbroad <- c(GRanges(H3K27acpeaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K27acpeaksbroad=H3K27acpeaksbroad[!(seqnames(H3K27acpeaksbroad) %in% c(crm))]
    H3K27acpeaksbroad <- H3K27acpeaksbroad%>%mutate(size=width(H3K27acpeaksbroad))
    H3K27acpeaksbroad$numGene <- countOverlaps(H3K27acpeaksbroad,grange[mcols(grange)$type=="protein_coding"])

    H4K16acpeaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H4K16ac_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(H4K16acpeaksbroad)[1:3] <- c("chrom","start","end")
    H4K16acpeaksbroad <- c(GRanges(H4K16acpeaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H4K16acpeaksbroad=H4K16acpeaksbroad[!(seqnames(H4K16acpeaksbroad) %in% c(crm))]
    H4K16acpeaksbroad <- H4K16acpeaksbroad%>%mutate(size=width(H4K16acpeaksbroad))
    H4K16acpeaksbroad$numGene <- countOverlaps(H4K16acpeaksbroad,grange[mcols(grange)$type=="protein_coding"])

    H3K9acpeaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K9ac_peaks.narrowPeak"),sep='\t',header = FALSE)
    colnames(H3K9acpeaks)[1:3] <- c("chrom","start","end")
    H3K9acpeaks <- c(GRanges(H3K9acpeaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K9acpeaks=H3K9acpeaks[!(seqnames(H3K9acpeaks) %in% c(crm))]
    H3K9acpeaks <- H3K9acpeaks%>%mutate(size=width(H3K9acpeaks))
    H3K9acpeaks$numGene <- countOverlaps(H3K9acpeaks,grange[mcols(grange)$type=="protein_coding"])

    H3K27acpeaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K27ac_peaks.narrowPeak"),sep='\t',header = FALSE)
    colnames(H3K27acpeaks)[1:3] <- c("chrom","start","end")
    H3K27acpeaks <- c(GRanges(H3K27acpeaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K27acpeaks=H3K27acpeaks[!(seqnames(H3K27acpeaks) %in% c(crm))]
    H3K27acpeaks <- H3K27acpeaks%>%mutate(size=width(H3K27acpeaks))
    H3K27acpeaks$numGene <- countOverlaps(H3K27acpeaks,grange[mcols(grange)$type=="protein_coding"])

    H4K16acpeaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H4K16ac_peaks.narrowPeak"),sep='\t',header = FALSE)
    colnames(H4K16acpeaks)[1:3] <- c("chrom","start","end")
    H4K16acpeaks <- c(GRanges(H4K16acpeaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H4K16acpeaks=H4K16acpeaks[!(seqnames(H4K16acpeaks) %in% c(crm))]
    H4K16acpeaks <- H4K16acpeaks%>%mutate(size=width(H4K16acpeaks))
    H4K16acpeaks$numGene <- countOverlaps(H4K16acpeaks,grange[mcols(grange)$type=="protein_coding"])

    H3K4me3peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K4me3_peaks.narrowPeak"),sep='\t',header = FALSE)
    colnames(H3K4me3peaks)[1:3] <- c("chrom","start","end")
    H3K4me3peaks <- c(GRanges(H3K4me3peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K4me3peaks=H3K4me3peaks[!(seqnames(H3K4me3peaks) %in% c(crm))]
    H3K4me3peaks <- H3K4me3peaks%>%mutate(size=width(H3K4me3peaks))
    H3K4me3peaks$numGene <- countOverlaps(H3K4me3peaks,grange[mcols(grange)$type=="protein_coding"])

    H3K4me2peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K4me2_peaks.narrowPeak"),sep='\t',header = FALSE)
    colnames(H3K4me2peaks)[1:3] <- c("chrom","start","end")
    H3K4me2peaks <- c(GRanges(H3K4me2peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K4me2peaks=H3K4me2peaks[!(seqnames(H3K4me2peaks) %in% c(crm))]
    H3K4me2peaks <- H3K4me2peaks%>%mutate(size=width(H3K4me2peaks))
    H3K4me2peaks$numGene <- countOverlaps(H3K4me2peaks,grange[mcols(grange)$type=="protein_coding"])

    H3K36me3peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K36me3_peaks.narrowPeak"),sep='\t',header = FALSE)
    colnames(H3K36me3peaks)[1:3] <- c("chrom","start","end")
    H3K36me3peaks <- c(GRanges(H3K36me3peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3K36me3peaks=H3K36me3peaks[!(seqnames(H3K36me3peaks) %in% c(crm))]
    H3K36me3peaks <- H3K36me3peaks%>%mutate(size=width(H3K36me3peaks))
    H3K36me3peaks$numGene <- countOverlaps(H3K36me3peaks,grange[mcols(grange)$type=="protein_coding"])

    if (species == "Atha" | species == "Ppat"){
        H3K9me3peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K9me2_peaks.narrowPeak"),sep='\t',header = FALSE)
        colnames(H3K9me3peaks)[1:3] <- c("chrom","start","end")
        H3K9me3peaks <- c(GRanges(H3K9me3peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K9me3peaks=H3K9me3peaks[!(seqnames(H3K9me3peaks) %in% c(crm))]
        H3K9me3peaks <- H3K9me3peaks%>%mutate(size=width(H3K9me3peaks))
        H3K9me3peaks$numGene <- countOverlaps(H3K9me3peaks,grange[mcols(grange)$type=="protein_coding"])
        H3K9me3peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K9me2_peaks.broadPeak"),sep='\t',header = FALSE)
        colnames(H3K9me3peaksbroad)[1:3] <- c("chrom","start","end")
        H3K9me3peaksbroad <- c(GRanges(H3K9me3peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K9me3peaksbroad=H3K9me3peaksbroad[!(seqnames(H3K9me3peaksbroad) %in% c(crm))]
        H3K9me3peaksbroad <- H3K9me3peaksbroad%>%mutate(size=width(H3K9me3peaksbroad))
        H3K9me3peaksbroad$numGene <- countOverlaps(H3K9me3peaksbroad,grange[mcols(grange)$type=="protein_coding"])
    } 
    if (species == "Nvec" | species == "Spun" | species == "Acas" | species == "Ddis" | species == "Tthe" | species == "Bnat" | species == "Ngru") {
        H3K9me3peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K9me3_peaks.narrowPeak"),sep='\t',header = FALSE)
        colnames(H3K9me3peaks)[1:3] <- c("chrom","start","end")
        H3K9me3peaks <- c(GRanges(H3K9me3peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K9me3peaks=H3K9me3peaks[!(seqnames(H3K9me3peaks) %in% c(crm))]
        H3K9me3peaks <- H3K9me3peaks%>%mutate(size=width(H3K9me3peaks))
        H3K9me3peaks$numGene <- countOverlaps(H3K9me3peaks,grange[mcols(grange)$type=="protein_coding"])
        H3K9me3peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K9me3_peaks.broadPeak"),sep='\t',header = FALSE)
        colnames(H3K9me3peaksbroad)[1:3] <- c("chrom","start","end")
        H3K9me3peaksbroad <- c(GRanges(H3K9me3peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K9me3peaksbroad=H3K9me3peaksbroad[!(seqnames(H3K9me3peaksbroad) %in% c(crm))]
        H3K9me3peaksbroad <- H3K9me3peaksbroad%>%mutate(size=width(H3K9me3peaksbroad))
        H3K9me3peaksbroad$numGene <- countOverlaps(H3K9me3peaksbroad,grange[mcols(grange)$type=="protein_coding"])
    }
    if (species == "Nvec" | species == "Spun" | species == "Acas" | species == "Ddis" | species == "Bnat" | species == "Atha" | species == "Ppat" | species == "Gthe" | species == "Ngru") {
        H3K9me1peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K9me1_peaks.narrowPeak"),sep='\t',header = FALSE)
        colnames(H3K9me1peaks)[1:3] <- c("chrom","start","end")
        H3K9me1peaks <- c(GRanges(H3K9me1peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K9me1peaks=H3K9me1peaks[!(seqnames(H3K9me1peaks) %in% c(crm))]
        H3K9me1peaks <- H3K9me1peaks%>%mutate(size=width(H3K9me1peaks))
        H3K9me1peaks$numGene <- countOverlaps(H3K9me1peaks,grange[mcols(grange)$type=="protein_coding"])
        H3K9me1peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K9me1_peaks.broadPeak"),sep='\t',header = FALSE)
        colnames(H3K9me1peaksbroad)[1:3] <- c("chrom","start","end")
        H3K9me1peaksbroad <- c(GRanges(H3K9me1peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K9me1peaksbroad=H3K9me1peaksbroad[!(seqnames(H3K9me1peaksbroad) %in% c(crm))]
        H3K9me1peaksbroad <- H3K9me1peaksbroad%>%mutate(size=width(H3K9me1peaksbroad))
        H3K9me1peaksbroad$numGene <- countOverlaps(H3K9me1peaksbroad,grange[mcols(grange)$type=="protein_coding"])
    }
    if (species == "Nvec" | species == "Spun" | species == "Acas" | species == "Ddis" | species == "Bnat" | species == "Atha" | species == "Ppat" | species == "Gthe" | species == "Tthe" | species == "Ngru") {
        H3K27me3peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K27me3_peaks.narrowPeak"),sep='\t',header = FALSE)
        colnames(H3K27me3peaks)[1:3] <- c("chrom","start","end")
        H3K27me3peaks <- c(GRanges(H3K27me3peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K27me3peaks=H3K27me3peaks[!(seqnames(H3K27me3peaks) %in% c(crm))]
        H3K27me3peaks <- H3K27me3peaks%>%mutate(size=width(H3K27me3peaks))
        H3K27me3peaks$numGene <- countOverlaps(H3K27me3peaks,grange[mcols(grange)$type=="protein_coding"])
        H3K27me3peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K27me3_peaks.broadPeak"),sep='\t',header = FALSE)
        colnames(H3K27me3peaksbroad)[1:3] <- c("chrom","start","end")
        H3K27me3peaksbroad <- c(GRanges(H3K27me3peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K27me3peaksbroad=H3K27me3peaksbroad[!(seqnames(H3K27me3peaksbroad) %in% c(crm))]
        H3K27me3peaksbroad <- H3K27me3peaksbroad%>%mutate(size=width(H3K27me3peaksbroad))
        H3K27me3peaksbroad$numGene <- countOverlaps(H3K27me3peaksbroad,grange[mcols(grange)$type=="protein_coding"])
    }
    tryCatch({
    if (species != "Atha" & species != "Ppat" & species != "Ngru") {
        H3K79me1peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K79me1_peaks.narrowPeak"),sep='\t',header = FALSE)
        colnames(H3K79me1peaks)[1:3] <- c("chrom","start","end")
        H3K79me1peaks <- c(GRanges(H3K79me1peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K79me1peaks=H3K79me1peaks[!(seqnames(H3K79me1peaks) %in% c(crm))]
        H3K79me1peaks <- H3K79me1peaks%>%mutate(size=width(H3K79me1peaks))
        H3K79me1peaks$numGene <- countOverlaps(H3K79me1peaks,grange[mcols(grange)$type=="protein_coding"])
        H3K79me1peaksbroad <- GRanges()
        H3K79me1peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K79me1_peaks.broadPeak"),sep='\t',header = FALSE)
        colnames(H3K79me1peaksbroad)[1:3] <- c("chrom","start","end")
        H3K79me1peaksbroad <- c(GRanges(H3K79me1peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K79me1peaksbroad=H3K79me1peaksbroad[!(seqnames(H3K79me1peaksbroad) %in% c(crm))]
        H3K79me1peaksbroad <- H3K79me1peaksbroad%>%mutate(size=width(H3K79me1peaksbroad))
        H3K79me1peaksbroad$numGene <- countOverlaps(H3K79me1peaksbroad,grange[mcols(grange)$type=="protein_coding"])
    }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    tryCatch({
    if (species != "Atha" & species != "Ppat" & species != "Ngru" & species != "Gthe" & species != "Tthe") {
        H3K79me2peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K79me2_peaks.narrowPeak"),sep='\t',header = FALSE)
        colnames(H3K79me2peaks)[1:3] <- c("chrom","start","end")
        H3K79me2peaks <- c(GRanges(H3K79me2peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K79me2peaks=H3K79me2peaks[!(seqnames(H3K79me2peaks) %in% c(crm))]
        H3K79me2peaks <- H3K79me2peaks%>%mutate(size=width(H3K79me2peaks))
        H3K79me2peaks$numGene <- countOverlaps(H3K79me2peaks,grange[mcols(grange)$type=="protein_coding"])
        H3K79me2peaksbroad <- GRanges()
        H3K79me2peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K79me2_peaks.broadPeak"),sep='\t',header = FALSE)
        colnames(H3K79me2peaksbroad)[1:3] <- c("chrom","start","end")
        H3K79me2peaksbroad <- c(GRanges(H3K79me2peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K79me2peaksbroad=H3K79me2peaksbroad[!(seqnames(H3K79me2peaksbroad) %in% c(crm))]
        H3K79me2peaksbroad <- H3K79me2peaksbroad%>%mutate(size=width(H3K79me2peaksbroad))
        H3K79me2peaksbroad$numGene <- countOverlaps(H3K79me2peaksbroad,grange[mcols(grange)$type=="protein_coding"])
    }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    tryCatch({
    if (species != "Atha" & species != "Ppat" & species != "Ngru" & species != "Gthe" & species != "Tthe") {
        H3K79me3peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K79me3_peaks.narrowPeak"),sep='\t',header = FALSE)
        colnames(H3K79me3peaks)[1:3] <- c("chrom","start","end")
        H3K79me3peaks <- c(GRanges(H3K79me3peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K79me3peaks=H3K79me3peaks[!(seqnames(H3K79me3peaks) %in% c(crm))]
        H3K79me3peaks <- H3K79me3peaks%>%mutate(size=width(H3K79me3peaks))
        H3K79me3peaks$numGene <- countOverlaps(H3K79me3peaks,grange[mcols(grange)$type=="protein_coding"])
        H3K79me3peaksbroad <- GRanges()
        H3K79me3peaksbroad <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3K79me3_peaks.broadPeak"),sep='\t',header = FALSE)
        colnames(H3K79me3peaksbroad)[1:3] <- c("chrom","start","end")
        H3K79me3peaksbroad <- c(GRanges(H3K79me3peaksbroad),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
        H3K79me3peaksbroad=H3K79me3peaksbroad[!(seqnames(H3K79me3peaksbroad) %in% c(crm))]
        H3K79me3peaksbroad <- H3K79me3peaksbroad%>%mutate(size=width(H3K79me3peaksbroad))
        H3K79me3peaksbroad$numGene <- countOverlaps(H3K79me3peaksbroad,grange[mcols(grange)$type=="protein_coding"])
    }}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    null_genome_perc[species,1] <- species
    null_genome_perc[species,2] <- (sum(width(union_ranges(H3K9acpeaksbroad,H3K9acpeaks)))/genome_size*100)
    null_genome_perc[species,3] <- (sum(width(union_ranges(H3K27acpeaksbroad,H3K27acpeaks)))/genome_size*100)
    null_genome_perc[species,4] <- (sum(width(union_ranges(H4K16acpeaksbroad,H4K16acpeaks)))/genome_size*100)
    null_genome_perc[species,5] <- (sum(width(union_ranges(H3K4me2peaksbroad,H3K4me2peaks)))/genome_size*100)
    null_genome_perc[species,6] <- (sum(width(union_ranges(H3K4me3peaksbroad,H3K4me3peaks)))/genome_size*100)
    null_genome_perc[species,7] <- (sum(width(union_ranges(H3K36me3peaksbroad,H3K36me3peaks)))/genome_size*100)
    if (species == "Nvec" | species == "Spun" | species == "Acas" | species == "Ddis" | species == "Bnat") {
        null_genome_perc[species,8] <- (sum(width(union_ranges(H3K79me1peaksbroad,H3K79me1peaks)))/genome_size*100)
        null_genome_perc[species,9] <- (sum(width(union_ranges(H3K79me2peaksbroad,H3K79me2peaks)))/genome_size*100)
        null_genome_perc[species,10] <- (sum(width(union_ranges(H3K79me3peaksbroad,H3K79me3peaks)))/genome_size*100)
        null_genome_perc[species,11] <- (sum(width(union_ranges(H3K27me3peaksbroad,H3K27me3peaks)))/genome_size*100)
        null_genome_perc[species,12] <- (sum(width(union_ranges(H3K9me1peaksbroad,H3K9me1peaks)))/genome_size*100)
        null_genome_perc[species,13] <- (sum(width(union_ranges(H3K9me3peaksbroad,H3K9me3peaks)))/genome_size*100)
        null_genome_perc[species,14] <- (genome_size-sum(width(
            union_ranges(H3K79me2peaksbroad,union_ranges(H3K79me3peaksbroad,union_ranges(H3K79me1peaksbroad,
            union_ranges(H3K79me1peaks,union_ranges(H3K79me2peaks,union_ranges(H3K79me3peaks,
            union_ranges(H3K9me1peaksbroad,union_ranges(H3K9me3peaksbroad,union_ranges(H3K27me3peaksbroad,
            union_ranges(H3K9me1peaks,union_ranges(H3K9me3peaks,union_ranges(H3K27me3peaks,
            union_ranges(H4K16acpeaksbroad,union_ranges(H3K9acpeaksbroad,union_ranges(H3K27acpeaksbroad,
            union_ranges(H4K16acpeaks,union_ranges(H3K9acpeaks,union_ranges(H3K27acpeaks,
            union_ranges(H3K36me3peaksbroad,union_ranges(H3K4me2peaksbroad,union_ranges(H3K4me3peaksbroad,
            union_ranges(H3K4me3peaks,union_ranges(H3K36me3peaks,H3K4me2peaks))))))))))))))))))))))))))/genome_size*100
        # pdf(paste0(plot_folder,"misc/peakWidths_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #     geom_density(data=as.data.frame(H3K9acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K27acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me2peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K79me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K79me2peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K79me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K9me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K9me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K27me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
        #     scale_colour_manual(values=cols) +
        #     theme_classic() +
        #     theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()
        # pdf(paste0(plot_folder,"misc/peakGeneOverlaps_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #   geom_histogram(data=as.data.frame(H3K9acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K27acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H4K16acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me2peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K36me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K79me1peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K79me2peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K79me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K9me1peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K9me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K27me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   scale_x_continuous(name="Number of Genes") +
        #   scale_colour_manual(values=cols) +
        #   theme_classic() +
        #   theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()
        pdf(paste0(plot_folder,"extData_Fig5_draft/peakWidth_GeneOverlaps_",species,".pdf"),width=8,height=8)
        print(ggplot() +
          geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K79me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K79me2peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K79me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K9me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K9me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K27me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=generegions,aes(log10(V11),col="H3")) +
          geom_text(aes(label=c("Mean","Max"),x=c(4.25,4.75),y=2),size=5) +
          geom_text(aes(label = c(round(mean(H4K16acpeaks$numGene),2),max(H4K16acpeaks$numGene)), x=c(4.25,4.75),y = 1.9,col="H4K16ac"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K36me3peaks$numGene),2),max(H3K36me3peaks$numGene)), x=c(4.25,4.75),y = 1.8,col="H3K36me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K79me1peaks$numGene),2),max(H3K79me1peaks$numGene)), x=c(4.25,4.75),y = 1.7,col="H3K79me1"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K79me2peaks$numGene),2),max(H3K79me2peaks$numGene)), x=c(4.25,4.75),y = 1.6,col="H3K79me2"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K79me3peaks$numGene),2),max(H3K79me3peaks$numGene)), x=c(4.25,4.75),y = 1.5,col="H3K79me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K27me3peaks$numGene),2),max(H3K27me3peaks$numGene)), x=c(4.25,4.75),y = 1.4,col="H3K27me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K9me1peaks$numGene),2),max(H3K9me1peaks$numGene)), x=c(4.25,4.75),y = 1.3,col="H3K9me1"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K9me3peaks$numGene),2),max(H3K9me3peaks$numGene)), x=c(4.25,4.75),y = 1.2,col="H3K9me3"),size=5,show.legend = FALSE) +
          scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
          scale_y_continuous(name="Density") +
          scale_colour_manual(values=cols,labels=c("H3"="Genes")) +
          theme_classic() +
          theme(axis.text.x = element_text(size=18),
                axis.text.y = element_text(size=18),
                axis.title.x = element_text(size=18),
                axis.title.y = element_text(size=18),
                legend.text = element_text(size=18),
                legend.title = element_blank(),
                strip.text.x = element_text(size=18)))
        dev.off()
        totPeakSizes <- rbind(totPeakSizes,as.data.frame(c(H3K9acpeaks,H3K27acpeaks,H4K16acpeaks,H3K4me2peaks,H3K4me3peaks,H3K36me3peaks,H3K79me1peaks,H3K79me2peaks,H3K79me3peaks,H3K27me3peaks,H3K9me1peaks,H3K9me3peaks))%>%mutate(mark=str_split(V4,"_",simplify = T)[,2],species=str_split(V4,"_",simplify = T)[,1]))
    } else if (species == "Atha" | species == "Ppat" | species == "Ngru"){
        null_genome_perc[species,8] <- 0
        null_genome_perc[species,9] <- 0
        null_genome_perc[species,10] <- 0
        null_genome_perc[species,11] <- (sum(width(union_ranges(H3K27me3peaksbroad,H3K27me3peaks)))/genome_size*100)
        null_genome_perc[species,12] <- (sum(width(union_ranges(H3K9me1peaksbroad,H3K9me1peaks)))/genome_size*100)
        null_genome_perc[species,13] <- (sum(width(union_ranges(H3K9me3peaksbroad,H3K9me3peaks)))/genome_size*100)
        null_genome_perc[species,14] <- (genome_size-sum(width(
            union_ranges(H3K9me1peaksbroad,union_ranges(H3K9me3peaksbroad,union_ranges(H3K27me3peaksbroad,
            union_ranges(H3K9me1peaks,union_ranges(H3K9me3peaks,union_ranges(H3K27me3peaks,
            union_ranges(H4K16acpeaksbroad,union_ranges(H3K9acpeaksbroad,union_ranges(H3K27acpeaksbroad,
            union_ranges(H4K16acpeaks,union_ranges(H3K9acpeaks,union_ranges(H3K27acpeaks,
            union_ranges(H3K36me3peaksbroad,union_ranges(H3K4me2peaksbroad,union_ranges(H3K4me3peaksbroad,
            union_ranges(H3K4me3peaks,union_ranges(H3K36me3peaks,H3K4me2peaks))))))))))))))))))))/genome_size*100      
        # pdf(paste0(plot_folder,"misc/peakWidths_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #     geom_density(data=as.data.frame(H3K9acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K27acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me2peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K9me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K9me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K27me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
        #     scale_colour_manual(values=cols) +
        #     theme_classic() +
        #     theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()
        # pdf(paste0(plot_folder,"misc/peakGeneOverlaps_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #   geom_histogram(data=as.data.frame(H3K9acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K27acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H4K16acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me2peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K36me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K9me1peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K9me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K27me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   scale_x_continuous(name="Number of Genes") +
        #   scale_colour_manual(values=cols) +
        #   theme_classic() +
        #   theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()
        pdf(paste0(plot_folder,"extData_Fig5_draft/peakWidth_GeneOverlaps_",species,".pdf"),width=8,height=8)
        print(ggplot() +
          geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K9me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K9me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K27me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=generegions,aes(log10(V11),col="H3")) +
          geom_text(aes(label=c("Mean","Max"),x=c(4.25,4.75),y=2),size=5) +
          geom_text(aes(label = c(round(mean(H4K16acpeaks$numGene),2),max(H4K16acpeaks$numGene)), x=c(4.25,4.75),y = 1.9,col="H4K16ac"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K36me3peaks$numGene),2),max(H3K36me3peaks$numGene)), x=c(4.25,4.75),y = 1.8,col="H3K36me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K27me3peaks$numGene),2),max(H3K27me3peaks$numGene)), x=c(4.25,4.75),y = 1.4,col="H3K27me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K9me1peaks$numGene),2),max(H3K9me1peaks$numGene)), x=c(4.25,4.75),y = 1.3,col="H3K9me1"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K9me3peaks$numGene),2),max(H3K9me3peaks$numGene)), x=c(4.25,4.75),y = 1.2,col="H3K9me3"),size=5,show.legend = FALSE) +
          scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
          scale_y_continuous(name="Density") +
          scale_colour_manual(values=cols,labels=c("H3"="Genes")) +
          theme_classic() +
          theme(axis.text.x = element_text(size=18),
                axis.text.y = element_text(size=18),
                axis.title.x = element_text(size=18),
                axis.title.y = element_text(size=18),
                legend.text = element_text(size=18),
                legend.title = element_blank(),
                strip.text.x = element_text(size=18)))
        dev.off()
        totPeakSizes <- rbind(totPeakSizes,as.data.frame(c(H3K9acpeaks,H3K27acpeaks,H4K16acpeaks,H3K4me2peaks,H3K4me3peaks,H3K36me3peaks,H3K27me3peaks,H3K9me1peaks,H3K9me3peaks))%>%mutate(mark=str_split(V4,"_",simplify = T)[,2],species=str_split(V4,"_",simplify = T)[,1]))
    } else if (species == "Scer" | species == "Cfra"){
        null_genome_perc[species,8] <- (sum(width(union_ranges(H3K79me1peaksbroad,H3K79me1peaks)))/genome_size*100)
        null_genome_perc[species,9] <- (sum(width(union_ranges(H3K79me2peaksbroad,H3K79me2peaks)))/genome_size*100)
        null_genome_perc[species,10] <- (sum(width(union_ranges(H3K79me3peaksbroad,H3K79me3peaks)))/genome_size*100)
        null_genome_perc[species,11] <- 0
        null_genome_perc[species,12] <- 0
        null_genome_perc[species,13] <- 0
        null_genome_perc[species,14] <- (genome_size-sum(width(
            union_ranges(H3K79me2peaksbroad,union_ranges(H3K79me3peaksbroad,union_ranges(H3K79me1peaksbroad,
            union_ranges(H3K79me1peaks,union_ranges(H3K79me2peaks,union_ranges(H3K79me3peaks,
            union_ranges(H4K16acpeaksbroad,union_ranges(H3K9acpeaksbroad,union_ranges(H3K27acpeaksbroad,
            union_ranges(H4K16acpeaks,union_ranges(H3K9acpeaks,union_ranges(H3K27acpeaks,
            union_ranges(H3K36me3peaksbroad,union_ranges(H3K4me2peaksbroad,union_ranges(H3K4me3peaksbroad,
            union_ranges(H3K4me3peaks,union_ranges(H3K36me3peaks,H3K4me2peaks))))))))))))))))))))/genome_size*100      
        # pdf(paste0(plot_folder,"misc/peakWidths_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #     geom_density(data=as.data.frame(H3K9acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K27acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me2peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K79me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K79me2peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K79me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
        #     scale_colour_manual(values=cols) +
        #     theme_classic() +
        #     theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()  
        # pdf(paste0(plot_folder,"misc/peakGeneOverlaps_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #   geom_histogram(data=as.data.frame(H3K9acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K27acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H4K16acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me2peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K36me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K79me1peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K79me2peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K79me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   scale_x_continuous(name="Number of Genes") +
        #   scale_colour_manual(values=cols) +
        #   theme_classic() +
        #   theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()
        pdf(paste0(plot_folder,"extData_Fig5_draft/peakWidth_GeneOverlaps_",species,".pdf"),width=8,height=8)
        print(ggplot() +
          geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K79me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K79me2peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K79me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=generegions,aes(log10(V11),col="H3")) +
          geom_text(aes(label=c("Mean","Max"),x=c(4.25,4.75),y=2),size=5) +
          geom_text(aes(label = c(round(mean(H4K16acpeaks$numGene),2),max(H4K16acpeaks$numGene)), x=c(4.25,4.75),y = 1.9,col="H4K16ac"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K36me3peaks$numGene),2),max(H3K36me3peaks$numGene)), x=c(4.25,4.75),y = 1.8,col="H3K36me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K79me1peaks$numGene),2),max(H3K79me1peaks$numGene)), x=c(4.25,4.75),y = 1.7,col="H3K79me1"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K79me2peaks$numGene),2),max(H3K79me2peaks$numGene)), x=c(4.25,4.75),y = 1.6,col="H3K79me2"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K79me3peaks$numGene),2),max(H3K79me3peaks$numGene)), x=c(4.25,4.75),y = 1.5,col="H3K79me3"),size=5,show.legend = FALSE) +
          scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
          scale_y_continuous(name="Density") +
          scale_colour_manual(values=cols,labels=c("H3"="Genes")) +
          theme_classic() +
          theme(axis.text.x = element_text(size=18),
                axis.text.y = element_text(size=18),
                axis.title.x = element_text(size=18),
                axis.title.y = element_text(size=18),
                legend.text = element_text(size=18),
                legend.title = element_blank(),
                strip.text.x = element_text(size=18)))
        dev.off()
        totPeakSizes <- rbind(totPeakSizes,as.data.frame(c(H3K9acpeaks,H3K27acpeaks,H4K16acpeaks,H3K4me2peaks,H3K4me3peaks,H3K36me3peaks,H3K79me1peaks,H3K79me2peaks,H3K79me3peaks))%>%mutate(mark=str_split(V4,"_",simplify = T)[,2],species=str_split(V4,"_",simplify = T)[,1]))
    } else if (species == "Gthe"){
        null_genome_perc[species,8] <- (sum(width(union_ranges(H3K79me1peaksbroad,H3K79me1peaks)))/genome_size*100)
        null_genome_perc[species,9] <- 0
        null_genome_perc[species,10] <- 0
        null_genome_perc[species,11] <- (sum(width(union_ranges(H3K27me3peaksbroad,H3K27me3peaks)))/genome_size*100)
        null_genome_perc[species,12] <- (sum(width(union_ranges(H3K9me1peaksbroad,H3K9me1peaks)))/genome_size*100)
        null_genome_perc[species,13] <- 0
        null_genome_perc[species,14] <- (genome_size-sum(width(
            union_ranges(H3K79me1peaksbroad,
            union_ranges(H3K79me1peaks,
            union_ranges(H3K9me1peaksbroad,union_ranges(H3K27me3peaksbroad,
            union_ranges(H3K9me1peaks,union_ranges(H3K27me3peaks,
            union_ranges(H4K16acpeaksbroad,union_ranges(H3K9acpeaksbroad,union_ranges(H3K27acpeaksbroad,
            union_ranges(H4K16acpeaks,union_ranges(H3K9acpeaks,union_ranges(H3K27acpeaks,
            union_ranges(H3K36me3peaksbroad,union_ranges(H3K4me2peaksbroad,union_ranges(H3K4me3peaksbroad,
            union_ranges(H3K4me3peaks,union_ranges(H3K36me3peaks,H3K4me2peaks))))))))))))))))))))/genome_size*100   
        # pdf(paste0(plot_folder,"misc/peakWidths_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #     geom_density(data=as.data.frame(H3K9acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K27acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me2peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K79me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K9me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K27me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
        #     scale_colour_manual(values=cols) +
        #     theme_classic() +
        #     theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()     
        # pdf(paste0(plot_folder,"misc/peakGeneOverlaps_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #   geom_histogram(data=as.data.frame(H3K9acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K27acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H4K16acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me2peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K36me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K79me1peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K9me1peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K27me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   scale_x_continuous(name="Number of Genes") +
        #   scale_colour_manual(values=cols) +
        #   theme_classic() +
        #   theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()
        pdf(paste0(plot_folder,"extData_Fig5_draft/peakWidth_GeneOverlaps_",species,".pdf"),width=8,height=8)
        print(ggplot() +
          geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K79me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K9me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K27me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=generegions,aes(log10(V11),col="H3")) +
          geom_text(aes(label=c("Mean","Max"),x=c(4.25,4.75),y=2),size=5) +
          geom_text(aes(label = c(round(mean(H4K16acpeaks$numGene),2),max(H4K16acpeaks$numGene)), x=c(4.25,4.75),y = 1.9,col="H4K16ac"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K36me3peaks$numGene),2),max(H3K36me3peaks$numGene)), x=c(4.25,4.75),y = 1.8,col="H3K36me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K79me1peaks$numGene),2),max(H3K79me1peaks$numGene)), x=c(4.25,4.75),y = 1.7,col="H3K79me1"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K27me3peaks$numGene),2),max(H3K27me3peaks$numGene)), x=c(4.25,4.75),y = 1.4,col="H3K27me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K9me1peaks$numGene),2),max(H3K9me1peaks$numGene)), x=c(4.25,4.75),y = 1.3,col="H3K9me1"),size=5,show.legend = FALSE) +
          scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
          scale_y_continuous(name="Density") +
          scale_colour_manual(values=cols,labels=c("H3"="Genes")) +
          theme_classic() +
          theme(axis.text.x = element_text(size=18),
                axis.text.y = element_text(size=18),
                axis.title.x = element_text(size=18),
                axis.title.y = element_text(size=18),
                legend.text = element_text(size=18),
                legend.title = element_blank(),
                strip.text.x = element_text(size=18)))
        dev.off()
        totPeakSizes <- rbind(totPeakSizes,as.data.frame(c(H3K9acpeaks,H3K27acpeaks,H4K16acpeaks,H3K4me2peaks,H3K4me3peaks,H3K36me3peaks,H3K79me1peaks,H3K27me3peaks,H3K9me1peaks))%>%mutate(mark=str_split(V4,"_",simplify = T)[,2],species=str_split(V4,"_",simplify = T)[,1]))
    } else if (species == "Tthe"){
        null_genome_perc[species,8] <- (sum(width(union_ranges(H3K79me1peaksbroad,H3K79me1peaks)))/genome_size*100)
        null_genome_perc[species,9] <- 0
        null_genome_perc[species,10] <- 0
        null_genome_perc[species,11] <- (sum(width(union_ranges(H3K27me3peaksbroad,H3K27me3peaks)))/genome_size*100)
        null_genome_perc[species,12] <- 0
        null_genome_perc[species,13] <- (sum(width(union_ranges(H3K9me3peaksbroad,H3K9me3peaks)))/genome_size*100)
        null_genome_perc[species,14] <- (genome_size-sum(width(
            union_ranges(H3K79me2peaksbroad,union_ranges(H3K79me3peaksbroad,union_ranges(H3K79me1peaksbroad,
            union_ranges(H3K79me1peaks,union_ranges(H3K79me2peaks,union_ranges(H3K79me3peaks,
            union_ranges(H3K9me3peaksbroad,union_ranges(H3K27me3peaksbroad,
            union_ranges(H3K9me3peaks,union_ranges(H3K27me3peaks,
            union_ranges(H4K16acpeaksbroad,union_ranges(H3K9acpeaksbroad,union_ranges(H3K27acpeaksbroad,
            union_ranges(H4K16acpeaks,union_ranges(H3K9acpeaks,union_ranges(H3K27acpeaks,
            union_ranges(H3K36me3peaksbroad,union_ranges(H3K4me2peaksbroad,union_ranges(H3K4me3peaksbroad,
            union_ranges(H3K4me3peaks,union_ranges(H3K36me3peaks,H3K4me2peaks))))))))))))))))))))))))/genome_size*100  
        # pdf(paste0(plot_folder,"misc/peakWidths_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #     geom_density(data=as.data.frame(H3K9acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K27acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me2peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K4me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K79me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K9me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     geom_density(data=as.data.frame(H3K27me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
        #     scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
        #     scale_colour_manual(values=cols) +
        #     theme_classic() +
        #     theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()      
        # pdf(paste0(plot_folder,"misc/peakGeneOverlaps_",species,".pdf"),width=8,height=8)
        # print(ggplot() +
        #   geom_histogram(data=as.data.frame(H3K9acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K27acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H4K16acpeaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me2peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K4me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K36me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K79me1peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K9me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   geom_histogram(data=as.data.frame(H3K27me3peaks),aes(numGene,col=strsplit(V4,split = "_")[[1]][2]),binwidth = 1,fill="NA") +
        #   scale_x_continuous(name="Number of Genes") +
        #   scale_colour_manual(values=cols) +
        #   theme_classic() +
        #   theme(axis.text.x = element_text(size=18),
        #         axis.text.y = element_text(size=18),
        #         axis.title.x = element_text(size=18),
        #         axis.title.y = element_text(size=18),
        #         legend.text = element_text(size=18),
        #         legend.title = element_blank(),
        #         strip.text.x = element_text(size=18)))
        # dev.off()
        pdf(paste0(plot_folder,"extData_Fig5_draft/peakWidth_GeneOverlaps_",species,".pdf"),width=8,height=8)
        print(ggplot() +
          geom_density(data=as.data.frame(H4K16acpeaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K36me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K79me1peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K9me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=as.data.frame(H3K27me3peaks),aes(log10(size),col=strsplit(V4,split = "_")[[1]][2])) +
          geom_density(data=generegions,aes(log10(V11),col="H3")) +
          geom_text(aes(label=c("Mean","Max"),x=c(4.25,4.75),y=2),size=5) +
          geom_text(aes(label = c(round(mean(H4K16acpeaks$numGene),2),max(H4K16acpeaks$numGene)), x=c(4.25,4.75),y = 1.9,col="H4K16ac"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K36me3peaks$numGene),2),max(H3K36me3peaks$numGene)), x=c(4.25,4.75),y = 1.8,col="H3K36me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K79me1peaks$numGene),2),max(H3K79me1peaks$numGene)), x=c(4.25,4.75),y = 1.7,col="H3K79me1"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K27me3peaks$numGene),2),max(H3K27me3peaks$numGene)), x=c(4.25,4.75),y = 1.4,col="H3K27me3"),size=5,show.legend = FALSE) +
          geom_text(aes(label = c(round(mean(H3K9me3peaks$numGene),2),max(H3K9me3peaks$numGene)), x=c(4.25,4.75),y = 1.2,col="H3K9me3"),size=5,show.legend = FALSE) +
          scale_x_continuous(name = "log10(Peak Width)",limits = c(2,5)) +
          scale_y_continuous(name="Density") +
          scale_colour_manual(values=cols,labels=c("H3"="Genes")) +
          theme_classic() +
          theme(axis.text.x = element_text(size=18),
                axis.text.y = element_text(size=18),
                axis.title.x = element_text(size=18),
                axis.title.y = element_text(size=18),
                legend.text = element_text(size=18),
                legend.title = element_blank(),
                strip.text.x = element_text(size=18)))
        dev.off()
        totPeakSizes <- rbind(totPeakSizes,as.data.frame(c(H3K9acpeaks,H3K27acpeaks,H4K16acpeaks,H3K4me2peaks,H3K4me3peaks,H3K36me3peaks,H3K79me1peaks,H3K27me3peaks,H3K9me3peaks))%>%mutate(mark=str_split(V4,"_",simplify = T)[,2],species=str_split(V4,"_",simplify = T)[,1]))
    }
}
colnames(null_genome_perc) <- c("species","H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","Null")
pdf(paste0(plot_folder,"misc/nullMark_genome_perc.pdf"),width=12,height=8)
ggplot(melt(null_genome_perc),aes(x=factor(variable, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","Null")),y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))),fill=value)) +
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", mid = "grey", high="#000000",name="% Genome") +
    theme_classic() + # minimal theme
    theme(axis.text.x = element_text(size=18,angle = 45,hjust = 1),
          axis.text.y = element_text(size=18),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          strip.text.x = element_text(size=18))
dev.off()

##By species
for (name in c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3")){
    pdf(paste0(plot_folder,"extData_Fig5_draft/peakGeneWidths.",name,".pdf"),width=12,height=12)
    print(totPeakSizes%>%mutate(mark=gsub("H3K9me2","H3K9me3",mark))%>%select(species,mark,size)%>%bind_rows(totGeneRegions%>%mutate(mark="gene",size=V11)%>%select(mark,species,size))%>%
        filter(mark==name | mark=="gene")%>%
        ggplot(aes(x=log10(size),y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))),fill=factor(mark))) +
            geom_density_ridges_gradient() +
            scale_fill_manual(values=c("gene"="white",cols[name]),guide=NULL) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            coord_cartesian(clip = "off",xlim = c(2,5)) +
            theme_minimal(base_size = 14) + 
            theme(axis.text.y = element_text(vjust = 0),
                  axis.title.y = element_blank()))
    dev.off()
}
totPeakSizes%>%mutate(mark=gsub("H3K9me2","H3K9me3",mark))%>%select(species,mark,size)%>%bind_rows(totGeneRegions%>%mutate(mark="gene",size=V11)%>%select(mark,species,size))%>%
    ggplot(aes(x=log10(size),y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))))) +
        geom_density_ridges_gradient() +
        facet_wrap(~mark) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        coord_cartesian(clip = "off",xlim = c(2,5)) +
        theme_minimal(base_size = 14) + 
        theme(axis.text.y = element_text(vjust = 0),
              axis.title.y = element_blank())

totPeakSizes%>%filter(mark=="H3K27me3")%>%mutate(numGeneBin=cut(numGene,breaks=c(-1,0,1,5,100)))%>%
    ggplot(aes(x=log10(size),y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))),fill=factor(numGeneBin))) +
        geom_density_ridges_gradient() +
        scale_fill_viridis_d(alpha = 0.5,name="Number of genes",labels=c("(-1,0]"=0,"(0,1]"=1,"(1,5]"="2-5","(5,100]"=">5")) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(name="H3K27me3",expand = c(0, 0)) +
        coord_cartesian(clip = "off") +
        theme_minimal(base_size = 14) +
        theme(axis.text.y = element_text(vjust = 0))
###############################################################################################################################################################################################################################################################
##Correlation of marks overall and per replicate
##DONE
final <- read.table(file="/home/smontgomery/data/chip/metadata_final.tsv",header = T,sep = "\t")
for (sp in specieslist$species){
    ##Working directory
    setwd(paste0(chip_folder,sp,"_all/"))
    ifelse(!dir.exists(file.path(paste0(chip_folder,sp,"_all/correlation/"))),
        dir.create(file.path(paste0(chip_folder,sp,"_all/correlation/"))),
        "Directory Exists")
    ##Merged files
    coveragefile <- as_tibble(read.table(paste0(chip_folder,sp,"_all/bed/",sp,".coverage.RPGC.txt.gz"), header = T))
    colnames(coveragefile) <- gsub("H3K9me2","H3K9me3",colnames(coveragefile))
    for (mark in c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")){
        if(length(grep(mark,colnames(coveragefile)))==0){
            coveragefile[,paste0(sp,"_",mark)] <- 0
        }
    }
    ##Plot coverage correlation
    hclust_plot <- hclust(as.dist(1-(cor(coveragefile[,4:min(c(grep("^0$",colSums(coveragefile[,4:ncol(coveragefile)]))+2,ncol(coveragefile)))],method="pearson"))),method = "complete")
    par(mar=c(0.1,4.1, 4.1, 2.1))
    plot(hclust_plot,xlab = "",main = "",sub = "",labels=FALSE)
    tmpplot <- recordPlot()
    pdf(paste0("correlation/",sp,".coverageCorrelation.pearson.pdf"),width=15,height=15)
    print(plot_grid(tmpplot,
              ggplot(melt(as.matrix(cor(coveragefile[,4:min(c(grep("^0$",colSums(coveragefile[,4:ncol(coveragefile)]))+2,ncol(coveragefile)))],method="pearson"))),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]),y=factor(Var1, levels=hclust_plot$labels[hclust_plot$order]), fill=value)) +
                geom_tile(color = "white")+
                scale_fill_gradient2(low = "steelblue", mid = "white", high = "red3", 
                                    midpoint = 0, limit = c(-1,1),
                                     space = "Lab", 
                                     name="Pearson\nCorrelation") +
                geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 5) +
                theme_classic()+ # minimal theme
                scale_x_discrete(position = "bottom") +
                labs(x = NULL, y = NULL) +
                theme(axis.text.x = element_text(size=18,angle = 45,hjust = 1),
                    axis.text.y = element_text(size=18),
                    plot.margin=unit(c(0,15,0,0),"mm")),
        ncol = 1,rel_heights = c(1,4),axis = "r"))
    dev.off()
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    hclust_plot <- hclust(as.dist(1-(cor(coveragefile[,4:min(c(grep("^0$",colSums(coveragefile[,4:ncol(coveragefile)]))+2,ncol(coveragefile)))],method="spearman"))),method = "complete")
    par(mar=c(0.1,4.1, 4.1, 2.1))
    plot(hclust_plot,xlab = "",main = "",sub = "",labels=FALSE)
    tmpplot <- recordPlot()
    pdf(paste0("correlation/",sp,".coverageCorrelation.spearman.pdf"),width=15,height=15)
    print(plot_grid(tmpplot,
              ggplot(melt(as.matrix(cor(coveragefile[,4:min(c(grep("^0$",colSums(coveragefile[,4:ncol(coveragefile)]))+2,ncol(coveragefile)))],method="spearman"))),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]),y=factor(Var1, levels=hclust_plot$labels[hclust_plot$order]), fill=value)) +
                geom_tile(color = "white")+
                scale_fill_gradient2(low = "steelblue", mid = "white", high = "red3", 
                                    midpoint = 0, limit = c(-1,1),
                                     space = "Lab", 
                                     name="Spearman\nCorrelation") +
                geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 5) +
                theme_classic()+ # minimal theme
                scale_x_discrete(position = "bottom") +
                labs(x = NULL, y = NULL) +
                theme(axis.text.x = element_text(size=18,angle = 45,hjust = 1),
                    axis.text.y = element_text(size=18),
                    plot.margin=unit(c(0,15,0,0),"mm")),
        ncol = 1,rel_heights = c(1,4),axis = "r"))
    dev.off()
    par(mar=c(5.1, 4.1, 4.1, 2.1))

    ##Per replicate
    coveragefile <- as_tibble(read.table(paste0("/home/smontgomery/no_backup/smontgomery/chip/",sp,"_all/bed_files/",sp,".coverage.200bp.txt.gz"), header = T,sep = "\t",check.names = FALSE)) %>% replace(is.na(.), 0)
    gencov <- as.data.frame(mcols(GRanges(coveragefile)),optional=TRUE)
    ##Pearson
    gencor <- cor(gencov[,subset(final,species==sp & final=="Y")$sampleName],method="pearson")
    colnames(gencor) <- unite(subset(final,species==sp & final=="Y")[c("species","antibody")],newcol)%>%mutate(newcol=gsub("_NA","_input",newcol))%>%mutate(counter=rep(1,nrow(.)))%>%group_by(newcol)%>%mutate(counter2=cumsum(counter))%>%select(-counter)%>%unite("newname",newcol:counter2,sep = "_rep")%>%pull(newname)
    rownames(gencor) <- unite(subset(final,species==sp & final=="Y")[c("species","antibody")],newcol)%>%mutate(newcol=gsub("_NA","_input",newcol))%>%mutate(counter=rep(1,nrow(.)))%>%group_by(newcol)%>%mutate(counter2=cumsum(counter))%>%select(-counter)%>%unite("newname",newcol:counter2,sep = "_rep")%>%pull(newname)
    ##Plot coverage correlation
    hclust_plot <- hclust(as.dist(1-gencor),method = "complete")
    par(mar=c(0.1,4.1, 4.1, 2.1))
    plot(hclust_plot,xlab = "",main = "",sub = "",labels=FALSE)
    tmpplot <- recordPlot()
    pdf(paste0("correlation/",sp,".coverageCorrelation.reps.pearson.renamed.pdf"),width=15,height=15)
    print(plot_grid(tmpplot,
              ggplot(melt(as.matrix(gencor)),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]),y=factor(Var1, levels=hclust_plot$labels[hclust_plot$order]), fill=value)) +
                geom_tile(color = "white")+
                scale_fill_gradient2(low = "steelblue", mid = "white", high = "red3", 
                                    midpoint = 0, limit = c(-1,1),
                                     space = "Lab", 
                                     name="Pearson\nCorrelation") +
                theme_classic()+ # minimal theme
                scale_x_discrete(position = "bottom") +
                labs(x = NULL, y = NULL) +
                theme(axis.text.x = element_text(size=18,angle = 45,hjust = 1),
                    axis.text.y = element_text(size=18),
                    plot.margin=unit(c(0,15,0,0),"mm")),
        ncol = 1,rel_heights = c(1,4),axis = "r"))
    dev.off()
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    ##Spearman
    gencor <- cor(gencov[,subset(final,species==sp & final=="Y")$sampleName],method="spearman")
    colnames(gencor) <- unite(subset(final,species==sp & final=="Y")[c("species","antibody")],newcol)%>%mutate(newcol=gsub("_NA","_input",newcol))%>%mutate(counter=rep(1,nrow(.)))%>%group_by(newcol)%>%mutate(counter2=cumsum(counter))%>%select(-counter)%>%unite("newname",newcol:counter2,sep = "_rep")%>%pull(newname)
    rownames(gencor) <- unite(subset(final,species==sp & final=="Y")[c("species","antibody")],newcol)%>%mutate(newcol=gsub("_NA","_input",newcol))%>%mutate(counter=rep(1,nrow(.)))%>%group_by(newcol)%>%mutate(counter2=cumsum(counter))%>%select(-counter)%>%unite("newname",newcol:counter2,sep = "_rep")%>%pull(newname)
    hclust_plot <- hclust(as.dist(1-gencor),method = "complete")
    par(mar=c(0.1,4.1, 4.1, 2.1))
    plot(hclust_plot,xlab = "",main = "",sub = "",labels=FALSE)
    tmpplot <- recordPlot()
    pdf(paste0("correlation/",sp,".coverageCorrelation.reps.spearman.renamed.pdf"),width=15,height=15)
    print(plot_grid(tmpplot,
              ggplot(melt(as.matrix(gencor)),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]),y=factor(Var1, levels=hclust_plot$labels[hclust_plot$order]), fill=value)) +
                geom_tile(color = "white")+
                scale_fill_gradient2(low = "steelblue", mid = "white", high = "red3", 
                                    midpoint = 0, limit = c(-1,1),
                                     space = "Lab", 
                                     name="Spearman\nCorrelation") +
                theme_classic()+ # minimal theme
                scale_x_discrete(position = "bottom") +
                labs(x = NULL, y = NULL) +
                theme(axis.text.x = element_text(size=18,angle = 45,hjust = 1),
                    axis.text.y = element_text(size=18),
                    plot.margin=unit(c(0,15,0,0),"mm")),
        ncol = 1,rel_heights = c(1,4),axis = "r"))
    dev.off()
    par(mar=c(5.1, 4.1, 4.1, 2.1))
}

###############################################################################################################################################################################################################################################################
##Make boxplots of log2 fold enrichment of mark over H3 per locus
k79cor <- data.frame()
for (species in specieslist$species){
    grange = get_our_grange(species=species)
    tes=grange[mcols(grange)$source == "EDTA"]
    if (species == "Acas"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Atha"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Bnat"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k5.reorder.txt"))
    } else if (species == "Cfra"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k5.txt"))
    } else if (species == "Ddis"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k5.reorder.txt"))
    } else if (species == "Gthe"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k8.reorder.txt"))
    } else if (species == "Ngru"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k2.txt"))
    } else if (species == "Nvec"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Ppat"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Scer"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k2.txt"))
    } else if (species == "Spun"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k3.reorder.txt"))
    } else if (species == "Tthe"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k1.txt"))
    }
    teregions$order <- row.names(teregions)
    teregions <- bed_intersect(teregions%>%rename(chrom=V1,start=V2,end=V3)%>%mutate(chrom=as.character(chrom)),
                                tes%>%as.data.frame()%>%mutate(start=start-1)%>%rename(chrom=seqnames))%>%
                    filter(width.y == .overlap & gsub("_r.*","",V4.x) == transcript_id.y)%>%
                    select(chrom,start.x,end.x,transcript_id.y,V5.x,V6.x,V7.x,V8.x,V9.x,V10.x,V11.x,V12.x,V13.x,order.x,class_id.y,order_id.y,family_id.y,type.y)%>%
                    mutate(V13.x=paste0("TE_",V13.x))
                    
    colnames(teregions) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","order","class_id","order_id","family_id","type")
    tmpdata=teregions%>%count(V13, order_id)%>%group_by(V13)%>%mutate(pct= prop.table(n) * 100)%>%summarize(total=sum(n))
    tes$order <- as.numeric(teregions[match(tes$transcript_id,teregions$V4),]$order)
    techisq=chi_sq_test_loop(teregions = teregions)
    write.table(techisq,file=paste0(plot_folder,"extData_Fig6_draft/",species,".chisqtest.tsv"),quote = FALSE,sep = "\t",row.names = TRUE,col.names = NA)
    pdf(paste0(plot_folder,"extData_Fig6_draft/TEannotations.stacked.",species,".pdf"),width=20,height=10)
    print(teregions%>%count(V13, order_id)%>%group_by(V13)%>%mutate(pct= prop.table(n) * 100)%>%ggplot(aes(x=gsub("TE_cluster_","TE Cluster ",V13),y=pct,fill=factor(order_id,levels=c("MITE","Helitron","TIR","TE_DNA_intact","TE_DNA_fragment","DNA","SINE","TE_RT_SINE","LINE","TE_RT_LINE","PLE","DIRS","DNAauto","DNAnona","TE_RT_LTR_intact","TE_RT_LTR_fragment","LTR","RT")))) +
              geom_bar(stat = "identity") +
              geom_text(data=tmpdata,aes(x=gsub("TE_cluster_","TE Cluster ",V13),label = total, y=101,fill=NULL)) +
              scale_fill_manual(values = TEcols,labels=c("DNA"="DNA","TIR"="TIR","Helitron"="Helitron","TE_DNA_fragment"="DNA fragment","TE_DNA_intact"="DNA intact","DNAauto"="DNA auto","DNAnona"="DNA nona","RT"="RT","LTR"="LTR","TE_RT_LTR_fragment"="LTR fragment","TE_RT_LTR_intact"="LTR intact","TE_RT_LINE"="LINE","TE_RT_SINE"="SINE","PLE"="Penelope","DIRS"="DIRS","MITE"="MITE","TE_RT_DIRS"="NA","NA"="NA")) +
              theme_classic() +
              theme(axis.text.x = element_text(color = "black",size = 12,hjust = 1,angle = 45),
                    axis.text.y = element_text(color = "black",size = 12),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.title = element_blank(),
                    legend.text = element_text(colour = "black",size = 12),
                    strip.text = element_text(color="black",size=12),
                    strip.background = element_rect(linetype = "blank")))
    dev.off()
    pdf(paste0(plot_folder,"extData_Fig6_draft/TEannotations.",species,".RPGC.pdf"))
        print(ggplot(teregions,aes(y=rev(order))) +
            geom_tile(aes(x=1,fill=factor(class_id,levels = c("RT","DNA","NA")))) +
            geom_tile(aes(x=2,fill=factor(order_id,levels = c("LTR","LINE","TE_RT_LINE","SINE","TE_RT_SINE","PLE","DIRS","TE_RT_DIRS","TIR","Helitron","DNAauto","DNAnona","MITE","NA")))) +
            geom_tile(aes(x=3,fill=factor(type,levels = c("TE_RT_LTR_intact","TE_RT_LTR_fragment","LINE","TE_RT_LINE","SINE","TE_RT_SINE","PLE","DIRS","TE_RT_DIRS","TE_DNA_intact","TE_DNA_fragment","NA")))) +
            facet_grid(rows = vars(V13), scales = "free", space = "free") +
            scale_fill_manual(values = TEcols,labels=c("DNA"="DNA","TIR"="TIR","Helitron"="Helitron","TE_DNA_fragment"="DNA fragment","TE_DNA_intact"="DNA intact","DNAauto"="DNA auto","DNAnona"="DNA nona","RT"="RT","LTR"="LTR","TE_RT_LTR_fragment"="LTR fragment","TE_RT_LTR_intact"="LTR intact","TE_RT_LINE"="NA","TE_RT_SINE"="NA","PLE"="Penelope","MITE"="MITE","TE_RT_DIRS"="NA","NA"="NA")) +
            theme_classic() +
            theme(legend.title = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.line = element_blank(),
                  axis.ticks = element_blank()))
    dev.off()
    if (species == "Atha" | species == "Ppat" | species == "Nvec"){
        png(filename=paste0(plot_folder,"extData_Fig6_draft/TEannotations.",species,".RPGC.png"), type="cairo",width = 1024, height = 1024, units = "px")
        print(ggplot(teregions,aes(y=rev(order))) +
          geom_tile(aes(x=1,fill=factor(class_id,levels = c("RT","DNA","NA")))) +
          geom_tile(aes(x=2,fill=factor(order_id,levels = c("LTR","LINE","TE_RT_LINE","SINE","TE_RT_SINE","PLE","DIRS","TE_RT_DIRS","TIR","Helitron","DNAauto","DNAnona","MITE","NA")))) +
          geom_tile(aes(x=3,fill=factor(type,levels = c("TE_RT_LTR_intact","TE_RT_LTR_fragment","LINE","TE_RT_LINE","SINE","TE_RT_SINE","PLE","DIRS","TE_RT_DIRS","TE_DNA_intact","TE_DNA_fragment","NA")))) +
          facet_grid(rows = vars(V13), scales = "free", space = "free") +
          scale_fill_manual(values = TEcols,labels=c("DNA"="DNA","TIR"="TIR","Helitron"="Helitron","TE_DNA_fragment"="DNA fragment","TE_DNA_intact"="DNA intact","DNAauto"="DNA auto","DNAnona"="DNA nona","RT"="RT","LTR"="LTR","TE_RT_LTR_fragment"="LTR fragment","TE_RT_LTR_intact"="LTR intact","TE_RT_LINE"="NA","TE_RT_SINE"="NA","PLE"="Penelope","MITE"="MITE","TE_RT_DIRS"="NA","NA"="NA"),guide="none") +
          theme_classic() +
          theme(legend.title = element_blank(),
                strip.background = element_blank(),
                strip.text.y = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank()))
        dev.off()
    }
    if (species == "Atha" | species == "Ppat" | species == "Nvec" | species == "Cfra" | species == "Spun" | species == "Acas" | species == "Ngru" | species == "Bnat" | species == "Gthe" | species == "Ddis"){
        tmptable <- as.matrix(table(teregions%>%select(V13,order_id)))
        pdf(paste0(plot_folder,"extData_Fig7_draft/TEenrichment.",species,".order.pdf"))
        print(mosaic(tmptable[,colnames(tmptable)[match(c("LTR","LINE","TE_RT_LINE","SINE","TE_RT_SINE","PLE","DIRS","TE_RT_DIRS","TIR","Helitron","DNAauto","DNAnona","MITE","NA"),colnames(tmptable))][!is.na(match(c("LTR","LINE","TE_RT_LINE","SINE","TE_RT_SINE","PLE","DIRS","TE_RT_DIRS","TIR","Helitron","DNAauto","DNAnona","MITE","NA"),colnames(tmptable)))]],
                     shade = TRUE,
                     labeling=labeling_border(labels=TRUE, rot_labels = c(90,0),varnames=FALSE,
                                              just_labels = c("center","right"),
                                              gp_labels = gpar(fontsize = 20)),
                     legend = legend_resbased(fontsize = 12),
                     gp = shading_Friendly2(h=c(0,260),lty=1:1),
                     margins = c(0,5,0,5)))
        dev.off()
    }
    ##Find top/bottom expressed gene cluster
    if (species == "Gthe"){
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.reorder.txt"))
    } else {
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.reorder.txt"))
    }
    RNAcountsfile <- read.table(paste0("/home/smontgomery/cluster/smontgomery/genomes/",species,"/",species,"_expr_sorted.tsv"),header=TRUE)
    if (species == "Cfra"){
        RNAcountsfile$gene_id <- gsub("Cfra_CFRG0","Cfra_CFRCFRG0",RNAcountsfile$gene_id)
    }
    generegions$V13 <- paste0("gene_",generegions$V13)
    generegions$RNAcounts <- RNAcountsfile$RNAcounts[match(generegions$V4,RNAcountsfile$gene_id)]
    ##Boxplots of UMI and length per gene cluster
    pdf(paste0(plot_folder,"Fig3_draft/RNAcounts.",species,".noregion.boxplot.pdf"),width=15,height=15)
    print(ggplot(subset(generegions,V13!="NA"),aes(y=factor(V13,levels = rev(unique(V13))),x=log10(RNAcounts + 1))) +
          geom_boxplot(outliers = FALSE) +
          theme_classic() +
          xlab(label = "log10(Counts + 1)") +
          scale_x_continuous(limits = c(0,6),breaks=c(0,2,4,6)) +
          theme(axis.text.x = element_text(size=18),
                axis.text.y = element_blank(),
                axis.title.x = element_text(size=18),
                axis.title.y = element_blank(),
                legend.text = element_text(size=18),
                legend.title = element_blank()))
    dev.off()
    pdf(paste0(plot_folder,"Fig3_draft/genelength.",species,".noregion.boxplot.pdf"),width=15,height=15)
    print(ggplot(subset(generegions,V13!="NA"),aes(y=factor(V13,levels = rev(unique(V13))),x=log10(V3-V2))) +
          geom_boxplot(outliers = FALSE) +
          xlab(label = "log10(Gene length)") +
          scale_x_continuous(limits = c(2,5),breaks=c(2,3,4)) +
          theme_classic() +
          theme(axis.text.x = element_text(size=18),
                axis.text.y = element_blank(),
                axis.title.x = element_text(size=18),
                axis.title.y = element_blank(),
                legend.text = element_text(size=18),
                legend.title = element_blank()))
    dev.off()
    allregions <- rbind(generegions[,1:13],teregions[,1:13])
    remove(RNAcountsfile,teregions,tes)
    ##Start processing
    grange_w_bins=add_bins(grange)
    grange_w_promoter=add_promoter(grange_w_bins)
    grange_wo_overlaps=merge_overlap_feature(grange_w_promoter)
    grange_w_intergenic=add_intergenic(grange_wo_overlaps)
    grange_w_expression=add_expression(grange_w_intergenic)
    feature=grange_to_feature(grange_w_expression)
    coveragefile <- as_tibble(read.table(paste0(chip_folder,species,"_all/bed/",species,".coverage.RPGC.txt.gz"), header = T))
    colnames(coveragefile) <- gsub("H3K9me2","H3K9me3",colnames(coveragefile))
    for (mark in c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")){
        if(length(grep(mark,colnames(coveragefile)))==0){
            coveragefile[,paste0(species,"_",mark)] <- 0
        }
    }
    fea_cov=bed_intersect(coveragefile,feature)%>%mutate(start=start.x,end=end.x)
    remove(coveragefile,feature,grange_w_bins,grange_w_promoter,grange_wo_overlaps,grange_w_intergenic,grange_w_expression)
    crm=c("Pt","Mt","Mito")
    H3peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(H3peaks)[1:3] <- c("chrom","start","end")
    H3peaks <- c(GRanges(H3peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3peaks=H3peaks[!(seqnames(H3peaks) %in% c(crm))]
    inputpeaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_input_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(inputpeaks)[1:3] <- c("chrom","start","end")
    inputpeaks <- c(GRanges(inputpeaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    inputpeaks=inputpeaks[!(seqnames(inputpeaks) %in% c(crm))]
    ctrlpeaks <- as_tibble(union_ranges(H3peaks,inputpeaks))%>%mutate(start=start-1)
    remove(H3peaks,inputpeaks)
    colnames(ctrlpeaks)[1] <- c("chrom")
    red_fea_cov <- bed_subtract(fea_cov,ctrlpeaks)%>%filter(chrom != "chrCH709182" | species != "Ddis")
    remove(ctrlpeaks,fea_cov)
    markCoverage_red=red_fea_cov%>%
        select(geneId.y,gene_id.y,feature.y,order_id.y,.overlap,!!as.name(paste0(quo_name(species),'_H3K27me3.x')),!!as.name(paste0(quo_name(species),'_H3K9me1.x')),!!as.name(paste0(quo_name(species),'_H3K9me3.x')),
               !!as.name(paste0(quo_name(species),'_H3K79me1.x')),!!as.name(paste0(quo_name(species),'_H3K79me2.x')),!!as.name(paste0(quo_name(species),'_H3K79me3.x'))
               ,!!as.name(paste0(quo_name(species),'_H3K9ac.x')),!!as.name(paste0(quo_name(species),'_H3K27ac.x')),!!as.name(paste0(quo_name(species),'_H4K16ac.x'))
               ,!!as.name(paste0(quo_name(species),'_H3K4me2.x')),!!as.name(paste0(quo_name(species),'_H3K4me3.x')),!!as.name(paste0(quo_name(species),'_H3K36me3.x'))
               ,!!as.name(paste0(quo_name(species),'_H3.x')),!!as.name(paste0(quo_name(species),'_input.x')))%>%
        mutate(feature.y=gsub("Bin 1 Low","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 2 Low","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 3 Low","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 4 Low","Gene",feature.y))%>%
        mutate(feature.y=gsub("Gene body Low","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 1 Med","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 2 Med","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 3 Med","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 4 Med","Gene",feature.y))%>%
        mutate(feature.y=gsub("Gene body Med","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 1 High","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 2 High","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 3 High","Gene",feature.y))%>%
        mutate(feature.y=gsub("Bin 4 High","Gene",feature.y))%>%
        mutate(feature.y=gsub("Gene body High","Gene",feature.y))%>%
        mutate(feature.y=gsub("TE RT LTR fragment","TE RT",feature.y))%>%
        mutate(feature.y=gsub("TE RT LINE","TE RT",feature.y))%>%
        mutate(feature.y=gsub("TE RT SINE","TE RT",feature.y))%>%
        mutate(feature.y=gsub("TE RT DIRS","TE RT",feature.y))%>%
        mutate(feature.y=gsub("TE DNA fragment","TE DNA",feature.y))%>%
        mutate(quanH3K27me3=(log2(!!as.name(paste0(quo_name(species),'_H3K27me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K27me3.x")]))))%>%
        mutate(quanH3K9me1=(log2(!!as.name(paste0(quo_name(species),'_H3K9me1.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K9me1.x")]))))%>%
        mutate(quanH3K9me3=(log2(!!as.name(paste0(quo_name(species),'_H3K9me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K9me3.x")]))))%>%
        mutate(quanH3K79me1=(log2(!!as.name(paste0(quo_name(species),'_H3K79me1.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K79me1.x")]))))%>%
        mutate(quanH3K79me2=(log2(!!as.name(paste0(quo_name(species),'_H3K79me2.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K79me2.x")]))))%>%
        mutate(quanH3K79me3=(log2(!!as.name(paste0(quo_name(species),'_H3K79me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K79me3.x")]))))%>%
        mutate(quanH3K9ac=(log2(!!as.name(paste0(quo_name(species),'_H3K9ac.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K9ac.x")]))))%>%
        mutate(quanH3K27ac=(log2(!!as.name(paste0(quo_name(species),'_H3K27ac.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K27ac.x")]))))%>%
        mutate(quanH4K16ac=(log2(!!as.name(paste0(quo_name(species),'_H4K16ac.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H4K16ac.x")]))))%>%
        mutate(quanH3K4me2=(log2(!!as.name(paste0(quo_name(species),'_H3K4me2.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K4me2.x")]))))%>%
        mutate(quanH3K4me3=(log2(!!as.name(paste0(quo_name(species),'_H3K4me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K4me3.x")]))))%>%
        mutate(quanH3K36me3=(log2(!!as.name(paste0(quo_name(species),'_H3K36me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K36me3.x")]))))%>%
        mutate(quanH3=(log2(!!as.name(paste0(quo_name(species),'_H3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3.x")]))))%>%
        mutate(quaninput=(log2(!!as.name(paste0(quo_name(species),'_input.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_input.x")]))))%>%
        group_by(geneId.y,gene_id.y,feature.y,order_id.y,.drop=TRUE)%>%
        summarize(H3K27me3=mean(quanH3K27me3[!is.infinite(quanH3K27me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K9me1=mean(quanH3K9me1[!is.infinite(quanH3K9me1)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K9me3=mean(quanH3K9me3[!is.infinite(quanH3K9me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K79me1=mean(quanH3K79me1[!is.infinite(quanH3K79me1)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K79me2=mean(quanH3K79me2[!is.infinite(quanH3K79me2)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K79me3=mean(quanH3K79me3[!is.infinite(quanH3K79me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K9ac=mean(quanH3K9ac[!is.infinite(quanH3K9ac)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K27ac=mean(quanH3K27ac[!is.infinite(quanH3K27ac)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H4K16ac=mean(quanH4K16ac[!is.infinite(quanH4K16ac)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K4me2=mean(quanH3K4me2[!is.infinite(quanH3K4me2)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K4me3=mean(quanH3K4me3[!is.infinite(quanH3K4me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3K36me3=mean(quanH3K36me3[!is.infinite(quanH3K36me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
                  H3=mean(quanH3[!is.infinite(quanH3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE))%>%
        ungroup()%>%filter(rowSums(select(., 5:17),na.rm=TRUE) != 0)%>%pivot_longer(cols = 5:17)%>%filter(!is.na(value))%>%
        mutate(genecluster=if_else(feature.y=="Gene",allregions$V13[match(gene_id.y, allregions$V4)],feature.y))%>%
        mutate(tecluster=if_else((feature.y=="TE RT" | feature.y=="TE RT LTR intact" | feature.y=="TE DNA intact" | feature.y=="TE DNA"),allregions$V13[match(geneId.y, allregions$V4)],feature.y))%>%
        mutate(genecluster=gsub("gene_cluster_","Gene ",genecluster))%>%
        mutate(RNAcounts=if_else(feature.y=="Gene",generegions$RNAcounts[match(gene_id.y, generegions$V4)],NA))%>%
        mutate(size=if_else(feature.y=="Gene",grange$size[match(gene_id.y, grange$gene_id)],NA))
    teCoverage_red=markCoverage_red%>%
        filter(!is.na(order_id.y) & !is.na(tecluster) & name != "H3")%>%
        group_by(name,tecluster,.drop=TRUE)%>%
        summarize(medVal=median(value,na.rm=TRUE))
    # pdf(paste0(plot_folder,"extData_Fig8_draft/TEmedEnrich.heatmap.",species,".pdf"))
    # print(ggplot(teCoverage_red,aes(x=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")),y=factor(tecluster,levels=rev(c("TE_cluster_1","TE_cluster_2","TE_cluster_3","TE_cluster_4","TE_cluster_5","TE_cluster_6","TE_cluster_7","TE_cluster_8"))),fill=medVal)) +
    #     geom_tile(color = "white")+
    #     scale_fill_gradient2(low = "steelblue", high = "red3", mid = "white",midpoint = 0, space = "Lab",name="Coverage",na.value = "grey",limits=c(-6,6)) +
    #     theme_classic()+ # minimal theme
    #     labs(x = NULL, y = NULL) +
    #     theme(axis.text.x = element_text(size=18,hjust=1,angle=45),
    #           axis.text.y = element_text(size=18),
    #           axis.title.x = element_blank(),
    #           axis.title.y = element_blank(),
    #           legend.text = element_text(size=18),
    #           legend.title = element_text(size=18),
    #           strip.text.x = element_text(size=18)))
    # dev.off()
    # pdf(paste0(plot_folder,"extData_Fig8_draft/boxplot.TEs.",species,".RPGC.pdf"),width=20,height=10)
    # print(subset(markCoverage_red, (tecluster=="TE_cluster_1" | tecluster=="TE_cluster_2" | tecluster=="TE_cluster_3" | tecluster=="TE_cluster_4" | tecluster=="TE_cluster_5" | tecluster=="TE_cluster_6" | tecluster=="TE_cluster_7" | tecluster=="TE_cluster_8") & 
    #              (name=="H3K27me3" | name=="H3K9me1" | name=="H3K9me3" |name=="H3K79me1" | name=="H3K79me2" | name=="H3K79me3")) %>%
    #       ggplot(aes(x=tecluster, y=value, color=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
    #           geom_hline(yintercept=0) +
    #           ggrastr::rasterise(geom_sina(alpha=0.2,size=2)) +
    #           scale_y_continuous(name = "log2(Input mean normalized coverage)",limits=c(-6.5,6.5)) +
    #           scale_colour_manual(name="Mark",values=cols) +
    #           theme_minimal() +
    #           theme(axis.text.x = element_text(color = "black",size = 12,hjust = 1,angle = 45),
    #                 axis.text.y = element_text(color = "black",size = 12),
    #                 axis.title.x = element_blank(),
    #                 axis.title.y = element_text(color = "black",size = 12),
    #                 legend.title = element_text(color = "black",size = 12),
    #                 legend.text = element_text(colour = "black",size = 12)))
    # dev.off()
    tmpcor <- NA
    tryCatch({
        tmpcor <- cor(subset(markCoverage_red,!is.na(RNAcounts) & name=="H3K4me3")$value,subset(markCoverage_red,!is.na(RNAcounts) & name=="H3K4me3")$RNAcounts,use = "pairwise.complete.obs",method = "spearman")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    k79cor[species,1] <- if_else(nrow(subset(markCoverage_red,name=="H3K4me3")) > 0, tmpcor,NA)
    tryCatch({
        tmpcor <- cor(subset(markCoverage_red,!is.na(RNAcounts) & name=="H3K79me1")$value,subset(markCoverage_red,!is.na(RNAcounts) & name=="H3K79me1")$RNAcounts,use = "pairwise.complete.obs",method = "spearman")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    k79cor[species,2] <- if_else(nrow(subset(markCoverage_red,name=="H3K79me1")) > 0, tmpcor,NA)
    tryCatch({
        tmpcor <- cor(subset(markCoverage_red,!is.na(RNAcounts) & name=="H3K79me2")$value,subset(markCoverage_red,!is.na(RNAcounts) & name=="H3K79me2")$RNAcounts,use = "pairwise.complete.obs",method = "spearman")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    k79cor[species,3] <- if_else(nrow(subset(markCoverage_red,name=="H3K79me2")) > 0, tmpcor,NA)
    tryCatch({
        tmpcor <- cor(subset(markCoverage_red,!is.na(RNAcounts) & name=="H3K79me3")$value,subset(markCoverage_red,!is.na(RNAcounts) & name=="H3K79me3")$RNAcounts,use = "pairwise.complete.obs",method = "spearman")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    k79cor[species,4] <- if_else(nrow(subset(markCoverage_red,name=="H3K79me3")) > 0, tmpcor,NA)
    tryCatch({
        k79cor[species,5] <- species
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    nam <- paste0("marCoverage_red_",species)
    assign(nam, markCoverage_red)
}
colnames(k79cor) <- c("H3K4me3","H3K79me1","H3K79me2","H3K79me3","species")
##Scatterplots of RNAcounts vs mark
pdf(paste0(plot_folder,"Fig3_draft/boxplot_quantile.Acas.H3K79me.RPGC.pdf"),width=5,height=10)
print(subset(marCoverage_red_Acas,!is.na(RNAcounts) & (name=="H3K79me1" | name=="H3K79me2" | name=="H3K79me3"))%>%mutate(quantile=as.factor(as.integer(cut(RNAcounts, quantile(RNAcounts, probs=0:5/5), include.lowest=TRUE))))%>%
          ggplot(aes(x=quantile,y=value,fill=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
          geom_hline(yintercept = 0) +
          geom_boxplot(outliers = FALSE) +
          facet_wrap(vars(name),nrow=3) +
          scale_fill_manual(name="Mark",values=cols) +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black",size = 12),
                axis.text.y = element_text(color = "black",size = 12),
                axis.title.x = element_text(color = "black",size = 12),
                axis.title.y = element_text(color = "black",size = 12),
                legend.title = element_text(color = "black",size = 12),
                legend.text = element_text(colour = "black",size = 12)))
dev.off()
pdf(paste0(plot_folder,"Fig3_draft/boxplot_quantile.Spun.H3K79me.RPGC.pdf"),width=5,height=10)
print(subset(marCoverage_red_Spun,!is.na(RNAcounts) & (name=="H3K79me1" | name=="H3K79me2" | name=="H3K79me3"))%>%mutate(quantile=as.factor(as.integer(cut(RNAcounts, quantile(RNAcounts, probs=0:5/5), include.lowest=TRUE))))%>%
          ggplot(aes(x=quantile,y=value,fill=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
          geom_hline(yintercept = 0) +
          geom_boxplot(outliers = FALSE) +
          facet_wrap(vars(name),nrow=3) +
          scale_fill_manual(name="Mark",values=cols) +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black",size = 12),
                axis.text.y = element_text(color = "black",size = 12),
                axis.title.x = element_text(color = "black",size = 12),
                axis.title.y = element_text(color = "black",size = 12),
                legend.title = element_text(color = "black",size = 12),
                legend.text = element_text(colour = "black",size = 12)))
dev.off()
pdf(paste0(plot_folder,"Fig3_draft/boxplot_quantile.Ddis.H3K79me.RPGC.pdf"),width=5,height=10)
print(subset(marCoverage_red_Ddis,!is.na(RNAcounts) & (name=="H3K79me1" | name=="H3K79me2" | name=="H3K79me3"))%>%mutate(quantile=as.factor(as.integer(cut(RNAcounts, quantile(RNAcounts, probs=0:5/5), include.lowest=TRUE))))%>%
          ggplot(aes(x=quantile,y=value,fill=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
          geom_hline(yintercept = 0) +
          geom_boxplot(outliers = FALSE) +
          facet_wrap(vars(name),nrow=3) +
          scale_fill_manual(name="Mark",values=cols) +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black",size = 12),
                axis.text.y = element_text(color = "black",size = 12),
                axis.title.x = element_text(color = "black",size = 12),
                axis.title.y = element_text(color = "black",size = 12),
                legend.title = element_text(color = "black",size = 12),
                legend.text = element_text(colour = "black",size = 12)))
dev.off()
# pdf(paste0(plot_folder,"Fig3_draft/scatter.enrich_RNAcounts.Acas.H3K79me2.RPGC.pdf"),width=15,height=15)
# print(subset(marCoverage_red_Acas,!is.na(RNAcounts) & name=="H3K79me2") %>%
#           ggplot(aes(x=value, y=log2(RNAcounts), color=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
#               geom_hline(yintercept=0) +
#               geom_vline(xintercept=0) +
#               ggrastr::rasterise(geom_point(alpha=0.2,size=2)) +
#               scale_x_continuous(name = "log2(Input mean normalized coverage)") +
#               scale_y_continuous(name = "log2(RNAcounts)") +
#               stat_cor(method = "pearson",show.legend = FALSE) +
#               scale_colour_manual(name="Mark",values=cols) +
#               theme_minimal() +
#               theme(axis.text.x = element_text(color = "black",size = 12),
#                     axis.text.y = element_text(color = "black",size = 12),
#                     axis.title.x = element_text(color = "black",size = 12),
#                     axis.title.y = element_text(color = "black",size = 12),
#                     legend.title = element_text(color = "black",size = 12),
#                     legend.text = element_text(colour = "black",size = 12)))
# dev.off()
# pdf(paste0(plot_folder,"Fig3_draft/scatter.enrich_RNAcounts.Acas.H3K79me3.RPGC.pdf"),width=15,height=15)
# print(subset(marCoverage_red_Acas,!is.na(RNAcounts) & name=="H3K79me3") %>%
#           ggplot(aes(x=value, y=log2(RNAcounts), color=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
#               geom_hline(yintercept=0) +
#               geom_vline(xintercept=0) +
#               ggrastr::rasterise(geom_point(alpha=0.2,size=2)) +
#               scale_x_continuous(name = "log2(Input mean normalized coverage)") +
#               scale_y_continuous(name = "log2(RNAcounts)") +
#               stat_cor(method = "pearson",show.legend = FALSE) +
#               scale_colour_manual(name="Mark",values=cols) +
#               theme_minimal() +
#               theme(axis.text.x = element_text(color = "black",size = 12),
#                     axis.text.y = element_text(color = "black",size = 12),
#                     axis.title.x = element_text(color = "black",size = 12),
#                     axis.title.y = element_text(color = "black",size = 12),
#                     legend.title = element_text(color = "black",size = 12),
#                     legend.text = element_text(colour = "black",size = 12)))
# dev.off()
# pdf(paste0(plot_folder,"Fig3_draft/scatter.enrich_RNAcounts.Spun.H3K79me3.RPGC.pdf"),width=15,height=15)
# print(subset(marCoverage_red_Spun,!is.na(RNAcounts) & name=="H3K79me3") %>%
#           ggplot(aes(x=value, y=log2(RNAcounts), color=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
#               geom_hline(yintercept=0) +
#               geom_vline(xintercept=0) +
#               ggrastr::rasterise(geom_point(alpha=0.2,size=2)) +
#               scale_x_continuous(name = "log2(Input mean normalized coverage)") +
#               scale_y_continuous(name = "log2(RNAcounts)") +
#               stat_cor(method = "pearson",show.legend = FALSE) +
#               scale_colour_manual(name="Mark",values=cols) +
#               theme_minimal() +
#               theme(axis.text.x = element_text(color = "black",size = 12),
#                     axis.text.y = element_text(color = "black",size = 12),
#                     axis.title.x = element_text(color = "black",size = 12),
#                     axis.title.y = element_text(color = "black",size = 12),
#                     legend.title = element_text(color = "black",size = 12),
#                     legend.text = element_text(colour = "black",size = 12)))
# dev.off()
# pdf(paste0(plot_folder,"Fig3_draft/scatter.enrich_RNAcounts.Ddis.H3K79me3.RPGC.pdf"),width=15,height=15)
# print(subset(marCoverage_red_Ddis,!is.na(RNAcounts) & name=="H3K79me3") %>%
#           ggplot(aes(x=value, y=log2(RNAcounts), color=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
#               geom_hline(yintercept=0) +
#               geom_vline(xintercept=0) +
#               ggrastr::rasterise(geom_point(alpha=0.2,size=2)) +
#               scale_x_continuous(name = "log2(Input mean normalized coverage)") +
#               scale_y_continuous(name = "log2(RNAcounts)") +
#               stat_cor(method = "pearson",show.legend = FALSE) +
#               scale_colour_manual(name="Mark",values=cols) +
#               theme_minimal() +
#               theme(axis.text.x = element_text(color = "black",size = 12),
#                     axis.text.y = element_text(color = "black",size = 12),
#                     axis.title.x = element_text(color = "black",size = 12),
#                     axis.title.y = element_text(color = "black",size = 12),
#                     legend.title = element_text(color = "black",size = 12),
#                     legend.text = element_text(colour = "black",size = 12)))
# dev.off()
pdf(paste0(plot_folder,"Fig3_draft/H3K79me_correlation.RPGC.pdf"),width=15,height=15)
p <- ggplot(melt(k79cor),aes(x=variable,y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))),fill=value)) +
    geom_tile(color="white") +
    scale_fill_gradient2(low = "steelblue", high = "red3", mid="white",midpoint = 0,limits=c(-1,1),na.value = "white",name="Spearman\nCorrelation") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 12),
          axis.text.y = element_text(color = "black",size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(color = "black",size = 12),
          legend.text = element_text(colour = "black",size = 12))
p + geom_segment(
    aes(x=xmin,xend=xmax,y=ymin,yend=ymax), 
    subset(ggplot_build(p)$data[[1]],fill=="white"),
    inherit.aes=F)
dev.off()
##Boxplots of coverage per TE and gene cluster
repStates <- data.frame()
for (species in specieslist$species){
    tmptable <- get(paste0("marCoverage_red_",species))
    geneclusterlist <- c("Gene 1","Gene 2","Gene 3","Gene 4","Gene 5","Gene 6")
    if (species == "Ngru" | species == "Scer"){
        teclusterlist <- c("TE_cluster_1","TE_cluster_2")
    } else if (species == "Tthe") {
        teclusterlist <- c("TE_cluster_1")
    } else if (species == "Spun") {
        teclusterlist <- c("TE_cluster_1","TE_cluster_2","TE_cluster_3")
    } else if (species == "Bnat" | species == "Ddis" | species == "Cfra"){
        teclusterlist <- c("TE_cluster_1","TE_cluster_2","TE_cluster_3","TE_cluster_4","TE_cluster_5")
    } else if (species == "Atha" | species == "Nvec" | species == "Ppat" | species == "Acas"){
        teclusterlist <- c("TE_cluster_1","TE_cluster_2","TE_cluster_3","TE_cluster_4","TE_cluster_5","TE_cluster_6")
    } else if (species == "Gthe"){
        teclusterlist <- c("TE_cluster_1","TE_cluster_2","TE_cluster_3","TE_cluster_4","TE_cluster_5","TE_cluster_6","TE_cluster_7","TE_cluster_8")
        geneclusterlist <- c("Gene 1","Gene 2","Gene 3","Gene 4","Gene 5","Gene 6","Gene 7","Gene 8")
    }
    for (cluster in teclusterlist){
        for (mark in c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3")){
            repStates[paste0(species,mark,cluster),1] <- species
            repStates[paste0(species,mark,cluster),2] <- mark
            repStates[paste0(species,mark,cluster),3] <- cluster
            if (tmptable%>%filter(tecluster==cluster & name==mark)%>%nrow()==0){
                repStates[paste0(species,mark,cluster),4:7] <- c(0,0,-1,paste0(species,cluster))
            } else {
                repStates[paste0(species,mark,cluster),4:5] <- tmptable%>%filter(tecluster==cluster & name==mark)%>%group_by(tecluster,name)%>%reframe(med=median(value,na.rm=T),q1=quantile(value,0.25))%>%select(med,q1)
                repStates[paste0(species,mark,cluster),6] <- if_else(repStates[paste0(species,mark,cluster),4] >= 1,2,if_else(repStates[paste0(species,mark,cluster),4] < 1 & repStates[paste0(species,mark,cluster),5] > 0,1,if_else(repStates[paste0(species,mark,cluster),5] <= 0,0,-1)))
                repStates[paste0(species,mark,cluster),7] <- paste0(species,cluster)
            }
        }
        numrow <- max(nrow(subset(tmptable, (tecluster==cluster) & (name=="H3K27me3"))),nrow(subset(tmptable, (tecluster==cluster) & (name=="H3K79me1"))))
        pdf(paste0(plot_folder,"extData_Fig8_draft/all_boxplots/boxplot.TEs.",species,".",cluster,".pdf"),width=10,height=10)
        print(subset(tmptable, (tecluster==cluster) & 
                 (name=="H3K27me3" | name=="H3K9me1" | name=="H3K9me3" |name=="H3K79me1" | name=="H3K79me2" | name=="H3K79me3")) %>%
          ggplot(aes(x=tecluster, y=value, fill=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
              geom_hline(yintercept=0) +
              geom_boxplot(outliers = FALSE) +
              annotate("text", x = 0.6, y = 5.5, label = paste0("n=",numrow),size=8) +
              scale_y_continuous(name = "log2(Input mean normalized coverage)",limits=c(-6.5,6.5),breaks=c(-4,4)) +
              scale_fill_manual(name="Mark",values=cols) +
              theme_minimal() +
              theme(axis.text.x = element_text(color = "black",size = 12,hjust = 1,angle = 45),
                    axis.text.y = element_text(color = "black",size = 12),
                    axis.title.x = element_blank(),
                    axis.title.y = element_text(color = "black",size = 12),
                    legend.title = element_text(color = "black",size = 12),
                    legend.text = element_text(colour = "black",size = 12),
                    panel.grid.minor.y = element_blank(),
                    panel.grid.major.x = element_blank(),
                    panel.grid.major.y = element_line(linetype = "dashed",colour = "black")))
        dev.off()
    }
    for (cluster in geneclusterlist){
        for (mark in c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3")){
            repStates[paste0(species,mark,cluster),1] <- species
            repStates[paste0(species,mark,cluster),2] <- mark
            repStates[paste0(species,mark,cluster),3] <- cluster
            if (tmptable%>%filter(genecluster==cluster & name==mark)%>%nrow()==0){
                repStates[paste0(species,mark,cluster),4:7] <- c(0,0,-1,paste0(species,cluster))
            } else {
                repStates[paste0(species,mark,cluster),4:5] <- tmptable%>%filter(genecluster==cluster & name==mark)%>%group_by(genecluster,name)%>%reframe(med=median(value,na.rm=T),q1=quantile(value,0.25))%>%select(med,q1)
                repStates[paste0(species,mark,cluster),6] <- if_else(repStates[paste0(species,mark,cluster),4] >= 1,2,if_else(repStates[paste0(species,mark,cluster),4] < 1 & repStates[paste0(species,mark,cluster),5] > 0,1,if_else(repStates[paste0(species,mark,cluster),5] <= 0,0,-1)))
                repStates[paste0(species,mark,cluster),7] <- paste0(species,cluster)
            }
        }
        numrow <- max(nrow(subset(tmptable, (genecluster==cluster) & (name=="H3K27me3"))),nrow(subset(tmptable, (genecluster==cluster) & (name=="H3K79me1"))))
        pdf(paste0(plot_folder,"extData_Fig8_draft/all_boxplots/boxplot.gene.",species,".",cluster,".pdf"),width=10,height=10)
        print(subset(tmptable, (genecluster==cluster) & 
                     (name=="H3K27me3" | name=="H3K9me1" | name=="H3K9me3" |name=="H3K79me1" | name=="H3K79me2" | name=="H3K79me3")) %>%
              ggplot(aes(x=genecluster, y=value, fill=factor(name, levels =c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")))) +
                  geom_hline(yintercept=0) +
                  geom_boxplot(outliers = FALSE) +
                  annotate("text", x = 0.6, y = 5.5, label = paste0("n=",numrow),size=8) +
                  scale_y_continuous(name = "log2(Input mean normalized coverage)",limits=c(-6.5,6.5),breaks=c(-4,4)) +
                  scale_fill_manual(name="Mark",values=cols) +
                  theme_minimal() +
                  theme(axis.text.x = element_text(color = "black",size = 12,hjust = 1,angle = 45),
                        axis.text.y = element_text(color = "black",size = 12),
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(color = "black",size = 12),
                        legend.title = element_text(color = "black",size = 12),
                        legend.text = element_text(colour = "black",size = 12),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.major.y = element_line(linetype = "dashed",colour = "black")))
        dev.off()
    }
}
colnames(repStates) <- c("species","mark","cluster","med","q1","lvl","spcluster")
pdf(paste0(plot_folder,"Fig4_draft/Repressive_State_Contributions.pdf"),width=10,height=30) #The police?
repStates%>%filter(mark=="H3K79me1" | mark=="H3K79me2" | mark=="H3K79me3" | mark=="H3K27me3" | mark=="H3K9me1" | mark=="H3K9me3")%>%
    ggplot(aes(x=factor(mark,levels=c("H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3")),y=factor(spcluster,levels=rev(unique(spcluster))),fill=factor(lvl,levels=c(-1,0,1,2)))) +
        geom_tile(color="grey") +
        scale_fill_manual(name="Contribution",values = c("-1"="yellow","0"="white","1"="grey","2"="black"),labels=c("-1"="Absent","0"="None","1"="Minor","2"="Major")) +
        theme_classic() +
        theme(axis.text.x = element_text(color = "black",size = 12,hjust=1,angle=45),
              axis.text.y = element_text(color = "black",size = 12),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.title = element_text(color = "black",size = 12),
              legend.text = element_text(colour = "black",size = 12))
dev.off()

###############################################################################################################################################################################################################################################################
##Play around with Fisher exact test on mark coverage and RNAcounts per 200bp bin in gene per mark per species
first_last_bin=function(grange){
    tes=grange[mcols(grange)$type != "protein_coding"]
    grange=grange[mcols(grange)$type == "protein_coding"]
    bin1=pintersect(grange%>%mutate(type=as.factor("bin1"))%>%resize(200),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    bin_n=pintersect(grange%>%mutate(type=as.factor("binN"))%>%resize(200,fix="end"),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    bin2=pintersect(grange%>%mutate(type=as.factor("bin2"))%>%resize(200)%>%shift_downstream(200),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    bin3=pintersect(grange%>%mutate(type=as.factor("bin3"))%>%resize(200)%>%shift_downstream(400),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    bin4=pintersect(grange%>%mutate(type=as.factor("bin4"))%>%resize(200)%>%shift_downstream(600),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    bin5=pintersect(grange%>%mutate(type=as.factor("bin5"))%>%resize(200)%>%shift_downstream(800),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    bin6=pintersect(grange%>%mutate(type=as.factor("bin6"))%>%resize(200)%>%shift_downstream(1000),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    bin7=pintersect(grange%>%mutate(type=as.factor("bin7"))%>%resize(200)%>%shift_downstream(1200),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    bin8=pintersect(grange%>%mutate(type=as.factor("bin8"))%>%resize(200)%>%shift_downstream(1400),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    bin9=pintersect(grange%>%mutate(type=as.factor("bin9"))%>%resize(200)%>%shift_downstream(1600),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    gene_body=pintersect(grange%>%mutate(type=as.factor("geneBody"))%>%shift_downstream(1800),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
    grangen=c(bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin_n,gene_body,tes)
    grangen
}
merge_overlap_feature=function(grange){
    tes=grange[mcols(grange)$source == "EDTA"]
    other=grange[mcols(grange)$type == "rRNA" | mcols(grange)$type == "tRNA" | mcols(grange)$type == "miRNA" | mcols(grange)$type == "lncRNA" | mcols(grange)$type == "ncRNA" | mcols(grange)$type == "snoRNA" | mcols(grange)$type == "snRNA" | mcols(grange)$type == "ribozyme" | mcols(grange)$type == "nontranslating_CDS"]
    grange=grange[mcols(grange)$type == "bin1" | mcols(grange)$type == "bin2" | mcols(grange)$type == "bin3" | mcols(grange)$type == "bin4" | mcols(grange)$type == "bin5" | mcols(grange)$type == "bin6" | mcols(grange)$type == "bin7" | mcols(grange)$type == "bin8" | mcols(grange)$type == "bin9" | mcols(grange)$type == "binN" | mcols(grange)$type == "promoter"]
    allgenes <- as.data.frame(grange)%>%
        mutate(start=start-1)
    colnames(allgenes)[1] <- "chrom"
    overlaponly <- as.data.frame(disjoin(grange,ignore.strand=TRUE)[countOverlaps(disjoin(grange,ignore.strand=TRUE),grange,ignore.strand=TRUE) > 1])%>%
        mutate(start=start-1)
    colnames(overlaponly)[1] <- "chrom"
    nooverlap <- as.data.frame(disjoin(grange,ignore.strand=TRUE)[countOverlaps(disjoin(grange,ignore.strand=TRUE),grange,ignore.strand=TRUE) == 1])%>%
        mutate(start=start-1)
    colnames(nooverlap)[1] <- "chrom" 
    overlaponlyfeature=bed_intersect(overlaponly,allgenes)%>%
        filter(.overlap>0)%>%
        mutate(start=start.x+1,end=end.x)%>%
        group_by(chrom,start,end)%>%
        summarise(type=paste(type.y, collapse = ""),gene_id=unique(gene_id.y),transcript_id=unique(transcript_id.y))%>%
        mutate(type=gsub(".*bin1.*","bin1",type))%>%
        mutate(type=gsub(".*binN.*","binN",type))%>%
        mutate(type=gsub("bin[1-9]","geneBody",type))%>%
        mutate(type=gsub(".*geneBodypromoter.*","promotergeneBody",type))%>%
        mutate(type=gsub("geneBody.*Body","geneBody",type))%>%
        mutate(type=gsub("prom.*moter","promoter",type))%>%
        ungroup()%>%
        GRanges()%>%
        group_by(type,gene_id,transcript_id)%>%reduce_ranges()%>%
        mutate(source="me")
    nooverlapfeature =bed_intersect(nooverlap,allgenes)%>%
        filter(.overlap>0)%>%
        mutate(start=start.x+1,end=end.x,strand=strand.y,source=source.y,type=type.y,score=score.y,phase=phase.y,transcript_id=transcript_id.y,gene_id=gene_id.y,class_id=class_id.y,order_id=order_id.y,family_id=family_id.y)%>%
        dplyr::select(chrom,start,end,strand,source,type,score,phase,transcript_id,gene_id,class_id,order_id,family_id)%>%
        GRanges()
    grangem=c(nooverlapfeature,overlaponlyfeature,tes,other)
    grangem
}
grange_to_feature=function(grange){
    feature=as_tibble(grange)[,c(1,6,7,2,3,8,5,9,10,11,12,13,14)]
    colnames(feature)=c("chrom","Source","feature","start","end","X6","strand","X8","geneId","gene_id","class_id","order_id","family_id")
    feature=feature%>%mutate(start=start-1)%>%
    mutate(feature=gsub("promoter","Promoter",feature))%>%
    mutate(feature=gsub("bin1","Bin 1",feature))%>%
    mutate(feature=gsub("bin2","Bin 2",feature))%>%
    mutate(feature=gsub("bin3","Bin 3",feature))%>%
    mutate(feature=gsub("bin4","Bin 4",feature))%>%
    mutate(feature=gsub("bin5","Bin 5",feature))%>%
    mutate(feature=gsub("bin6","Bin 6",feature))%>%
    mutate(feature=gsub("bin7","Bin 7",feature))%>%
    mutate(feature=gsub("bin8","Bin 8",feature))%>%
    mutate(feature=gsub("bin9","Bin 9",feature))%>%
    mutate(feature=gsub("binN","Last bin",feature))%>%
    mutate(feature=gsub("geneBody","Gene body",feature))%>%
    mutate(feature=gsub("TE_RT_LTR_intact","TE RT LTR intact",feature))%>%
    mutate(feature=gsub("TE_RT_LTR_fragment","TE RT LTR fragment",feature))%>%
    mutate(feature=gsub("TE_RT_LINE","TE RT LINE",feature))%>%
    mutate(feature=gsub("TE_RT_SINE","TE RT SINE",feature))%>%
    mutate(feature=gsub("TE_RT_DIRS","TE RT DIRS",feature))%>%
    mutate(feature=gsub("TE_DNA_intact","TE DNA intact",feature))%>%
    mutate(feature=gsub("TE_DNA_fragment","TE DNA fragment",feature))%>%
    mutate(feature=gsub("intergenic","Intergenic",feature))%>%
    mutate(feature=gsub("interte","Inter-TE",feature))%>%
    mutate(feature=gsub("promotergeneBody","Promoter + Gene body",feature))%>%
    mutate(feature=gsub("promoter","Misc Promoter",feature))%>%
    mutate(feature=gsub("geneBody","Misc Gene body",feature))%>%
    mutate(feature=gsub("protein_coding","Gene",feature))%>%
    mutate(feature=factor(feature,levels=c("Gene","Promoter","Bin 1","Bin 2","Bin 3","Bin 4","Bin 5","Bin 6","Bin 7","Bin 8","Bin 9","Last bin","Gene body",
    "TE RT LTR intact", "TE RT LTR fragment","TE RT LINE","TE RT SINE","TE RT DIRS","TE DNA intact","TE DNA fragment",
    "Intergenic","Inter-TE","Promoter + Gene body","Misc Promoter","Misc Gene body",
    "rRNA","tRNA","miRNA","lncRNA","ncRNA","sRNA","snoRNA","snRNA","ribozyme","SRP_RNA","RNase_MRP_RNA","sense_intronic","telomerase_RNA","nontranslating_CDS","pseudogene","misc_RNA")))
}
for (species in specieslist$species){
    grange = get_our_grange(species=species)
    tes=grange[mcols(grange)$source == "EDTA"]
    if (species == "Acas"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Atha"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Bnat"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k5.reorder.txt"))
    } else if (species == "Cfra"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k5.txt"))
    } else if (species == "Ddis"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k5.reorder.txt"))
    } else if (species == "Gthe"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k8.reorder.txt"))
    } else if (species == "Ngru"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k2.txt"))
    } else if (species == "Nvec"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Ppat"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k6.reorder.txt"))
    } else if (species == "Scer"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k2.txt"))
    } else if (species == "Spun"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k3.reorder.txt"))
    } else if (species == "Tthe"){
        teregions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".TEnogene.region.heatmap.k1.txt"))
    }
    teregions$order <- row.names(teregions)
    tes$order <- as.numeric(teregions[match(tes$transcript_id,teregions$V4),]$order)
    teregions$class_id <- tes$class_id[match(teregions$V4, tes$transcript_id)]
    teregions$order_id <- tes$order_id[match(teregions$V4, tes$transcript_id)]
    teregions$family_id <- tes$family_id[match(teregions$V4, tes$transcript_id)]
    teregions$type <- tes$type[match(teregions$V4, tes$transcript_id)]
    teregions$V13 <- paste0("TE_",teregions$V13)
    ##Find top/bottom expressed gene cluster
    if (species == "Gthe"){
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.reorder.txt"))
    } else {
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.reorder.txt"))
    }
    RNAcountsfile <- read.table(paste0("/home/smontgomery/cluster/smontgomery/genomes/",species,"/",species,"_expr_sorted.tsv"),header=TRUE)
    if (species == "Cfra"){
        RNAcountsfile$gene_id <- gsub("Cfra_CFRG0","Cfra_CFRCFRG0",RNAcountsfile$gene_id)
    }
    generegions$V13 <- paste0("gene_",generegions$V13)
    generegions$RNAcounts <- RNAcountsfile$RNAcounts[match(generegions$V4,RNAcountsfile$gene_id)]
    allregions <- rbind(generegions[,1:13],teregions[,1:13])
    remove(RNAcountsfile,teregions,tes)
    ##Start processing
    grange_w_bins=first_last_bin(grange)
    grange_w_promoter=add_promoter(grange_w_bins)
    grange_wo_overlaps=merge_overlap_feature(grange_w_promoter)
    grange_w_intergenic=add_intergenic(grange_wo_overlaps)
    feature=grange_to_feature(grange_w_intergenic)
    coveragefile <- as_tibble(read.table(paste0(chip_folder,species,"_all/bed/",species,".coverage.RPGC.txt.gz"), header = T))
    colnames(coveragefile) <- gsub("H3K9me2","H3K9me3",colnames(coveragefile))
    for (mark in c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3","input")){
        if(length(grep(mark,colnames(coveragefile)))==0){
            coveragefile[,paste0(species,"_",mark)] <- 0
        }
    }
    fea_cov=bed_intersect(coveragefile,feature)%>%mutate(start=start.x,end=end.x)
    remove(coveragefile,feature,grange_w_bins,grange_w_promoter,grange_wo_overlaps,grange_w_intergenic)
    crm=c("Pt","Mt","Mito")
    H3peaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_H3_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(H3peaks)[1:3] <- c("chrom","start","end")
    H3peaks <- c(GRanges(H3peaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    H3peaks=H3peaks[!(seqnames(H3peaks) %in% c(crm))]
    inputpeaks <- read.table(paste0(chip_folder,species,"_all/macs2/",species,"_input_peaks.broadPeak"),sep='\t',header = FALSE)
    colnames(inputpeaks)[1:3] <- c("chrom","start","end")
    inputpeaks <- c(GRanges(inputpeaks),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,species,"/",species,".chromsizes.txt"))$V1)))
    inputpeaks=inputpeaks[!(seqnames(inputpeaks) %in% c(crm))]
    ctrlpeaks <- as_tibble(union_ranges(H3peaks,inputpeaks))%>%mutate(start=start-1)
    remove(H3peaks,inputpeaks)
    colnames(ctrlpeaks)[1] <- c("chrom")
    red_fea_cov <- bed_subtract(fea_cov,ctrlpeaks)%>%filter(chrom != "chrCH709182" | species != "Ddis")
    remove(ctrlpeaks,fea_cov)
    markCoverage_red=red_fea_cov%>%
        select(geneId.y,gene_id.y,feature.y,.overlap,!!as.name(paste0(quo_name(species),'_H3K27me3.x')),!!as.name(paste0(quo_name(species),'_H3K9me1.x')),!!as.name(paste0(quo_name(species),'_H3K9me3.x')),
           !!as.name(paste0(quo_name(species),'_H3K79me1.x')),!!as.name(paste0(quo_name(species),'_H3K79me2.x')),!!as.name(paste0(quo_name(species),'_H3K79me3.x'))
           ,!!as.name(paste0(quo_name(species),'_H3K9ac.x')),!!as.name(paste0(quo_name(species),'_H3K27ac.x')),!!as.name(paste0(quo_name(species),'_H4K16ac.x'))
           ,!!as.name(paste0(quo_name(species),'_H3K4me2.x')),!!as.name(paste0(quo_name(species),'_H3K4me3.x')),!!as.name(paste0(quo_name(species),'_H3K36me3.x'))
           ,!!as.name(paste0(quo_name(species),'_H3.x')),!!as.name(paste0(quo_name(species),'_input.x')))%>%
        mutate(feature.y=gsub("TE RT LTR fragment","TE RT",feature.y))%>%
        mutate(feature.y=gsub("TE RT LINE","TE RT",feature.y))%>%
        mutate(feature.y=gsub("TE RT SINE","TE RT",feature.y))%>%
        mutate(feature.y=gsub("TE RT DIRS","TE RT",feature.y))%>%
        mutate(feature.y=gsub("TE DNA fragment","TE DNA",feature.y))%>%
        mutate(quanH3K27me3=(log2(!!as.name(paste0(quo_name(species),'_H3K27me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K27me3.x")]))))%>%
        mutate(quanH3K9me1=(log2(!!as.name(paste0(quo_name(species),'_H3K9me1.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K9me1.x")]))))%>%
        mutate(quanH3K9me3=(log2(!!as.name(paste0(quo_name(species),'_H3K9me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K9me3.x")]))))%>%
        mutate(quanH3K79me1=(log2(!!as.name(paste0(quo_name(species),'_H3K79me1.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K79me1.x")]))))%>%
        mutate(quanH3K79me2=(log2(!!as.name(paste0(quo_name(species),'_H3K79me2.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K79me2.x")]))))%>%
        mutate(quanH3K79me3=(log2(!!as.name(paste0(quo_name(species),'_H3K79me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K79me3.x")]))))%>%
        mutate(quanH3K9ac=(log2(!!as.name(paste0(quo_name(species),'_H3K9ac.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K9ac.x")]))))%>%
        mutate(quanH3K27ac=(log2(!!as.name(paste0(quo_name(species),'_H3K27ac.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K27ac.x")]))))%>%
        mutate(quanH4K16ac=(log2(!!as.name(paste0(quo_name(species),'_H4K16ac.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H4K16ac.x")]))))%>%
        mutate(quanH3K4me2=(log2(!!as.name(paste0(quo_name(species),'_H3K4me2.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K4me2.x")]))))%>%
        mutate(quanH3K4me3=(log2(!!as.name(paste0(quo_name(species),'_H3K4me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K4me3.x")]))))%>%
        mutate(quanH3K36me3=(log2(!!as.name(paste0(quo_name(species),'_H3K36me3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3K36me3.x")]))))%>%
        mutate(quanH3=(log2(!!as.name(paste0(quo_name(species),'_H3.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_H3.x")]))))%>%
        mutate(quaninput=(log2(!!as.name(paste0(quo_name(species),'_input.x')))-log2(colMeans(red_fea_cov[,paste0(species,"_input.x")]))))%>%
        group_by(geneId.y,gene_id.y,feature.y,.drop=TRUE)%>%
        summarize(H3K27me3=mean(quanH3K27me3[!is.infinite(quanH3K27me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K9me1=mean(quanH3K9me1[!is.infinite(quanH3K9me1)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K9me3=mean(quanH3K9me3[!is.infinite(quanH3K9me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K79me1=mean(quanH3K79me1[!is.infinite(quanH3K79me1)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K79me2=mean(quanH3K79me2[!is.infinite(quanH3K79me2)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K79me3=mean(quanH3K79me3[!is.infinite(quanH3K79me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K9ac=mean(quanH3K9ac[!is.infinite(quanH3K9ac)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K27ac=mean(quanH3K27ac[!is.infinite(quanH3K27ac)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H4K16ac=mean(quanH4K16ac[!is.infinite(quanH4K16ac)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K4me2=mean(quanH3K4me2[!is.infinite(quanH3K4me2)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K4me3=mean(quanH3K4me3[!is.infinite(quanH3K4me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3K36me3=mean(quanH3K36me3[!is.infinite(quanH3K36me3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE),
              H3=mean(quanH3[!is.infinite(quanH3)],na.rm=TRUE)-mean(quaninput[!is.infinite(quaninput)],na.rm=TRUE))%>%
        ungroup()%>%filter(rowSums(select(., 4:16),na.rm=TRUE) != 0)%>%pivot_longer(cols = 4:16)%>%filter(!is.na(value))%>%
        mutate(genecluster=if_else((feature.y=="Gene" | feature.y=="Promoter" | feature.y=="Bin 1" | feature.y=="Bin 2" | feature.y=="Bin 3" | feature.y=="Bin 4" | feature.y=="Bin 5" | feature.y=="Bin 6" | feature.y=="Bin 7" | feature.y=="Bin 8" | feature.y=="Bin 9" | feature.y=="Gene body" | feature.y=="Last bin"),allregions$V13[match(gene_id.y, allregions$V4)],feature.y))%>%
        mutate(tecluster=if_else((feature.y=="TE RT" | feature.y=="TE RT LTR intact" | feature.y=="TE DNA intact" | feature.y=="TE DNA"),allregions$V13[match(geneId.y, allregions$V4)],feature.y))%>%
        mutate(genecluster=gsub("gene_cluster_","Gene ",genecluster))%>%
        mutate(RNAcounts=if_else((feature.y=="Gene" | feature.y=="Promoter" | feature.y=="Bin 1" | feature.y=="Bin 2" | feature.y=="Bin 3" | feature.y=="Bin 4" | feature.y=="Bin 5" | feature.y=="Bin 6" | feature.y=="Bin 7" | feature.y=="Bin 8" | feature.y=="Bin 9" | feature.y=="Gene body" | feature.y=="Last bin"),generegions$RNAcounts[match(gene_id.y, generegions$V4)],NA))%>%
        mutate(size=if_else((feature.y=="Gene" | feature.y=="Promoter" | feature.y=="Bin 1" | feature.y=="Bin 2" | feature.y=="Bin 3" | feature.y=="Bin 4" | feature.y=="Bin 5" | feature.y=="Bin 6" | feature.y=="Bin 7" | feature.y=="Bin 8" | feature.y=="Bin 9" | feature.y=="Gene body" | feature.y=="Last bin"),grange$size[match(gene_id.y, grange$gene_id)],NA))
    nam <- paste0("marCoverage_red_",species)
    assign(nam, markCoverage_red)
}
#Show coverage only
fish_out <- data.frame()
for (species in specieslist$species){
    tmptable <- get(paste0("marCoverage_red_",species))
    for (mark in c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3")){
        for (feature in c("Promoter","Bin 1","Bin 2","Bin 3","Bin 4","Bin 5","Bin 6","Bin 7","Bin 8","Bin 9","Gene body","Last bin")){
            for (cluster in c("Gene 1","Gene 2","Gene 3","Gene 4","Gene 5","Gene 6")){
                tmpsamp <- subset(tmptable,name==mark & feature.y==feature & genecluster==cluster)
                fish_out[paste0(species,mark,feature,cluster),1] <- species
                fish_out[paste0(species,mark,feature,cluster),2] <- mark
                fish_out[paste0(species,mark,feature,cluster),3] <- feature
                fish_out[paste0(species,mark,feature,cluster),4] <- cluster
                fish_out[paste0(species,mark,feature,cluster),5] <- if_else(nrow(tmpsamp) > nrow(subset(tmptable,name==mark & feature.y=="Bin 1" & genecluster==cluster))/5, mean(tmpsamp$value,rm.na=TRUE),NA)
                fish_out[paste0(species,mark,"RNAcounts",cluster),1] <- species
                fish_out[paste0(species,mark,"RNAcounts",cluster),2] <- mark
                fish_out[paste0(species,mark,"RNAcounts",cluster),3] <- "RNAcounts"
                fish_out[paste0(species,mark,"RNAcounts",cluster),4] <- cluster
                fish_out[paste0(species,mark,"RNAcounts",cluster),5] <- median(asinh(subset(tmptable,name==mark & feature.y=="Bin 1" & genecluster==cluster)$RNAcounts),na.rm=TRUE)
            }
        }
    }
}
colnames(fish_out) <- c("species","mark","feature","cluster","coverage")
pdf(paste0(plot_folder,"Fig3_draft/Acetylation_coverage.RPGC.pdf"),width=15,height=10)
fish_out%>%filter(feature != "Gene body" & feature != "RNAcounts" & (mark=="H3K9ac" | mark=="H3K27ac" | mark=="H4K16ac" | mark=="H3K27me3" | mark=="H3K9me3") & (species=="Atha" | species=="Spun" | species=="Cfra" | species=="Ngru"))%>%
              ggplot(aes(x=factor(feature,levels=c("Promoter","Bin 1","Bin 2","Bin 3","Bin 4","Bin 5","Bin 6","Bin 7","Bin 8","Bin 9","Last bin")),y=factor(cluster,levels=rev(c("Gene 1","Gene 2","Gene 3","Gene 4","Gene 5","Gene 6"))),fill=scales::oob_squish(coverage,range = c(-2,2)))) +
              geom_tile(color = "white")+
              scale_fill_gradient2(low = "steelblue", high = "red3", mid = "white",midpoint = 0, space = "Lab",name="Coverage",na.value = "grey") +
              facet_wrap(vars(factor(species,levels = c("Cfra","Spun","Atha","Ngru")),factor(mark,levels=c("H3K9ac","H3K27ac","H4K16ac","H3K27me3","H3K9me3"))), strip.position = c("bottom"),nrow = 4) +
              theme_classic()+ # minimal theme
              labs(x = NULL, y = NULL) +
              theme(axis.text.x = element_text(size=18,hjust=1,angle=45),
                    axis.text.y = element_text(size=18),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.text = element_text(size=18),
                    legend.title = element_text(size=18),
                    strip.text.x = element_text(size=18))
dev.off()



###############################################################################################################################################################################################################################################################
##Pearson correlation and hierarchical tree clustering
##Combine overlap tables
##Combine emission tables
i=1
for (species in specieslist$species){
    statenumber <- specieslist$statenumber[grep(species,specieslist$species)]
    tmpfile <- t(read.table(paste0(chip_folder,species,"_all/state_",statenumber,"/emissions_",statenumber,".txt"), sep = '\t', header = TRUE, row.names = 1))
    row.names(tmpfile) <- gsub("H3K9me2","H3K9me3",gsub(".*_","",row.names(tmpfile)))
    colnames(tmpfile) <- paste0(species,"_State",statenumber,"_",colnames(tmpfile))
    tmpfile <- tmpfile[order(row.names(tmpfile)),]
    if (i == 1){
        i=0
        emission_table <- tmpfile
    } else {
        emission_table <- merge.data.frame(emission_table,tmpfile,all = TRUE,by = 0)
        emission_table[is.na(emission_table)] <- 0
        row.names(emission_table) <- emission_table[,1]
        emission_table <- emission_table[,-1]
    }
}
i=1
for (species in specieslist$species){
    statenumber <- specieslist$statenumber[grep(species,specieslist$species)]
    cs=c(brewer.pal(7,"Reds"),brewer.pal(7,"Blues"),brewer.pal(7,"Greys"),brewer.pal(7,"Purples"),brewer.pal(7,"Oranges"),brewer.pal(7,"Greens"),brewer.pal(8,"Set3"))[1:statenumber]
    map_to_paper = data.frame(
      from=c(as.character(1:statenumber)),
      to=factor(c(paste0("State",1:statenumber)),
                levels = c(paste0("State",1:statenumber))),
      cols=cs)
    large_col=read_and_fix_states(paste0(chip_folder,species,"_all/state_",statenumber,"/",species,"_",statenumber,"_segments.bed")) # large WT 
    # Change state names to fit paper (full model) using map_to_paper defined in settings.R
    large_col=mutate(large_col,state=str_remove(state,"^E"))%>%
        left_join(map_to_paper,by=c("state"="from"))
    large_col=dplyr::select(large_col,-state)%>%mutate(state=to)
    large_col=large_col%>%mutate(state=factor(state,levels = c(paste0("State",1:statenumber))))# annotations
    ##Get genome coverage percentages
    genome_dist = large_col%>%mutate(length=end-start)%>%group_by(state)%>%summarise(tot=sum(length))%>%
          ungroup()%>%mutate(all=sum(tot),proc=tot/all*100)%>%
          mutate(state=factor(paste0(species,"_State",statenumber,"_",gsub("State","",state)),levels=c(paste0(species,"_State",statenumber,"_",1:statenumber))))%>%
          column_to_rownames(var="state")
    ##Do it with granges
    grange = get_our_grange(species=species)
    grange_w_bins=add_bins(grange)
    grange_w_promoter=add_promoter(grange_w_bins)
    grange_wo_overlaps=merge_overlap_feature(grange_w_promoter)
    grange_w_intergenic=add_intergenic(grange_wo_overlaps)
    grange_w_expression=add_expression(grange_w_intergenic)
    feature=grange_to_feature(grange_w_expression)
    ##Find top/bottom expressed gene cluster
    if (species == "Gthe"){
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.reorder.txt"))
    } else {
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.reorder.txt"))
    }
    RNAcountsfile <- read.table(paste0("/home/smontgomery/cluster/smontgomery/genomes/",species,"/",species,"_expr_sorted.tsv"),header=TRUE)
    if (species == "Cfra"){
        RNAcountsfile$gene_id <- gsub("Cfra_CFRG0","Cfra_CFRCFRG0",RNAcountsfile$gene_id)
    }
    generegions$V13 <- paste0("gene_",generegions$V13)
    generegions$RNAcounts <- RNAcountsfile$RNAcounts[match(generegions$V4,RNAcountsfile$gene_id)]
    remove(RNAcountsfile)
    ##Reduce number of features to be plotted
    feature_red=feature%>%mutate(feature=gsub("TE DNA fragment","TE DNA",feature))%>%
                mutate(feature=gsub("TE DNA intact","TE DNA",feature))%>%
                mutate(feature=gsub("TE RT SINE","TE RT nonLTR",feature))%>%
                mutate(feature=gsub("TE RT DIRS","TE RT nonLTR",feature))%>%
                mutate(feature=gsub("TE RT LINE","TE RT nonLTR",feature))%>%
                mutate(feature=gsub("TE RT LTR fragment","TE RT LTR",feature))%>%
                mutate(feature=gsub("TE RT LTR intact","TE RT LTR",feature))%>%
                mutate(feature=gsub("Misc Gene body","Med",feature))%>%
                mutate(feature=gsub("Misc Promoter","Med",feature))%>%
                mutate(feature=gsub("Promoter \\+ Gene body","Med",feature))%>%
                mutate(feature=gsub("Inter-TE","Med",feature))%>%
                mutate(feature=gsub("rRNA","Med",feature))%>%
                mutate(feature=gsub("tRNA","Med",feature))%>%
                mutate(feature=gsub("miRNA","Med",feature))%>%
                mutate(feature=gsub("lncRNA","Med",feature))%>%
                mutate(feature=gsub("ncRNA","Med",feature))%>%
                mutate(feature=gsub("sRNA","Med",feature))%>%
                mutate(feature=gsub("snoRNA","Med",feature))%>%
                mutate(feature=gsub("snRNA","Med",feature))%>%
                mutate(feature=gsub("ribozyme","Med",feature))%>%
                mutate(feature=gsub("nontranslating_CDS","Med",feature))%>%
                mutate(feature=gsub("pseudogene","Med",feature))%>%
                mutate(feature=gsub("SRP_RNA","Med",feature))%>%
                mutate(feature=gsub("telomerase_RNA","Med",feature))%>%
                mutate(feature=factor(feature,levels=c("Gene","Promoter High","Bin 1 High","Bin 2 High","Bin 3 High","Bin 4 High","Gene body High",
                    "Promoter Med","Bin 1 Med","Bin 2 Med","Bin 3 Med","Bin 4 Med","Gene body Med",
                    "Promoter Low","Bin 1 Low","Bin 2 Low","Bin 3 Low","Bin 4 Low","Gene body Low",
                    "TE RT LTR","TE RT nonLTR","TE DNA",
                    "Intergenic","Med")))%>%
                mutate(RNAcounts=if_else((feature=="Gene" | feature=="Promoter High" | feature=="Bin 1 High" | feature=="Bin 2 High" | feature=="Bin 3 High" | feature=="Bin 4 High" | feature=="Gene body High" | 
                    feature=="Promoter Med" | feature=="Bin 1 Med" | feature=="Bin 2 Med" | feature=="Bin 3 Med" | feature=="Bin 4 Med" | feature=="Gene body Med" | 
                    feature=="Promoter Low" | feature=="Bin 1 Low" | feature=="Bin 2 Low" | feature=="Bin 3 Low" | feature=="Bin 4 Low" | feature=="Gene body Low"),generegions$RNAcounts[match(gene_id, generegions$V4)],NA))
    f_col=bed_intersect(large_col,feature_red)%>%mutate(state.x=factor(paste0(species,"_State",statenumber,"_",gsub("State","",state.x)),levels=c(paste0(species,"_State",statenumber,"_",1:statenumber))))
    tmpfile=f_col%>%group_by(state.x,feature.y)%>%
            summarize(ov=sum(.overlap))%>%
            mutate(n=sum(ov),proc=(ov/n)*100)%>%dplyr::select(-ov,-n)%>%
            pivot_wider(names_from = state.x, values_from = proc,values_fill = 0)%>%
            column_to_rownames(var = "feature.y")
    tmp_enrich=f_col%>%group_by(state.x,feature.y)%>%
            summarize(C=sum(.overlap))%>%mutate(A=sum(C),proc=(C/A)*100)%>%
            ungroup()%>%group_by(feature.y)%>%mutate(B=sum(C))%>%ungroup()%>%
            mutate(D=sum(C),prox=((C/A)/(B/D)))%>%dplyr::select(-A,-B,-C,-D,-proc)%>%
            pivot_wider(names_from = state.x, values_from = prox,values_fill = 0)%>%
            column_to_rownames(var = "feature.y")
    tmp_RNAcounts=f_col%>%group_by(state.x)%>%
            summarize(RNAcounts_med=median(RNAcounts.y,na.rm=TRUE),RNAcounts_mean=mean(RNAcounts.y,na.rm=TRUE))%>%
            column_to_rownames("state.x")
    if (i == 1){
        i=0
        overlap_table <- tmpfile
        enrichment_table <- tmp_enrich
        genome_dist_table <- as.data.frame(genome_dist)
        RNAcounts_table <- tmp_RNAcounts
    } else {
        overlap_table <- merge.data.frame(overlap_table,tmpfile,all = TRUE,by=0)
        overlap_table[is.na(overlap_table)] <- 0
        row.names(overlap_table) <- overlap_table[,1]
        overlap_table <- overlap_table[,-1]
        enrichment_table <- merge.data.frame(enrichment_table,tmp_enrich,all = TRUE,by=0)
        enrichment_table[is.na(enrichment_table)] <- 0
        row.names(enrichment_table) <- enrichment_table[,1]
        enrichment_table <- enrichment_table[,-1]
        genome_dist_table <- rbind(genome_dist_table,as.data.frame(genome_dist))
        RNAcounts_table <- rbind(RNAcounts_table,tmp_RNAcounts)
    }
}

##Calculate pearson correlation and plot hierarchical clustering
##Pull out null states
null_emission_table <- emission_table[,names(emission_table) %in% nullList]
nonull_emission_table <- emission_table[-1,!names(emission_table) %in% nullList]
null_enrichment_table <- enrichment_table[,names(enrichment_table) %in% nullList]
nonull_enrichment_table <- enrichment_table[,!names(enrichment_table) %in% nullList]
null_genome_dist_table <- genome_dist_table[row.names(genome_dist_table) %in% nullList,]
nonull_genome_dist_table <- genome_dist_table[!row.names(genome_dist_table) %in% nullList,]
null_RNAcounts_table <- t(RNAcounts_table[row.names(RNAcounts_table) %in% nullList,])
nonull_RNAcounts_table <- t(RNAcounts_table[!row.names(RNAcounts_table) %in% nullList,])

##Cut tree by height
dissEmission <- 1-(cor(nonull_emission_table,method="pearson"))
hclust_plot <- hclust(as.dist(dissEmission),method = "complete")
cut_avg <- cutree(hclust_plot, h = 0.55)

# ##Plot full tree and heatmaps with boxes around metastates
# par(mar=c(0.1,4.1, 4.1, 2.1))
# plot(hclust_plot,xlab = "",main = "",sub = "")
# rect.hclust(hclust_plot , h = 0.55, border = 2:6)
# tmpplot <- recordPlot()
# pdf(paste0(plot_folder,"metastate_full.pdf"),width=15,height=15)
# plot_grid(tmpplot,
#             ggplot(melt(as.matrix(nonull_emission_table)),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]), y=factor(Var1, levels = rev(c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3"))), fill=value)) +
#                 geom_tile(color = "white")+
#                 scale_fill_scico(palette = 'bilbao',direction = -1,na.value="white",limits=c(0,1),name="Emission") +
#                 theme_classic()+ # minimal theme
#                 scale_x_discrete(position = "bottom") +
#                 labs(x = NULL, y = NULL) +theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7),axis.text.y = element_text(size=7),plot.margin=unit(c(0,22,0,17),"mm")),
#             ggplot(melt(as.matrix(nonull_enrichment_table[-grep("Med",row.names(nonull_enrichment_table)),])),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]), 
#                                                                                                                 y=factor(Var1, levels = c("Intergenic","TE DNA","TE RT nonLTR","TE RT LTR",
#                                                                                                                                           "Gene body Low","Bin 4 Low","Bin 3 Low","Bin 2 Low","Bin 1 Low","Promoter Low",
#                                                                                                                                           "Gene body High","Bin 4 High","Bin 3 High","Bin 2 High","Bin 1 High","Promoter High")), 
#                                                                                                                 fill=log2(value))) +
#                 geom_tile(color = "white")+
#                 scale_fill_scico(palette = 'tokyo',limits=c(log2(1),log2(max(nonull_enrichment_table))),direction = -1,na.value="white",space="Lab",name="Enrichment") +
#                 theme_classic()+ # minimal theme
#                 scale_x_discrete(position = "bottom") +
#                 labs(x = NULL, y = NULL) +theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7),axis.text.y = element_text(size=7),plot.margin=unit(c(0,22,0,7),"mm")),
#             ggplot(melt(as.matrix(nonull_RNAcounts_table))%>%filter(Var1=="RNAcounts_med"),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]), y=factor(Var1), fill=value)) +
#                 geom_tile(color = "white")+
#                 scale_fill_scico(palette = 'grayC',direction = -1,na.value="white",name="RNAcounts") +
#                 theme_classic()+ # minimal theme
#                 scale_x_discrete(position = "bottom") +
#                 labs(x = NULL, y = NULL) +theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7),axis.text.y = element_text(size=7),plot.margin=unit(c(0,22,0,17),"mm")),
#         ncol = 1,rel_heights = c(2,1,1,0.3),axis = "r")
# dev.off()

# ##Calculate metastate emissions
# metastate_emission_tibble <- cbind(as.data.frame(t(nonull_emission_table)),cut_avg)%>%
#     group_by(cut_avg)%>%
#     filter(n() >= 3)%>%
#     summarize(H3K27ac=median(H3K27ac),
#                 H3K27me3=median(H3K27me3),
#                 H3K36me3=median(H3K36me3),
#                 H3K4me2=median(H3K4me2),
#                 H3K4me3=median(H3K4me3),
#                 H3K79me1=median(H3K79me1),
#                 H3K79me2=median(H3K79me2),
#                 H3K79me3=median(H3K79me3),
#                 H3K9ac=median(H3K9ac),
#                 H3K9me1=median(H3K9me1),
#                 H3K9me3=median(H3K9me3),
#                 H4K16ac=median(H4K16ac))%>%
#     column_to_rownames(var = "cut_avg")
# ##Plot metastate heatmap
# pdf(paste0(plot_folder,"metastate_emission.pdf"),width=15,height=15)
# ggplot(melt(as.matrix(metastate_emission_tibble)),aes(x=factor(Var1, levels=c(7,9,14,4,15,6,10,3,5,8,11,18,12,16,13,19,1,2,17)), y=factor(Var2, levels = rev(c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3"))), fill=value)) +
#     geom_tile(color = "white",width=0.8,height=0.8)+
#     scale_fill_scico(palette = 'bilbao',direction = -1,na.value="white",limits=c(0.25,1),breaks=c(0.25,0.5,0.75,1),name="Emission") +
#     theme_classic() +
#     scale_x_discrete(position = "bottom",name="MetaState",labels=c(1:length(unique(cut_avg)))) +
#     theme(
#         axis.text.x = element_text(color = "black",size = 12),
#         axis.text.y = element_text(color = "black",size = 12),
#         axis.title.x = element_text(color = "black",size = 12),
#         axis.title.y = element_blank(),
#         legend.title = element_text(color = "black",size = 12),
#         legend.text = element_text(colour = "black",size = 12))
# dev.off()

# ##Calculate metastate genome coverage percentages
# metastate_genome_dist_tibble <- cbind(nonull_genome_dist_table,cut_avg)%>%
#     mutate(species=gsub("_.*","",row.names(nonull_genome_dist_table)))%>%
#     group_by(cut_avg)%>%
#     filter(n() >= 3)%>%
#     ungroup()%>%
#     group_by(species,cut_avg)%>%
#     summarize(proc=sum(proc))%>%
#     rename("metastate"="cut_avg")
# ##Show heatmap of percentage of genome in metastates per species
# pdf(paste0(plot_folder,"metastate_genome_perc.pdf"),width=15,height=15)
# ggplot(metastate_genome_dist_tibble,aes(x=factor(metastate,levels=c(7,9,14,4,15,6,10,3,5,8,11,18,12,16,13,19,1,2,17)),y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))))) +
#     geom_tile(color = "black",aes(fill=proc),width=0.8,height=0.8,size=1) +
#     theme_classic() +
#     scale_x_discrete(position = "bottom",name="MetaState",labels=c(1:length(unique(cut_avg)))) +
#     scale_fill_gradient2(low = "white", mid = "#000000", high="#000000",midpoint = 33,limits=c(0,68),name=NULL) +
#     geom_text(aes(label = round(proc,0)), color = "black", size = 5) +
#     theme(axis.text.x = element_text(color = "black",size = 12),
#         axis.text.y = element_text(color = "black",size = 12),
#         axis.title.x = element_text(color = "black",size = 12),
#         axis.title.y = element_blank(),
#         legend.title = element_text(color = "black",size = 12),
#         legend.text = element_text(colour = "black",size = 12))
# dev.off()

##Plot null state feature enrichments
pdf(paste0(plot_folder,"extData_Fig3_draft/metastate_null_full_feature.pdf"),width=15,height=15)
ggplot(melt(as.matrix(null_enrichment_table[-grep("Med",row.names(null_enrichment_table)),])),aes(x=Var2, 
                                                                                                  y=factor(Var1, levels = c("Intergenic","TE DNA","TE RT nonLTR","TE RT LTR",
                                                                                                                            "Gene body Low","Bin 4 Low","Bin 3 Low","Bin 2 Low","Bin 1 Low","Promoter Low",
                                                                                                                            "Gene body High","Bin 4 High","Bin 3 High","Bin 2 High","Bin 1 High","Promoter High")), 
                                                                                                  fill=log2(value))) +
    geom_tile(color = "white")+
    scale_fill_scico(palette = 'tokyo',limits=c(log2(1),log2(max(nonull_enrichment_table))),direction = -1,na.value="white",space="Lab",name="Enrichment") +
    theme_classic()+ # minimal theme
    scale_x_discrete(position = "bottom") +
    labs(x = NULL, y = NULL) +theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7),axis.text.y = element_text(size=7),plot.margin=unit(c(0,22,0,7),"mm"))
dev.off()
##Plot null state percentages and emissions
pdf(paste0(plot_folder,"extData_Fig3_draft/metastate_null_full_emission.pdf"),width=15,height=15)
ggplot(melt(as.matrix(null_emission_table)),aes(x=Var2, y=factor(Var1, levels = rev(c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3"))), fill=value)) +
    geom_tile(color = "white")+
    scale_fill_scico(palette = 'bilbao',direction = -1,na.value="white",limits=c(0,1),name="Emission") +
    theme_classic()+ # minimal theme
    scale_x_discrete(position = "bottom") +
    labs(x = NULL, y = NULL) +theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7),axis.text.y = element_text(size=7))
dev.off()
##Calculate null metastate emissions
metastate_null_emission_tibble <- as.data.frame(t(null_emission_table))%>%
    summarize(H3K27ac=median(H3K27ac),
                H3K27me3=median(H3K27me3),
                H3K36me3=median(H3K36me3),
                H3K4me2=median(H3K4me2),
                H3K4me3=median(H3K4me3),
                H3K79me1=median(H3K79me1),
                H3K79me2=median(H3K79me2),
                H3K79me3=median(H3K79me3),
                H3K9ac=median(H3K9ac),
                H3K9me1=median(H3K9me1),
                H3K9me3=median(H3K9me3),
                H4K16ac=median(H4K16ac))%>%
    mutate(tmprowname="null")%>%
    column_to_rownames("tmprowname")
##Plot metastate heatmap
pdf(paste0(plot_folder,"Fig2_draft/metastate_null_emission.pdf"),width=15,height=15)
ggplot(melt(as.matrix(metastate_null_emission_tibble)),aes(x=Var1, y=factor(Var2, levels = rev(c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3","H3"))), fill=value)) +
    geom_tile(color = "white",width=0.8,height=0.8)+
    scale_fill_scico(palette = 'bilbao',limits=c(0.25,1),breaks=c(0.25,0.5,0.75,1),direction = -1,na.value="white",name="Emission") +
    theme_classic() +
    scale_x_discrete(position = "bottom",name="MetaState",labels=c(1:length(unique(cut_avg)))) +
    theme(axis.text.x = element_text(color = "black",size = 12),
          axis.text.y = element_text(color = "black",size = 12),
          axis.title.x = element_text(color = "black",size = 12),
          axis.title.y = element_blank(),
          legend.title = element_text(color = "black",size = 12),
          legend.text = element_text(colour = "black",size = 12))
dev.off()
##Calculate metastate genome coverage percentages
metastate_null_genome_dist_tibble <- null_genome_dist_table%>%
    mutate(species=gsub("_.*","",row.names(null_genome_dist_table)))%>%
    group_by(species)%>%
    summarize(proc=sum(proc))%>%
    mutate(metastate="null")
##Show heatmap of percentage of genome in metastates per species
pdf(paste0(plot_folder,"Fig2_draft/metastate_null_genome_perc.pdf"),width=15,height=15)
ggplot(metastate_null_genome_dist_tibble,aes(x=metastate,y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))))) +
    geom_tile(color = "black",aes(fill=proc),width=0.8,height=0.8,size=1) +
    theme_classic() +
    scale_x_discrete(position = "bottom",name="MetaState",labels=c(1:length(unique(cut_avg)))) +
    scale_fill_gradient2(low = "white", mid = "darkgrey", high="#000000",midpoint = 35,limits=c(0,70),name=NULL) +
    geom_text(aes(label = round(proc,0)), color = "white", size = 5) +
    theme(
        axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12),
        axis.title.x = element_text(color = "black",size = 12),
        axis.title.y = element_blank(),
        legend.title = element_text(color = "black",size = 12),
        legend.text = element_text(colour = "black",size = 12))
dev.off()

##See how clusters change
k=22
hclust_reorder <- hclust(as.dist(dissEmission[hclust_plot$order,hclust_plot$order]),method = "complete")
cut_avg <- cutree(hclust_reorder, k=k)
tmp <- NULL
for (j in 1:30){
    tmp[j] <- list(cutree(hclust_reorder, k=j))
}
df <- data.frame(tmp)
# add a prefix to the column names
colnames(df) <- seq(1:30)
colnames(df) <- paste0("k",colnames(df))
# get individual PCA
df.pca <- prcomp(df, center = TRUE, scale. = FALSE)
ind.coord <- df.pca$x
ind.coord <- ind.coord[,1:2]
df <- bind_cols(as.data.frame(df), as.data.frame(ind.coord))
# pdf(paste0(plot_folder,"cluster_tree.pdf"),width=15,height=15)
#    clustree(df, prefix = "k")
# dev.off()
##Testing cutting tree by max dissimilarity within a cluster
##Plot heatmap of dissimilarity scores per number of hierarchical clusters
dissim <- data.frame()
for (j in 1:30){
    cut_avg <- cutree(hclust_reorder, k=j)
    for (i in 1:j){
        dissim[j,i] <- cbind(as.data.frame(dissEmission[hclust_plot$order,hclust_plot$order]),cut_avg)%>%
            filter(cut_avg==i)%>%
            select(-cut_avg)%>%
            select(rownames(.))%>%
            daisy()%>%
            max()%>%
            round(2)
    }
}
pdf(paste0(plot_folder,"dissimilarity_hierarchical.pdf"),width=15,height=15)
    ggplot(melt(as.matrix(dissim[1:29,1:29])),aes(x=Var2,y=rev(factor(Var1)),fill=log2(value))) +
        geom_tile() +
        scale_fill_gradient2(low="steelblue",high="firebrick",mid="white",midpoint=0,na.value = "white") +
        coord_fixed() +
        scale_y_discrete(labels=paste0("k",29:1),name=NULL) +
        scale_x_discrete(labels=paste0("cl",1:29),name=NULL) +
        theme_minimal()
dev.off()
##Plot cluster tree with dissimilarity scores
dissim_frame <- df
for (j in 1:30){
    cut_avg <- cutree(hclust_reorder, k=j)
    dissim_frame[,j] <- paste(dissim_frame[,j],unlist(rep(dissim[j,!is.na(dissim[j,])],times=c(table(cut_avg)))),sep = "\n")
}
dissim_frame <- dissim_frame%>%mutate(species=gsub("_.*","",rownames(.)))
pdf(paste0(plot_folder,"cluster_tree_dissimilarity.pdf"),width=15,height=15)
    clustree(dissim_frame, prefix = "k", node_label = "species",
             node_label_aggr = "label_position",node_label_nudge=-0.5,node_label_size=2.5)
dev.off()
##Plot data with ideal number of clusters "k"
library(dendextend)
k=22
cut_avg <- cutree(hclust_reorder, k=k)
timesx=c(4,4,3,5,4,7,15,8,11,2,13,2,7,3,6,3,5,6,5,5,6,4) #MANUALLY SET FROM CLUSTER DISSIMILARITY < 1
mslabel=c(1,2,6,7,3,4,5,10,9,0,8,0,16,20,19,17,18,11,12,13,14,15)
cut_avg[1:128] <- rep(c(1:22),times=timesx) 
col_set <- rep(c("steelblue","firebrick"),11)
hclust_reorder_den <- as.dendrogram(hclust_reorder)
hResDen <- color_branches(hclust_reorder_den,clusters=as.numeric(as.factor(cut_avg)),col=col_set)
lableCol<- col_set
names(lableCol)<- unique(cut_avg)
hResDen <- color_labels(hResDen,col=lableCol[as.character(cut_avg)])
par(mar=c(6.1, 2.1, 2.1, 0.1))
pdf(paste0(plot_folder,"hclust_reorder.pdf"),width=15,height=15)
    plot(hResDen)
    text(unique(cut_avg[hclust_reorder$order]),y = 0.5,x=match(unique(cut_avg[hclust_reorder$order]), cut_avg[hclust_reorder$order]))
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))
##Plot hierarchical by number of clusters
metastate_emission_tibble <- cbind(as.data.frame(t(nonull_emission_table)[hclust_plot$order,]),cut_avg)%>%
    group_by(cut_avg)%>%
    filter(n() >= 3)%>%
    summarize(H3K27ac=median(H3K27ac),
              H3K27me3=median(H3K27me3),
              H3K36me3=median(H3K36me3),
              H3K4me2=median(H3K4me2),
              H3K4me3=median(H3K4me3),
              H3K79me1=median(H3K79me1),
              H3K79me2=median(H3K79me2),
              H3K79me3=median(H3K79me3),
              H3K9ac=median(H3K9ac),
              H3K9me1=median(H3K9me1),
              H3K9me3=median(H3K9me3),
              H4K16ac=median(H4K16ac))%>%
    column_to_rownames(var = "cut_avg")
msorder=c(1,2,5,6,7,3,4,11,9,8,18,19,20,21,22,13,16,17,15,14) #MANUALLY REORDER TREE FOR EASE OF READING
pdf(paste0(plot_folder,"Fig2_draft/hclust_reorder_metastate.pdf"),width=8,height=8)
ggplot(melt(as.matrix(metastate_emission_tibble)),aes(x=factor(Var1, levels=msorder), y=factor(Var2, levels = rev(c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3"))), fill=value)) +
    geom_tile(color = "white",width=0.8,height=0.8)+
    scale_fill_scico(palette = 'bilbao',direction = -1,na.value="white",limits=c(0.25,1),breaks=c(0.25,0.5,0.75,1),name="Emission") +
    theme_classic() +
    scale_x_discrete(position = "bottom",name="MetaState",labels=c(1:length(unique(cut_avg)))) +
    theme(
        axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12),
        axis.title.x = element_text(color = "black",size = 12),
        axis.title.y = element_blank(),
        legend.title = element_text(color = "black",size = 12),
        legend.text = element_text(colour = "black",size = 12))
dev.off()

##Plot final labelled tree
pdf(paste0(plot_folder,"Fig2_draft/hclust_reorder_final.pdf"),width=15,height=5)
    plot(hResDen)
    text(mslabel,y = 0.5,x=match(unique(cut_avg[hclust_reorder$order]), cut_avg[hclust_reorder$order])+(cumsum(timesx)-match(unique(cut_avg[hclust_reorder$order]), cut_avg[hclust_reorder$order]))/2)
    # text(dissim,y = c(0.55,0.6),x=match(unique(cut_avg[hclust_reorder$order]), cut_avg[hclust_reorder$order]),col="darkgreen")
dev.off()

##Plot RNA counts of metastates
metastate_RNAcounts_tibble <- cbind(as.data.frame(t(nonull_RNAcounts_table)[hclust_plot$order,]),cut_avg)%>%
    group_by(cut_avg)%>%
    filter(n() >= 3)%>%
    summarize(RNAcounts_med=median(RNAcounts_med))%>%
    column_to_rownames(var = "cut_avg")
pdf(paste0(plot_folder,"Fig2_draft/hclust_reorder_RNAcounts.pdf"),width=15,height=5)
ggplot(melt(as.matrix(metastate_RNAcounts_tibble)),aes(x=factor(Var1, levels=msorder), y=factor(Var2), fill=value)) +
    geom_tile(color = "white",width=0.8,height=0.8)+
    scale_fill_scico(palette = 'grayC',direction = -1,na.value="white",name="RNAcounts",limits=c(0,max(metastate_RNAcounts_tibble))) +
    theme_classic() +
    scale_x_discrete(position = "bottom",name="MetaState",labels=c(1:length(unique(cut_avg)))) +
    theme(
        axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12),
        axis.title.x = element_text(color = "black",size = 12),
        axis.title.y = element_blank(),
        legend.title = element_text(color = "black",size = 12),
        legend.text = element_text(colour = "black",size = 12))
dev.off()

##Calculate metastate genome coverage percentages
metastate_genome_dist_tibble <- cbind(nonull_genome_dist_table,cutavg=cut_avg[match(rownames(nonull_genome_dist_table),names(cut_avg))])%>%
    mutate(species=gsub("_.*","",row.names(nonull_genome_dist_table)))%>%
    group_by(cutavg)%>%
    filter(n() >= 3)%>%
    ungroup()%>%
    group_by(species,cutavg)%>%
    summarize(proc=sum(proc))%>%
    rename("metastate"="cutavg")
##Show heatmap of percentage of genome in metastates per species
pdf(paste0(plot_folder,"Fig2_draft/metastate_genome_perc.pdf"),width=15,height=15)
ggplot(metastate_genome_dist_tibble,aes(x=factor(metastate,levels=msorder),y=factor(species,levels=rev(c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru"))))) +
    geom_tile(color = "black",aes(fill=proc),width=0.8,height=0.8,size=1) +
    theme_classic() +
    scale_x_discrete(position = "bottom",name="MetaState",labels=c(1:length(unique(cut_avg)))) +
    scale_fill_gradient2(low = "white", mid = "#000000", high="#000000",midpoint = 30,limits=c(0,60),name=NULL) +
    geom_text(aes(label = round(proc,0)), color = "black", size = 5) +
    theme(axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12),
        axis.title.x = element_text(color = "black",size = 12),
        axis.title.y = element_blank(),
        legend.title = element_text(color = "black",size = 12),
        legend.text = element_text(colour = "black",size = 12))
dev.off()

##Plot full metastate
par(mar=c(6.1,2.1, 2.1, 0.1))
plot(hResDen)
text(mslabel,y = 0.5,x=match(unique(cut_avg[hclust_reorder$order]), cut_avg[hclust_reorder$order])+(cumsum(timesx)-match(unique(cut_avg[hclust_reorder$order]), cut_avg[hclust_reorder$order]))/2)
tmpplot <- recordPlot()
pdf(paste0(plot_folder,"extData_Fig3_draft/metastate_full.pdf"),width=15,height=15)
plot_grid(tmpplot,
            ggplot(melt(as.matrix(nonull_emission_table)),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]), y=factor(Var1, levels = rev(c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me3"))), fill=value)) +
                geom_tile(color = "white")+
                scale_fill_scico(palette = 'bilbao',direction = -1,na.value="white",limits=c(0,1),name="Emission") +
                theme_classic()+ # minimal theme
                scale_x_discrete(position = "bottom") +
                labs(x = NULL, y = NULL) +theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7),axis.text.y = element_text(size=7),plot.margin=unit(c(0,22,0,17),"mm")),
            ggplot(melt(as.matrix(nonull_enrichment_table[-grep("Med",row.names(nonull_enrichment_table)),])),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]), 
                                                                                                                y=factor(Var1, levels = c("Intergenic","TE DNA","TE RT nonLTR","TE RT LTR",
                                                                                                                                          "Gene body Low","Bin 4 Low","Bin 3 Low","Bin 2 Low","Bin 1 Low","Promoter Low",
                                                                                                                                          "Gene body High","Bin 4 High","Bin 3 High","Bin 2 High","Bin 1 High","Promoter High")), 
                                                                                                                fill=log2(value))) +
                geom_tile(color = "white")+
                scale_fill_scico(palette = 'tokyo',limits=c(log2(1),log2(max(nonull_enrichment_table))),direction = -1,na.value="white",space="Lab",name="Enrichment") +
                theme_classic()+ # minimal theme
                scale_x_discrete(position = "bottom") +
                labs(x = NULL, y = NULL) +theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7),axis.text.y = element_text(size=7),plot.margin=unit(c(0,22,0,7),"mm")),
            ggplot(melt(as.matrix(nonull_RNAcounts_table))%>%filter(Var1=="RNAcounts_med"),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]), y=factor(Var1), fill=value)) +
                geom_tile(color = "white")+
                scale_fill_scico(palette = 'grayC',direction = -1,na.value="white",name="RNAcounts") +
                theme_classic()+ # minimal theme
                scale_x_discrete(position = "bottom") +
                labs(x = NULL, y = NULL) +theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7),axis.text.y = element_text(size=7),plot.margin=unit(c(0,22,0,17),"mm")),
        ncol = 1,rel_heights = c(2,1,1,0.3),axis = "r")
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))

##Plot correlation matrix of all chromatin states
pdf(paste0(plot_folder,"extData_Fig3_draft/metastate_correlation_all.pdf"),width=15,height=15)
ggplot(melt(as.matrix(1 - dissEmission)),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]), y=factor(Var1, levels=hclust_plot$labels[hclust_plot$order]), fill=value)) +
    geom_tile(color = "white")+
    scale_fill_scico(palette = 'tokyo',direction = -1,na.value="white",limits=c(-1,1),name="Correlation") +
    theme_classic()+ # minimal theme
    scale_x_discrete(position = "bottom") +
    labs(x = NULL, y = NULL) +theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7),axis.text.y = element_text(size=7),plot.margin=unit(c(0,22,0,17),"mm"))
dev.off()


#########################################################################################################################################################################################################################################################################
##Genome browser figures
options(ucscChromosomeNames=FALSE)
trackrange=c(0,2)
for (species in specieslist$species){
    ##Prep per species
    grange = get_our_grange(species=species)
    tenogene <- read.table(paste0(genome_folder,species,"/",species,".TE.nogene.bed"))
    colnames(tenogene) <- c("chromosome","start","end","transcript","score","strand","source","type","dot")
    if (species == "Scer"){
        genetrack <- data.frame(rtracklayer::import(paste0(raw_folder,species,"_long.annot.gtf")))%>%filter(type=="CDS")%>%mutate(type="protein_coding")
    } else if (species == "Ddis") {
        genetrack <- data.frame(rtracklayer::import(paste0(raw_folder,species,"_long.annot.gtf")))%>%filter(type=="exon")%>%mutate(type="protein_coding")
        genetrack[,1] <- paste0("chr",genetrack[,1])
    } else {
        genetrack <- data.frame(rtracklayer::import(paste0(raw_folder,species,"_long.annot.gtf")))%>%filter(type=="exon")%>%mutate(type="protein_coding")
    }
    colnames(genetrack) <- c("chromosome","start","end","width","strand","source","feature","score","phase","transcript","gene")
    if (species=="Nvec" | species=="Ddis" | species=="Ppat" | species=="Gthe" | species=="Bnat" | species=="Ngru"){
        colnames(genetrack) <- c("chromosome","start","end","width","strand","source","feature","score","phase","transcript","gene","locus")
    }
    f <- paste0(chip_folder,species,"/bw/",species,"_geneExp.bw")
    if (file.exists(f)){
        tmpRNAold <- read_bigwig(f)
    }
    f <- paste0(chip_folder,species,"_all/bw/",species,"_merged.sorted.plus.bw")
    if (file.exists(f)){
        tmpRNAfwd <- read_bigwig(f)
    }
    f <- paste0(chip_folder,species,"_all/bw/",species,"_merged.sorted.minus.bw")
    if (file.exists(f)){
        tmpRNArev <- read_bigwig(f)
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9ac.log2r.input.bw")
    if (file.exists(f)){
        tmpK9ac <- read_bigwig(f)
        tmpK9ac[[4]] <- ifelse(tmpK9ac[[4]] < 0, 0, tmpK9ac[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K27ac.log2r.input.bw")
    if (file.exists(f)){
        tmpK27ac <- read_bigwig(f)
        tmpK27ac[[4]] <- ifelse(tmpK27ac[[4]] < 0, 0, tmpK27ac[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H4K16ac.log2r.input.bw")
    if (file.exists(f)){
        tmpK16ac <- read_bigwig(f)
        tmpK16ac[[4]] <- ifelse(tmpK16ac[[4]] < 0, 0, tmpK16ac[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K4me2.log2r.input.bw")
    if (file.exists(f)){
        tmpK4me2 <- read_bigwig(f)
        tmpK4me2[[4]] <- ifelse(tmpK4me2[[4]] < 0, 0, tmpK4me2[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K4me3.log2r.input.bw")
    if (file.exists(f)){
        tmpK4me3 <- read_bigwig(f)
        tmpK4me3[[4]] <- ifelse(tmpK4me3[[4]] < 0, 0, tmpK4me3[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K36me3.log2r.input.bw")
    if (file.exists(f)){
        tmpK36me3 <- read_bigwig(f)
        tmpK36me3[[4]] <- ifelse(tmpK36me3[[4]] < 0, 0, tmpK36me3[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K79me1.log2r.input.bw")
    if (file.exists(f)){
        tmpK79me1 <- read_bigwig(f)
        tmpK79me1[[4]] <- ifelse(tmpK79me1[[4]] < 0, 0, tmpK79me1[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K79me2.log2r.input.bw")
    if (file.exists(f)){
        tmpK79me2 <- read_bigwig(f)
        tmpK79me2[[4]] <- ifelse(tmpK79me2[[4]] < 0, 0, tmpK79me2[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K79me3.log2r.input.bw")
    if (file.exists(f)){
        tmpK79me3 <- read_bigwig(f)
        tmpK79me3[[4]] <- ifelse(tmpK79me3[[4]] < 0, 0, tmpK79me3[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K27me3.log2r.input.bw")
    if (file.exists(f)){
        tmpK27me3 <- read_bigwig(f)
        tmpK27me3[[4]] <- ifelse(tmpK27me3[[4]] < 0, 0, tmpK27me3[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9me1.log2r.input.bw")
    if (file.exists(f)){
        tmpK9me1 <- read_bigwig(f)
        tmpK9me1[[4]] <- ifelse(tmpK9me1[[4]] < 0, 0, tmpK9me1[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9me3.log2r.input.bw")
    if (file.exists(f)){
        tmpK9me3 <- read_bigwig(f)
        tmpK9me3[[4]] <- ifelse(tmpK9me3[[4]] < 0, 0, tmpK9me3[[4]])
    }
    f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9me2.log2r.input.bw")
    if (file.exists(f)){
        tmpK9me3 <- read_bigwig(f)
        tmpK9me3[[4]] <- ifelse(tmpK9me3[[4]] < 0, 0, tmpK9me3[[4]])
    }
    ##Prep per region
    if (species == "Ppat"){
        #Fig4 HC Ppat
        figname="RS3"
        chr="GWHFIHF00000015.1"
        start=9172276
        end=9196580
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene%>%filter(chromosome==chr),name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()  
        #Fig4 HC Ppat
        figname="RS2.1"
        chr="GWHFIHF00000016.1"
        start=13993700
        end=14044732
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene%>%filter(chromosome==chr),name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()      
        #Fig4 HC Ppat
        figname="RS2.2"
        chr="GWHFIHF00000014.1"
        start=5133344
        end=5180551
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene%>%filter(chromosome==chr),name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()   
        #Fig4 HC Ppat
        figname="RS1"
        chr="GWHFIHF00000024.1"
        start=6795223
        end=6799629
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene%>%filter(chromosome==chr),name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()  
        #Fig4 HC Ppat
        figname="RS2.3"
        chr="GWHFIHF00000014.1"
        start=8005440
        end=8041530
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene%>%filter(chromosome==chr),name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()        
    }
    if (species == "Gthe"){
        #Fig4 HC Gthe
        figname="RS1"
        chr="GUITHscaffold_137"
        start=63197
        end=72781
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        # #Fig4 HC Gthe
        # figname="EDF8"
        # chr="GUITHscaffold_99"
        # start=19472
        # end=24442
        # ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        # pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        # plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
        #                 AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
        #             type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        # dev.off()
        #Fig4 HC Gthe
        figname="RS3"
        chr="GUITHscaffold_142"
        start=84724
        end=87497
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig4 HC Gthe
        figname="RS1.2"
        chr="GUITHscaffold_325"
        start=1459
        end=6885
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig4 HC Gthe
        figname="RS4"
        chr="GUITHscaffold_13"
        start=44430
        end=70180
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
    } else if (species == "Ddis"){
        #Fig3 K79 Ddis
        figname="Fig3.a"
        chr="chr1"
        start=954395
        end=965571
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig3 K79 Ddis
        figname="Fig3.b"
        chr="chr2"
        start=946052
        end=952395
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig4 HC Ddis
        figname="RS2"
        chr="chr5"
        start=3760528
        end=3766440
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        # #Fig4 HC Ddis
        # figname="EDF8.2"
        # chr="chr1"
        # start=2067602
        # end=2078633
        # ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        # pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        # plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
        #                 AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
        #             type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        # dev.off()
        # #Fig4 HC Ddis
        # figname="EDF8.3"
        # chr="chr6"
        # start=532556
        # end=560879
        # ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        # pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        # plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
        #                 AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
        #             type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        # dev.off()
        #Fig4 HC Ddis
        figname="RS1"
        chr="chr1"
        start=98239
        end=100452  
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()      
    } else if (species == "Acas"){
        #Fig2 Acas
        figname="Fig2.1"
        chr="scaffold_1"
        start=854332
        end=883478
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig2_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK36me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#B82300","#B82300"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig3 K79 Acas
        figname="Fig3.b"
        chr="scaffold_3"
        start=817194
        end=825916
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig4 HC Acas
        figname="RS2"
        chr="scaffold_7"
        start=1319568
        end=1327585
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()     
        #Fig4 HC Acas
        figname="RS1"
        chr="scaffold_6"
        start=749366
        end=810620
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()          
    } else if (species == "Spun"){
        #Fig2.1 Spun
        figname="Fig2.1"
        chr="KQ257458"
        start=877396
        end=955197        
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig2_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK36me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#B82300","#B82300"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig2.2 Spun
        figname="Fig2.2"
        chr="KQ257458"
        start=268495
        end=343564
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig2_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK36me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#B82300","#B82300"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Spun Fig3.1
        figname="Fig3.1"
        chr="KQ257452"
        start=803080
        end=813235
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Spun Fig3.2
        figname="Fig3.2"
        chr="KQ257461"
        start=305790
        end=314364
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Spun Fig3.3
        figname="Fig3.3"
        chr="KQ257458"
        start=720683
        end=725552
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Spun Fig3.b
        figname="Fig3.b"
        chr="KQ257457"
        start=76924
        end=82495
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        ##Spun Fig4.1
        figname="RS1"
        chr="KQ257453"
        start=351228
        end=376000
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off() 
        #Spun EDF5a
        figname="EDF5a.1"
        chr="KQ257451"
        start=1198965
        end=1214722
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig5_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
    } else if (species == "Cfra"){
        #Fig2.2 Cfra
        figname="Fig2.1"
        chr="scaffold07"
        start=45651
        end=86313
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig2_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK36me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#B82300","#B82300"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Cfra Fig3.1
        figname="Fig3.1"
        chr="scaffold03"
        start=1366344
        end=1369822
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Cfra Fig3.2
        figname="Fig3.2"
        chr="scaffold07"
        start=58080
        end=71973
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Cfra EDF5a.1
        figname="EDF5a.1"
        chr="scaffold28"
        start=291418
        end=300732
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig5_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Cfra EDF5a.2
        figname="EDF5a.2"
        chr="scaffold09"
        start=532395
        end=545911
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig5_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
    } else if (species == "Ngru"){
        #Fig4 HC Ngru
        figname="RS1"
        chr="NAEGRscaffold_29"
        start=162479
        end=208602
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig4 HC Ngru
        figname="RS2"
        chr="NAEGRscaffold_23"
        start=369274
        end=389183
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()    
        #Fig3 ac Ngru
        figname="Fig3.1"
        chr="NAEGRscaffold_3"
        start=436097
        end=452976
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()     
    } else if (species == "Nvec"){
        #Fig2 Nvec
        figname="Fig2.1"
        chr="NC_064037.1"
        start=7951374
        end=8032618   
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig2_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK36me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#B82300","#B82300"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()  
        # #Fig3 Nvec ac
        # figname="Fig3.1"
        # chr="NC_064035.1"
        # start=3436822
        # end=3453909
        # ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        # pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        # plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
        #                 AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
        #             type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        # dev.off()    
        # #Fig3 Nvec ac
        # figname="Fig3.2"
        # chr="NC_064040.1"
        # start=10583679
        # end=10609552
        # ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        # pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        # plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
        #                 AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
        #             type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        # dev.off()   
        #Fig4 Nvec HC
        figname="RS1"
        chr="NC_064035.1"
        start=7389095
        end=7414284   
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()    
        #Fig4 Nvec HC
        figname="RS2"
        chr="NC_064043.1"
        start=13403228
        end=13414177   
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()    
        # #Fig4 Nvec HC
        # figname="EDF8"
        # chr="NC_064043.1"
        # start=2320747
        # end=2356575   
        # ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        # pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        # plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
        #                 AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
        #             type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        # dev.off()  
        #Fig4 Nvec HC
        figname="RS3"
        chr="NC_064037.1"
        start=12326974
        end=12334753  
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()   
        #Nvec SuppFig4.y.1
        figname="EDF5a.1"
        chr="NC_064035.1"
        start=4441330
        end=4469207  
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig5_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()  
    } else if (species == "Atha"){
        #Atha ac
        figname="Fig3.2"
        chr="Chr5"
        start=309182
        end=318055
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig3_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig4 Atha HC
        figname="RS4"
        chr="Chr1"
        start=986250
        end=997584
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig4 Atha HC
        figname="RS1"
        chr="Chr2"
        start=7252739
        end=7268897
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()  
        #Fig4 Atha HC
        figname="RS3"
        chr="Chr4"
        start=12664591
        end=12682125
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()  
        #Fig4 Atha HC
        figname="RS2"
        chr="Chr1"
        start=11903664
        end=11931837
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()      
    } else if (species == "Bnat"){
        #Fig2 Bnat
        figname="Fig2.1"
        chr="scaffold_3"
        start=1512250
        end=1554671
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"Fig2_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK36me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#B82300","#B82300"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
        #Fig4 Bnat HCscaffold_4:1,011,562-1,021,988
        figname="RS2"
        chr="scaffold_3"
        start=1543361
        end=1553885
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()   
        #Fig4 Bnat HC
        figname="RS1.1"
        chr="scaffold_13"
        start=694358
        end=702894
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off() 
        #Fig4 Bnat HC
        figname="RS3.1"
        chr="scaffold_7"
        start=1549198
        end=1552883
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()  
        #Fig4 Bnat HC
        figname="RS3.2"
        chr="scaffold_90"
        start=193476
        end=198829
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()  
        #Bnat EDF5a.1
        figname="RS1.2"
        chr="scaffold_1"
        start=726140
        end=741461
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig5_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAold%>%select("chrom","start","end","value")%>%GRanges(),name="RNAold",col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()    
        #Bnat H3K4me2 silence
        figname="RS1.3"
        chr="scaffold_41"
        start=258581
        end=273987 
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#94FA7E","#94FA7E"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27ac%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#257604","#257604"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK16ac%>%select("chrom","start","end","value")%>%GRanges(),name="H4K16ac",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#112B0A","#112B0A"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#FFE04D","#FFE04D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK4me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K4me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F46D43","#F46D43"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK36me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#B82300","#B82300"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me2",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#D448BA","#D448BA"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#70015D","#70015D"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#41B6C4","#41B6C4"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
    } else if (species == "Tthe"){
        # #Tthe Fig4.1
        # figname="EDF8.1"
        # chr="chr_135"
        # start=737672
        # end=743030
        # ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        # pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        # plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
        #                 AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
        #             type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        # dev.off()     
        # #Tthe Fig4.2
        # figname="EDF8.2"
        # chr="chr_002"
        # start=148941
        # end=155838
        # ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        # pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        # plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
        #                 DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
        #                 AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
        #                 GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
        #             type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        # dev.off()      
        #Tthe Fig4.3
        figname="RS1"
        chr="chr_004"
        start=227565
        end=239093
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig8_draft/",species,".",figname,".pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK79me1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K79me1",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#F7BEEE","#F7BEEE"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK27me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K27me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#4120A9","#4120A9"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpK9me3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K9me3",ylim=trackrange,col.histogram=NULL,fill.histogram=c("#3262AB","#3262AB"),showTitle=FALSE,showAxis=FALSE),
                        AnnotationTrack(tenogene,name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()        
    } 
}   


#########################################################################################################################################################################################################################################################################
# ##Whole chromosome plotting
# for (species in specieslist$species){
#     ##Prep per species
#     grange = get_our_grange(species=species)
#     tenogene <- read.table(paste0(genome_folder,species,"/",species,".TE.nogene.bed"))
#     colnames(tenogene) <- c("chromosome","start","end","transcript","score","strand","source","type","dot")
#     if (species == "Scer"){
#         genetrack <- data.frame(rtracklayer::import(paste0(raw_folder,species,"_long.annot.gtf")))%>%filter(type=="CDS")%>%mutate(type="protein_coding")
#     } else if (species == "Ddis") {
#         genetrack <- data.frame(rtracklayer::import(paste0(raw_folder,species,"_long.annot.gtf")))%>%filter(type=="exon")%>%mutate(type="protein_coding")
#         genetrack[,1] <- paste0("chr",genetrack[,1])
#     } else {
#         genetrack <- data.frame(rtracklayer::import(paste0(raw_folder,species,"_long.annot.gtf")))%>%filter(type=="exon")%>%mutate(type="protein_coding")
#     }
#     colnames(genetrack) <- c("chromosome","start","end","width","strand","source","feature","score","phase","transcript","gene")
#     if (species=="Nvec" | species=="Ddis" | species=="Ppat" | species=="Gthe" | species=="Bnat" | species=="Ngru"){
#         colnames(genetrack) <- c("chromosome","start","end","width","strand","source","feature","score","phase","transcript","gene","locus")
#     }
#     f <- paste0(chip_folder,species,"_all/bw/",species,"_merged.sorted.plus.bw")
#     if (file.exists(f)){
#         tmpRNAfwd <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/bw/",species,"_merged.sorted.minus.bw")
#     if (file.exists(f)){
#         tmpRNArev <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9ac.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK9ac <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K27ac.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK27ac <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H4K16ac.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK16ac <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K4me2.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK4me2 <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K4me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK4me3 <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K36me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK36me3 <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K79me1.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK79me1 <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K79me2.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK79me2 <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K79me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK79me3 <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K27me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK27me3 <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9me1.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK9me1 <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK9me3 <- read_bigwig(f)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9me2.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK9me3 <- read_bigwig(f)
#     }
#     ##Plot
#     if (species == "Atha" | species=="Ppat" | species=="Ngru"){
#         #Biggest chromsome
#         chr=names(sort(seqlengths(grange),decreasing = T)[1])
#         pdf(paste0(plot_folder,"misc/",species,".",chr,".pdf"),width=15,height=15)
#             par(mfrow=c(1,1),mai=c(0.1,0.1,0.2,0.1))
#             par(fig=c(0.1,1,0.8,0.9))
#             plot(density(subset(genetrack,chromosome==chr)$end),axes=F,main="Gene density",xlab=NA,ylab=NA,col="black")
#             par(fig=c(0.1,1,0.65,0.75), new=TRUE)
#             plot(density(subset(tenogene,chromosome==chr)$end),axes=F,main="TE density",xlab=NA,ylab=NA,col="#808080")
#             par(fig=c(0.1,1,0.1,0.6), new=TRUE)
#             plot(NULL,xlim=c(1,seqlengths(grange)[chr]),ylim=c(-2,2), col="white")
#                 abline(h = 0,col="grey")
#                 lines(smooth.spline(subset(tmpK16ac, chrom==chr)[,c("end","value")]), col="#112B0A")
#                 lines(smooth.spline(subset(tmpK9me1, chrom==chr)[,c("end","value")]), col="#41B6C4")
#                 lines(smooth.spline(subset(tmpK9me3, chrom==chr)[,c("end","value")]), col="#3262AB")
#                 lines(smooth.spline(subset(tmpK27me3, chrom==chr)[,c("end","value")]), col="#4120A9")
#             par(fig=c(0,1,0,1),mai=c(1.02,0.82,0.82,0.42),new=FALSE)
#         dev.off()
#     } else if (species=="Cfra" | species=="Scer"){
#         #Biggest chromsome
#         chr=names(sort(seqlengths(grange),decreasing = T)[1])
#         pdf(paste0(plot_folder,"misc/",species,".",chr,".pdf"),width=15,height=15)
#             par(mfrow=c(1,1),mai=c(0.1,0.1,0.2,0.1))
#             par(fig=c(0.1,1,0.8,0.9))
#             plot(density(subset(genetrack,chromosome==chr)$end),axes=F,main="Gene density",xlab=NA,ylab=NA,col="black")
#             par(fig=c(0.1,1,0.65,0.75), new=TRUE)
#             plot(density(subset(tenogene,chromosome==chr)$end),axes=F,main="TE density",xlab=NA,ylab=NA,col="#808080")
#             par(fig=c(0.1,1,0.1,0.6), new=TRUE)
#             plot(NULL,xlim=c(1,seqlengths(grange)[chr]),ylim=c(-2,2), col="white")
#                 abline(h = 0,col="grey")
#                 lines(smooth.spline(subset(tmpK16ac, chrom==chr)[,c("end","value")]), col="#112B0A")
#                 lines(smooth.spline(subset(tmpK79me1, chrom==chr)[,c("end","value")]), col="#F7BEEE")
#                 lines(smooth.spline(subset(tmpK79me2, chrom==chr)[,c("end","value")]), col="#D448BA")
#                 lines(smooth.spline(subset(tmpK79me3, chrom==chr)[,c("end","value")]), col="#70015D")
#             par(fig=c(0,1,0,1),mai=c(1.02,0.82,0.82,0.42),new=FALSE)
#         dev.off()
#     } else if (species=="Gthe"){
#         #Biggest chromsome
#         chr=names(sort(seqlengths(grange),decreasing = T)[1])
#         pdf(paste0(plot_folder,"misc/",species,".",chr,".pdf"),width=15,height=15)
#             par(mfrow=c(1,1),mai=c(0.1,0.1,0.2,0.1))
#             par(fig=c(0.1,1,0.8,0.9))
#             plot(density(subset(genetrack,chromosome==chr)$end),axes=F,main="Gene density",xlab=NA,ylab=NA,col="black")
#             par(fig=c(0.1,1,0.65,0.75), new=TRUE)
#             plot(density(subset(tenogene,chromosome==chr)$end),axes=F,main="TE density",xlab=NA,ylab=NA,col="#808080")
#             par(fig=c(0.1,1,0.1,0.6), new=TRUE)
#             plot(NULL,xlim=c(1,seqlengths(grange)[chr]),ylim=c(-2,2), col="white")
#                 abline(h = 0,col="grey")
#                 lines(smooth.spline(subset(tmpK16ac, chrom==chr)[,c("end","value")]), col="#112B0A")
#                 lines(smooth.spline(subset(tmpK79me1, chrom==chr)[,c("end","value")]), col="#F7BEEE")
#                 lines(smooth.spline(subset(tmpK9me1, chrom==chr)[,c("end","value")]), col="#41B6C4")
#                 lines(smooth.spline(subset(tmpK27me3, chrom==chr)[,c("end","value")]), col="#4120A9")
#             par(fig=c(0,1,0,1),mai=c(1.02,0.82,0.82,0.42),new=FALSE)
#         dev.off()
#     } else if (species=="Tthe"){
#         #Biggest chromsome
#         chr=names(sort(seqlengths(grange),decreasing = T)[1])
#         pdf(paste0(plot_folder,"misc/",species,".",chr,".pdf"),width=15,height=15)
#             par(mfrow=c(1,1),mai=c(0.1,0.1,0.2,0.1))
#             par(fig=c(0.1,1,0.8,0.9))
#             plot(density(subset(genetrack,chromosome==chr)$end),axes=F,main="Gene density",xlab=NA,ylab=NA,col="black")
#             par(fig=c(0.1,1,0.1,0.6), new=TRUE)
#             plot(NULL,xlim=c(1,seqlengths(grange)[chr]),ylim=c(-2,2), col="white")
#                 abline(h = 0,col="grey")
#                 lines(smooth.spline(subset(tmpK16ac, chrom==chr)[,c("end","value")]), col="#112B0A")
#                 lines(smooth.spline(subset(tmpK79me1, chrom==chr)[,c("end","value")]), col="#F7BEEE")
#                 lines(smooth.spline(subset(tmpK9me3, chrom==chr)[,c("end","value")]), col="#3262AB")
#                 lines(smooth.spline(subset(tmpK27me3, chrom==chr)[,c("end","value")]), col="#4120A9")
#             par(fig=c(0,1,0,1),mai=c(1.02,0.82,0.82,0.42),new=FALSE)
#         dev.off()
#     } else {
#         #Biggest chromsome
#         chr=names(sort(seqlengths(grange),decreasing = T)[1])
#         pdf(paste0(plot_folder,"misc/",species,".",chr,".pdf"),width=15,height=15)
#             par(mfrow=c(1,1),mai=c(0.1,0.1,0.2,0.1))
#             par(fig=c(0.1,1,0.8,0.9))
#             plot(density(subset(genetrack,chromosome==chr)$end),axes=F,main="Gene density",xlab=NA,ylab=NA,col="black")
#             par(fig=c(0.1,1,0.65,0.75), new=TRUE)
#             plot(density(subset(tenogene,chromosome==chr)$end),axes=F,main="TE density",xlab=NA,ylab=NA,col="#808080")
#             par(fig=c(0.1,1,0.1,0.6), new=TRUE)
#             plot(NULL,xlim=c(1,seqlengths(grange)[chr]),ylim=c(-2,2), col="white")
#                 abline(h = 0,col="grey")
#                 lines(smooth.spline(subset(tmpK16ac, chrom==chr)[,c("end","value")]), col="#112B0A")
#                 lines(smooth.spline(subset(tmpK79me1, chrom==chr)[,c("end","value")]), col="#F7BEEE")
#                 lines(smooth.spline(subset(tmpK79me2, chrom==chr)[,c("end","value")]), col="#D448BA")
#                 lines(smooth.spline(subset(tmpK79me3, chrom==chr)[,c("end","value")]), col="#70015D")
#                 lines(smooth.spline(subset(tmpK9me1, chrom==chr)[,c("end","value")]), col="#41B6C4")
#                 lines(smooth.spline(subset(tmpK9me3, chrom==chr)[,c("end","value")]), col="#3262AB")
#                 lines(smooth.spline(subset(tmpK27me3, chrom==chr)[,c("end","value")]), col="#4120A9")
#             par(fig=c(0,1,0,1),mai=c(1.02,0.82,0.82,0.42),new=FALSE)
#         dev.off()
#     }
# }





# pp <- getDefaultPlotParams(plot.type=4)
# pp$leftmargin <- 0.15
# pp$topmargin <- 15
# pp$bottommargin <- 15
# pp$ideogramheight <- 5
# pp$data1inmargin <- 10
# pp$data1outmargin <- 0

# options(ucscChromosomeNames=FALSE)
# binsize=250
# for (species in specieslist$species){
#     ##Prep per species
#     grange = get_our_grange(species=species)
#     # f <- paste0(chip_folder,species,"_all/bw/",species,"_merged.sorted.plus.bw")
#     # if (file.exists(f)){
#     #     tmpRNAfwd <- read_bigwig(f)
#     # }
#     # f <- paste0(chip_folder,species,"_all/bw/",species,"_merged.sorted.minus.bw")
#     # if (file.exists(f)){
#     #     tmpRNArev <- read_bigwig(f)
#     # }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K27ac.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK27ac <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H4K16ac.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK16ac <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K4me2.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK4me2 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K4me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK4me3 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K36me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK36me3 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K79me1.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK79me1 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K79me2.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK79me2 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K79me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK79me3 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K27me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK27me3 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9me1.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK9me1 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9me3.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK9me3 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
#     f <- paste0(chip_folder,species,"_all/log2r/",species,"_H3K9me2.log2r.input.bw")
#     if (file.exists(f)){
#         tmpK9me3 <- unlist(summary(BigWigFile(f),which=seqinfo(BigWigFile(f))[names(sort(seqlengths(BigWigFile(f)),decreasing = T)[1:2])], size = c(binsize,binsize)))%>%mutate(value=score)
#     }
    
#     ##Plot
#     if (species == "Atha" | species=="Ppat" | species=="Ngru"){
#         histone.marks <- c(H3K9me3=tmpK9me3,H3K9me1=tmpK9me1,H3K27me3=tmpK27me3,H4K16ac=tmpK16ac)
#     } else if (species=="Cfra" | species=="Scer"){
#         histone.marks <- c(H3K79me3=tmpK79me3,H3K79me2=tmpK79me2,H3K79me1=tmpK79me1,H4K16ac=tmpK16ac)
#     } else if (species=="Gthe"){
#         histone.marks <- c(H3K9me1=tmpK9me1,H3K27me3=tmpK27me3,H3K79me1=tmpK79me1,H4K16ac=tmpK16ac)
#     } else if (species=="Tthe"){
#         histone.marks <- c(H3K9me3=tmpK9me3,H3K27me3=tmpK27me3,H3K79me1=tmpK79me1,H4K16ac=tmpK16ac)
#     } else {
#         histone.marks <- c(H3K9me3=tmpK9me3,H3K9me1=tmpK9me1,H3K27me3=tmpK27me3,H3K79me3=tmpK79me3,H3K79me2=tmpK79me2,H3K79me1=tmpK79me1,H4K16ac=tmpK16ac)
#     }
#     custom.genome <- toGRanges(read.table(paste0(genome_folder,species,"/",species,".circlize.txt"),col.names = c("chr","start","end")))
#     pp$data1min <- min(mcols(unlist(as(histone.marks, "GRangesList")))$score,na.rm=T)
#     pp$data1max <- max(mcols(unlist(as(histone.marks, "GRangesList")))$score,na.rm=T)
#     pdf(paste0(plot_folder,"misc/",species,".chromosomeCov.pdf"),width=15,height=15)
#         ##Plot chromsome
#         kp <- plotKaryotype(genome=custom.genome, plot.type=4, ideogram.plotter = NULL, labels.plotter = NULL,plot.params = pp,chromosomes=names(sort(seqlengths(grange),decreasing = T)[1:2]))
#         kpAddCytobandsAsLine(kp)
#         kpAddChromosomeNames(kp, srt=45)
#         ##Set margins
#         total.tracks <- length(histone.marks)+2
#         out.at <- autotrack(1:length(histone.marks), total.tracks, margin = 0.1, r0=min(mcols(unlist(as(histone.marks, "GRangesList")))$score,na.rm=T), r1=max(mcols(unlist(as(histone.marks, "GRangesList")))$score,na.rm=T))
#         ##Plot genes
#         at <- autotrack(1, total.tracks, r0=out.at$r0, r1=out.at$r1, margin = 0.1)
#         kpPlotDensity(kp, grange[mcols(grange)$type == "protein_coding"], r0=at$r0, r1=at$r1,window.size = mean(sort(seqlengths(grange),decreasing = T)[1:2])/250)
#         kpAddLabels(kp, labels = "Genes", r0=at$r0, r1=at$r1, cex=1.6, label.margin = 0.035)
#         ##Plot TEs
#         at <- autotrack(2, total.tracks, r0=out.at$r0, r1=out.at$r1, margin = 0.1)
#         kpPlotDensity(kp, grange[mcols(grange)$source == "EDTA"], r0=at$r0, r1=at$r1,window.size = mean(sort(seqlengths(grange),decreasing = T)[1:2])/250)
#         kpAddLabels(kp, labels = "TEs", r0=at$r0, r1=at$r1, cex=1.6, label.margin = 0.035)
#         ##Plot marks
#         for(i in seq_len(length(histone.marks))) {
#             at <- autotrack(i+2, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
#             kp <- kpLines(kp, data=unlist(as(histone.marks[i], "GRangesList")),col=cols[names(histone.marks)[i]], r0=at$r0, r1=at$r1,lwd=3)
#             kpAxis(kp,r0=at$r0, r1=at$r1, cex=1.6)
#             kpAddLabels(kp, labels = names(histone.marks)[i], r0=at$r0, r1=at$r1, 
#                         cex=2.2, label.margin = 0.035)
#             kpAbline(kp,r0=at$r0,r1=at$r1, h=0,col="grey")
#         }
#     dev.off()
# }

#########################################################################################################################################################################################################################################################################
##Bnat H3K4me2 peak investigation
# ##R
# clust2 <- read.table("no_backup/smontgomery/chip/Bnat_all/macs2/Bnat_H3K4me2_peaks.narrowPeak.filt100.regions.k10.RNAcounts.cluster_2.txt")
# clust9 <- read.table("no_backup/smontgomery/chip/Bnat_all/macs2/Bnat_H3K4me2_peaks.narrowPeak.filt100.regions.k10.RNAcounts.cluster_9.txt")
# clust2$up <- as.numeric(gsub(",.*","",clust2$V4))
# clust2$down <- as.numeric(gsub(".*,","",clust2$V4))
# clust9$down <- as.numeric(gsub(",.*","",clust9$V4))
# clust9$up <- as.numeric(gsub(".*,","",clust9$V4))

# rbind(clust2,clust9) %>%
#     select(up,down) %>%melt()%>%
#     mutate(paired = rep(1:(n()/2),2))%>%
#     ggplot(aes(variable, value, fill = variable)) +
#     geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA) +
#     geom_line(aes(group = paired), color = "grey70", alpha = 0.7, position = position_dodge(0.2)) +
#     geom_point(aes(group = paired), shape = 21, size = 2, position = position_dodge(0.2)) +
#     # Add the summary line for the mean
#     stat_summary(
#         aes(group = 1), # Group all points together for the summary
#         fun = "mean",
#         geom = "line",
#         color = "red",
#         linewidth = 1.2,
#         linetype = "dashed"
#     ) +
#     stat_compare_means(
#         method = "t.test",
#         paired = TRUE,
#         label = "p.format",
#         size=5,
#         label.x=1.5
#     ) +
#     annotate("text",x=1.5,y=35000,label=paste0("Cohen's d = ",round(cohen.d(rbind(clust2,clust9)$up, rbind(clust2,clust9)$down)$estimate,3)),size=5) +
#     scale_x_discrete(labels=c("Upstream","Downstream")) +
#     scale_y_continuous(name = "RNAcounts") +
#     theme_minimal() +
#     theme(
#         axis.text.x = element_text(color = "black",size = 12),
#         axis.text.y = element_text(color = "black",size = 12),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(color = "black",size = 12),
#         legend.position="none")



##FOR REAL NOW
mindis=200
maxdis=2000
peaksummarytable=data.frame()
orisummarytable=data.frame()
geneRNAsummarytable=data.frame()
peakstats=data.frame()
for (sp in specieslist$species){
    RNAcountsfile <- read.table(paste0("/home/smontgomery/cluster/smontgomery/genomes/",sp,"/",sp,"_expr_sorted.tsv"),header=TRUE)
    if (sp == "Cfra"){
        RNAcountsfile$gene_id <- gsub("Cfra_CFRG0","Cfra_CFRCFRG0",RNAcountsfile$gene_id)
    }
    ##Orientation of gene pairs
    grange = get_our_grange(species=sp)%>%filter(type=="protein_coding")
    intergen=add_intergenic(grange)%>%mutate(size=width(.))%>%filter(type=="intergenic" & size > mindis & size < maxdis)
    genepairs=bed_closest(as.data.frame(intergen)%>%mutate(start=start-1)%>%mutate(chrom=seqnames),
        as.data.frame(grange)%>%mutate(start=start-1)%>%mutate(chrom=seqnames))%>%
        filter(!is.na(.dist))%>%
        group_by(seqnames.x,start.x,end.x)%>%
        summarise(genepair=paste(gene_id.y, collapse = "-pair-"))%>%
        mutate(gene1=strsplit(genepair,split = "-pair-")[[1]][1])%>%
        mutate(gene2=strsplit(genepair,split = "-pair-")[[1]][2])%>%
        filter(!is.na(gene2))%>%ungroup()
    genepairs[,"gene1strand"]=grange[match(genepairs$gene1,mcols(grange)$gene_id)]%>%as.data.frame()%>%select(strand)%>%as.vector()
    genepairs[,"gene2strand"]=grange[match(genepairs$gene2,mcols(grange)$gene_id)]%>%as.data.frame()%>%select(strand)
    genepairs$gene1RNA=RNAcountsfile[match(genepairs$gene1,RNAcountsfile$gene_id),"RNAcounts"]
    genepairs$gene2RNA=RNAcountsfile[match(genepairs$gene2,RNAcountsfile$gene_id),"RNAcounts"]
    genepairs$orientation=paste0(genepairs$gene1strand,genepairs$gene2strand)
    ##Distribution of peaks in the genome
    grange = get_our_grange(species=sp)
    grange_w_bins=c(grange%>%filter(type=="protein_coding")%>%mutate(type="geneBody"),grange%>%filter(source=="EDTA")%>%mutate(type="TE"))
    grange_w_promoter=c(grange_w_bins,grange_w_bins%>%filter(type=="geneBody")%>%promoters(upstream=500,downstream=0)%>%mutate(type="promoter")%>%trim())
    grange_wo_overlaps=merge_overlap_feature(grange_w_promoter)%>%mutate(type=gsub("promotergeneBody","promoter",type))
    grange_w_intergenic=add_intergenic(grange_wo_overlaps)%>%mutate(type=gsub("interte","intergenic",type))

    for (mark in c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3")){
        tmppeaks <- read.table(paste0("no_backup/smontgomery/chip/",sp,"_all/macs2/",sp,"_",mark,"_peaks.narrowPeak.filt100"))
        peakfeature=bed_intersect(tmppeaks%>%mutate(chrom=V1,start=V2,end=V3),as.data.frame(grange_w_intergenic)%>%mutate(start=start-1)%>%mutate(chrom=seqnames))%>%
            group_by(chrom,start.x,end.x,V4.x)%>%
            filter(.overlap==max(.overlap))%>%ungroup()
        ##Distribution of peaks in the genome
        peaksumtmptable=peakfeature%>%count(type.y)%>%mutate(pct=prop.table(n)*100)%>%
            mutate(species=sp,name=mark)
        peaksummarytable <- rbind(peaksummarytable,data.frame(peaksumtmptable))
        ##Orientation of gene pairs
        intergenpeak=peakfeature%>%filter(type.y=="intergenic")
        tmptable=bed_intersect(genepairs%>%rename(chrom=seqnames.x,start=start.x,end=end.x),
                    intergenpeak%>%rename(start=start.x,end=end.x))%>%
                    group_by(V4.x.y)%>%filter(start.x==min(start.x))%>%ungroup()
        orisumtmptable=tmptable%>%count(orientation.x)%>%mutate(pct=prop.table(n)*100)%>%
            mutate(species=sp,name=mark)
        orisummarytable <- rbind(orisummarytable,data.frame(orisumtmptable))
        geneRNAtmptable=tmptable%>%
            mutate(lfc=log2((gene1RNA.x+1)/(gene2RNA.x+1)))%>%
            mutate(lfc=if_else(orientation.x=="--",-lfc,lfc))%>%
            mutate(species=sp,name=mark)
        peakfeature[match(geneRNAtmptable$V4.x.y,peakfeature$V4.x),"orientation"] <- geneRNAtmptable$orientation.x
        peakfeature[match(geneRNAtmptable$V4.x.y,peakfeature$V4.x),"lfc"] <- geneRNAtmptable$lfc
        peakfeature=peakfeature%>%select(V4.x,type.y,orientation,lfc)%>%mutate(species=sp,name=mark)
        peakstats <- rbind(peakstats,peakfeature)
        geneRNAtmptable=geneRNAtmptable%>%
            filter(orientation.x=="--" | orientation.x=="++")%>%
            select(species,name,lfc)
        geneRNAsummarytable <- rbind(geneRNAsummarytable,data.frame(geneRNAtmptable))
    }
    peaksumtmptable=grange_w_intergenic%>%as.data.frame()%>%group_by(type)%>%summarize(n=sum(width))%>%
        mutate(pct=prop.table(n)*100)%>%rename(type.y=type)%>%mutate(species=sp,name="genomic")
    peaksummarytable <- rbind(peaksummarytable,data.frame(peaksumtmptable))
    orisumtmptable=genepairs%>%count(orientation)%>%mutate(pct=prop.table(n)*100)%>%
        mutate(species=sp,name="genomic")%>%rename(orientation.x=orientation)
    orisummarytable <- rbind(orisummarytable,data.frame(orisumtmptable))
    geneRNAtmptable=genepairs%>%select(gene1RNA,gene2RNA,orientation)%>%
        mutate(lfc=log2((gene1RNA+1)/(gene2RNA+1)))%>%
        filter(orientation=="--" | orientation=="++")%>%
        mutate(lfc=if_else(orientation=="--",-lfc,lfc))%>%
        select(lfc)%>%mutate(species=sp,name="genomic")
        geneRNAsummarytable <- rbind(geneRNAsummarytable,data.frame(geneRNAtmptable))
}
##Distribution of peaks in the genome
tmpdata=peaksummarytable%>%group_by(species,name)%>%summarize(total=sum(n))%>%mutate(total=ifelse(name=="genomic",NA,total))
pdf(paste0(plot_folder,"misc/PeakDistributionSummary.pdf"),width=15,height=15)
ggplot(peaksummarytable,aes(x=factor(name,levels=c("genomic","H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3")),y=pct,fill=factor(type.y,levels=c("TE","geneBody","promoter","intergenic")))) +
    geom_bar(stat = "identity") +
    geom_text(data=tmpdata,aes(x=name,label = total, y=101,fill=NULL)) +
    facet_wrap(vars(species),nrow=3) +
    # scale_fill_manual(values = cols) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 12,hjust = 1,angle = 45),
          axis.text.y = element_text(color = "black",size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(colour = "black",size = 12),
          strip.text = element_text(color="black",size=12),
          strip.background = element_rect(linetype = "blank"))
dev.off()
##Orientation of gene pairs
tmpdata=orisummarytable%>%group_by(species,name)%>%summarize(total=sum(n))
pdf(paste0(plot_folder,"misc/GeneOrientationSummary.pdf"),width=15,height=15)
ggplot(orisummarytable,aes(x=factor(name,levels=c("genomic","H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3")),y=pct,fill=factor(orientation.x,levels=c("+-","-+","--","++")))) +
    geom_bar(stat = "identity") +
    geom_text(data=tmpdata,aes(x=name,label = total, y=101,fill=NULL)) +
    facet_wrap(vars(species),nrow=3) +
    scale_fill_manual(values=c("+-"="#257604","-+"="#F46D43","--"="#D448BA","++"="#70015D"),labels=c("+-"="Convergent","-+"="Divergent","--"="Head-to-tail","++"="Tail-to-head")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 12,hjust = 1,angle = 45),
          axis.text.y = element_text(color = "black",size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(colour = "black",size = 12),
          strip.text = element_text(color="black",size=12),
          strip.background = element_rect(linetype = "blank"))
dev.off()
##LFC of same-strand genes
tmpdata=geneRNAsummarytable%>%group_by(species,name)%>%count()
pdf(paste0(plot_folder,"misc/LFCbyOrientationSummary.pdf"),width=15,height=15)
ggplot(geneRNAsummarytable,aes(x=factor(name,levels=c("genomic","H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3")),y=lfc)) +
    geom_hline(yintercept=0,col="grey") +
    geom_boxplot() +
    geom_text(data=tmpdata,aes(x=name,label = n, y=19,fill=NULL)) +
    facet_wrap(vars(species),nrow=3) +
    stat_compare_means(method="t.test",ref.group = "genomic",label = "p.signif") +
    scale_y_continuous(name="LFC(Up/Down)",limits=c(-20,20),breaks=c(-20,-10,0,10,20)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 12,hjust = 1,angle = 45),
          axis.text.y = element_text(color = "black",size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black",size = 12),
          legend.title = element_blank(),
          legend.text = element_text(colour = "black",size = 12),
          strip.text = element_text(color="black",size=12),
          strip.background = element_rect(linetype = "blank"))
dev.off()



chi_sq_test_loop=function(teregions=teregions){
    tmptable=orisummarytable%>%filter(species==sp & (name==mark | name=="genomic"))
    tmptable$n <- tmptable$n - c(0,0,0,0,subset(tmptable, name==mark)$n)
    df=tmptable%>%mutate(orientation=if_else(orientation.x=="++" | orientation.x=="--","SameStrand",if_else(orientation.x=="+-","Convergent",ifelse(orientation.x=="-+","Divergent",NA))))%>%
        group_by(orientation,name)%>%
        summarize(total=sum(n))%>%
        pivot_wider(names_from = orientation,values_from = total)%>%
        column_to_rownames(var="name")

    df=as.matrix(table(teregions[,c("V13","order_id")]))
    mat_ps_chisq=as.data.frame(matrix(ncol=ncol(df),nrow=nrow(df)-1));colnames(mat_ps_chisq)=colnames(df);rownames(mat_ps_chisq)=rownames(head(df,-1))

    for(marky in rownames(df)){
        for (direc in colnames(df)){ 
            niche_values=c(df[marky,direc],sum(df[marky,which(colnames(df)!=direc)]))
            other_sn_values=c(sum(df[which(rownames(df)!=marky),direc]),sum(df[which(rownames(df)!=marky),which(colnames(df)!=direc)]))
            matrix=rbind(niche_values,other_sn_values);colnames(matrix)=c(direc,"other direc");rownames(matrix)=c(marky,"other marky")      
            #print(matrix)
            mat_ps_chisq[marky,direc]=fisher.test(matrix,alternative="greater")$p.value
        }
    }   
    mat_ps_chisq=as.matrix(mat_ps_chisq)
    mat_ps_chisq[is.nan(mat_ps_chisq)] <- 1
    a1 = mat_ps_chisq %>% melt() 
    colnames(a1)=c("marky","direc","qval")
    a= a1 %>% mutate(qval = p.adjust(qval,method = "BH")) %>% dcast(marky ~ direc, value.var='qval')
    mat_ps_chisq_qv=a[,-1];rownames(mat_ps_chisq_qv)=a[,1]
    return(mat_ps_chisq_qv)
}



##Mess around with TE PFAMs
totCatRegions <- data.frame()
for (species in specieslist$species){
    if (species == "Gthe"){
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.reorder.txt"))
    } else {
        generegions <- read.table(paste0(chip_folder,species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.reorder.txt"))
    }
    pfam_table = read_tsv("/home/smontgomery/OneDrive/Proteomics/MassSpec/annotations/pfam_GOs_table.tsv")
    pfam_TEs = read.table("/home/smontgomery/Downloads/TE_Pfam_domains")%>%mutate(category="TE")%>%rename("pfam_family"=V1)

    pfam_annot = read_tsv(paste0("/home/smontgomery/OneDrive/Proteomics/MassSpec/SPACE/",species,"/annotations/",species,"_long.pep.annotations.tsv"), col_names = c("transcript_id", "gene_id", "pfam_arch"))
    pfam_scan = read.table(paste0("/home/smontgomery/OneDrive/Proteomics/MassSpec/SPACE/",species,"/annotations/",species,"_long.pep.pfamscan.tsv"),skip = 1,col.names = c("transcript_id","alignment_start","alignment_end","envelope_start","envelope_end","pfam_id","pfam_family","type", "hmm_start","hmm_end","hmm_length","bit_score","E_value","significance","clan"))%>%as_tibble()%>%mutate(pfam_id=gsub("\\..*","",pfam_id))
    if (species == "Spun"){
        prot2gene <- read_tsv(paste0("/home/smontgomery/OneDrive/Proteomics/MassSpec/SPACE/",species,"/annotations/",species,"_NCBIprot2gene.tsv"))%>%
            mutate(transcript_id=paste0(species,"_",gsub("\\..*","",`Protein accession`)),prot_id=`Locus tag`)%>%
            select(transcript_id,prot_id)%>%mutate(prot_id=paste0(species,"_",prot_id))

        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            left_join(prot2gene,by="transcript_id",relationship = "many-to-many")%>%
            mutate(prot_id=ifelse(is.na(prot_id),transcript_id,prot_id))%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))
    } else if (species == "Acas"){
        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            mutate(prot_id=transcript_id)%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))%>%
            mutate(prot_id=gsub(".mRNA.*","",prot_id))
    } else if (species == "Cfra"){
        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            mutate(prot_id=transcript_id)%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))%>%
            mutate(prot_id=gsub("Cfra_CFRG0","Cfra_CFRCFRG0",gsub("T1","",prot_id)))
    } else if (species == "Scer"){
        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            mutate(prot_id=transcript_id)%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))%>%
            mutate(prot_id=gsub("_mRNA","",gsub("_id.*","",prot_id)))
    } else if (species == "Ppat"){
        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            mutate(prot_id=transcript_id)%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))%>%
            mutate(prot_id=gsub("-.*","",prot_id))
    } else if (species == "Tthe"){
        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            mutate(prot_id=transcript_id)%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))%>%
            mutate(prot_id=gsub(".t1.*","",prot_id))
    } else if (species == "Ddis"){
        prot2gene <- read_tsv(paste0("/home/smontgomery/OneDrive/Proteomics/MassSpec/SPACE/",species,"/annotations/",species,"_gene_length.tsv"),col_names = c("transcript_id", "prot_id", "gene_length"))%>%
            mutate(prot_id=paste0(species,"_",prot_id))
        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            left_join(prot2gene,by="transcript_id",relationship = "many-to-many")%>%
            mutate(prot_id=ifelse(is.na(prot_id),transcript_id,prot_id))%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))%>%
            mutate(prot_id=gsub("_mRNA","",prot_id))
    } else if (species == "Bnat"){
        prot2gene <- read_tsv(paste0("/home/smontgomery/OneDrive/Proteomics/MassSpec/SPACE/",species,"/annotations/",species,"_expr_sorted.tsv"))%>%
            rename("prot_id"=gene_id)
        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            left_join(prot2gene,by="transcript_id",relationship = "many-to-many")%>%
            mutate(prot_id=ifelse(is.na(prot_id),transcript_id,prot_id))%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))%>%
            mutate(prot_id=gsub("_mRNA","",prot_id))
    } else if (species == "Ngru" | species=="Gthe"){
        prot2gene <- data.frame(rtracklayer::import(paste0(raw_folder,species,"_long.annot.gtf")))%>%
            filter(type=="transcript")%>%rename("prot_id"=gene_id)
        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            left_join(prot2gene,by="transcript_id",relationship = "many-to-many")%>%
            mutate(prot_id=ifelse(is.na(prot_id),transcript_id,prot_id))%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))%>%
            mutate(prot_id=gsub("_mRNA","",prot_id))
    } else {
        sp_annot_final=pfam_annot %>% 
            select(-pfam_arch) %>% 
            left_join(pfam_scan) %>% 
            select(transcript_id,gene_id,pfam_id) %>%
            left_join(pfam_table,by="pfam_id") %>% 
            select(-clan_id,-clan_name) %>% 
            left_join(pfam_TEs,by="pfam_family")%>%
            mutate(prot_id=transcript_id)%>%
            mutate(category=ifelse(is.na(category),ifelse(is.na(pfam_id),"Domainless","nonTE"),category))%>%
            group_by(prot_id)%>%summarize(category=paste(unique(category), collapse = ";"))%>%
            mutate(category=gsub(";nonTE","",gsub("nonTE;","",gsub(";Domainless","",gsub("Domainless;","",category)))))
    }

    generegions$category <- sp_annot_final$category[match(generegions$V4,sp_annot_final$prot_id)]
    generegions$species <- species

    totCatRegions <- rbind(totCatRegions,generegions)
}
print(totCatRegions%>%count(V13, category,species)%>%group_by(V13,species)%>%mutate(pct= prop.table(n) * 100)%>%ggplot(aes(x=V13,y=pct,fill=factor(category))) +
          geom_bar(stat = "identity") +
          facet_wrap(~species) +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black",size = 12,hjust = 1,angle = 45),
                axis.text.y = element_text(color = "black",size = 12),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.title = element_blank(),
                legend.text = element_text(colour = "black",size = 12),
                strip.text = element_text(color="black",size=12),
                strip.background = element_rect(linetype = "blank")))