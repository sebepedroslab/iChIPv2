# install.packages("idr")
# install.packages("remotes")
# install.packages("~/Downloads/spp_1.16.0.tar.gz", repos = NULL, type = "source")
# remotes::install_github("imbforge/encodeChIPqc")
library(readr)
library(GenomicRanges)
library(idr)
library(ggplot2)
library(readxl)
library(encodeChIPqc)
library(spp)
library("tidyverse")
library("rtracklayer")
library(MatrixGenerics)
library("Gviz")
library("plyranges")
library("valr")

##Set functions
crm=c("Pt","Mt","Mito","ChrM","ChrmC","chrmt")
get_our_grange=function(crm,species){
            grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".TEgene.all.gtf"))
            grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
            grange$size <- width(grange)
            grange = grange[!(seqnames(grange) %in% c(crm))]
        }
##Set directories
raw_folder="/home/smontgomery/cluster/xgraubove/genomes/data/"
genome_folder="/home/smontgomery/cluster/smontgomery/genomes/"
chip_folder="/home/smontgomery/no_backup/smontgomery/chip/"
plot_folder="/home/smontgomery/OneDrive/iChIPv2_SAM_plots/00Revision/"

# ##Sample table prep for iChIP2 pipeline
# for (species in c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru")){
#     write.table(data.frame(sample_id=gsub("_S.*fastq.gz","",list.files(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/"),pattern = paste0(species,".*R1"))),
#                             fastq_1= list.files(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/"),pattern = paste0(species,".*R1")),
#                             fastq_2= list.files(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/"),pattern = paste0(species,".*R2")),
#                             sps_id=species,
#                             is.control=ifelse(sapply(strsplit(list.files(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/"),pattern = paste0(species,".*R1")),"_"),"[[",2)=="H3" | 
#                                 sapply(strsplit(list.files(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/"),pattern = paste0(species,".*R1")),"_"),"[[",2)=="input","Y","N"),
#                             control=ifelse(sapply(strsplit(list.files(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/"),pattern = paste0(species,".*R1")),"_"),"[[",2)=="H3" | 
#                                 sapply(strsplit(list.files(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/"),pattern = paste0(species,".*R1")),"_"),"[[",2)=="input","",
#                                     gsub("_S.*fastq.gz","",list.files(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/"),pattern = paste0(species,".*R1")))[length(gsub("_S.*fastq.gz","",list.files(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/"),pattern = paste0(species,".*R1"))))])),
#         file=paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",species,"/sample_table.tsv"),quote = F,col.names = T,row.names = F,na = "",sep = "\t")
# }

##Run on command line to prepare files
# mkdir /home/smontgomery/data/chip/ssp/
# for file in /home/smontgomery/no_backup/smontgomery/chip/*_all/bam_files/*/*.bam; do NAME=`basename ${file} | cut -d '_' -f1`; FILENAME=`basename ${file} | sed 's/.q30.rmdup.bam//'`; /opt/R/4.5.0/bin/Rscript /home/smontgomery/Documents/phantompeakqualtools/run_spp.R -c=${file} -out=/home/smontgomery/data/chip/ssp/${FILENAME}.out; done
# cat /home/smontgomery/data/chip/ssp/* > /home/smontgomery/data/chip/ssp_all.tsv
# cat /home/smontgomery/no_backup/smontgomery/chip/*_all/QC/gathered_QC_FRiP.tsv | grep -v sampleName | cat <(cat /home/smontgomery/no_backup/smontgomery/chip/*_all/QC/gathered_QC_FRiP.tsv | head -1) - > /home/smontgomery/data/chip/gathered_QC_all.tsv

##Read in tables
# ##Run once to prep
# metadata <- read.table("/home/smontgomery/data/chip/metadata.tsv",header = T,sep = "\t")
# metadata=metadata%>%mutate(present=ifelse(((species=="Atha" | species=="Ppat" | species=="Ngru") & (antibody=="H3K79me1" | antibody=="H3K79me2" | antibody=="H3K79me3")) |
#                                         ((species == "Scer" | species == "Cfra") & (antibody=="H3K27me3" | antibody=="H3K9me1" | antibody=="H3K9me3")) |
#                                         ((species=="Gthe") & (antibody=="H3K9me3" | antibody=="H3K79me2" | antibody=="H3K79me3")) |
#                                         ((species=="Tthe") & (antibody=="H3K9me1" | antibody=="H3K79me2" | antibody=="H3K79me3")),"N","Y"))
# bamlist <- c()
# for (species in c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru")){
#     bamlist <- c(paste0(species,"_all/bam_files/",species,"/",list.files(path=paste0("/home/smontgomery/no_backup/smontgomery/chip/",species,"_all/bam_files/",species,"/"),pattern="*.bam$",recursive=TRUE)),bamlist) 
# }
# for (file in metadata$sampleName){
#     tryCatch({pbc <- NA
#         bamfile <- grep(file, bamlist, value = TRUE, fixed = TRUE)
#         pbc <- round(PBC(paste0("/home/smontgomery/no_backup/smontgomery/chip/",bamfile)),5)
#         metadata[grep(file,metadata$sampleName),"PBC"] <- as.numeric(pbc)
#     }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# write.table(metadata,file="/home/smontgomery/data/chip/metadata_pbc.tsv",quote = F,sep = "\t",col.names = T,row.names = F)

metadata <- read.table(file="/home/smontgomery/data/chip/metadata_pbc.tsv",header = T,sep = "\t")
qctable <- read.table("/home/smontgomery/data/chip/gathered_QC_all.tsv",header=T,sep = "\t",fill=NA)
qctable$alignedRate <- as.numeric(gsub("%","",qctable$alignedRate))/100
qctable$dupRate <- as.numeric(gsub("%","",qctable$dupRate))/100
ssptable <- read.table("/home/smontgomery/data/chip/ssp_all.tsv",sep = "\t",col.names = c("Filename","numReads","estFragLen","corr_estFragLen","phatomPeak","corr_PhantomPeak","argmin_corr","min_corr","NSC","RSC","QualityTag"))
ssptable$sampleName <- gsub(".q30.rmdup.bam","",ssptable$Filename)

qctable$maxPeak <- pmax(qctable$N_narrowPeak,qctable$N_broadPeak)
qctable$maxFRiP <- pmax(qctable$narrowFRiP,qctable$broadFRiP)
qctable$maxAllFRiP <- pmax(qctable$narrowAllFRiP,qctable$broadAllFRiP)
metadata <- left_join(metadata,qctable)
metadata$NSC <- as.numeric(ssptable[match(metadata$sampleName,ssptable$sampleName),]$NSC)
metadata$RSC <- as.numeric(ssptable[match(metadata$sampleName,ssptable$sampleName),]$RSC)
metadata$QT <- as.factor(ssptable[match(metadata$sampleName,ssptable$sampleName),]$QualityTag)

##Correlation of all replicates of mark in species
for (sp in c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru")){
    grange = get_our_grange(crm,sp)
    genome_size <- sum(read.table(paste0(genome_folder,sp,"/",sp,".chromsizes.txt"))$V2)
    coveragefile <- as_tibble(read.table(paste0("/home/smontgomery/no_backup/smontgomery/chip/",sp,"_all/bed_files/",sp,".coverage.200bp.txt.gz"), header = T,sep = "\t",check.names = FALSE)) %>% replace(is.na(.), 0)
    ctrlsamp <- read.table(paste0("/home/smontgomery/no_backup/smontgomery/raw_reads/fastq_newgenome/",sp,"/sample_table.tsv"),sep="\t",header = T)%>%select(control)%>%table()%>%sort(decreasing = T)%>%names()%>%.[1]
    for (ab in c("H3K9ac","H3K27ac","H4K16ac","H3K4me2","H3K4me3","H3K36me3","H3K79me1","H3K79me2","H3K79me3","H3K27me3","H3K9me1","H3K9me2","H3K9me3")){
        tryCatch({
            tmptable <- metadata%>%filter(species==sp & antibody==ab)
            tmptable=tmptable%>%mutate(rankFRiP=rank(-maxFRiP,ties.method = "max"),rankPeakN=rank(-maxPeak,ties.method = "max")) %>%
                mutate(topQuad=if_else(maxFRiP > max(tmptable$maxFRiP)/2 & maxPeak > max(100,max(tmptable$maxPeak)/10),"Y","N")) %>%
                mutate(topQuad=if_else(table(.$topQuad)["Y"]==1 & rankFRiP == 2 & !is.na(table(.$topQuad)["Y"]),"lowfrip",topQuad))
            metadata[match(tmptable$sampleName,metadata$sampleName),"rankFRiP"] <- tmptable[,"rankFRiP"]
            metadata[match(tmptable$sampleName,metadata$sampleName),"rankPeakN"] <- tmptable[,"rankPeakN"]
            ##Spearman correlation with union of all peaks
            tmppeaks <- GRanges(seqlengths=setNames(read.table(paste0(genome_folder,sp,"/",sp,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,sp,"/",sp,".chromsizes.txt"))$V1))
            for (sample in tmptable$sampleName){
                tryCatch({
                    narrowpeak <- read.table(paste0("/home/smontgomery/no_backup/smontgomery/chip/",sp,"_all/peakCalling/",sp,"/",sample,"_peaks.narrowPeak"))
                    colnames(narrowpeak)[1:3] <- c("chrom","start","end")
                    narrowpeak <- c(GRanges(narrowpeak),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,sp,"/",sp,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,sp,"/",sp,".chromsizes.txt"))$V1)))
                    narrowpeak=narrowpeak[!(seqnames(narrowpeak) %in% c(crm))]
                    broadpeak <- read.table(paste0("/home/smontgomery/no_backup/smontgomery/chip/",sp,"_all/peakCalling/",sp,"/",sample,"_peaks.broadPeak"))
                    colnames(broadpeak)[1:3] <- c("chrom","start","end")
                    broadpeak <- c(GRanges(broadpeak),GRanges(seqlengths=setNames(read.table(paste0(genome_folder,sp,"/",sp,".chromsizes.txt"))$V2,read.table(paste0(genome_folder,sp,"/",sp,".chromsizes.txt"))$V1)))
                    broadpeak=broadpeak[!(seqnames(broadpeak) %in% c(crm))]
                    peaks <- union_ranges(narrowpeak,broadpeak)
                    tmppeaks <- union_ranges(tmppeaks,peaks)
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            }
            peakcov <- as.data.frame(mcols(find_overlaps(GRanges(coveragefile),tmppeaks)),optional=TRUE)
            peakcor <- cor(peakcov[,tmptable$sampleName],method="spearman")
            peakctrlcor <- cor(peakcov[,c(tmptable$sampleName,ctrlsamp)],method="spearman")[ctrlsamp,1:nrow(peakcor)]
            for (i in c(1:nrow(peakcor))){
                peakcor[i,i] <- 0
            }
            if (tmptable%>%filter(sampleName==names(sort(peakcor[as.character(tmptable%>%filter(rankFRiP==1)%>%select(sampleName)),],decreasing = T)[1]))%>%select(topQuad)=="N"){
                print(tmptable%>%filter(sampleName==names(sort(peakcor[as.character(tmptable%>%filter(rankFRiP==1)%>%select(sampleName)),],decreasing = T)[1]))%>%select(sampleName))
                if (tmptable[order(tmptable$rankFRiP),][2,"topQuad"]=="Y"){
                    tmptable=tmptable%>%mutate(topPeakQuad=if_else(maxFRiP > (tmptable[order(tmptable$rankFRiP),][2,"maxFRiP"])/2 & maxPeak > (tmptable[order(tmptable$rankFRiP),][2,"maxPeak"])/10,
                                        if_else(topQuad=="Y","Y","Yy"),"N"))
                } else if (tmptable[order(tmptable$rankFRiP),][2,"topQuad"]=="lowfrip"){
                    tmptable=tmptable%>%mutate(topPeakQuad=if_else(maxFRiP > (tmptable[order(tmptable$rankFRiP),][2,"maxFRiP"])/2 & maxPeak > (tmptable[order(tmptable$rankFRiP),][2,"maxPeak"])/10,"lowfripY","N"))
                }
                tmptable=tmptable%>%mutate(topPeakQuad=ifelse(rankFRiP==1,"lowcor",topPeakQuad))
            } else if (tmptable%>%filter(sampleName==names(sort(peakcor[as.character(tmptable%>%filter(rankFRiP==1)%>%select(sampleName)),],decreasing = T)[1]))%>%select(topQuad)=="lowfrip"){
                print(tmptable%>%filter(sampleName==names(sort(peakcor[as.character(tmptable%>%filter(rankFRiP==1)%>%select(sampleName)),],decreasing = T)[1]))%>%select(sampleName))
                tmptable=tmptable%>%mutate(topPeakQuad=if_else(maxFRiP > (tmptable[order(tmptable$rankFRiP),][2,"maxFRiP"])/2 & maxPeak > (tmptable[order(tmptable$rankFRiP),][2,"maxPeak"])/10,"lowfripY","N"))%>%mutate(topPeakQuad=ifelse(rankFRiP==1,"lowcor",topPeakQuad))
            } else {
                tmptable$topPeakQuad=tmptable$topQuad
            }
            ##Check if top quadrant has more than two reps
            if (table(tmptable$topPeakQuad)["Y"] > 2 & !is.na(table(tmptable$topPeakQuad)["Y"])){
                tmptable=tmptable%>%mutate(topPeakQuad=if_else(topPeakQuad=="Y",if_else(peakcor[subset(tmptable,rankFRiP==1)$sampleName,sampleName] >= 0.8 | rankFRiP==1,"Y","dropcor"),topPeakQuad))
            }
            metadata[match(tmptable$sampleName,metadata$sampleName),"topPeakQuad"] <- tmptable[,"topPeakQuad"]
            metadata[match(tmptable$sampleName,metadata$sampleName),"maxPeakcor"] <- rowMaxs(peakcor,useNames = F)
            metadata[match(tmptable$sampleName,metadata$sampleName),"bestPeakcor"] <- c(rowMaxs(peakcor,useNames = F)==max(rowMaxs(peakcor,useNames = F)))
            metadata[match(tmptable$sampleName,metadata$sampleName),"peakctrlcor"] <- peakctrlcor
            hclust_plot <- hclust(as.dist(1-cor(peakcov[,tmptable$sampleName],method="pearson")),method = "complete")
            par(mar=c(0.1,4.1, 4.1, 2.1))
            if (nrow(peakcor) <= 2){
                plot(as.dendrogram(hclust_plot),xlab = "",main = "",sub = "")
            } else {
                plot(hclust_plot,xlab = "",main = "",sub = "",labels=FALSE)
            }
            tmpplot <- recordPlot()
            pdf(paste0("/home/smontgomery/data/chip/QCplots/correlation/",sp,".",ab,".peakCoverageCorrelation.spearman.pdf"),width=15,height=15)
            print(plot_grid(tmpplot,
                      ggplot(reshape2::melt(as.matrix(cor(peakcov[,tmptable$sampleName],method="spearman"))),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]),y=factor(Var1, levels=hclust_plot$labels[hclust_plot$order]), fill=value)) +
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

            ##Pearson correlation over whole genome
            gencov <- as.data.frame(mcols(GRanges(coveragefile)),optional=TRUE)
            gencor <- cor(gencov[,tmptable$sampleName],method="pearson")
            genctrlcor <- cor(gencov[,c(tmptable$sampleName,ctrlsamp)],method="pearson")[ctrlsamp,1:nrow(gencor)]
            for (i in c(1:nrow(gencor))){
                gencor[i,i] <- 0
            }
            ##Check if top correlate of top FRiP in top quadrant
            if (tmptable%>%filter(sampleName==names(sort(gencor[as.character(tmptable%>%filter(rankFRiP==1)%>%select(sampleName)),],decreasing = T)[1]))%>%select(topQuad)=="N"){
                print(tmptable%>%filter(sampleName==names(sort(gencor[as.character(tmptable%>%filter(rankFRiP==1)%>%select(sampleName)),],decreasing = T)[1]))%>%select(sampleName))
                if (tmptable[order(tmptable$rankFRiP),][2,"topQuad"]=="Y"){
                    tmptable=tmptable%>%mutate(topGenQuad=if_else(maxFRiP > (tmptable[order(tmptable$rankFRiP),][2,"maxFRiP"])/2 & maxPeak > (tmptable[order(tmptable$rankFRiP),][2,"maxPeak"])/10,
                                        if_else(topQuad=="Y","Pass","Pass"),"Fail"))
                } else if (tmptable[order(tmptable$rankFRiP),][2,"topQuad"]=="lowfrip"){
                    tmptable=tmptable%>%mutate(topGenQuad=if_else(maxFRiP > (tmptable[order(tmptable$rankFRiP),][2,"maxFRiP"])/2 & maxPeak > (tmptable[order(tmptable$rankFRiP),][2,"maxPeak"])/10,"Pass","Fail"))
                } else if (tmptable[order(tmptable$rankFRiP),][1,"topQuad"]=="N" | tmptable[order(tmptable$rankFRiP),][1,"topQuad"]=="lowcor"){
                    tmptable=tmptable%>%mutate(topGenQuad="Fail")
                }
                tmptable=tmptable%>%mutate(topGenQuad=ifelse(rankFRiP==1,"Check",topGenQuad))
            } else if (tmptable%>%filter(sampleName==names(sort(gencor[as.character(tmptable%>%filter(rankFRiP==1)%>%select(sampleName)),],decreasing = T)[1]))%>%select(topQuad)=="lowfrip"){
                print(tmptable%>%filter(sampleName==names(sort(gencor[as.character(tmptable%>%filter(rankFRiP==1)%>%select(sampleName)),],decreasing = T)[1]))%>%select(sampleName))
                tmptable=tmptable%>%mutate(topGenQuad=if_else(maxFRiP > (tmptable[order(tmptable$rankFRiP),][2,"maxFRiP"])/2 & maxPeak > (tmptable[order(tmptable$rankFRiP),][2,"maxPeak"])/10,"Pass","Fail"))%>%mutate(topGenQuad=ifelse(rankFRiP==1,"Check",topGenQuad))
            } else {
                tmptable$topGenQuad=gsub("Y","Pass",gsub("lowfrip","Y",gsub("N","Fail",tmptable$topQuad)))
            }
            ##Check if top quadrant has more than two reps
            if (table(tmptable$topGenQuad)["Pass"] > 2 & tmptable%>%filter(rankFRiP==1)%>%select(topGenQuad)=="Pass"){
                tmptable=tmptable%>%mutate(topGenQuad=if_else(topGenQuad=="Pass",if_else(gencor[subset(tmptable,rankFRiP==1)$sampleName,sampleName] >= 0.8 | rankFRiP==1,"Pass","Drop"),topGenQuad))
            }
            ##Add correlation to top frip
            tmptable=tmptable%>%mutate(topfripGencor=if_else(rankFRiP==1,1,gencor[subset(tmptable,rankFRiP==1)$sampleName,sampleName]))
            tmptable=tmptable%>%mutate(topGenQuad=if_else(present=="N","Fail",topGenQuad))
            metadata[match(tmptable$sampleName,metadata$sampleName),"topGenQuad"] <- tmptable[,"topGenQuad"]
            metadata[match(tmptable$sampleName,metadata$sampleName),"maxGencor"] <- rowMaxs(gencor,useNames = F)
            metadata[match(tmptable$sampleName,metadata$sampleName),"bestGencor"] <- c(rowMaxs(gencor,useNames = F)==max(rowMaxs(gencor,useNames = F)))
            metadata[match(tmptable$sampleName,metadata$sampleName),"genctrlcor"] <- genctrlcor
            metadata[match(tmptable$sampleName,metadata$sampleName),"topfripGencor"] <- tmptable[,"topfripGencor"]
            hclust_plot <- hclust(as.dist(1-cor(gencov[,tmptable$sampleName],method="pearson")),method = "complete")
            par(mar=c(0.1,4.1, 4.1, 2.1))
            if (nrow(gencor) <= 2){
                plot(as.dendrogram(hclust_plot),xlab = "",main = "",sub = "")
            } else {
                plot(hclust_plot,xlab = "",main = "",sub = "",labels=FALSE)
            }
            tmpplot <- recordPlot()
            pdf(paste0("/home/smontgomery/data/chip/QCplots/correlation/",sp,".",ab,".genomeCoverageCorrelation.pearson.pdf"),width=15,height=15)
            print(plot_grid(tmpplot,
                      ggplot(reshape2::melt(as.matrix(cor(gencov[,tmptable$sampleName],method="pearson"))),aes(x=factor(Var2, levels=hclust_plot$labels[hclust_plot$order]),y=factor(Var1, levels=hclust_plot$labels[hclust_plot$order]), fill=value)) +
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
            pdf(paste0("/home/smontgomery/data/chip/QCplots/frip/",sp,".",ab,".fripvspeak.pdf"),width=10,height=10)
            print(tmptable%>%ggplot(aes(x=maxFRiP,y=maxPeak,size=topfripGencor)) +
                      geom_vline(linetype = 2,xintercept = if_else("Check" %in% tmptable$topGenQuad,tmptable[order(tmptable$rankFRiP),][2,"maxFRiP"],tmptable[order(tmptable$rankFRiP),][1,"maxFRiP"])/2) +
                      geom_hline(linetype = 2,yintercept = max(if_else("Check" %in% tmptable$topGenQuad,tmptable[order(tmptable$rankFRiP),][2,"maxPeak"],tmptable[order(tmptable$rankFRiP),][1,"maxPeak"])/10,100)) +
                      geom_point(color="black",shape=21,aes(fill=factor(topGenQuad))) +
                      scale_x_continuous(name="FRiP",limits = c(0,1)) +
                      scale_y_continuous(name="N peaks",limits = c(0,max(max(tmptable$maxPeak),100))) +
                      scale_fill_manual(name="QC outcome",values = c("Pass"="#748E54","Fail"="#F93943","Check"="#E0BAD7","Drop"="#FFC43D")) +
                      scale_size(name="Correlation \nto top FRiP",limits = c(0,1)) +
                      theme_classic() +
                      theme(axis.text.x = element_text(color = "black",size = 36),
                            axis.text.y = element_text(color = "black",size = 36),
                            axis.title.x = element_text(color = "black",size = 36),
                            axis.title.y = element_text(color = "black",size = 36),
                            legend.title = element_text(color = "black",size = 36),
                            legend.text = element_text(colour = "black",size = 36)))
            dev.off()
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
}

# write.table(metadata,file="/home/smontgomery/data/chip/metadata_full.tsv",quote = F,sep = "\t",col.names = T,row.names = F)

##Read in file saved as checkpoint
metadata <- read.table(file="/home/smontgomery/data/chip/metadata_full.tsv",header = T,sep = "\t")
final <- read.table(file="/home/smontgomery/data/chip/metadata_final.tsv",header = T,sep = "\t")

##Antibody table
pdf(paste0("/home/smontgomery/data/chip/QCplots/antibody_summary.pdf"),width=10,height=10)
final%>%filter(!is.na(antibody))%>%group_by(species,antibody,refAb)%>%summarise(Passes=any(topGenQuad=="Pass"),Fails=any(topGenQuad=="Fail"),Checks=any(topGenQuad=="Check"),Drops=any(topGenQuad=="Drop"))%>%mutate(Status=if_else(antibody=="H3" | Passes==TRUE,"Pass",if_else(Checks==TRUE,"Check",if_else(Drops==TRUE,"Drop",if_else(Fails==TRUE,"Fail","Absent")))))%>%
    ggplot(aes(x=factor(species,levels=c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru")),y=factor(refAb,levels=rev(c("ab1791","ab32356","ab176878","07-473","17-658","MA5-33384","ab9045","MA5-33385","ab1220","ab176916","C15410056","713008","ab4729","17-683","39136","07-449","17-622","ab192985","9733","C15410069","MA5-11198","ab9050","ab177183","ab3594","C15410068","07-329"))),fill=Status)) +
        geom_tile() +
        scale_fill_manual(name="QC outcome",values = c("Pass"="#748E54","Fail"="#F93943","Check"="#E0BAD7","Drop"="#FFC43D","Absent"="grey"))  +
        theme_classic() +
        theme(axis.text.x = element_text(color = "black",size = 12),
              axis.text.y = element_text(color = "black",size = 12),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.title = element_text(color = "black",size = 12),
              legend.text = element_text(colour = "black",size = 12))
dev.off()

##Write table for GEO submission
unite(subset(final,final=="Y")[c("species","antibody")],newcol)%>%
    mutate(newcol=gsub("_NA","_input",newcol))%>%
    mutate(counter=rep(1,nrow(.)))%>%
    group_by(newcol)%>%
    mutate(counter2=cumsum(counter))%>%
    select(-counter)%>%
    unite("newname",newcol:counter2,sep = "_rep")%>%
    cbind(.,subset(final,final=="Y"))%>%
    mutate(procfile=paste0(species,"_",antibody,".log2r.input.bw"))%>%
    mutate(fastq1=paste0(newname,"_R1_001.fastq.gz"))%>%
    mutate(fastq2=paste0(newname,"_R2_001.fastq.gz"))%>%
    mutate(newname=gsub("Acas","Acanthamoeba castellani",gsub("Atha","Arabidopsis thaliana",gsub("Bnat","Bigelowiella natans",gsub("Cfra","Creolimax fragrantissima",gsub("Ddis","Dictyostelium discoideum",gsub("Gthe","Guillardia theta",
        gsub("Ngru","Naegleria gruberi",gsub("Nvec","Nematostella vectensis",gsub("Ppat","Physcomitrium patens",gsub("Scer","Saccharomyces cerevisiae",gsub("Spun","Spizellomyces punctatus",gsub("Tthe","Tetrahymena thermophila",newname)))))))))))))%>%
    mutate(species=gsub("Acas","Acanthamoeba castellani",gsub("Atha","Arabidopsis thaliana",gsub("Bnat","Bigelowiella natans",gsub("Cfra","Creolimax fragrantissima",gsub("Ddis","Dictyostelium discoideum",gsub("Gthe","Guillardia theta",
        gsub("Ngru","Naegleria gruberi",gsub("Nvec","Nematostella vectensis",gsub("Ppat","Physcomitrium patens",gsub("Scer","Saccharomyces cerevisiae",gsub("Spun","Spizellomyces punctatus",gsub("Tthe","Tetrahymena thermophila",species)))))))))))))%>%
    mutate(refAb=gsub("ab4729","H3K27ac (Abcam, ab4729)",gsub("07-449","H3K27me3 (Millipore, #07-449)",gsub("ab9050","H3K36me3 (Abcam, ab9050)",gsub("07-473","H3K4me3 (Millipore, #07-473)",gsub("ab32356","H3K4me2 (Abcam, ab32356)",
        gsub("ab177183","H3K79me1 (Abcam, ab177183)",gsub("ab3594","H3K79me2 (Abcam, ab3594)",gsub("C15410068","H3K79me3 (Diagenode, C15410068)",gsub("17-658","H3K9ac (Millipore, #17-658)",gsub("ab9045","H3K9me1 (Abcam, ab9045)",
            gsub("ab176916","H3K9me3 (Abcam, ab176916)",gsub("07-329","H4K16ac (Millipore, #07-329)",gsub("ab1791","H3 (Abcam, ab1791)",gsub("17-683","H3K27ac (Millipore, #17-683)",gsub("17-622","H3K27me3 (Millipore, #17-622)",
                gsub("ab176878","H3K4me2 (Abcam, ab176878)",gsub("713008","H3K9me3 (Invitrogen, #713008)",gsub("ab1220","H3K9me2 (Abcam, ab1220)",gsub("ab192985","H3K27me3 (Abcam, ab192985)",gsub("MA5-11198","H3K27me3 (Invitrogen, MA5-11198)",
                    gsub("39136","H3K27ac (Active Motif, #39136)",gsub("C15410069","H3K27me3 (Diagenode, C15410069)",refAb)))))))))))))))))))))))%>%
    select(newname,species,refAb,procfile,fastq1,fastq2)%>%
    mutate(refAb=if_else(is.na(refAb),"none",refAb))%>%
    write.table(file="Downloads/temptable.tsv",quote = F,sep = "\t",row.names = F,col.names = F)





options(ucscChromosomeNames=FALSE)

for (sp in c("Acas","Ppat","Ngru")){
    if (sp=="Acas"){
        mark="H3K36me3"
    } else if (sp=="Ppat"){
        mark="H3K4me3"
    } else if (sp=="Ngru"){
        mark="H3K9me3"
    }
    grange = get_our_grange(crm,sp)
    tenogene <- read.table(paste0(genome_folder,sp,"/",sp,".TE.nogene.bed"))
    colnames(tenogene) <- c("chromosome","start","end","transcript","score","strand","source","type","dot")
    if (sp == "Scer"){
        genetrack <- data.frame(rtracklayer::import(paste0(raw_folder,sp,"_long.annot.gtf")))%>%filter(type=="CDS")%>%mutate(type="protein_coding")
    } else if (sp == "Ddis") {
        genetrack <- data.frame(rtracklayer::import(paste0(raw_folder,sp,"_long.annot.gtf")))%>%filter(type=="exon")%>%mutate(type="protein_coding")
        genetrack[,1] <- paste0("chr",genetrack[,1])
    } else {
        genetrack <- data.frame(rtracklayer::import(paste0(raw_folder,sp,"_long.annot.gtf")))%>%filter(type=="exon")%>%mutate(type="protein_coding")
    }
    colnames(genetrack) <- c("chromosome","start","end","width","strand","source","feature","score","phase","transcript","gene")
    if (sp=="Nvec" | sp=="Ddis" | sp=="Ppat" | sp=="Gthe" | sp=="Bnat" | sp=="Ngru"){
        colnames(genetrack) <- c("chromosome","start","end","width","strand","source","feature","score","phase","transcript","gene","locus")
    }
    f <- paste0(chip_folder,sp,"/bw/",sp,"_geneExp.bw")
    if (file.exists(f)){
        tmpRNAold <- read_bigwig(f)
    }
    f <- paste0(chip_folder,sp,"_all/bw/",sp,"_merged.sorted.plus.bw")
    if (file.exists(f)){
        tmpRNAfwd <- read_bigwig(f)
    }
    f <- paste0(chip_folder,sp,"_all/bw/",sp,"_merged.sorted.minus.bw")
    if (file.exists(f)){
        tmpRNArev <- read_bigwig(f)
    }
    filelist=list.files(path=paste0(chip_folder,sp,"_all/bw_files/",sp),pattern = paste0(".*",mark,".*"))
    namlist=c()
    for (i in 1:length(filelist)){
        f=paste0(chip_folder,sp,"_all/bw_files/",sp,"/",filelist[i])
        nam=paste(gsub(".q30.rmdup.CPM.bw","",filelist[i]),i,sep="_")
        namlist=c(namlist,nam)
        assign(nam,read_bigwig(f))
    }
    if (sp=="Acas" & mark=="H3K36me3"){
        figname="EDF2.1"
        chr="scaffold_1"
        start=1594689
        end=1623069
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig2_draft/",sp,".",mark,".",figname,"_checkFAIL.pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(Acas_H3K36me3_191121_6%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#E0BAD7","#E0BAD7"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Acas_H3K36me3_170821_2_3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#748E54","#748E54"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Acas_H3K36me3_170821_4%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#748E54","#748E54"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Acas_H3K36me3_070722_1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#F93943","#F93943"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Acas_H3K36me3_091122_2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#F93943","#F93943"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Acas_H3K36me3_190123_5%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#F93943","#F93943"),showTitle=FALSE,showAxis=TRUE),
                        AnnotationTrack(tenogene%>%filter(chromosome==chr),name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
    } else if (sp=="Ppat" & mark=="H3K4me3"){
        figname="EDF2.1"
        chr="GWHFIHF00000003.1"
        start=19402790
        end=20245202
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig2_draft/",sp,".",mark,".",figname,"_checkPASS.pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(Ppat_H3K4me3_090622_2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#E0BAD7","#E0BAD7"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ppat_H3K4me3_A_091122_3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#748E54","#748E54"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ppat_H3K4me3_B_091122_4%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#748E54","#748E54"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ppat_H3K4me3_010223_1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#F93943","#F93943"),showTitle=FALSE,showAxis=TRUE),
                        AnnotationTrack(tenogene%>%filter(chromosome==chr),name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
    } else if (sp=="Ngru" & mark=="H3K9me3"){
        figname="EDF2.1"
        chr="NAEGRscaffold_1"
        start=1166556
        end=1538193
        ylimRNA=range(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range(),tmpRNArev%>%select("chrom","start","end","value")%>%GRanges()%>%subsetByOverlaps(.,GRanges(seqnames = chr,ranges = c(start:end)))%>%as.data.frame()%>%select(value)%>%range())
        pdf(paste0(plot_folder,"extData_Fig2_draft/",sp,".",mark,".",figname,"_checkPASS.pdf"),width=15,height=15)
        plotTracks(list(DataTrack(tmpRNAfwd%>%select("chrom","start","end","value")%>%GRanges(),name="RNAplus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(tmpRNArev%>%select("chrom","start","end","value")%>%GRanges(),name="RNAminus",ylim=ylimRNA,col.histogram=NULL,fill.histogram=c("#232323","#232323"),showTitle=FALSE,showAxis=FALSE),
                        DataTrack(Ngru_H3K9me3_010823_2%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#E0BAD7","#E0BAD7"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_080224_5%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#748E54","#748E54"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_250124_8%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#748E54","#748E54"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_ab176916_160323_10%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#748E54","#748E54"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_005SDS_281122_1%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#748E54","#748E54"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_01SDS_281122_3%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#748E54","#748E54"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_210721_6%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#F93943","#F93943"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_ab8898_160323_11%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#F93943","#F93943"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_713008_160323_9%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#F93943","#F93943"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_210722_7%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#F93943","#F93943"),showTitle=FALSE,showAxis=TRUE),
                        DataTrack(Ngru_H3K9me3_020322_4%>%select("chrom","start","end","value")%>%GRanges(),name="H3K36me3",col.histogram=NULL,fill.histogram=c("#F93943","#F93943"),showTitle=FALSE,showAxis=TRUE),
                        AnnotationTrack(tenogene%>%filter(chromosome==chr),name = "TEs",shape = "box",col.title="black",fill="#808080",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="box",showTitle=FALSE),
                        GeneRegionTrack(genetrack, feature="protein_coding",genome = genome(grange),name = "Genes",col.title="black",fill="#232323",col.line="#232323",shape="arrow",showTitle=FALSE)),
                    type="histogram",legend=TRUE,chromosome = chr,from = start,to=end,col="transparent",col.axis="black",background.title="transparent")
        dev.off()
    }
}
