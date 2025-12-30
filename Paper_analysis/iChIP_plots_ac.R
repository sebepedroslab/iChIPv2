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
##Acetylation enrichment stuff
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

