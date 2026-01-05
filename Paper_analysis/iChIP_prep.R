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

##Source functions
source("~/OneDrive/Scripts/iChIP2_paper/functions.R")

##Species names
specieslist=data.frame(species=c("Nvec","Cfra","Scer","Spun","Acas","Ddis","Atha","Ppat","Gthe","Tthe","Bnat","Ngru") ,statenumber=c(16,9,11,16,14,9,12,10,9,11,14,12))

##Set colours
cols <- c("H3K9ac"="#94FA7E","H3K27ac"="#257604","H4K16ac"="#112B0A","H3K4me2"="#FFE04D","H3K4me3"="#F46D43","H3K36me3"="#B82300","H3K79me1"="#F7BEEE","H3K79me2"="#D448BA","H3K79me3"="#70015D","H3K27me3"="#4120A9","H3K9me1"="#41B6C4","H3K9me2"="#3262AB","H3K9me3"="#3262AB","H3"="grey","input"="black")
##Probably need to add more colours and classes
TEcols=c("DNA"="steelblue","TIR"="steelblue1","Helitron"="steelblue4","TE_DNA_fragment"="cyan","TE_DNA_intact"="darkblue","DNAauto"="aquamarine","DNAnona"="lightgreen",
       "RT"="firebrick3","LTR"="firebrick","TE_RT_LTR_fragment"="pink","TE_RT_LTR_intact"="darkred","TE_RT_LINE"="orange","TE_RT_SINE"="orange3","LINE"="orange","SINE"="orange3","PLE"="yellow4",
       "MITE"="purple","TE_RT_DIRS"="grey","DIRS"="violet","NA"="grey")

#########################################################################################################################################################################################################################################################################
##Write RNAcounts levels plus gene length files
for (species in specieslist$species){
    if (species == "Nvec"){
        rnaseq = read_tsv(paste0("/home/smontgomery/cluster/asebe/proj/RNAseq_multisps_NG_revision/Gene_expression_counts/",species,"_RNA_counts.txt"), skip = 1, col_select = c(1,2), col_names = c("gene_id","RNAcounts"))
        RNAcountsfile <- read_tsv(paste0("/home/smontgomery/cluster/xgraubove/genomes/data/",species,"_long.annot.gtf"),col_select = c(1,3,4,5,7,9),col_names=c("chrom","2","feature","start","end","6","strand","8","id"))
        RNAcountsfile=RNAcountsfile%>%
            filter(feature=="transcript")%>%
            mutate(transcript_id=str_split_fixed(id,'\"',5)[,4], gene_id=str_split_fixed(id,'\"',5)[,2], gene_length=end-start)%>%
            select(-id)%>%
            left_join(rnaseq, by = "gene_id") %>%
            arrange(desc(RNAcounts))
    } else if (species=="Cfra"){
        rnaseq = read_tsv(paste0("/home/smontgomery/cluster/asebe/proj/RNAseq_multisps_NG_revision/Gene_expression_counts/",species,"_RNA_counts.txt"), col_select = c(1,2), col_names = c("gene_id","RNAcounts"))
        RNAcountsfile <- read_tsv(paste0("/home/smontgomery/cluster/xgraubove/genomes/data/",species,"_long.annot.gtf"),col_select = c(1,3,4,5,7,9),col_names=c("chrom","2","feature","start","end","6","strand","8","id"))
        RNAcountsfile=RNAcountsfile%>%
            filter(feature=="transcript")%>%
            mutate(transcript_id=str_split_fixed(id,'\"',5)[,2], gene_id=str_split_fixed(id,'\"',5)[,4], gene_length=end-start)%>%
            select(-id)%>%
            left_join(rnaseq, by = "gene_id") %>%
            arrange(desc(RNAcounts))
    } else {
        rnaseq = read_tsv(paste0("/home/smontgomery/cluster/asebe/proj/RNAseq_multisps_NG_revision/Gene_expression_counts/",species,"_RNA_counts.txt"), skip = 1, col_select = c(1,2,3), col_names = c("gene_id","rep1","rep2"))
        RNAcountsfile <- read_tsv(paste0("/home/smontgomery/cluster/xgraubove/genomes/data/",species,"_long.annot.gtf"),col_select = c(1,3,4,5,7,9),col_names=c("chrom","2","feature","start","end","6","strand","8","id"))
        RNAcountsfile=RNAcountsfile%>%
            filter(feature=="transcript")%>%
            mutate(transcript_id=str_split_fixed(id,'\"',5)[,2], gene_id=str_split_fixed(id,'\"',5)[,4], gene_length=end-start)%>%
            select(-id)%>%
            left_join(rnaseq, by = "gene_id") %>%
            mutate(RNAcounts = (rep1 + rep2)) %>%
            arrange(desc(RNAcounts))
    }
    write_tsv(RNAcountsfile, paste0("/home/smontgomery/cluster/smontgomery/genomes/",species,"/",species,"_expr_sorted.tsv"))

    ##Reorder gene clusters by expression
    if (species == "Gthe"){
        generegions <- read.table(paste0("/home/smontgomery/no_backup/smontgomery/chip/",species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.txt"))
    } else {
        generegions <- read.table(paste0("/home/smontgomery/no_backup/smontgomery/chip/",species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.txt"))
    }
    generegions$V13 <- paste0("gene_",generegions$V13)
    if (species == "Cfra"){
        RNAcountsfile$gene_id <- gsub("Cfra_CFRG0","Cfra_CFRCFRG0",RNAcountsfile$gene_id)
    }
    generegions$RNAcounts <- RNAcountsfile$RNAcounts[match(generegions$V4,RNAcountsfile$gene_id)]
    print(species)
    print(generegions%>%group_by(V13)%>%summarise(mean=mean(RNAcounts,na.rm=T),med=median(RNAcounts,na.rm=T)))
}

#########################################################################################################################################################################################################################################################################
##Remove genes with TSS in 5kb upstream after reordering gene clusters by expression
for (species in specieslist$species){
    if (species == "Gthe"){
        generegions <- read.table(paste0("/home/smontgomery/no_backup/smontgomery/chip/",species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.reorder.txt"))
    } else {
        generegions <- read.table(paste0("/home/smontgomery/no_backup/smontgomery/chip/",species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.reorder.txt"))
    }
    generegions$V13 <- paste0("gene_",generegions$V13)
    RNAcountsfile <- read.table(paste0("/home/smontgomery/cluster/smontgomery/genomes/",species,"/",species,"_expr_sorted.tsv"),header=TRUE)
    if (species == "Cfra"){
        RNAcountsfile$gene_id <- gsub("Cfra_CFRG0","Cfra_CFRCFRG0",RNAcountsfile$gene_id)
    }
    generegions$RNAcounts <- RNAcountsfile$RNAcounts[match(generegions$V4,RNAcountsfile$gene_id)]
    if (species == "Nvec"){
        generegions$RNAcounts <- RNAcountsfile$RNAcounts[match(generegions$V4,RNAcountsfile$transcript_id)]
    }
    grange = get_our_grange(species=species)
    TSSgrange <- add_TSS(grange)
    nohead2head <- bed_intersect(as.data.frame(promoters(TSSgrange,upstream=5000,downstream=0)%>%
                                                   mutate(type="promoter")%>%
                                                   trim())%>%
                                     mutate(end=end-1)%>%
                                     mutate(chrom=seqnames),
                                as.data.frame(TSSgrange)%>%
                                     mutate(start=start-1)%>%
                                     mutate(chrom=seqnames),
                                invert = TRUE)%>%
        select(gene_id)
    head2head <- bed_intersect(as.data.frame(promoters(TSSgrange,upstream=5000,downstream=0)%>%
                                     mutate(type="promoter")%>%
                                     trim())%>%
                       mutate(end=end-1)%>%
                       mutate(chrom=seqnames),
                   as.data.frame(TSSgrange)%>%
                       mutate(start=start-1)%>%
                       mutate(chrom=seqnames),
                   invert = FALSE)%>%
         select(gene_id.x)%>%unique()
    TSSoverlap <- bed_intersect(as.data.frame(TSSgrange%>%shift_upstream(2))%>%
                      mutate(start=start-1)%>%
                      mutate(chrom=seqnames),
                  as.data.frame(grange[mcols(grange)$type == "protein_coding"])%>%
                      mutate(start=start-1)%>%
                      mutate(chrom=seqnames),
                  invert = FALSE)%>%
        select(gene_id.x)%>%unique()
    if (species == "Gthe"){
	    subset(generegions, V4 %in% subset(nohead2head,!(gene_id %in% TSSoverlap$gene_id.x))$gene_id)%>%
	        write.table(file=paste0("/home/smontgomery/no_backup/smontgomery/chip/",species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.reorder.no5kb.txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
	    subset(generegions, V4 %in% head2head$gene_id.x)%>%
	        write.table(file=paste0("/home/smontgomery/no_backup/smontgomery/chip/",species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k8.reorder.only5kb.txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
    } else {
	    subset(generegions, V4 %in% subset(nohead2head,!(gene_id %in% TSSoverlap$gene_id.x))$gene_id)%>%
	        write.table(file=paste0("/home/smontgomery/no_backup/smontgomery/chip/",species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.reorder.no5kb.txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
	    subset(generegions, V4 %in% head2head$gene_id.x)%>%
	        write.table(file=paste0("/home/smontgomery/no_backup/smontgomery/chip/",species,"_all/profiles/log2r.",species,".gene.noregion.heatmap.k6.reorder.only5kb.txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
    }
}
