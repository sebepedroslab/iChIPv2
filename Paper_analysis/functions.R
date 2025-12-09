# install.packages("devtools")
# install.packages("versions-package")
# install.packages("remotes")
# install.packages("BiocManager")
# library(remotes)
# library(BiocManager)
# BiocManager::install(version = "3.14")
# BiocManager::install("ComplexHeatmap", version = "3.14")
# install_version("ggalluvial", "0.12.5")
# install_version("ggpubr", "0.4.0")
# install_version("RColorBrewer", "1.1-2")
# install_version("testthat","3.1.4")
# install_version("tidyverse","1.3.1")
# BiocManager::install("rtracklayer")
# install_version("valr","0.6.6")
# install_version("vdiffr","1.0.5")
# install.packages("plotly")
library("ComplexHeatmap")
library("ggalluvial")
library("ggpubr")
library("RColorBrewer")
library("testthat")
library("tidyverse")
library("valr")
library("vdiffr")
library("reshape2")
library("rtracklayer")
library("GenomicRanges")
library("plyranges")

##From support_script/functions.R
		#####################################
		## functions used in several figures
		#####################################

		## function takes emission matrix as read in (raw) and rearanges the columns according to sample order (samo)
		## outputs a matrix that can be plotted 
		raw_to_matrix=function(raw,samo){
		  mat = raw%>%dplyr::select(-`State (Emission order)`)
		  mat = as.matrix(mat)
		  rownames(mat)=raw$`State (Emission order)`
		  mat=mat[,samo]
		}

		## reads in ChromHMM output bed files, renames column X4 to state and excludes mitochondria and chloroplast
		## outputs tibble with genomic coordinates and state assignment
		read_and_fix_states=function(file,crm=c("Pt","Mt","Mito")){
		  states = valr::read_bed(file,n_fields = 4)%>%dplyr::rename(state=name)
		  states = filter(states,!chrom%in%crm)
		  states
		  }

		## reads in the gff file used in this paper (included in the data subfolder), set useful colnames, excludes mitochondria and chloroplast
		## extracts gene id
		# get_our_gff=function(crm=c("Pt","Mt")){
		#   gff = as_tibble(as.data.frame(import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".TEgene.all.gtf")))[,c(1,6,7,2,3,8,5,9,10,11,12,13,14)])
		#   colnames(gff)=c("chrom","Source","feature","start","end","X6","strand","X8","geneId","gene_id","class_id","order_id","family_id")
		#   gff = filter(gff,!chrom%in%c(crm))
		#   # files from chromHMM uses 1:5, Mt, Pt rather than Chr1:5, ChrC, ChrM,
		#   # gff = gff%>%mutate(chrom=str_match(chrom,"^Chr(.+)")[,2])
		#   # gff= gff%>%mutate(geneId=str_match(info,"^ID=([A-Za-z0-9]+)")[,2])
		#   gff
		# }
		get_our_grange=function(crm=c("Pt","Mt","Mito"),species=species){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".TEgene.all.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange$size <- width(grange)
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		## splits gene bodies into 200bp bins starting from TSS
		## output tibble with these features
		add_bins=function(grange){
			tes=grange[mcols(grange)$type != "protein_coding"]
			grange=grange[mcols(grange)$type == "protein_coding"]
			bin1=pintersect(grange%>%mutate(type=as.factor("bin1"))%>%resize(200),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
			bin2=pintersect(grange%>%mutate(type=as.factor("bin2"))%>%resize(200)%>%shift_downstream(200),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
			bin3=pintersect(grange%>%mutate(type=as.factor("bin3"))%>%resize(200)%>%shift_downstream(400),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
			bin4=pintersect(grange%>%mutate(type=as.factor("bin4"))%>%resize(200)%>%shift_downstream(600),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
			gene_body=pintersect(grange%>%mutate(type=as.factor("geneBody"))%>%shift_downstream(800),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
			grangen=c(bin1,bin2,bin3,bin4,gene_body,tes)
			grangen
		}
		## adds introns to gff. introns= regions overlapping mRNA but not exons, assigns a geneId to each intro
		## output tibble with extra introns 
		# add_intron=function(gff){
		#   intron_strand=list()
		#   g_size = gff%>%group_by(chrom)%>%summarize(size=max(end)) # chrom size
		#   for (strand in c("+","-")){
		#     intergenic=bed_complement(gff%>%filter(feature=="transcript" & strand==!!strand),g_size)
		#     exon = gff%>%filter(feature=="exon" & strand==!!strand)%>%dplyr::select(chrom,start,end) # exons on + strand
		#     both = rbind(intergenic,exon) 
		#     intron_strand[[strand]]= bed_complement(both,g_size)%>%mutate(feature="intron",strand=strand)
		#   }
		#   introns = do.call("rbind",intron_strand)
		#   genes = filter(gff,feature=="transcript")%>%dplyr::select(chrom,start,end,feature,strand,geneId,gene_id)%>%distinct()
		#   an_introns = bed_intersect(introns,genes)%>%filter(.overlap>0)%>%dplyr::select(-start.y,-end.y)%>%distinct()%>%
		#     filter(strand.x==strand.y)%>%dplyr::select(-feature.y,-strand.y,-.overlap)%>%mutate(Source="-",X6=as.double("."),X8=as.integer("."),class_id="",order_id="",family_id="")
		#   names(an_introns) <- c("chrom","start","end","feature","strand","geneId","gene_id","Source","X6","X8","class_id","order_id","family_id")
		#   gfft = full_join(gff,an_introns)
		#   gfft
		# }
		add_intron=function(grange){
			introns <- setdiff_ranges(grange[mcols(grange)$type == "protein_coding"],grange[mcols(grange)$type == "exon"])
			mcols(introns) <- mcols(grange[mcols(grange)$type == "protein_coding"][findOverlaps(introns,grange[mcols(grange)$type == "protein_coding"],select="arbitrary")]%>%mutate(type="intron"))
			granget=c(grange,introns)
			granget
		}
		## adds promoter to gff. promoter= regions 500bp upstream of gene but not into upstream genes, though allows overlap of promoters with each other (problematic?)
		## output tibble with promoters
		# add_promoter=function(gff){
		#   promoter_strand_full=as_tibble(as.data.frame(promoters(import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".TEgene.all.gtf")),upstream = 500, downstream = 0))[,c(1,6,7,2,3,8,5,9,10,11,12,13,14)])%>%filter(type=="transcript")%>%mutate(type="promoter")
		#   colnames(promoter_strand_full)=c("chrom","Source","feature","start","end","X6","strand","X8","geneId","gene_id","class_id","order_id","family_id")
		#   promoter_strand=bed_subtract(promoter_strand_full,gff)
		#   gffp = full_join(gff,promoter_strand)
		#   gffp
		# }
		add_promoter=function(grange){
			promoter=filter(grange,type=="bin1")%>%promoters(upstream=500,downstream=0)%>%mutate(type="promoter")%>%trim()
			grangep=c(grange,promoter)
		}
		merge_overlap_feature=function(grange){
			tenooverlap <- grange[mcols(grange)$source == "EDTA"][countOverlaps(grange[mcols(grange)$source == "EDTA"],grange[mcols(grange)$type == "protein_coding"],ignore.strand=TRUE) == 0]
			other=grange[mcols(grange)$type == "rRNA" | mcols(grange)$type == "tRNA" | mcols(grange)$type == "miRNA" | mcols(grange)$type == "lncRNA" | mcols(grange)$type == "ncRNA" | mcols(grange)$type == "snoRNA" | mcols(grange)$type == "snRNA" | mcols(grange)$type == "ribozyme" | mcols(grange)$type == "nontranslating_CDS"]
			grange=grange[mcols(grange)$type == "bin1" | mcols(grange)$type == "bin2" | mcols(grange)$type == "bin3" | mcols(grange)$type == "bin4" | mcols(grange)$type == "geneBody" | mcols(grange)$type == "promoter"]
			allgenes <- as.data.frame(grange)%>%
				mutate(start=start-1)
			colnames(allgenes)[1] <- "chrom"
			overlaponly <- as.data.frame(disjoin(grange,ignore.strand=TRUE)[countOverlaps(disjoin(grange,ignore.strand=TRUE),grange,ignore.strand=TRUE) > 1])%>%
				mutate(start=start-1)
			colnames(overlaponly)[1] <- "chrom"
			nooverlap <- as.data.frame(disjoin(grange,ignore.strand=TRUE)[countOverlaps(disjoin(grange,ignore.strand=TRUE),grange,ignore.strand=TRUE) == 1])%>%
				mutate(start=start-1)
			colnames(nooverlap)[1] <- "chrom" 
			overlaponlyfeature =bed_intersect(overlaponly,allgenes)%>%
				filter(.overlap>0)%>%
				mutate(start=start.x+1,end=end.x)%>%
				group_by(chrom,start,end)%>%
				summarise(type=paste(type.y, collapse = ""))%>%
				mutate(type=gsub("bin[1-4]","geneBody",type))%>%
				mutate(type=gsub(".*geneBodypromoter.*","promotergeneBody",type))%>%
				mutate(type=gsub("geneBody.*Body","geneBody",type))%>%
				mutate(type=gsub("prom.*moter","promoter",type))%>%
				mutate(source="me")%>%
				GRanges()
			nooverlapfeature =bed_intersect(nooverlap,allgenes)%>%
			    filter(.overlap>0)%>%
			    mutate(start=start.x+1,end=end.x,strand=strand.y,source=source.y,type=type.y,score=score.y,phase=phase.y,transcript_id=transcript_id.y,gene_id=gene_id.y,class_id=class_id.y,order_id=order_id.y,family_id=family_id.y)%>%
			    dplyr::select(chrom,start,end,strand,source,type,score,phase,transcript_id,gene_id,class_id,order_id,family_id)%>%
			    GRanges()
			
			grangem=c(nooverlapfeature,overlaponlyfeature,tenooverlap,other)
			grangem
		}
		remove_te_in_gene=function(grange){
			tenooverlap <- grange[mcols(grange)$source == "EDTA"][countOverlaps(grange[mcols(grange)$source == "EDTA"],grange[mcols(grange)$type == "protein_coding"],ignore.strand=TRUE) == 0]
			tenooverlap%>%as.data.frame()%>%
		        mutate(start=start-1)%>%
		        select(c(seqnames,start,end,transcript_id,score,strand,source,type,phase))%>%
		        mutate(score=".",phase=".")%>%
		        write.table(file=paste0("~/cluster/smontgomery/genomes/",species,"/",species,".TE.nogene.bed"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
		}
		## adds intergenic to gff. intergenic= regions that are not in input gff (best to input gff with all genes, promoters, and TEs included)
		## output tibble with promoters
		# add_intergenic=function(gff){
		#   intergenic=bed_complement(gff,chromsizes)%>%mutate(feature="intergenic",strand="*",length=end-start)%>%filter(length > 1)%>%dplyr::select(-length)
		#   gffi = full_join(gff,intergenic)
		#   gffi
		# }
		add_intergenic=function(grange){
			tes=grange[mcols(grange)$source == "EDTA"]
			grange=grange[mcols(grange)$source != "EDTA"]
			interte=tes%>%mutate(strand="*")%>%gaps()%>%filter(strand=="*")%>%filter(mcols(distanceToNearest(.,grange))[,1]>0)
			tryCatch({
				mcols(interte)$type <- as.factor("interte")
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
			grange=c(grange,tes,interte)
			intergenic=grange%>%mutate(strand="*")%>%gaps()%>%filter(strand=="*")
			mcols(intergenic)$type <- as.factor("intergenic")
			grangei=c(grange,intergenic)
			grangei
		}
		## add gene expression level by top and bottom quintile
		add_expression=function(grange){
			tmpfile <- read.table(paste0("/home/smontgomery/cluster/smontgomery/genomes/",species,"/",species,"_expr_sorted.tsv"),header=TRUE)
			grangehigh <- split(grange, mcols(grange)$gene_id %in% tmpfile[1:round(nrow(tmpfile)/10),]$gene_id)
			grangehigh <- c(grangehigh$`FALSE`,grangehigh$`TRUE`%>%mutate(type=paste0(type,"High")))
			grangemed <- split(grangehigh, mcols(grangehigh)$gene_id %in% tmpfile[(round(nrow(tmpfile)/10)+1):(nrow(tmpfile)-round(nrow(tmpfile)/10)-1),]$gene_id)
			grangemed <- c(grangemed$`FALSE`,grangemed$`TRUE`%>%mutate(type=paste0(type,"Med")))
			grangelow <- split(grangemed, mcols(grangemed)$gene_id %in% tmpfile[(nrow(tmpfile)-round(nrow(tmpfile)/10)):nrow(tmpfile),]$gene_id)
			grangelow <- c(grangelow$`FALSE`,grangelow$`TRUE`%>%mutate(type=paste0(type,"Low")))
			grangelow
		}
		add_TSS=function(grange){
		    grange=grange[mcols(grange)$type == "protein_coding"]
		    TSS=pintersect(grange%>%mutate(type=as.factor("TSS"))%>%resize(1),grange)%>%filter(hit=="TRUE")%>%dplyr::select(-hit)
		    TSS
		}
		get_our_PCG=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.PCG.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_rRNA=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.rRNA.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_tRNA=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.tRNA.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_miRNA=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.miRNA.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_lncRNA=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.lncRNA.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_ncRNA=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.ncRNA.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_nontranslatingRNA=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.nontranslating.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_snoRNA=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.snoRNA.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_snRNA=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.snRNA.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_ribozyme=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.ribozyme.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		get_our_tRNA_pseudogene=function(crm=c("Pt","Mt")){
			grange = import.gff(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".gene.tRNA_pseudogene.gtf"))
			grange = c(grange,GRanges(seqlengths=setNames(read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V2,read.table(paste0("~/cluster/smontgomery/genomes/",species,"/",species,".chromsizes.txt"))$V1)))
			grange = grange[!(seqnames(grange) %in% c(crm))]
		}
		## extract protein coding genes from gff
		get_gene_pos=function(gff){
		  gene_pos=gff%>%filter(feature=="transcript")%>%
		    dplyr::select(chrom,start,end,strand,geneId)
		  gene_pos
		}
		## extract TE genes from gff
		get_te_gene_pos=function(gff){
		  te_gene_pos = gff%>%filter(Source=="EDTA")%>%
		    dplyr::select(chrom,start,end,strand,geneId)
		  te_gene_pos
		}

		## keep only features of interest, remove te gene exons and intro, rename the features and relevel them to give preferred order
		gff_to_feature=function(gff){
		    feature=gff%>%filter(feature%in%c("promoter","bin1","bin2","bin3","bin4","geneBody",
		                                      "TE_RT_LTR_intact","TE_RT_LTR_fragment","TE_RT_LINE",,"TE_RT_SINE",
		                                      "TE_DNA_intact","TE_DNA_fragment","intergenic"))%>%
		        mutate(feature=gsub("promoter","Promoter",feature))%>%
		        mutate(feature=gsub("bin1","Bin 1",feature))%>%
		        mutate(feature=gsub("bin2","Bin 2",feature))%>%
		        mutate(feature=gsub("bin3","Bin 3",feature))%>%
		        mutate(feature=gsub("bin4","Bin 4",feature))%>%
		        mutate(feature=gsub("geneBody","Gene body",feature))%>%
		        mutate(feature=gsub("TE_RT_LTR_intact","TE RT LTR intact",feature))%>%
		        mutate(feature=gsub("TE_RT_LTR_fragment","TE RT LTR fragment",feature))%>%
		        mutate(feature=gsub("TE_RT_LINE","TE RT LINE",feature))%>%
		        mutate(feature=gsub("TE_RT_SINE","TE RT SINE",feature))%>%
		        mutate(feature=gsub("TE_DNA_intact","TE DNA intact",feature))%>%
		        mutate(feature=gsub("TE_DNA_fragment","TE DNA fragment",feature))%>%
		        mutate(feature=gsub("intergenic","Intergenic",feature))%>%
		        mutate(feature=factor(feature,levels=c("Promoter","Bin 1","Bin 2","Bin 3","Bin 4","Gene body","TE RT LTR intact", "TE RT LTR fragment",
		                                               "TE RT LINE","TE RT SINE","TE DNA intact","TE DNA fragment","Intergenic")))
		    te_gene_id=filter(feature,feature=="TE gene")$geneId
		    feature = feature%>%filter((feature%in%c("exon","intro"))+geneId%in%te_gene_id<2)
		    feature
		}
		# grange_to_feature=function(grange){
		#     feature=as_tibble(grange)[,c(1,6,7,2,3,8,5,9,10,11,12,13,14)]
		#     colnames(feature)=c("chrom","Source","feature","start","end","X6","strand","X8","geneId","gene_id","class_id","order_id","family_id")
		#     feature=feature%>%mutate(start=start-1)%>%
		#         mutate(feature=gsub("bin1","Bin 1",feature))%>%
		#         mutate(feature=gsub("bin2","Bin 2",feature))%>%
		#         mutate(feature=gsub("bin3","Bin 3",feature))%>%
		#         mutate(feature=gsub("bin4","Bin 4",feature))%>%
		        # mutate(feature=gsub("TE_RT_LTR_intact","TE RT LTR intact",feature))%>%
		        # mutate(feature=gsub("TE_RT_LTR_fragment","TE RT LTR fragment",feature))%>%
		        # mutate(feature=gsub("TE_RT_LINE","TE RT LINE",feature))%>%
		        # mutate(feature=gsub("TE_RT_SINE","TE RT SINE",feature))%>%
		        # mutate(feature=gsub("TE_DNA_intact","TE DNA intact",feature))%>%
		        # mutate(feature=gsub("TE_DNA_fragment","TE DNA fragment",feature))%>%
		        # mutate(feature=gsub("intergenic","Intergenic",feature))%>%
		        # mutate(feature=gsub("interte","Inter-TE",feature))%>%
		        # mutate(feature=gsub("promotergeneBody","Promoter + Gene body",feature))%>%
		        # mutate(feature=gsub("promoter","Promoter",feature))%>%
		        # mutate(feature=gsub("geneBody","Gene body",feature))%>%
		        # mutate(feature=gsub("promoter","Misc Promoter",feature))%>%
		        # mutate(feature=gsub("geneBody","Misc Gene body",feature))%>%
		        # mutate(feature=factor(feature,levels=c("Promoter","Bin 1","Bin 2","Bin 3","Bin 4","Gene body","TE RT LTR intact", "TE RT LTR fragment",
		        #                                        "TE RT LINE","TE RT SINE","TE DNA intact","TE DNA fragment","Intergenic","Inter-TE","Promoter + Gene body","Misc Promoter","Misc Gene body")))
		#         feature
		# }
		# gff_to_feature=function(gff){
		#     feature=gff%>%filter(feature%in%c("promoterHigh","bin1High","bin2High","bin3High","bin4High","geneBodyHigh",
		#                                       "promoterMed","bin1Med","bin2Med","bin3Med","bin4Med","geneBodyMed",
		#                                       "promoterLow","bin1Low","bin2Low","bin3Low","bin4Low","geneBodyLow",
		#                                       "TE_RT_LTR_intact","TE_RT_LTR_fragment","TE_RT_LINE",,"TE_RT_SINE",
		                                      # "TE_DNA_intact","TE_DNA_fragment","intergenic"))%>%
		#         mutate(feature=gsub("promoterHigh","Promoter High",feature))%>%
		#         mutate(feature=gsub("bin1High","Bin 1 High",feature))%>%
		#         mutate(feature=gsub("bin2High","Bin 2 High",feature))%>%
		#         mutate(feature=gsub("bin3High","Bin 3 High",feature))%>%
		#         mutate(feature=gsub("bin4High","Bin 4 High",feature))%>%
		#         mutate(feature=gsub("geneBodyHigh","Gene body High",feature))%>%
		#         mutate(feature=gsub("promoterMed","Promoter Med",feature))%>%
		#         mutate(feature=gsub("bin1Med","Bin 1 Med",feature))%>%
		#         mutate(feature=gsub("bin2Med","Bin 2 Med",feature))%>%
		#         mutate(feature=gsub("bin3Med","Bin 3 Med",feature))%>%
		#         mutate(feature=gsub("bin4Med","Bin 4 Med",feature))%>%
		#         mutate(feature=gsub("geneBodyMed","Gene body Med",feature))%>%
		#         mutate(feature=gsub("promoterLow","Promoter Low",feature))%>%
		#         mutate(feature=gsub("bin1Low","Bin 1 Low",feature))%>%
		#         mutate(feature=gsub("bin2Low","Bin 2 Low",feature))%>%
		#         mutate(feature=gsub("bin3Low","Bin 3 Low",feature))%>%
		#         mutate(feature=gsub("bin4Low","Bin 4 Low",feature))%>%
		#         mutate(feature=gsub("geneBodyLow","Gene body Low",feature))%>%
		        # mutate(feature=gsub("TE_RT_LTR_intact","TE RT LTR intact",feature))%>%
		        # mutate(feature=gsub("TE_RT_LTR_fragment","TE RT LTR fragment",feature))%>%
		        # mutate(feature=gsub("TE_RT_LINE","TE RT LINE",feature))%>%
		        # mutate(feature=gsub("TE_RT_SINE","TE RT SINE",feature))%>%
		        # mutate(feature=gsub("TE_DNA_intact","TE DNA intact",feature))%>%
		        # mutate(feature=gsub("TE_DNA_fragment","TE DNA fragment",feature))%>%
		#         mutate(feature=gsub("intergenic","Intergenic",feature))%>%
		#         mutate(feature=factor(feature,levels=c("Promoter High","Bin 1 High","Bin 2 High","Bin 3 High","Bin 4 High","Gene body High",
		#                                                "Promoter Med","Bin 1 Med","Bin 2 Med","Bin 3 Med","Bin 4 Med","Gene body Med",
		#                                                "Promoter Low","Bin 1 Low","Bin 2 Low","Bin 3 Low","Bin 4 Low","Gene body Low",
		# 												 "TE RT LTR intact", "TE RT LTR fragment","TE RT LINE","TE RT SINE","TE DNA intact","TE DNA fragment","Intergenic")))
		#     te_gene_id=filter(feature,feature=="TE gene")$geneId
		#     feature = feature%>%filter((feature%in%c("exon","intro"))+geneId%in%te_gene_id<2)
		#     feature
		# }
		grange_to_feature=function(grange){
		    feature=as_tibble(grange)[,c(1,6,7,2,3,8,5,9,10,11,12,13,14)]
		    colnames(feature)=c("chrom","Source","feature","start","end","X6","strand","X8","geneId","gene_id","class_id","order_id","family_id")
		    feature=feature%>%mutate(start=start-1)%>%
		        mutate(feature=gsub(".*([0-9]+)$","geneBody",feature))%>%
		        mutate(feature=gsub("promoterHigh","Promoter High",feature))%>%
		        mutate(feature=gsub("bin1High","Bin 1 High",feature))%>%
		        mutate(feature=gsub("bin2High","Bin 2 High",feature))%>%
		        mutate(feature=gsub("bin3High","Bin 3 High",feature))%>%
		        mutate(feature=gsub("bin4High","Bin 4 High",feature))%>%
		        mutate(feature=gsub("geneBodyHigh","Gene body High",feature))%>%
		        mutate(feature=gsub("promoterMed","Promoter Med",feature))%>%
		        mutate(feature=gsub("bin1Med","Bin 1 Med",feature))%>%
		        mutate(feature=gsub("bin2Med","Bin 2 Med",feature))%>%
		        mutate(feature=gsub("bin3Med","Bin 3 Med",feature))%>%
		        mutate(feature=gsub("bin4Med","Bin 4 Med",feature))%>%
		        mutate(feature=gsub("geneBodyMed","Gene body Med",feature))%>%
		        mutate(feature=gsub("promoterLow","Promoter Low",feature))%>%
		        mutate(feature=gsub("bin1Low","Bin 1 Low",feature))%>%
		        mutate(feature=gsub("bin2Low","Bin 2 Low",feature))%>%
		        mutate(feature=gsub("bin3Low","Bin 3 Low",feature))%>%
		        mutate(feature=gsub("bin4Low","Bin 4 Low",feature))%>%
		        mutate(feature=gsub("geneBodyLow","Gene body Low",feature))%>%
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
		        mutate(feature=factor(feature,levels=c("Gene","Promoter High","Bin 1 High","Bin 2 High","Bin 3 High","Bin 4 High","Gene body High",
		                                               "Promoter Med","Bin 1 Med","Bin 2 Med","Bin 3 Med","Bin 4 Med","Gene body Med",
		                                               "Promoter Low","Bin 1 Low","Bin 2 Low","Bin 3 Low","Bin 4 Low","Gene body Low",
		                                               "TE RT LTR intact", "TE RT LTR fragment","TE RT LINE","TE RT SINE","TE RT DIRS","TE DNA intact","TE DNA fragment",
		                                               "Intergenic","Inter-TE","Promoter + Gene body","Misc Promoter","Misc Gene body",
		                                               "rRNA","tRNA","miRNA","lncRNA","ncRNA","sRNA","snoRNA","snRNA","ribozyme","SRP_RNA","RNase_MRP_RNA","sense_intronic","telomerase_RNA","nontranslating_CDS","pseudogene","misc_RNA")))
		}
		## keep only rna features of interest, relevel them to give preferred order
		gff_to_rnatype=function(gff){
		  feature=gff%>%filter(feature%in%c("mRNA","snoRNA","ncRNA","snRNA","rRNA",
		                                    "lnc_RNA","miRNA","antisense_lncRNA",
		                                    "antisense_RNA","tRNA",
		                                    "pseudogenic_tRNA"))%>%
		    mutate(feature=factor(feature,levels=c("mRNA","snoRNA","ncRNA","snRNA",
		                                           "rRNA","lnc_RNA","miRNA",
		                                           "antisense_lncRNA","antisense_RNA",
		                                           "tRNA","pseudogenic_tRNA")))
		  feature
		}  

		## how to set range when plotting without outliers  
		find_outlier_range <- function(x) {
		  return( quantile(x, .75,na.rm=T) + 1.5*IQR(x,na.rm=T))
		}

		## plot distance from TSS to each state within genes
		find_gene_states <- function(grange) {
			tmpgrange=as_tibble(grange)[,c(1,6,7,2,3,8,5,9,10,11,12,13,14)]
		    colnames(tmpgrange)=c("chrom","Source","feature","start","end","X6","strand","X8","geneId","gene_id","class_id","order_id","family_id")
		    tmpgrange=tmpgrange%>%mutate(start=start-1)
			genestates=bed_intersect(filter(tmpgrange,feature=="protein_coding"),bed_subtract(large_col,bed_subtract(large_col,filter(tmpgrange,feature=="protein_coding"))))
			genestate_strand=list()
			  for (strand in c("+","-")){
			  	if (strand == "-"){
				    genestate_strand[[strand]]=filter(genestates,strand.x==!!strand)%>%mutate(distance=end.x-end.y)%>%dplyr::select(state.y,cols.y,distance)%>%mutate(distance = if_else(distance < 0, 0, distance))
				}
				else if (strand == "+"){
					genestate_strand[[strand]]=filter(genestates,strand.x==!!strand)%>%mutate(distance=start.y-start.x)%>%dplyr::select(state.y,cols.y,distance)%>%mutate(distance = if_else(distance < 0, 0, distance))
				}
			  }
			  genestate_distances = do.call("rbind",genestate_strand)
			  colnames(genestate_distances)=c("state","cols","distance")
			  genestate_distances
		}
		# find_gene_states <- function(gff) {
		# 	genestates=bed_intersect(filter(gff,feature=="transcript"),bed_subtract(large_col,bed_subtract(large_col,filter(gff,feature=="transcript"))))
		# 	genestate_strand=list()
		# 	  for (strand in c("+","-")){
		# 	  	if (strand == "-"){
		# 		    genestate_strand[[strand]]=filter(genestates,strand.x==!!strand)%>%mutate(distance=end.x-end.y)%>%dplyr::select(state.y,cols.y,distance)%>%mutate(distance = if_else(distance < 0, 0, distance))
		# 		}
		# 		else if (strand == "+"){
		# 			genestate_strand[[strand]]=filter(genestates,strand.x==!!strand)%>%mutate(distance=start.y-start.x)%>%dplyr::select(state.y,cols.y,distance)%>%mutate(distance = if_else(distance < 0, 0, distance))
		# 		}
		# 	  }
		# 	  genestate_distances = do.call("rbind",genestate_strand)
		# 	  colnames(genestate_distances)=c("state","cols","distance")
		# 	  genestate_distances
		# }
		# Get lower triangle of the correlation matrix
		get_lower_tri<-function(tmp){
		    tmp[upper.tri(tmp)] <- NA
		    return(tmp)
		}
		# Get upper triangle of the correlation matrix
		get_upper_tri <- function(tmp){
		    tmp[lower.tri(tmp)]<- NA
		    return(tmp)
		}


chi_sq_test_loop=function(teregions=teregions){
	df=as.matrix(table(teregions[,c("V13","order_id")]))
	mat_ps_chisq=as.data.frame(matrix(ncol=ncol(df),nrow=nrow(df)-1));colnames(mat_ps_chisq)=colnames(df);rownames(mat_ps_chisq)=rownames(head(df,-1))

	for(tecluster in rownames(df)){
		for (teorder in colnames(df)){ 
			niche_values=c(df[tecluster,teorder],sum(df[tecluster,which(colnames(df)!=teorder)]))
			other_sn_values=c(sum(df[which(rownames(df)!=tecluster),teorder]),sum(df[which(rownames(df)!=tecluster),which(colnames(df)!=teorder)]))
			matrix=rbind(niche_values,other_sn_values);colnames(matrix)=c(teorder,"other teorder");rownames(matrix)=c(tecluster,"other tecluster")      
			#print(matrix)
			mat_ps_chisq[tecluster,teorder]=fisher.test(matrix,alternative="greater")$p.value
		}
	}   
	mat_ps_chisq=as.matrix(mat_ps_chisq)
	mat_ps_chisq[is.nan(mat_ps_chisq)] <- 1
	a1 = mat_ps_chisq %>% melt() 
	colnames(a1)=c("tecluster","teorder","qval")
	a= a1 %>% mutate(qval = p.adjust(qval,method = "BH")) %>% dcast(tecluster ~ teorder, value.var='qval')
	mat_ps_chisq_qv=a[,-1];rownames(mat_ps_chisq_qv)=a[,1]
	return(mat_ps_chisq_qv)
}
 
 