 
chi_sq_test_loop=function(teregions=teregions){
	# footprint=mc_object@mc_fp
	# if(!is.null(niche_order)){
	# 	footprint=footprint[,niche_order]
	# }
	# #genes=names(which(rowSums(as.matrix(mat_object@mat[intersect(rownames(footprint),rownames(mat_object@mat)),])) > T_totumis_seen))
	# #f_genes_fc=names(which(apply(footprint[genes,],1,max) > fc_threshold))
	# PS_table=read.table(PS_table_file,h=FALSE,row.names=1,stringsAsFactors=F);colnames(PS_table)="PS" 
	# genes=intersect(rownames(footprint),rownames(PS_table))
	# footprint=footprint[genes,]
	# PS_table=PS_table[genes,,drop=F];colnames(PS_table)="PS" 
	# message("N seen genes: ",length(genes))

	# df=matrix(ncol=length(unique(PS_table$PS)),nrow=ncol(footprint))
	# colnames(df)=PS_order;rownames(df)=colnames(footprint)
	# Niche_genes=c()  #we compile here genes that are considered as part of a niche gene module (used later to define backgroud)
	# Niche_genes_list=list() ##we want to return this list to then explore which genes correspond to each enrichment
	# for (niche in colnames(footprint)){
	# 	if(ct_unique_genes==TRUE) { 
	# 		Genes=unique(names(which(footprint[,niche]>fc_threshold & apply(footprint[,setdiff(colnames(footprint),niche)],1,max) < fc_threshold)))
	# 	}else{
	# 		Genes=unique(names(which(footprint[,niche]>fc_threshold)))
	# 	}
	# 	Niche_genes_list[[niche]]=Genes
	# 	df[niche,]=table(factor(PS_table[Genes,"PS"],levels=PS_order))
	# 	Niche_genes=unique(c(Niche_genes,Genes))
	# }  
	# IDs_OTHER=setdiff(genes,unique(Niche_genes))
	# df=rbind(df,table(factor(PS_table[IDs_OTHER,"PS"],levels=PS_order)))
	# rownames(df)=c(colnames(footprint),"other_genes")
	# df=df[,PS_order] #we chose a custom order of phylostrata ("Eukaryota", etc)
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
	# df_freq=apply(df,2,function(x) x/rowSums(df))
	# df_freq=df_freq[!grepl("Unclass*", rownames(df_freq)),]
	# bck_freq=table(factor(PS_table$PS,levels=PS_order))/nrow(PS_table)  ##phylostratigraphy table previously REDUCED to only genes present in the footprint
	# mat=apply(df_freq,1,function(x) log2((x)/(bck_freq)))
	# mat[which(!is.finite(mat))]=-2
	# return(list(fg_freq=t(df_freq),bg_freq=bck_freq,fg_counts=df,all_genes=genes,plot_matrix=mat,niche_gene_list=Niche_genes_list))
}
 
 