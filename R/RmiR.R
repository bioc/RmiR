RmiR <-
function(mirna=NULL,genes=NULL,annotation=NULL,dbname="targetscan",id="probes",id.out="symbol",verbose=FALSE)
{
	if(is.null(c(mirna,genes))) stop("missing mir and genes input")
	if(ncol(genes)!=2 | ncol(mirna)!=2) stop("Both files must have two colums!")
	if(is.null(annotation)) stop("You have to specify an annotation package for the microarray platform or the organism you are working with.")
	if(id.out=="probes" & id!="probes") stop("If you want probe IDs in the output you should have probe IDs also in the input.")
	names(mirna)=c("mature_miRNA","mirExpr")
	names(genes)[2]="geneExpr"
	if(id.out=="probes") (names(genes)[1]="probe_id")
	### Annotate the Genes File ###
	if (!is.null(annotation)){
   	   require (annotation,character.only=TRUE) || stop(paste(annotation, "package must be installed first"))
	con_comm <- paste(as.character(strsplit(annotation,split=".db")),"dbconn",sep="_")
	usedDB <- eval(as.name(con_comm))
	con <- usedDB()
	}
	if (id=="genes")
		{
		names(genes)[1] <- "gene_id"
		annot.match ="gene_id"
		}
	if (id=="probes")
		{
		annot.match ="probe_id"
		}
	if (id=="ensembl")
		{
		annot.match ="ensembl_id"
		}
	if (id=="unigene")
		{
		annot.match ="unigene_id"
		}
	if (id=="alias")
		{
		annot.match ="alias_symbol"
		}

	annot <- dbReadTable(con,id)
	if (id!="genes")
		{
		annot.e <- dbReadTable(con,"genes")
		} else  annot.e <- annot

	if (id!="genes")
		{
		annot <- merge(annot,annot.e)[,c(annot.match,"gene_id")]
		genes <- na.exclude(merge(x=genes,y=annot,by.x=names(genes)[1],by.y=annot.match,all.x=T,all.y=F))
		}

	if (id.out=="symbol")
		{		
		annot.s <- dbReadTable(con,"gene_info")[,c("X_id","symbol")]
		annot <- merge(annot.e,annot.s)
#		genes <- na.exclude(merge(x=genes,y=annot,by.x=names(genes)[1],by.y=annot.match,all.x=T,all.y=F))
		genes <- na.exclude(merge(annot,genes))
		}
	if(nrow(genes)==0) stop ("May be you have not a file with the annotation in the first column, or the id variable is wrong. Please adjust and try again!")
	### mean of replicated or duplicated genes

	if (id.out!="probes")
		{
		if (id.out=="symbol")
			{
			genes <- genes[,c("symbol","gene_id","geneExpr")]
			} else genes <- genes[,c("gene_id","geneExpr")]
		genes <- genes[order(genes$gene_id),]
		dup <- which(duplicated(genes$gene_id) | duplicated(genes$gene_id,fromLast=T)== TRUE)
		reps <- table(as.vector(genes$gene_id[dup]))
		rN <- names(reps)
		if (id.out=="symbol")
			{
			symb.reps <- table(as.vector(genes$symbol[dup]))
			name.symb <- names(symb.reps)
			}
		Lreps <- length(rN)
		if (Lreps != 0) 
			{
			MEDmat <- as.data.frame(rN)
			colnames(MEDmat)<-"gene_id"
			if (id.out=="symbol")
                        	{
				MEDmat$symbol <- name.symb
                        	}
			MEDmat$geneExpr <- 0 
			MEDmat$geneCV <- NA
			genes$geneCV <- NA
 			for (i in 1:Lreps) 
				{
				index = which(genes$gene_id== rN[i])
				MED = apply(as.data.frame(genes$geneExpr[index]), 2, median)
 			#	tvalue <- function(x)
 			#		{
 			#		nx <- length(x)
			#		tx <- mean(x)/sqrt(1/(nx*(nx-1))*sum((x-mean(x))^2))
			#		return(tx)
  			#		}
			#	pvalue <- function (x)
			#		{
			#		2*pt(-abs(tvalue(x)),df=(length(x)-1))
			#		}
			#	PV = apply(as.data.frame(genes$geneExpr[index]), 2, pvalue)
			cvres <- function(x)
				{
				abs(sd(x)/mean(x))
				}
				CV = apply(as.data.frame(genes$geneExpr[index]), 2, cvres)

				MEDmat[i,"geneExpr"] = MED
				MEDmat[i,"geneCV"] = CV
       				}
			genes <- genes[-dup,]
			genes <- rbind(genes,MEDmat) 
			}
		} else genes <- genes[,c("probe_id","gene_id","geneExpr")]
	###  mean of replicated or duplicated miRNAs
	mirna <- mirna[order(mirna$mature_miRNA),]
        mdup <- which(duplicated(mirna$mature_miRNA) | duplicated(mirna$mature_miRNA,fromLast=T)== TRUE)
        mreps <- table(as.vector(mirna$mature_miRNA[mdup]))
        mrN <- names(mreps)
        mLreps <- length(mrN)
        if (mLreps != 0) 
		{
        	mMEDmat <- as.data.frame(mrN)
        	colnames(mMEDmat)<-"mature_miRNA"
        	mMEDmat$mirExpr <- 0
        	mMEDmat$mirCV <- NA
        	mirna$mirCV <- NA
        	for (m in 1:mLreps)
                	{
                	mindex = which(mirna$mature_miRNA== mrN[m])
                	mMED = apply(as.data.frame(mirna$mirExpr[mindex]), 2, median)
                #	tvalue <- function(x)
                #		{
                #		nx <- length(x)
                #		tx <- mean(x)/sqrt(1/(nx*(nx-1))*sum((x-mean(x))^2))
                #		return(tx)
                #       	}
                #	pvalue <- function (x)
                #        	{
                #        	2*pt(-abs(tvalue(x)),df=(length(x)-1))
                #        	}
               	#	mCV = apply(as.data.frame(mirna$mirExpr[mindex]), 2, pvalue)
			cvres <- function(x)
                                {
                                abs(sd(x)/mean(x))
                                }
                                mCV = apply(as.data.frame(mirna$mirExpr[mindex]), 2, cvres)

                	mMEDmat[m,"mirExpr"] = mMED
                	mMEDmat[m,"mirCV"] = mCV
		       		}
          		mCV = apply(as.data.frame(mirna$mirExpr[mindex]), 2, cvres)

        	mirna <- mirna[-mdup,]
		mirna <- rbind(mirna,mMEDmat)
        	}
	###
	### Selecting the miRNA database and reduce the data as needed ### 	
	miRDBs <- RmiR_dbconn()
	mirs_targets <- dbReadTable(miRDBs,dbname)
	mirs_targets <- mirs_targets[mirs_targets$gene_id%in%genes$gene_id & mirs_targets$mature_miRNA%in%mirna$mature_miRNA, 1:2]
	found_gene <- length(unique(na.exclude(mirs_targets$gene_id)))
	found_mirs <- length(unique(na.exclude(mirs_targets$mature_miRNA)))
	if (verbose==TRUE){
	warn1 <- paste("In ",dbname," database there are ",found_gene," genes and ",found_mirs," microRNA which are in your files.",sep="")
#	cat("          --------------------------------------------")
	cat("\n")	
	cat(warn1)
	cat("\n")
	cat("          --------------------------------------------")
	cat("\n")
	}
	###
	### Remove not present miRNA and Genes from input files ###
	genes <- genes[genes$gene_id%in%mirs_targets$gene_id,]
	mirna <- mirna[mirna$mature_miRNA%in%mirs_targets$mature_miRNA,]
	###
	### Merging results
	tot <- merge(mirna,mirs_targets)
	tot <- merge(tot,genes)
	unique(tot)
}

