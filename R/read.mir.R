read.mir <- 
function(mirna=NULL, genes=NULL, annotation=NULL, id="probes", dbname=c("targetscan","pictar"), org="hsa", at.least=2,id.out="symbol", verbose=FALSE)
	{
	tot <- RmiR(genes=genes,mirna=mirna,annotation=annotation,id=id,dbname=dbname[1],id.out=id.out,verbose=verbose,org=org)
	if (length(dbname)>1)
		{	
		if (is.null(tot$mirCV)) tot$mirCV = NA
		if (is.null(tot$geneCV)) tot$geneCV = NA

		for(i in 2:length(dbname))
			{
			tmp <- RmiR(genes=genes,mirna=mirna,annotation=annotation,id=id,dbname=dbname[i],id.out=id.out,verbose=verbose,org=org)
			if (nrow(tmp)!=0)
				{
				if (is.null(tmp$mirCV)) tmp$mirCV = NA
        			if (is.null(tmp$geneCV)) tmp$geneCV = NA
				tot <- rbind(tot,tmp)
				}
			}
		tot_vect<- paste(tot$gene_id,tot$mature_miRNA,sep="_")
		reps <- table(as.vector(tot_vect))>=at.least
		idx <- tot_vect%in%names(reps[reps])
		tot <- unique(tot[idx,])
		} 
	return(tot)
	}
