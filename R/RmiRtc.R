RmiRtc <- function(timeline=NULL,timevalue=NULL,method="pearson")
{
	if(length(timeline)!=length(timeline)) stop("timeline and timevalue must have the same length")
	if(!is.numeric(timevalue)) stop("timevalue must be a numeric vector")
	TCtab <- eval(as.name(timeline[1]))[,c("gene_id","mature_miRNA")]
	for(i in 2:length(timeline))
		{
		TCtab <- unique(rbind(TCtab,eval(as.name(timeline[i]))[,c("gene_id","mature_miRNA")]))
		}

	TCtab<- TCtab[order(TCtab$gene_id),]
	TCtab<- TCtab[order(TCtab$mature_miRNA),]
	mirTC <- new("miRtcList")

	mirTC$couples <- TCtab
	mirTC$mirExpr <- matrix(ncol=length(timeline),nrow=nrow(mirTC$couples))
	colnames(mirTC$mirExpr) <- timevalue
	mirTC$geneExpr <-matrix(ncol=length(timeline),nrow=nrow(mirTC$couples))
	colnames(mirTC$geneExpr) <- timevalue
	mirTC$mirCV  <-matrix(ncol=length(timeline),nrow=nrow(mirTC$couples))
	colnames(mirTC$mirCV) <- timevalue
	mirTC$geneCV <-matrix(ncol=length(timeline),nrow=nrow(mirTC$couples))
	colnames(mirTC$geneCV) <- timevalue
	mirTC$correlation <- vector(length=nrow(mirTC$couples))
	for (n in 1:length(timeline))
		{
		tmp <- eval(as.name(timeline[n]))
		tmp <- tmp[order(tmp$gene_id),]
		tmp <- tmp[order(tmp$mature_miRNA),]
		isok <- which(mirTC$couples$gene_id%in%tmp$gene_id &  mirTC$couples$mature_miRNA%in%tmp$mature_miRNA)
	
		mirTC$mirExpr[isok,n] <- tmp$mirExpr
		mirTC$geneExpr[isok,n] <- tmp$geneExpr
		if(!is.null(tmp$mirCV))
			{
			mirTC$mirCV[isok,n] <- tmp$mirCV
			} else mirTC$mirCV[isok,n] <- 0
		if(!is.null(tmp$geneCV))
			{
			mirTC$geneCV[isok,n] <- tmp$geneCV
			} else mirTC$gene[isok,n] <- 0
		}
	
	corr_mir.gene <- function(x,type=NULL)
		{
		cor(mirTC$mirExpr[x,],mirTC$geneExpr[x,],method=type)
		}

	mirTC$correlation<- apply(t(c(1:length(mirTC$correlation))),2,corr_mir.gene,type=method)



	return(mirTC)
}
