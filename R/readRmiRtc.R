readRmiRtc <- 
function(miRtcObj,correlation=-0.75,exprLev=1,annotation=NULL,fileName="miRNA_for_genes")
	{
        if(is.null(annotation)) stop("You have to specify an annotation package for the microarray platform or the organism you are working with.")
	if (class(miRtcObj)!="miRtcList") stop("The input must be an object of 'miRtcList' class.")
	if (correlation >=0)
		{
		filt.list <- which(miRtcObj$correlation >= correlation)
		} 
	if (correlation <=0)
		{
		filt.list <- which(miRtcObj$correlation <= correlation)
		}

	if(length(filt.list)==0) stop("There are not couples miRNA/targets matching your criteria.")

	miRtcObj$correlation <- miRtcObj$correlation[filt.list]
	miRtcObj$couples <- miRtcObj$couples[filt.list,]
	miRtcObj$mirExpr <- as.data.frame(miRtcObj$mirExpr)[filt.list,]
	miRtcObj$geneExpr <- as.data.frame(miRtcObj$geneExpr)[filt.list,]
	miRtcObj$mirCV <- as.data.frame(miRtcObj$mirCV)[filt.list,]
	miRtcObj$geneCV <- as.data.frame(miRtcObj$geneCV)[filt.list,]
	
        ### Filter modulated genes
	gene.mod <- which(abs(miRtcObj$geneExpr[,ncol(miRtcObj$geneExpr)]) >= exprLev)
	
	if(length(gene.mod)==0) stop("There are not couples miRNA/targets matching your criteria.")

	miRtcObj$correlation <- miRtcObj$correlation[gene.mod]
        miRtcObj$couples <- miRtcObj$couples[gene.mod,]
        miRtcObj$mirExpr <- miRtcObj$mirExpr[gene.mod,]
        miRtcObj$geneExpr <- miRtcObj$geneExpr[gene.mod,]
        miRtcObj$mirCV <- miRtcObj$mirCV[gene.mod,]
        miRtcObj$geneCV <- miRtcObj$geneCV[gene.mod,]

	### Order by gene_id
	
	ord.gene <- order(miRtcObj$couples$gene_id)

	miRtcObj$correlation <- miRtcObj$correlation[ord.gene]
        miRtcObj$couples <- miRtcObj$couples[ord.gene,]
        miRtcObj$mirExpr <- miRtcObj$mirExpr[ord.gene,]
        miRtcObj$geneExpr <- miRtcObj$geneExpr[ord.gene,]
        miRtcObj$mirCV <- miRtcObj$mirCV[ord.gene,]
        miRtcObj$geneCV <- miRtcObj$geneCV[ord.gene,]

	### Count the genes with more miRNA and write the file
        if (!is.null(annotation)){
   	   require (annotation,character.only=TRUE) || stop(paste(annotation, "package must be installed first"))
	con_comm <- paste(as.character(strsplit(annotation,split=".db")),"dbconn",sep="_")
        usedDB <- eval(as.name(con_comm))
        con <- usedDB()
	}
	annot.e <- dbReadTable(con,"genes")
	annot <- dbReadTable(con,"gene_info")[,c("X_id","symbol")]
	annot <- merge(annot.e,annot)
	
	rep.genes <- table(miRtcObj$couples$gene_id)
	rep.data <- as.data.frame(names(rep.genes))
	rep.data$reps <- rep.genes
	names(rep.data) <- c("gene_id","miRNAs")

	rep.data<- merge(rep.data,annot)
	rep.data <- rep.data[order(rep.data$miRNAs,decreasing=T),c("symbol","miRNAs","gene_id")]
	if (!is.null(fileName))
		{
		write.table(rep.data,col.names=T,row.names=F,sep="\t",file=as.character(paste(fileName,"txt",sep=".")))
		}
	miRtcObj$reps <- rep.data
	miRtcObj
	}

plotRmiRtc <-
function(miRtcObj,gene_id=NULL,timeunit="Time",legend.x=NULL,legend.y=NULL,svgTips=FALSE,svgname=NULL,height=10,width=15)
	{
	if(class(miRtcObj)=="miRtcList")
		{
		gname <- miRtcObj$reps$symbol[miRtcObj$reps$gene_id==gene_id]
		mname <- as.vector(miRtcObj$couples$mature_miRNA[miRtcObj$couples$gene_id==gene_id])
		title <- paste(gname,"and its miRNAs expression trends",sep=" ")
		if(svgTips==TRUE) 
			{
			require(RSVGTipsDevice)
			if (is.null(svgname)) 
				{
				svgname <- paste("Gene",gname,"mirnas.svg",sep="_")
				}
			}
		timevalue <- as.vector(colnames(miRtcObj$geneExpr))
		geneE <- as.vector(unique(miRtcObj$geneExpr[miRtcObj$couples$gene_id==gene_id,]))
		mirE <- as.matrix(miRtcObj$mirExpr[miRtcObj$couples$gene_id==gene_id,])
		ymin <- min(min(mirE),min(geneE))	
		ymax <- max(max(mirE),max(geneE))
		if(svgTips==TRUE) devSVGTips(as.character(svgname),height=height,width=width,toolTipMode=1,title=title)	
		plot(x=timevalue,y=geneE,type="l",col="blue",ylab="Expression",xlab=timeunit,ylim=c(ymin,ymax),main=title)
		if(svgTips==TRUE)
			{
			for (g in 1:length(geneE))
				{
				setSVGShapeToolTip(title=as.name(miRtcObj$reps$symbol[miRtcObj$reps$gene_id == gene_id]),desc=paste("Entrez ID:",gene_id,sep=" "))
		       		geneURL=paste("http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch=",gene_id,sep="")
        			setSVGShapeURL(url=geneURL,target="_blank")
			       	points(x=timevalue[g],y=geneE[g],col="blue")
				}
			}
		for (i in 1:nrow(mirE))
			{
			mirEi <- mirE[i,]
			lines(x=timevalue,y=mirEi,col = "red")
			if(svgTips==TRUE)
		      	        {
				for (m in 1:length(mirEi))
					{
					mirURL <- paste("http://microrna.sanger.ac.uk/cgi-bin/sequences/query.pl?terms=",as.vector(miRtcObj$couples$mature_miRNA[miRtcObj$couples$gene_id == gene_id])[i],sep="")
	        	       		setSVGShapeToolTip(title=as.name(as.vector(miRtcObj$couples$mature_miRNA[miRtcObj$couples$gene_id == gene_id])[i]),desc="click to go on miRBase")
	              			setSVGShapeURL(url=mirURL,target="_blank")
	               			points(x=timevalue[m],y=mirEi[m],col = "red",type="b")
					}
				}
			}
		if(!is.null(legend.x)&!is.null(legend.y))
			{
			col <- c("blue",rep("red",length(mname)))
			legend(x=legend.x,y=legend.y,c(gname,mname),text.col=col)
			}
		 if(svgTips==TRUE) dev.off()
		} else {
			if (svgTips==FALSE) stop("plotRmiR for not 'miRtcObj' must have the 'svgTips' flag enabled")
			if (is.null(svgname)) stop("svgname must be specified")
			require(RSVGTipsDevice)
			devSVGTips(as.name(svgname),height=height,width=width,toolTipMode=2,title="microRNA and relative targets expressions")
			plot(miRtcObj$geneExpr,miRtcObj$mirExpr,ylab="miRNA expression",xlab="Gene Target Expression")
			for (l in 1:nrow(miRtcObj))
				{
				if (is.vector(miRtcObj$symbol))
					{
					geneame <- miRtcObj[l,"symbol"]
					} else	geneame <- miRtcObj[l,"gene_id"]
				setSVGShapeToolTip(title=miRtcObj[l,"mature_miRNA"], desc=paste("Target=",geneame))
				points(y=miRtcObj[l,"mirExpr"],x=miRtcObj[l,"geneExpr"],pch=18,cex=2,col="black")
				}
			dev.off()
			}
	}
