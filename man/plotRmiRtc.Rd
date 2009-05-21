\name{plotRmiRtc}
\alias{plotRmiRtc}
\title{Plot object from read.mir or a selected gene and respective miRNAs from a miRtcList object}

\description{
  Ploting function for object coming from read.mir or a selected gene and respective miRNAs from a miRtcList object
} 

\usage{
  plotRmiRtc(miRtcObj,gene_id=NULL,timeunit="Time",legend.x=NULL,legend.y=NULL,svgTips=FALSE,svgname=NULL,height=10,width=15)
}

\arguments{
  \item{miRtcObj}{A data.frame resulting from \code{read.mir} or \code{RmiR} functions or a \code{miRtcList-class} object.}

  \item{gene_id}{Selected \code{gene_id} contained in a \code{miRtcList-class} object.}

  \item{timeunit}{Name for the abscissae axes, normally a time unit like "Hours", "PD" etc.}

  \item{legend.x}{Position of the legend in the x-axes.}

  \item{legend.y}{Position of the legend in the y-axes.}

  \item{svgTips}{TRUE if you want to use the \code{RSVGTipsDevice}, default is FALSE.}

  \item{svgname}{Name for the SVG image output.}

  \item{height}{Height of the graphs.}

  \item{width}{Width of the graphs.}
}

\details{
  The function plots the trends of a gene target with the specified \code{gene_id} and respective miRNA contained in a \code{miRtcList-class} object.
  
  If the \code{miRtcObj} argument is a dataframe coming from \code{read.mir} function, the resulting plot will be a point graph in SVG format. 
  Each couple miRNA/Target is a point, the x value is the gene target expression value and the y value is the microRNA expression value.
  To decrease the size of the graph is possible to select just the desired miRNAs or gene targets in the data.frame 

}

\seealso{
  \code{readRmiRtc},\code{miRtcList}
}

\examples{
	
        data(RmiR)
        res1 <- read.mir(genes=gene1, mirna=mir1, annotation="hgug4112a.db")
        res2 <- read.mir(genes=gene2, mirna=mir2, annotation="hgug4112a.db")
        res3 <- read.mir(genes=gene3, mirna=mir3, annotation="hgug4112a.db")

        res_tc <- RmiRtc(timeline=c("res1", "res2", "res3"),
			 timevalue=c(12, 48, 72))
        res <- readRmiRtc(res_tc, correlation=-0.9, exprLev=1,
			  annotation="hgug4112a.db")

	## List of genes with anti-correlated miRNAs:
        
	res$reps

	## Plot of the first gene of the list:
	plotRmiRtc (res, gene_id=351, timeunit="Hours")

	## Setting the position of the legend:
	
	plotRmiRtc (res,gene_id=351, legend.x=50, legend.y=0, timeunit="Hours")

	## Plot with RSVGTipsDevice:
	
	plotRmiRtc (res,gene_id=351, legend.x=50, legend.y=0, timeunit="Hours",
		    svgTips=TRUE)

	## Plot of a read.mir results:
	
	plotRmiRtc (res1, svgname="gene1.svg", svgTips=TRUE)
   
}

