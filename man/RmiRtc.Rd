\name{RmiRtc}
\alias{RmiRtc}
\alias{readRmiRtc}
\alias{miRtcList-class}
\title{Time Course relationship between microRNA and Genes}

\description{
Given a timeline of experiments resulting from \code{RmiR} or \code{read.mir}, it calculates the correlation between the trend of miRNA and corresponding gene targets.
} 

\usage{
	RmiRtc(timeline = NULL, timevalue = NULL, method = "pearson")
	readRmiRtc(miRtcObj, correlation = -0.75, exprLev = 1, annotation= NULL, fileName = "miRNA_for_genes")
}

\arguments{
  \item{timeline}{A vector with the names of the experiments resulting from \code{RmiR} or \code{read.mir}, in chronological order.}

  \item{timevalue}{A vector of numbers with the unity of time correspondig to \code{timeline}.}

  \item{method}{Method to use to calculate the correlation between miRNA and gene expression, default is "pearson". For other see \code{cor} from \code{stats} package.}

  \item{miRtcObj}{An object resulting from \code{RmiRtc}.}

  \item{annotation}{The annotation package to retrive the corresponding symbol given the \code{gene_id} . eg: Agilent 44k annotation="hgug4112a.db" or 
annotation="org.Hs.eg.db".}
 
  \item{correlation}{The correlation level desired to filter the \code{miRtcList} object created with the \code{RmiRtc} function.}
  
  \item{exprLev}{The absolute value of gene expression as cut-off to filter the \code{miRtcList} object created with the \code{RmiRtc} function.}

  \item{fileName}{The file name to print the file with the gene targets with the number of miRNAs matching the \code{correlation} criteria. If nothing is specified, no file will be created.}
}

\details{
  \code{RmiRtc} creates an \code{miRtcList} wich includes all the information of the time course experiment: couples of miRNA and gene target, expression of gene and miRNA in the time, the correlation between the miRNA and the gene expression trends.
  
  \code{readRmiRtc} subsets the \code{miRtcList} created with \code{RmiRtc}. We can select a correlation level, if positive we select the correlated genes and miRNas, if negative the anti-correlated couples. Also we can decrease the data by setting a log ratio cut off for the gene expression, to select only the case which the a gene is op or down regulated.

}
\value{
  \item{couples}{The couples of \code{mature_miRNA} and targets in entrez gene annotation.}
  
  \item{mirExpr}{A matrix with the expression of miRNA in order by \code{timeline}.}

  \item{geneExpr}{A matrix with the expression of miRNA in order by \code{timeline}.}

  \item{mirCV}{A matrix with the coefficents of variation of the miRNAs from \code{RmiR} or \code{read.mir}.}

  \item{geneCV}{A matrix with the coefficents of variation of the genes resulting from \code{RmiR} or \code{read.mir}.}

  \item{correlation}{A vector with the correlation value between miRNAs and gene targets.}

  \item{reps}{With \code{readRmiRtc} we list all the gene targets ordered by the number of miRNAs matching the correlation criteria.}
}

\seealso{
  \code{RmiR}, \code{read.mir}, \code{plotRmiRtc}
}

\examples{
	##An example without the data
	data(RmiR)
	res1 <- read.mir(genes=gene1, mirna=mir1, annotation="hgug4112a.db")
	res2 <- read.mir(genes=gene2, mirna=mir2, annotation="hgug4112a.db")
	res3 <- read.mir(genes=gene3, mirna=mir3, annotation="hgug4112a.db")
	
	res_tc <- RmiRtc(timeline=c("res1", "res2", "res3"),
			 timevalue=c(12,48,72))
	res <- readRmiRtc(res_tc, correlation=-0.9, exprLev=1, 
			  annotation="hgug4112a.db")
	res$reps
}

