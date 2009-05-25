\name{RmiR}
\alias{RmiR}
\alias{read.mir}
\title{Coupling miRNA and Gene expression results}

\description{
  Coupling miRNA and Gene expression results for a selected target database.  
} 

\usage{
  RmiR(mirna=NULL,genes=NULL,annotation=NULL,dbname="targetscan",org="hsa",id="probes",id.out="symbol",verbose=FALSE)
  read.mir(mirna=NULL,genes=NULL,annotation=NULL,id="probes",dbname=c("targetscan","pictar"),org="hsa",at.least=2,id.out="symbol",verbose=FALSE)
}

\arguments{
  \item{mirna}{A data.frame with two columns, the first with the microRNA names, the second with the expression values.}

  \item{genes}{A data.frame with two columns, the first with gene ID (probes, symbols, ensembl, entrez...), the second with the expression values.}

  \item{annotation}{The annotation package to annotate the genes file with entrez gene ID, eg: Agilent 44k annotation="hgug4112a.db" or annotation="org.Hs.eg.db" for human not using microarrays probes.}

  \item{dbname}{A selected database of miRNA target. See 
\code{RmiR.hsa_dbconn}, default is "targetscan". If using \code{read.mir} it can be a vector of databases, 
default are "targetscan" and "pictar".}

  \item{id}{The type of annotation of the genes input file. An accepted value is one of: "genes" for entrez gene id, "probes" for microarray probes id, "ensembl" for ensembl gene id, "unigene" for unigene gene id and "alias" for official gene symbols and aliases.}

  \item{id.out}{The annotation of the genes in the output. The default it is "symbol", to have the HGNC symbols, it can be also "probes" if the input \code{id} is "probes" or "gene" to leave just the entrez gene annotation.}

  \item{at.least}{Minimum number of databases that should yeld the result, when the search is performed in multiple databases with \code{read.mir}. If it is 1 it is basically an union between databases. Default is 2.}
 
  \item{org}{Define the targets database package of the desired organism. Default is "hsa"}

  \item{verbose}{If it is desired or not to have some verbose output while analysing the data. Default is FALSE}
}

\details{
  \code{RmiR} couples the gene expression and microRNA expression. It uses the AnnotationDbi package to annotate the gene expression file. 
  We intend to put already filtered and significant values in the input file, so in case of duplicate probes or different sequences identifying the same gene or more than one values for a miRNA, the function will take just the mean of the different results and give the corresponding coefficent of variation.
  Each input file must have two columns. The first one for annotation, the second for expression value. The name of the columns does not matter.

  \code{read.mir} uses \code{RmiR} but performs the search in one or more databases and returns only the object present in \code{at.least} databases. If \code{at.least} is equal to 1 we basically do an union between the results from the databases of choice, if we specify just a database in \code{dbname} it is exactly the same of the \code{RmiR} function. 

}
\value{
  \item{mature_miRNA}{The resulting miRNAs present in the input file with at least one target in the selected database.}
  \item{gene_id}{The resulting entrez gene ids present in the input file that are also targets in the selected database.}
  \item{mirExpr}{microRNA expression value}
  \item{geneExpr}{Gene expression Value}
  \item{mirCV}{miRNA expression coefficent of variation in case of duplication otherwise is NA}
  \item{geneCV}{Gene expression coefficent of variation in case of duplication otherwise is NA}
  \item{symbol}{If the \code{id.out} is "symbol".}
  \item{probe_id}{If the \code{id.out} is "probes".}
}

\seealso{
  \code{RmiR.hsa_dbconn},
}

\examples{
  ## Merge gene expression and mirna expression for agilent IDs

  genes <- data.frame(genes=c("A_23_P171258", "A_23_P150053", "A_23_P150053", 
			      "A_23_P150053", "A_23_P202435", "A_24_P90097",
			      "A_23_P127948"))
  genes$expr <- c(1.21, -1.50, -1.34, -1.45, -2.41, -2.32, -3.03)

  mirna <- data.frame(mirna=c("hsa-miR-148b", "hsa-miR-27b", "hsa-miR-25",
			      "hsa-miR-181a", "hsa-miR-27a", "hsa-miR-7",
			      "hsa-miR-32", "hsa-miR-32", "hsa-miR-7"))
  mirna$expr <- c(1.23, 3.52, -2.42, 5.2, 2.2, -1.42, -1.23, -1.20, -1.37)

  RmiR(genes=genes, mirna=mirna, annotation="hgug4112a.db", id="probes")

  ## Search in pictar

  RmiR(genes=genes, mirna=mirna, annotation="hgug4112a.db", id="probes", 
       dbname="pictar")
  ## or

  read.mir(genes=genes, mirna=mirna, annotation="hgug4112a.db", id="probes",
           dbname="pictar", at.least=1)

  ## Search in miranda, pictar and targetscan, present in each database:

  read.mir(genes=genes, mirna=mirna, annotation="hgug4112a.db", id="probes",
	   dbname=c("miranda", "pictar", "targetscan"), at.least=3)

  ## Search in miranda, pictar and targetscan, present in at least 2 database:

  read.mir(genes=genes, mirna=mirna, annotation="hgug4112a.db", id="probes",
	   dbname=c("miranda", "pictar", "targetscan"), at.least=2)
	
}


