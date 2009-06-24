\name{gene.similarity}
\alias{gene.similarity}
\title{Calculate adjacency matrix for gene-gene interaction}
\description{
Calculate an adjacency matrix for gene-gene interaction (using correlation/mutual information metric). For gene expression data with M genes and N experiments, the adjacency matrix is in size of MxM. It is optional to get a trimmed adjacency matrix according to the net.trim argument mrnet(),clrnet() or aracnenet() (from the minet package).
}
\usage{
gene.similarity(EXP, measure,net.trim,na.replica = 50)
}
\arguments{
  \item{EXP}{ Gene expression data in form of a matrix. Row stands for genes and column for experiments.}
  \item{measure}{ Metric used to calculate similarity between genes: "corr" for correlation, "MI" for mutual information.}
  \item{net.trim}{Method used to trim the adjacency matrix: "mrnet", "clr", "aracne" and "none" for no trim. 
  Option "mrnet" infers a network using the maximum relevance/minimum redundancy feature selection method; 
  option "clr" use the CLR algorithm; option "aracne" applies the data processing inequality to all triplets of nodes in order to remove the least significant edge in each triplet, options from minet package.
  As the above methods are only for mutual information, option "none" should be chosen if correlation is set to be the metric.}
  \item{na.replica}{ Impute method. Times of replicate for filling NANs; default value is 50. The (smooth) bootstrapping approach is used to give an estimation to missing value in the data.}
}
\value{
An adjacency matrix in size of MxM with rows and columns both standing for genes. Element in row i and column j indicates the similarity between gene i and gene j.
}
\author{ Yin Jin, Hesen Peng, Lei Wang and Christine Nardini. }
\examples{
data(copasi)
mat=as.matrix(copasi)[1:5,] 
gene.similarity(mat,measure="corr",net.trim="none")
}
\keyword{arith}


