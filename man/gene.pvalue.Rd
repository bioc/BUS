\name{gene.pvalue}
\alias{gene.pvalue}
\title{Calculates p-value for gene-gene interaction}
\description{
Calculates p-value for the null hypothesis that there is no gene-gene interaction. For gene data with M genes, a p-value matrix under MxM single null hypotheses (each two genes have no interaction) is computed; besides, matrices with correct the p-values are output: corrected permutation method using a distribution of MxMxP (P number of permutations) null hypotheses tests (multi.perm.p.value) and using a meta-analysis based method (MMcorrected.p.value). p-values are calculated based for the adjacency matrix for gene-gene interaction computed by function "gene.similarity".
}
\usage{
gene.pvalue(EXP,measure,net.trim,n.replica=400)
}
\arguments{
  \item{EXP}{Gene expression data in form of a matrix. Row stands for genes and column for experiments.}
  \item{measure}{Metric used to calculate similarity between genes: "corr" for correlation, "MI" for mutual information.}
  \item{net.trim}{Method used to trim the adjacency matrix: "mrnet", "clr", "aracne" and "none" for no trim. 
  Option "mrnet" infers a network using the maximum relevance/minimum redundancy feature selection method; 
  option "clr" use the CLR algorithm; option "aracne" applies the data processing inequality to all triplets of nodes in order to remove the least significant edge in each triplet, these options come from the package minet.
  As the above methods are only for mutual information, option "none" should be chosen if correlation is set to be the metric.}
  \item{n.replica}{Number of permutations; default value is 400.}
}
\value{
\item{single.perm.p.value}{A matrix of single p-values obtained with permutation method + beta distribution for extreme values (for MI) or obtained with the exact distribution computed directly by cor.test (for correlation)}
\item{multi.perm.p.value}{A matrix of corrected p-values obtained with permutation method}
\item{MMcorrected.p.value}{A matrix of multiple hypothesis corrected p-values obtained with MM-correction}
}
\details{
Normally, in a permutation method, we use the empirical distribution of some statistics to estimate the p-value. To get a simple p-value for no interaction between gene i and j, empirical distribution of a vector with length of P (number of replicates) is used; to correct for multiple hypothesis with permutations, an empirical distribution of a vector with length of PxM (M being the number of hypotheses tested) is used. 
}
\author{Yin Jin, Hesen Peng, Lei Wang and Christine Nardini}
\seealso{ \code{\link{gene.similarity}}}
\examples{
data(copasi)
mat=as.matrix(copasi)[1:5,]
gene.pvalue(mat,measure="MI",net.trim="mrnet")
}
\keyword{htest}
