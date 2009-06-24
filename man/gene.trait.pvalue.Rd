\name{gene.trait.pvalue}
\alias{gene.trait.pvalue}
\title{ Calculate p-value for gene-trait interaction}
\description{
Calculates p-value for null hypothesis that there is no interaction between gene and trait. There are MxT interactions between M genes and T traits. Results are given with 4 possibilities 1 for single p-value, and 3 for different types of correction.
p-values are calculated based for the adjacency matrix for gene-gene interaction computed by function "gene.trait.similarity". 
}
\usage{
gene.trait.pvalue(EXP,trait,measure,method.permut=2,n.replica=400)
}
\arguments{
  \item{EXP}{Gene expression data in form of a matrix. Row stands for genes and column for experiments.}
  \item{trait}{Trait data in form of matrix. The row stands for traits and column for experiments. }
  \item{measure}{Metric used to calculate similarity: "corr" for correlation, "MI" for mutual information.}
  \item{method.permut}{A flag to indicate which dimension is considered when a correcting for multiple hypotheses. 1 for traits correction, 2 for genes and 3 for both genes and traits. The default value is 2.}
  \item{n.replica}{Number of permutations; default value is 400.}
}
\value{
\item{single.perm.p.value}{A matrix of single p-values obtained with permutation method + beta distribution for extreme values (for MI) or obtained with the exact distribution computed directly by cor.test (for correlation)}
\item{multi.perm.p.value}{A matrix of corrected p-values obtained with permutation method}
\item{MMcorrected.p.value}{A matrix of multiple hypothesis corrected p-values obtained with MM-correction}
}
\details{
According to a permutation method, we use the empirical distribution of some statistics to estimate the p-value. For single p-value the empirical distribution is a vector of P (number of random replicates for each test) test values. It is then possible to correct p-value in different ways: method.permut = 1, it is the empirical distribution of a vector with length of TxP, corrects for the multiple traits tested; method.permut = 2, it is the empirical distribution of a vector with length of MxP, corrects for the multiple genes tested; method.permut = 3, it is empirical distribution of a vector with length of MxTxP, corrects for the multiple traits and genes tested. 
}
\author{Yin Jin, Hesen Peng, Lei Wang and Christine Nardini }
\seealso{ \code{\link{gene.trait.similarity}}}
\examples{
data(copasi)
mat=as.matrix(copasi)[1:5,]
gt<-matrix(rnorm(200),2,100)
gene.trait.pvalue(EXP=mat,trait=gt,measure="corr")
}
\keyword{htest}
