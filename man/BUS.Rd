\name{BUS}
\alias{BUS}
\title{A wrapper function for matrices of p-value and predicted network}
\description{
A wrapper function to calculate the computation of two types of similarities
(correlation and mutual information) with two different
 goals: (i) identification of the statistically significant similarities among the activity of molecules sampled
 across different experiments (option Unsupervised, U), (ii) identification
 of the statistically significant similarities between such molecules and other types of information (clinical etc., option supervised, S) .
}
\usage{
BUS(EXP,trait=NULL,measure,method.permut=2,method.correct,n.replica=400,net.trim=NULL,thresh=NULL,nflag)
}
\arguments{
  \item{EXP}{Gene expression data in form of a matrix. Row stands for genes and column for experiments.}
  \item{trait}{Trait data in form of a matrix. The row stands for traits and column for experiments. }
  \item{measure}{Metric used to calculate similarity: "corr" for correlation, "MI" for mutual information.}
  \item{method.permut}{A flag to indicate which method is used to correct permutation p-values, default as 2. See gene.pvalue for details.}
  \item{method.correct}{Option for method to calculate p-value corrected for multiple hypothesis: "MMcorrection" for MM-correction, "permutation" for corrected permutation, "both" for both MM-correction and permutation methods.}
  \item{n.replica}{Number of permutations; default value is 400.}
  \item{net.trim}{Method used to trim the network: "mrnet", "clr", "aracne" and "none" for no trim. Option "mrnet" infers a network using the maximum relevance/minimum redundancy feature selection method; option "clr" use the CLR algorithm; option "aracne" applies the data processing inequality to all triplets of nodes in order to remove the least significant edge in each triplet. These options come from the package minet. As these methods are only for mutual information, option "none" should be chosen if correlation is set to be the .}
  \item{thresh}{Threshold for significance of the corrected p-value, in option Unsupervised is used to trim the adjacency matrix (contains the results of the gene-gene association based on the chosen metric) and obtain a predicted gene interaction network. In the Supervised option, the network is not output, hence default is NULL.}
  \item{nflag}{A flag to indicate a gene-gene interaction case (Unsupervised) or a gene-trait interaction case (Supervised); 1 for Unsupervised and 2 for Supervised.}
}

\value{
\item{single.perm.p.value}{A matrix of single p-values obtained with permutation method}
\item{multi.perm.p.value}{A matrix of corrected p-values obtained with permutation method}
\item{MMcorrected.p.value}{A matrix of multiple hypothesis corrected p-values obtained with MM-correction}
\item{net.pred.MMcorrected}{Predicted network obtained trimming non significant values, as they are computed with the MM-correction.}
\item{net.pred.permut}{Predicted network obtained trimming non significant values, as they are computed with the permutation based multiple-hypothesis correction.}
}
\author{Yin Jin, Hesen Peng, Lei Wang and Christine Nardini }
\seealso{\code{\link{gene.pvalue}},\code{\link{gene.trait.pvalue}},\code{\link{pred.network}}
}
\examples{
data(copasi)
mat=as.matrix(copasi)[1:5,] 
BUS(mat,measure="corr",method.correct="both",net.trim="none",thresh=0.05,nflag=1)
}
\keyword{multivariate}
\keyword{htest}
