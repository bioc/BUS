\name{BUS}
\alias{BUS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A wrapper function for matrices of p-value and predicted network }
\description{
A wrapper function to calculate the computation of two types of similarities
(correlation and mutual information) with two different
 goals: (i) identification of the statistically significant similarities among the activity of molecules sampled
 across different experiments (option Unsupervised, U), (ii) identification
 of the statistically significant similarities between such molecules and other types of information (clinical etc., option supervised, S) .
}
\usage{
BUS(EXP, trait = NULL, measure, method.permut = 2, n.replica = 400, net.trim = NULL, thresh = NULL, nflag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EXP}{Gene expression data in form of a matrix. Row stands for genes and column for experiments.}
  \item{trait}{Trait data in form of a matrix. The row stands for traits and column for experiments. }
  \item{measure}{Metric used to calculate similarity: "corr" for correlation, "MI" for mutual information.}
  \item{method.permut}{A flag to indicate which method is used to correct permutation p-values, default as 2. See gene.pvalue for details.}
  \item{n.replica}{Number of permutations used for the correction of multiple hypothesis testing; default value is 400.}
  \item{net.trim}{Method used to trim the network: "mrnet", "clr", "aracne" and "none" . "mrnet" infers a network using the maximum relevance/minimum redundancy feature selection method; "clr" use the CLR algorithm; "aracne" applies the data processing inequality to all triplets of nodes in order to remove the least significant edge in each triplet. These options come from the package \emph{minet}, and they are used only for mutual information. "none" indicates no trim operation. It should be chosen when correlation is considered.}
  \item{thresh}{Threshold for significance of the corrected p-value. It is used, in the Unsupervised case, to trim the adjacency matrix (contains the results of the gene-gene association based on the chosen metric) and obtain a predicted gene interaction network. In the Supervised case, since no network is predicted, it is set as NULL.}
  \item{nflag}{A flag to indicate a gene-gene interaction case (Unsupervised) or a gene-trait interaction case (Supervised); 1 for Unsupervised and 2 for Supervised.}
}

\value{
\item{similarity}{A matrix of similarity, which could be correlation or mutual information}
\item{single.perm.p.value}{A matrix of single p-values}
\item{multi.perm.p.value}{A matrix of corrected p-values}
\item{net.pred.permut}{Predicted network obtained trimming non-significant values}
}

\author{Yin Jin, Hesen Peng, Lei Wang, Raffaele Fronza, Yuanhua Liu and Christine Nardini}
\seealso{\code{\link{gene.pvalue}},\code{\link{gene.trait.pvalue}},\code{\link{pred.network}}
}
\examples{

data(copasi)
mat<-as.matrix(copasi[1:10,])
rownames(mat)<-paste("G",1:nrow(mat), sep="")
BUS(EXP=mat,measure="corr",net.trim="none",thresh=0.05,nflag=1)

}
 
\keyword{multivariate}
\keyword{htest}
