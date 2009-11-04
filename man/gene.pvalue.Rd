\name{gene.pvalue}
\alias{gene.pvalue}
\title{ Calculates p-value for gene-gene interaction }
\description{
 To calculate p-value for the null hypothesis that there is no gene-gene interaction. For gene expression data with M genes, a p-value matrix under MxM single null hypotheses (each two genes have no interaction) is computed; besides, matrices with correct p-values are output: corrected permutation method using a distribution of MxMxP (P number of permutations) null hypotheses tests (multi.perm.p.value). p-values are calculated based on the adjacency matrix for gene-gene interaction computed by function \emph{gene.similarity}.
 }
\usage{
gene.pvalue(EXP, measure, net.trim, n.replica = 400)
}
 
\arguments{
  \item{EXP}{ Gene expression data in form of a matrix. Row stands for genes and column for experiments.}
  \item{measure}{ Metric used to calculate similarity between genes: "corr" for correlation, "MI" for mutual information.}
  \item{net.trim}{Method used to trim the network: "mrnet", "clr", "aracne" and "none" . "mrnet" infers a network using the maximum relevance/minimum redundancy feature selection method; "clr" use the CLR algorithm; "aracne" applies the data processing inequality to all triplets of nodes in order to remove the least significant edge in each triplet. These options come from the package \emph{minet}, and they are used only for mutual information. "none" indicates no trim operation. It should be chosen when correlation is considered.}
  \item{n.replica}{Number of permutations used for the correction of multiple hypothesis testing; default value is 400.}
 }
 
\value{
\item{single.perm.p.value}{A matrix of single p-values obtained with permutation method + beta distribution for extreme values (for MI) or obtained with the exact distribution computed directly by \emph{cor.test} (for correlation)}
\item{multi.perm.p.value}{A matrix of corrected p-values obtained with permutation method}
 }
 
\details{
Normally, in a permutation method, we use the empirical distribution of some statistics to estimate the p-value. To get a simple p-value for no interaction between gene i and j, empirical distribution of a vector with length of P (number of replicates) is used; to correct for multiple hypothesis with permutations, an empirical distribution of a vector with length of PxM (M being the number of hypotheses tested) is used. 
}
\author{ Yin Jin, Hesen Peng, Lei Wang, Raffaele Fronza, Yuanhua Liu and Christine Nardini }

\seealso{ \code{\link{gene.similarity}}}
\examples{
data(copasi)
mat=as.matrix(copasi)[1:10,]
rownames(mat)<-paste("G",1:nrow(mat), sep="")
gene.pvalue(mat,measure="MI",net.trim="mrnet")
}
\keyword{htest}