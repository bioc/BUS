\name{gene.similarity}
\alias{gene.similarity}
 \title{ Calculate adjacency matrix for gene-gene interaction }
\description{
To calculate an adjacency matrix for gene-gene interaction (using correlation/mutual information metric). For gene expression data with M genes and N experiments, the adjacency matrix is in size of MxM. It is optional to get a trimmed adjacency matrix according to the argument \emph{net.trim}, i.e. \emph{mrnet}, \emph{clr} and\emph{aracne} (from the package \emph{minet}).
}
\usage{
gene.similarity(EXP, measure, net.trim, na.replica = 50)
}

\arguments{
  \item{EXP}{ Gene expression data in form of a matrix. Row stands for genes and column for experiments.}
  \item{measure}{ Metric used to calculate similarity between genes: "corr" for correlation, "MI" for mutual information.}
  \item{net.trim}{Method used to trim the adjacency matrix: "mrnet", "clr", "aracne" and "none". "mrnet" infers a network using the maximum relevance/minimum redundancy feature selection method; "clr" use the CLR algorithm; "aracne" applies the data processing inequality to all triplets of nodes in order to remove the least significant edge in each triplet. These options come from the package \emph{minet}, and they are used only for mutual information. "none" indicates no trim operation. It should be chosen when correlation is considered.}
  \item{na.replica}{ Times of replication for filling NANs in the impute method; default value is 50. The (smooth) bootstrapping approach is used to give an estimation to missing value in the data.}

}

\value{
  An adjacency matrix in size of MxM with rows and columns both standing for genes. Element in row i and column j indicates the similarity between gene i and gene j.
  }
\author{ Yin Jin, Hesen Peng, Lei Wang, Raffaele Fronza, Yuanhua Liu and Christine Nardini}
\examples{
data(copasi)
mat=as.matrix(copasi)[1:10,] 
rownames(mat)<-paste("G",1:nrow(mat), sep="")
res<-gene.similarity(mat,measure="corr",net.trim="none")
}
\keyword{arith}  }
