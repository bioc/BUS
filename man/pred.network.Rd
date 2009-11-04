\name{pred.network}
\alias{pred.network}
\title{Predict the network}
\description{
To predict the matrix of gene network, based on the similarity matrix and filtered according to a corrected p-value matrix.
}
\usage{pred.network(pM,similarity,thresh)
}
\arguments{
  \item{pM}{A corrected p-value matrix, a MxM matrix for significance of similarity among M genes.}
  \item{similarity}{A MxM matrix for similarity between genes.}
  \item{thresh}{Threshold for significance of the p-value.}
}
\value{
  A MxM matrix of the predicted network, where cell emph{ij} infers a link between gene i and j and set 0 when the p-value is not significant (no link).
  }
\author{Yin Jin, Hesen Peng, Lei Wang, Raffaele Fronza, Yuanhua Liu and Christine Nardini}
\examples{
data(copasi)
mat<-as.matrix(copasi[1:10,])
rownames(mat)<-paste("G",1:nrow(mat), sep="")
similarity=gene.similarity(mat,measure="MI",net.trim="mrnet")
pM=gene.pvalue(mat,measure="MI",net.trim="mrnet")$single.perm.p.value
pred.network(pM,similarity,thresh=0.05)
}
\keyword{arith}