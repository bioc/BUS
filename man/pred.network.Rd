\name{pred.network}
\alias{pred.network}
\title{Predict the network}
\description{
Predict the matrix of gene network, based on the similarity matrix and filtered according to a corrected p-value matrix.
}
\usage{pred.network(pM,similarity,thresh)
}
\arguments{
  \item{pM}{A corrected p-value matrix, a MxM matrix for significance of similarity among M genes.}
  \item{similarity}{A MxM matrix for similarity between genes.}
  \item{thresh}{Threshold for significance of the p-value.}
}
\value{
  A MxM matrix of the predicted network, cell ij infers a link between gene i and j and set 0 when the p-value is not significant (no link).
  }
\author{Yin Jin, Hesen Peng, Lei Wang and Christine Nardini }
\examples{
data(copasi)
mat=as.matrix(copasi)[1:5,]
similarity=gene.similarity(mat,measure="MI",net.trim="mrnet")
pM=gene.pvalue(mat,measure="MI",net.trim="mrnet")$single.perm.p.value
pred.network(pM,similarity,thresh=0.05)
}
\keyword{arith}
