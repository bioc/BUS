\name{gene.trait.similarity}
\alias{gene.trait.similarity}
\title{Calculate similarity for gene-trait interaction}
\description{
Calculate similarity for gene-trait interaction (using correlation/mutual information metric). 
}
\usage{
gene.trait.similarity(EXP, trait, measure, data.type=c("continuous","continuous"), na.replica = 50)
}
\arguments{
  \item{EXP}{ Gene expression data in form of a matrix. Row stands for genes and column for experiments.}
  \item{trait}{Trait data in form of matrix. The row stands for traits and column for experiments. }
  \item{measure}{Metric used to calculate similarity: "corr" for correlation, "MI" for mutual information.}
  \item{data.type}{A vector for the type of gene expression data and trait data. Each element in the vector has two options, "discrete" and "countinuous". The default value is c("continuous","continuous"). }
  \item{na.replica}{Impute method. Times of replicate for filling NANs; default value is 50. The (smooth) bootstrapping approach is used to give an estimation to missing value in the data.}
}
\value{
A matrix, row stands for gene and column for trait. Element in row i and column j stands for the link between the gene i and trait j.
}
\author{Yin Jin, Hesen Peng, Lei Wang and Christine Nardini }
\examples{
data(copasi)
mat=as.matrix(copasi)[1:5,]
gt<-matrix(rnorm(200),2,100)
gene.trait.similarity(EXP=mat,trait=gt,measure="MI")
}
\keyword{arith}
