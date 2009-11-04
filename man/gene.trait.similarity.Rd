\name{gene.trait.similarity}
\alias{gene.trait.similarity}
 \title{Calculate similarity for gene-trait interaction }
\description{
 To calculate similarity for gene-trait interaction (using correlation/mutual information metric). 
}
\usage{
gene.trait.similarity(EXP, trait, measure, na.replica = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EXP}{ Gene expression data in form of a matrix. Row stands for genes and column for experiments.}
  \item{trait}{Trait data in form of matrix. Row stands for traits and column for experiments. }
  \item{measure}{Metric used to calculate similarity: "corr" for correlation, "MI" for mutual information.}
  \item{na.replica}{Times of replicates for filling NANs in impute method; default value is 50. The (smooth) bootstrapping approach is used to give an estimation to missing value in the data.}

}
\value{
A matrix, row stands for gene and column for trait. Element in row i and column j stands for the association between the gene i and trait j.
}
 \author{ Yin Jin, Hesen Peng, Lei Wang, Raffaele Fronza, Yuanhua Liu and Christine Nardini  }
\examples{
data(tumors.mRNA)
data(tumors.miRNA)
exp<-tumors.mRNA
trait<-tumors.miRNA
gene.trait.similarity(EXP= exp[1:10, ],trait= trait[1:5, ],measure="MI")
}
\keyword{arith}