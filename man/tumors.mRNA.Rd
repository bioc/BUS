\name{tumors.mRNA}
\alias{tumors.mRNA}
\docType{data}
\title{ Gene expression data from Human brain tumors}
\description{
 Gene expression data obtained by microarray from human brain tumors. 12 brain tumors at different levels are analyzed for both mRNA and miRNA levels to study the correlation of any mRNA-miRNA pair in the \emph{reference}  .
}
\usage{data(tumors.mRNA)}
\format{
  \code{tumors.mRNA} is a matrix with mRNA probe IDs as rows and tumor type as columns.
}
\references{
Liu T, Papagiannakopoulos T, Puskar K, Qi S,  Santiago F, Clay W, Lao K, Lee Y, Nelson SF, Kornblum HI, Doyle F, Petzold L, Shraiman B, Kosik KS. Detection of a microRNA signal in an in vivo expression set of mRNAs. Plos One. 2007; 2(8):e804.
}
\examples{
data(tumors.mRNA)
tumors.mRNA[1:10,]
}
\keyword{datasets}
