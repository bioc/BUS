\name{MM.correction}
\alias{MM.correction}
\title{MM-correction method}
\description{
Correct single p-value for Multiple Hypothesis Testing 
}
\usage{
MM.correction(data)
}
\arguments{
  \item{data}{ A square matrix of single p-values; each element is a p-value to be corrected for multi-hypothesis testing.}
}
\value{
A square matrix with the same size of data, elements are corrected p-values.
}
\references{
C. Nardini, H. Peng, L. Wang, L. Benini, M.D. Kuo, MM-Correction: Meta-analysis-Based Multiple Hypotheses Correction in Omic Studies, Springer CCIS, Vol 25, pp 242-255, 2008.
}
\author{Yin Jin, Hesen Peng, Lei Wang and Christine Nardini}
\examples{
pmat=matrix(c(0.0000,0.6675,0.5800,0.7375,0.3800,0.6675,0.0000,0.1975,0.3075,0.1400,0.5800,0.1975,0.0000,0.1625,0.1125,0.7375,0.3075,0.1625,0.0000,0.0000,0.3800,0.1400,0.1125,0.0000,0.0000),nr=5,nc=5,byrow=TRUE)
MM.correction(pmat)
}
\keyword{htest}
