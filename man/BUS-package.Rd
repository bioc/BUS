\name{BUS-package}
\alias{BUS-package}
\docType{package}
\title{For network reconstruction.}
\description{
This package can be used to compute associations among genes (gene-networks) or between genes and some external traits (i.e. clinical). [Function: BUS]

Both associations can be computed via correlation or mutual information (MI). [Functions: gene.similarity (for adjacency matrix for genes) and gene.trait.similarity (for association between genes and traits)]

Statistical significance of the association is computed for single and multiple hypotheses testing, using (i) random permutations method and (ii) MM-correction. [Functions: gene.pvalue,  gene.trait.pvalue and MM.correction]

The package can handle data with missing values using bootstrapping methods to fill NAs. [Arguments: na.replica]
}
\details{
\tabular{ll}{
Package: \tab BUS \cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2008-12-09\cr
License: \tab GPL-3\cr
}
}
\author{
Yin Jin, Hesen Peng, Lei Wang and Christine Nardini

Maintainer: Yin Jin<golden.jinyin@gmail.com>
}
\references{
C. Nardini, H. Peng, L. Wang, L. Benini, M.D. Kuo, MM-Correction: Meta-analysis-Based Multiple Hypotheses Correction in Omic Studies, Springer CCIS, Vol 25, pp 242-255, 2008.
}
\keyword{ multivariate }


