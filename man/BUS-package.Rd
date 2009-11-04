\name{BUS-package}
\alias{BUS-package}
\docType{package}
\title{For network reconstruction.}
\description{
This package can be used to compute associations among genes (gene-networks) or between genes and some external traits (i.e. clinical). [Function: BUS]

Both associations can be computed via correlation or mutual information (MI). [Functions: gene.similarity (gene-gene associations) and gene.trait.similarity (gene-trait associations)]

Statistical significance of the association is computed for single and multiple hypotheses testing, using random permutations method [Functions: gene.pvalue,  gene.trait.pvalue]

The package can handle data with missing values using bootstrapping methods to fill NAs. [Arguments: na.replica]
}
\details{
\tabular{ll}{
Package: \tab BUS \cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2009-10-31\cr
License: \tab GPL-3\cr
}
}
\author{
Yin Jin, Hesen Peng, Lei Wang, Raffaele Fronza, Yuanhua Liu and Christine Nardini

Maintainer: Yuanhua Liu<liuyuanhua@picb.ac.cn.com>
}
\keyword{ multivariate }
