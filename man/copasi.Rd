\name{copasi}
\alias{copasi}
\title{copasi data}
\description{
This dataset is taken from Copasi2 (Complex Pathway Simulator), a software for simulation and analysis of biochemical networks. The system generates random 
artificial gene networks according to well-defined topological and kinetic properties. These are used to run in silico experiments simulating real laboratory micro-array experiments. Noise with controlled properties is added to the simulation results several times emulating measurement replicates, before expression ratios
are calculated. This series consists of 150 artificial gene networks. Each network consists of 100 genes with a total of 200 gene interactions (on average each gene has 2 modulators). 
}
\format{
A data frame is size of 100x200, the 100 rows represent 100 genes and 200 columns for 200 experiments.
}
\references{
See http://www.comp-sys-bio.org/AGN/data.html for detailed information.
}
\keyword{datasets}