
\name{REBET}
\alias{REBET}
\docType{package}
\title{
 The REBET package
}
\description{
 An R package for the subREgion-based BurdEn Test (REBET).
}
\details{
In rare-variant association studies, aggregating rare and/or low frequency variants,
may increase statistical power for detection of the underlying susceptibility gene or region.
However, it is unclear which variants, or class of them, in a gene contribute most to the association.
This subregion-based burden test (REBET) simultaneously selects susceptibility genes and
identifies important underlying sub-regions.
The sub-regions are predefined based on shared common biologic characteristics,
such as the protein domain or possible functional impact.
Based on a subset-based approach considering local correlations between combinations
of test statistics of sub-regions, REBET is able to properly control the type I error rate
while adjusting for multiple comparisons in a computationally efficient manner. 
See the reference for the details of this test.
The main function in this package is \code{\link{rebet}}, which performs the REBET test.
}
\references{
Zhu, B., Mirabello, L., Chatterjee, N. (2018) A Subregion-based Burden Test for 
Simultaneous Identification of Susceptibility Loci and Sub-regions within
Genetic Epidemiology. 
In press. https://doi.org/10.1002/gepi.22134
}
\author{Bin Zhu <bin.zhu@nih.gov>, Lisa Mirabello and Nilanjan Chatterjee}
\keyword{package}