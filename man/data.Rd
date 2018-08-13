
\name{data}
\alias{data}
\alias{response}
\alias{genotypes}
\alias{subRegions}


\title{
  Data for the example
}
\description{
 Data for the example.
}
\details{ The data contains a binary phenotype vector \code{response},
 a genotype matrix \code{genotypes} consisting of 20 rare-variant SNPs, 
 and the sub-region annotation vector \code{subRegions}  
  for the \code{\link{rebet}} example.
}
\examples{

 data(data, package="REBET")

 # Display some of the data
 table(response)
 dim(genotypes)
 subRegions
}
\keyword{data}
\seealso{\code{\link{rebet}}}



