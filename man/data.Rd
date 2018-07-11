
\name{data}
\alias{data}
\alias{Y}
\alias{G}
\alias{E}


\title{
  Data for the example
}
\description{
 Data for the example.
}
\details{ The data contains a binary phenotype vector \code{Y},
 a genotype matrix \code{G} consisting of 20 rare-variant SNPs, 
 and the sub-region annotation vector \code{E}  
  for the \code{\link{rebet}} example.
}
\examples{

 data(data, package="REBET")

 # Display some of the data
 table(Y)
 dim(G)
 E
}
\keyword{data}
\seealso{\code{\link{rebet}}}



