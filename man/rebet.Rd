
\name{rebet}
\alias{rebet}
\title{The subREgion-based BurdEn Test (REBET)}
\description{
 A Subregion-based Burden Test for Simultaneous Identification of Susceptibility Loci and Sub-regions within                              
}
\usage{ rebet(response, genotypes, subRegions, responseType=NULL, 
        covariates=NULL, shape1=1, shape2=1, saveMem=FALSE)
}
\arguments{
  \item{response}{Numerical vector of phenotypes. A binary phenotype must be coded as 0 and 1.}
  \item{genotypes}{Matrix of genotypes with each column as a locus. }
  \item{subRegions}{Sub-region annotation vector with length equal to the number of columns of \code{genotypes}.
           In the returned object, these regions will appear as 
           \code{paste("Region_", subRegions, sep="")}. }
  \item{responseType}{NULL, "continuous" or "binary". 
        If NULL, then "continuous" or "binary" will be chosen based on \code{Y}.
        The default is NULL.} 
  \item{covariates}{NULL or matrix of covariates. The default is NULL.}
  \item{shape1}{The \code{shape1} parameter in the beta distribution. The default is 1.}
  \item{shape2}{The \code{shape2} parameter in the beta distribution. The default is 1.}
  \item{saveMem}{TRUE or FALSE to conserve memory (see details). The default is FALSE.}
}
\details{
 See the reference for details of this method.

 Missing values in any of \code{response}, \code{genotypes} or \code{covariates} will be removed
 before the analysis. Setting \code{saveMem} to TRUE will allow for the
 analysis of a much larger number of subjects, but will take more time to compute.
 When \code{saveMem} is FALSE, there needs to be enough memory available to hold
 two or three NxN matrices, where N is the number of subjects. 

 This function calls the \code{\link[ASSET]{h.traits}} function in the 
  \code{\link[ASSET]{ASSET}} package.
}
\value{
  The object returned from \code{\link[ASSET]{h.traits}} in the 
  \code{\link[ASSET]{ASSET}} package.
}
\references{
Zhu, B., Mirabello, L., Chatterjee, N. (2018) A Subregion-based Burden Test for 
Simultaneous Identification of Susceptibility Loci and Sub-regions within
Genetic Epidemiology. 
In press. https://doi.org/10.1002/gepi.22134

}
\author{Bin Zhu <bin.zhu@nih.gov>, Lisa Mirabello and Nilanjan Chatterjee}
\examples{
  data(data, package="REBET")

  res <- rebet(response, genotypes, subRegions)
  h.summary(res)

}
\keyword{ SNP, gene, rare variant } 
%\seealso{\code{\link{REBET}}}
