test_rebet <- function() {

  data(data, package="REBET")

  res         <- rebet(response, genotypes, subRegions)
  meta.beta   <- res$Meta$beta
  side1.beta  <- res$Subset.1sided$beta
  side1.se    <- res$Subset.1sided$sd
  side2.beta1 <- res$Subset.2sided$beta.1
  side2.se1   <- res$Subset.2sided$sd.1
  side2.p2    <- res$Subset.2sided$pval.2
  
  checkEqualsNumeric(meta.beta,   4.398939, tolerance=1.0e-4)
  checkEqualsNumeric(side1.beta,  7.037888, tolerance=1.0e-4)
  checkEqualsNumeric(side1.se,    1.339332, tolerance=1.0e-4)
  checkEqualsNumeric(side2.beta1, 7.037888, tolerance=1.0e-4)
  checkEqualsNumeric(side2.se1,   1.348216, tolerance=1.0e-4)
  checkEqualsNumeric(side2.p2,    0.644174, tolerance=1.0e-4)
  
}
