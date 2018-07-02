#' Log Mean Devisia Index method for energy decomposition analysis
#'
#' @usage LMDI(C0, CT, X0, XT)
#' \method{print}{LMDI}(x, ...)
#' 
#' @aliases LMDI print.LMDI
#' 
#' @param C0 A numeric vector or a data.frame of carbon emission or energy consumption 
#' in the initial year.
#' @param CT A numeric vector or a data.frame of carbon emission or energy consumption 
#' in the year T.
#' @param X0 A numeric vector or a data.frame of an impact factor in the initial year.
#' @param XT A numeric vector or a data.frame of an impact factor in the year T.
#' @param x A list of LMDI result.
#' @param ... Ignore
#'
#' @author Yongze Song \email{yongze.song@postgrad.curtin.edu.au}
#' and Peng Wu \email{peng.wu@curtin.edu.au}.
#' 
#' @references Ang, B. W. (2005). The LMDI approach to decomposition 
#' analysis: a practical guide. Energy policy, 33(7), 867-871.
#' 
#' @seealso \code{\link{EDA}}
#' 
#' @examples
#' library(EDA)
#' data(carbon)
#' data(factordata)
#' ## set parameters
#' cdata <- carbon[,-c(1,2)]
#' C0 <- cdata[1,]
#' CT <- cdata[2,]
#' X0 <- factordata[[2]][1,]
#' XT <- factordata[[2]][2,]
#' ## run LMDI model
#' ed1 <- LMDI(C0, CT, X0, XT)
#' ed1
#'
#' @export

LMDI <- function(C0, CT, X0, XT){
  Dx <- (CT - C0)/(log(CT) - log(C0)) * log(XT/X0)
  Dx.Fuel <- colSums(Dx, na.rm = TRUE)
  Dx.sum <- sum(Dx, na.rm = TRUE)
  result <- list("Dx"=Dx, "Dx.Fuel"=Dx.Fuel, "Dx.sum"=Dx.sum)
  class(result) <- "LMDI"
  return(result)
}

print.LMDI <- function(x, ...){
  cat("Total change: ", x$Dx.sum)
  cat("\nChange by fuel types:\n")
  print(x$Dx.Fuel)
}
