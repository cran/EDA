#' Log Mean Devisia Index method for energy decomposition analysis
#'
#' @param C0 A data.frame
#' @param CT A data.frame
#' @param X0 A data.frame
#' @param XT A data.frame 
#'
#' @examples
#' data(CarbonEmission)
#' carbon <- CarbonEmission$carbon
#' cdata <- carbon[,-c(1,2)]
#' xdata <- CarbonEmission$xdata
#' C0 <- cdata[1,]
#' CT <- cdata[2,]
#' X0 <- xdata[[1]][1,]
#' XT <- xdata[[2]][2,]
#' LMDI(C0, CT, X0, XT)
#'
#' @export

LMDI <- function(C0, CT, X0, XT){
  Dx <- (CT - C0)/(log(CT) - log(C0)) * log(XT/X0)
  Dx.category <- rowSums(Dx, na.rm = TRUE)
  Dx.component <- colSums(Dx, na.rm = TRUE)
  Dx.sum <- sum(Dx, na.rm = TRUE)
  return(list("Dx"=Dx, "Dx.category"=Dx.category, 
              "Dx.component"=Dx.component, "Dx.sum"=Dx.sum))
}
