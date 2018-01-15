#' Energy consumption analysis for calculating carbon emission changes
#'
#' @param cdata A data.frame
#' @param xdata A list
#' @param Year A numeric vector
#' @param Category A vector
#' @param Factor A vector
#' @param Component A vector
#' @param method A character chosen from 
#' "LMDI" or "Laspeyres" or "Paasche" or "Marshall-Edgeworth" or "Walsh"
#' 
#' @importFrom stats aggregate
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw scale_fill_hue 
#' geom_line geom_point 
#'
#' @examples
#' data(CarbonEmission)
#' carbon <- CarbonEmission$carbon
#' cdata <- carbon[,-c(1,2)]
#' xdata <- CarbonEmission$xdata
#' Year <- carbon$year
#' Category <- carbon$building
#' Component <- colnames(carbon[,-c(1:2)])
#' Factor <- c(1:length(xdata))
#' EDA(cdata, xdata, Year = Year, Category = Category, Factor = Factor, 
#' Component = Component, method = "LMDI")
#'
#' @export

EDA <- function(cdata, xdata, Year = Year, Category = Category, 
                Factor = Factor, Component = Component, method = "LMDI"){
  NYear <- levels(factor(Year))
  NCategory <- levels(factor(Category))
  Dx <-  expand.grid(Year = NYear[-1], Category = NCategory, Factor = Factor)
  Dxc <- matrix(NA, nrow(Dx), length(Component))
  Dxc <- as.data.frame(Dxc)
  names(Dxc) <- Component
  Dx <- cbind(Dx, Dxc)
  
  LF <- c(1:length(Factor))
  
  if (method == "LMDI"){
    for (u in 1:(length(NYear)-1)){
      C0 <- cdata[which(Year == NYear[u]),]
      CT <- cdata[which(Year == NYear[u+1]),]
      for (v in LF){
        X0 <- xdata[[v]][which(Year == NYear[u]),]
        XT <- xdata[[v]][which(Year == NYear[u+1]),]
        Dx[which(Dx$Year == NYear[u+1] & Dx$Factor == Factor[v]),-c(1:3)] <- LMDI(C0, CT, X0, XT)$Dx
      }
    }
  } else if (method == "Laspeyres"){
    for (u in 1:(length(NYear)-1)){
      X0 <- list(); XT <- list()
      for (v in LF){
        X0[[v]] <- xdata[[v]][which(Year == NYear[u]),]
        XT[[v]] <- xdata[[v]][which(Year == NYear[u+1]),]
      }
      
      for (w in LF){
        Vx1 <- XT[[w]] - X0[[w]]
        for (r in LF[-w]){
          Vx1 <- Vx1 * X0[[r]]
        }
        Dx[which(Dx$Year == NYear[u+1] & Dx$Factor == w),-c(1:3)] <- Vx1
      }
    }
  } else if (method == "Paasche"){
    for (u in 1:(length(NYear)-1)){
      X0 <- list(); XT <- list()
      for (v in LF){
        X0[[v]] <- xdata[[v]][which(Year == NYear[u]),]
        XT[[v]] <- xdata[[v]][which(Year == NYear[u+1]),]
      }
      
      for (w in LF){
        Vx1 <- XT[[w]] - X0[[w]]
        for (r in LF[-w]){
          Vx1 <- Vx1 * XT[[r]]
        }
        Dx[which(Dx$Year == NYear[u+1] & Dx$Factor == w),-c(1:3)] <- Vx1
      }
    }
  } else if (method == "Marshall-Edgeworth"){
    for (u in 1:(length(NYear)-1)){
      X0 <- list(); XT <- list()
      for (v in LF){
        X0[[v]] <- xdata[[v]][which(Year == NYear[u]),]
        XT[[v]] <- xdata[[v]][which(Year == NYear[u+1]),]
      }
      
      for (w in LF){
        Vx1 <- XT[[w]] - X0[[w]]
        for (r in LF[-w]){
          Vx1 <- Vx1 * (X0[[r]] + XT[[r]])/2
        }
        Dx[which(Dx$Year == NYear[u+1] & Dx$Factor == w),-c(1:3)] <- Vx1
      }
    }
  } else if (method == "Walsh"){
    for (u in 1:(length(NYear)-1)){
      X0 <- list(); XT <- list()
      for (v in LF){
        X0[[v]] <- xdata[[v]][which(Year == NYear[u]),]
        XT[[v]] <- xdata[[v]][which(Year == NYear[u+1]),]
      }
      
      for (w in LF){
        Vx1 <- XT[[w]] - X0[[w]]
        for (r in LF[-w]){
          Vx1 <- Vx1 * sqrt(X0[[r]] * XT[[r]])
        }
        Dx[which(Dx$Year == NYear[u+1] & Dx$Factor == w),-c(1:3)] <- Vx1
      }
    }
  } else {
    stop("Input method is not correct")
  }
  
  Dx.Year <- rowSums(Dx[,-c(1:3)], na.rm = TRUE)
  Dx.Year <- cbind(Dx[,c(1:3)], Dx.Year)
  
  Dx.Category <- aggregate(Dx.Year$Dx.Year, 
                              by=list(Dx.Year$Year, Dx.Year$Category), FUN=sum, na.rm = TRUE)
  names(Dx.Category) <- c("Year", "Category", "CEC")
  
  Dx.Factor <- aggregate(Dx.Year$Dx.Year, 
                        by=list(Dx.Year$Year, as.character(Dx.Year$Factor)), FUN=sum, na.rm = TRUE)
  names(Dx.Factor) <- c("Year", "Factor", "CEC")
  Dx.Factor0 <- Dx.Factor
  Dx.Factor <- as.data.frame(matrix(Dx.Factor$CEC, length(NYear[-1]), length(Factor)))
  Dx.Factor <- cbind(NYear[-1], Dx.Factor)
  names(Dx.Factor) <- c("Year", as.character(Factor))
  
  Dx.Factor0$CECsum <- aggregate(Dx.Factor0$CEC, by=list(Dx.Factor0$Year), FUN=sum, na.rm=TRUE)$x
  CEC <- NULL; CECsum <- NULL
  plotDxFactor <- ggplot(Dx.Factor0, aes(Year, CEC)) + 
    geom_bar(aes(fill = Factor), stat = "identity", width = 0.6) +
    scale_fill_hue(l=70) +
    geom_line(data = subset(Dx.Factor0, variable=Factor),   
              aes(x=Year, y=CECsum, group = 1), size = 1) +
    geom_point(data = subset(Dx.Factor0, variable=Factor),   
               aes(x=Year, y=CECsum, group = 1), size = 2, shape = 15) + 
    theme_bw()
  
  Dxresult <- list("DX" = Dx, "Dx.Factor" = Dx.Factor, "Dx.Category" = Dx.Category, 
                   "Dx.Year" = Dx.Year, "plotDxFactor" = plotDxFactor)
  
  return(Dxresult)
}
  
