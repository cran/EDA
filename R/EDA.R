#' Energy consumption analysis for calculating carbon emission changes
#'
#' @usage EDA(cdata, factordata, Year = Year, Factor = Factor, 
#'     Fuel = 1, Sector = 1, method = "LMDI")
#' \method{print}{EDA}(x, ...)
#' \method{plot}{EDA}(x, ...)
#' 
#' @aliases EDA print.EDA plot.EDA
#' 
#' @param cdata A data.frame of annual carbon emission or energy consumption 
#' data, which can include multiple Fuels stored by columns.
#' @param factordata A list of factors' data.frame.
#' @param Year A numeric vector of year.
#' @param Sector A vector of carbon emission or energy consumption 
#' sector names or number. If only one sector of carbon emission or 
#' energy consumption, set \code{Sector = 1}. 
#' @param Factor A vector of factor names.
#' @param Fuel A vector of fuel names.
#' @param method A character of energy consumption analysis method's name. 
#' One of "\code{\link{LMDI}}", "\code{Laspeyres}", "\code{Paasche}", 
#' "\code{Marshall-Edgeworth}" or "\code{Walsh}".
#' @param x A list of EDA result.
#' @param ... Ignore
#' 
#' @importFrom stats aggregate
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw scale_fill_hue 
#' geom_line geom_point labs
#'
#' @author Yongze Song \email{yongze.song@postgrad.curtin.edu.au}
#' and Peng Wu \email{peng.wu@curtin.edu.au}.
#' 
#' @references Ang, B. W. (2005). The LMDI approach to decomposition 
#' analysis: a practical guide. Energy policy, 33(7), 867-871.
#' @references Marlay, R. C. (1984). Trends in industrial use of energy. 
#' Science, 226, 1277-1284.
#' @references Paasche, H. Uber die Preisentwicklung der letzten Jahre. 
#' Jahrbiicher fur Nationalokonomie und Statistik, 23(1874), 168.
#' @references Marshall, A. (1887). Remedies for fluctuations of 
#' general prices.
#' @references Edgeworth, F. Y. (1925). Papers relating to political economy.
#' @references Walsh, C. M. (1921). The Problem of Estimation, a 
#' Seventeenth-century Controversy and Its Bearing on Modern Statistical 
#' Questions, Especially Index-numbers, by Correa Moylan Walsh.
#'
#' @seealso \code{\link{LMDI}}
#' 
#' @examples
#' library(EDA)
#' data(carbon)
#' data(factordata)
#' ## set parameters
#' cdata <- carbon[,-c(1,2)]
#' Year <- 2001:2005
#' Sector <- c("b1", "b2", "b3")
#' Fuel <- colnames(cdata)
#' Factor <- names(factordata)
#' ## run EDA model
#' eda1 <- EDA(cdata, factordata, Year = Year, Factor = Factor, 
#'     Fuel = Fuel, Sector = Sector, method = "LMDI")
#' eda1
#' plot(eda1)
#' 
#' @export

EDA <- function(cdata, factordata, Year = Year, Factor = Factor, 
                Fuel = 1, Sector = 1, method = "LMDI"){
  Dx <-  expand.grid(Year = Year[-1], Sector = Sector, Factor = Factor)
  Dx$Sector <- factor(Dx$Sector, as.character(Sector))
  Dx$Factor <- factor(Dx$Factor, as.character(Factor))
  Dxc <- matrix(NA, nrow(Dx), length(Fuel))
  Dxc <- as.data.frame(Dxc)
  names(Dxc) <- Fuel
  Dx <- cbind(Dx, Dxc)
  
  LF <- c(1:length(Factor))
  
  if (method == "LMDI"){
    for (u in 1:(length(Year)-1)){
      C0 <- cdata[which(Year == Year[u]),]
      CT <- cdata[which(Year == Year[u+1]),]
      for (v in LF){
        X0 <- factordata[[v]][which(Year == Year[u]),]
        XT <- factordata[[v]][which(Year == Year[u+1]),]
        Dx[which(Dx$Year == Year[u+1] & Dx$Factor == Factor[v]),-c(1:3)] <- LMDI(C0, CT, X0, XT)$Dx
      }
    }
  } else if (method == "Laspeyres"){
    for (u in 1:(length(Year)-1)){
      X0 <- list(); XT <- list()
      for (v in LF){
        X0[[v]] <- factordata[[v]][which(Year == Year[u]),]
        XT[[v]] <- factordata[[v]][which(Year == Year[u+1]),]
      }
      
      for (w in LF){
        Vx1 <- XT[[w]] - X0[[w]]
        for (r in LF[-w]){
          Vx1 <- Vx1 * X0[[r]]
        }
        Dx[which(Dx$Year == Year[u+1] & Dx$Factor == w),-c(1:3)] <- Vx1
      }
    }
  } else if (method == "Paasche"){
    for (u in 1:(length(Year)-1)){
      X0 <- list(); XT <- list()
      for (v in LF){
        X0[[v]] <- factordata[[v]][which(Year == Year[u]),]
        XT[[v]] <- factordata[[v]][which(Year == Year[u+1]),]
      }
      
      for (w in LF){
        Vx1 <- XT[[w]] - X0[[w]]
        for (r in LF[-w]){
          Vx1 <- Vx1 * XT[[r]]
        }
        Dx[which(Dx$Year == Year[u+1] & Dx$Factor == w),-c(1:3)] <- Vx1
      }
    }
  } else if (method == "Marshall-Edgeworth"){
    for (u in 1:(length(Year)-1)){
      X0 <- list(); XT <- list()
      for (v in LF){
        X0[[v]] <- factordata[[v]][which(Year == Year[u]),]
        XT[[v]] <- factordata[[v]][which(Year == Year[u+1]),]
      }
      
      for (w in LF){
        Vx1 <- XT[[w]] - X0[[w]]
        for (r in LF[-w]){
          Vx1 <- Vx1 * (X0[[r]] + XT[[r]])/2
        }
        Dx[which(Dx$Year == Year[u+1] & Dx$Factor == w),-c(1:3)] <- Vx1
      }
    }
  } else if (method == "Walsh"){
    for (u in 1:(length(Year)-1)){
      X0 <- list(); XT <- list()
      for (v in LF){
        X0[[v]] <- factordata[[v]][which(Year == Year[u]),]
        XT[[v]] <- factordata[[v]][which(Year == Year[u+1]),]
      }
      
      for (w in LF){
        Vx1 <- XT[[w]] - X0[[w]]
        for (r in LF[-w]){
          Vx1 <- Vx1 * sqrt(X0[[r]] * XT[[r]])
        }
        Dx[which(Dx$Year == Year[u+1] & Dx$Factor == w),-c(1:3)] <- Vx1
      }
    }
  } else {
    stop("Input method is not correct")
  }
  
  Dx.Year <- rowSums(Dx[,-c(1:3)], na.rm = TRUE)
  Dx.Year <- cbind(Dx[,c(1:3)], Dx.Year)
  
  Dx.Sector <- aggregate(Dx.Year$Dx.Year, 
                              by=list(Dx.Year$Year, Dx.Year$Sector), FUN=sum, na.rm = TRUE)
  names(Dx.Sector) <- c("Year", "Sector", "CEC")
  
  Dx.Factor <- aggregate(Dx.Year$Dx.Year, 
                        by=list(Dx.Year$Year, Dx.Year$Factor), FUN=sum, na.rm = TRUE)
  names(Dx.Factor) <- c("Year", "Factor", "CEC")
  
  Dx.Factor <- as.data.frame(matrix(Dx.Factor$CEC, length(Year[-1]), length(Factor)))
  Dx.Factor <- cbind(Year[-1], Dx.Factor)
  names(Dx.Factor) <- c("Year", as.character(Factor))

  Dxresult <- list("Dx.Year" = Dx.Year, "Dx.Factor" = Dx.Factor, 
                   "Dx.Sector" = Dx.Sector, "DX" = Dx)
  class(Dxresult) <- "EDA"
  return(Dxresult)
}
  
print.EDA <- function(x, ...){
  cat("Annual changes:\n")
  Dx.Year <- x$Dx.Year
  ac <- aggregate(Dx.Year$Dx.Year, by = list(Dx.Year$Year), FUN = sum, na.rm = TRUE)
  names(ac) <- c("Year", "Change")
  print(ac)
  cat("Total change: ", sum(ac$Change))
}

plot.EDA <- function(x, ...){
  Dx.Year <- x$Dx.Year
  Dx.Sector <- x$Dx.Sector
  Year <- NULL; Sector <- NULL; Factor <- NULL
  # summary by sector
  plotDxSector <- ggplot(Dx.Sector, aes(Year, CEC)) + 
    geom_bar(aes(fill = Sector), stat = "identity", width = 0.6) +
    scale_fill_hue(l=70) +
    theme_bw() + 
    labs(y = "Carbon Emission Change")
  print(plotDxSector)
  # summary by factor
  Dx.Factor0 <- aggregate(Dx.Year$Dx.Year, 
                         by=list(Dx.Year$Year, Dx.Year$Factor), FUN=sum, na.rm = TRUE)
  names(Dx.Factor0) <- c("Year", "Factor", "CEC")
  Dx.Factor0$CECsum <- aggregate(Dx.Factor0$CEC, by=list(Dx.Factor0$Year), FUN=sum, na.rm=TRUE)$x
  CEC <- NULL; CECsum <- NULL
  plotDxFactor <- ggplot(Dx.Factor0, aes(Year, CEC)) + 
    geom_bar(aes(fill = Factor), stat = "identity", width = 0.6) +
    scale_fill_hue(l=70) +
    geom_line(data = subset(Dx.Factor0, variable=Factor),   
              aes(x=Year, y=CECsum, group = 1), size = 1) +
    geom_point(data = subset(Dx.Factor0, variable=Factor),   
               aes(x=Year, y=CECsum, group = 1), size = 2, shape = 15) + 
    theme_bw() + 
    labs(y = "Carbon Emission Change")
  print(plotDxFactor)
}