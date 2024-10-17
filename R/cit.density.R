#' cit.density
#'
#' @param x ...
#' @param doplot ...
#' @param pc ...
#' @param ... ...
#'
#' @import graphics
#'
#' @return


cit.density <-
  function(x,doplot=FALSE,pc=.05,...){
    dx <- density(x,na.rm=TRUE,...)
    ymax <- diff(range(dx$y))
    n <- length(dx$y)
    croissance <- as.numeric((dx$y[-1] - dx$y[-n])>0)
    wB <- 1+which(diff(croissance)== 1)
    wH <- 1+which(diff(croissance)== -1)
    pointsBas   <- dx$x[wB]
    pointsHauts <- dx$x[intersect(wH,which(dx$y>pc*ymax))]

    if(length(pointsHauts)>0)pointsBas   <- sapply(split(pointsBas,cit.discretize(pointsBas,pointsHauts)),median)

    if(doplot){
      plot(dx,...)
      abline(v=pointsBas,lty=3,col="blue")
      abline(v=pointsHauts,lty=3,col="red")
    }
    L <- c(dx,list("down"=pointsBas,"top"=pointsHauts))
    attr(L,"class") <- "density"
    L
  }
