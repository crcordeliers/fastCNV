#' cit.discretize
#'
#' @param x ...
#' @param lim ...
#' @param quant ...
#' @param addlevels ...
#'
#' @import stats


cit.discretize <-
  function(x, lim, quant = FALSE, addlevels = FALSE){

    lim <- sort(lim)

    if( quant & ( any(lim>1) | any(lim<0) ) )
      stop("lim must be [0;1] as quant=T\n")

    res <- rep(NA, length(x))

    if(quant) lim  <- stats::quantile(x, probs=lim, na.rm=TRUE)
    n <- length(lim)
    for(i in n:1) res[which(x<lim[i])] <- i
    res[which(x>=lim[n])] <- n+1

    if (addlevels) {
      res <- as.factor(res)

      if(quant)
        lim <- gsub(" ", "", prettyNum(lim, format="g", digits=1))

      if(length(lim)==1)
        lev <- c(paste("<", lim[1], sep = ""), paste(">=", lim[1], sep = ""))
      else
        lev <- c(paste("<", lim[1], sep = ""), paste("[",lim[-length(lim)],";",lim[-1],"[",sep=""), paste(">=", lim[length(lim)], sep = ""))
      levels(res) <- lev[as.numeric(levels(res))]
    }
    res
  }
