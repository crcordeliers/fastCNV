#' cit.peaks
#'
#' @param x ...
#' @param percentHighestPeak ...
#' @param maxNbPeaks ...
#' @param minDeltaBetweenPeaks ...
#' @param deltaApproach ...
#' @param doplot ...
#' @param bw ...
#' @param ... ...
#'
#' @import graphics
#'
#' @return ...
#'

cit.peaks <-
  function(x,
           percentHighestPeak=.2,
           maxNbPeaks=NULL,
           minDeltaBetweenPeaks=.03,
           deltaApproach=1,
           doplot=FALSE, bw = "nrd0",
           ...){

    if(is.na( minDeltaBetweenPeaks ))minDeltaBetweenPeaks <- NULL
    v <- cit.density(x,pc =percentHighestPeak,bw=bw)
    m <- v$y[which(v$x %in% v$top)]
    w <- which(m/max(m) > percentHighestPeak)

    xw <- v$top[w]

    which.eq <- function(z){ order(abs(v$x-z))[1]}

    # umin[which(apply(t(sapply(umin,function(z)c( max(which(umax <=z )),min(which(umax>=z))))),1,function(z)!any(is.infinite(z))))]

    if(deltaApproach==1 & !is.null( minDeltaBetweenPeaks ) & length(xw)>1){
      for(i in 1:(length(xw)-1)){
        if(xw[i+1]-xw[i] <  minDeltaBetweenPeaks) {
          xw[i:(i+1)] <- xw[c(i:(i+1))][which.max(m[i:(i+1)])]
        }
      }

      xw <- unique(xw)

    }


    if(!is.null(maxNbPeaks)){
      if(length(xw)> maxNbPeaks)
        for(i in 1:(length(xw)-maxNbPeaks)){
          oneamong <- xw[c(0,1)+which.min(diff(xw))]
          w1 <-  which.eq(oneamong[1])
          w2 <-  which.eq(oneamong[2])
          out <- oneamong[which.min(v$y[c(w1,w2)])]
          xw <- setdiff(xw,out)

        }
    }


    umin <- NULL
    if(length(xw)>1){
      for(i in 1:(length(xw)-1)){
        possiblemin <- v$down[which(v$down >= xw[i] & v$down <= xw[i+1])]
        if(length(possiblemin)>1){
          w<- sapply(possiblemin,which.eq)
          possiblemin <- possiblemin[which.min(v$y[w])]
        }
        umin <- c(umin,possiblemin)
      }

    }


    #  if(doplot){
    #            plot(v,...)
    #            abline(v=xw,col="red")
    #            if(length(umin)>0)abline(v=umin,col="green",lty=3)
    #        }

    L <- list("x abciss big peaks"=xw,
              "x abciss inter-peaks"=umin,
              "nb big peaks"=length(xw))

    peaks <- L[[1]]
    interpeaks <- L[[2]]
    temp <- as.data.frame(t(sapply(peaks,function(pic) {
      wleft <- which(interpeaks < pic)
      if(length(wleft)>0){
        wleft <- max(interpeaks[wleft])
      }else{
        wleft <- min(x,na.rm=TRUE)
      }
      wright <- which(interpeaks > pic)
      if(length(wright)>0){
        wright <- min(interpeaks[wright])
      }else{
        wright <- max(x,na.rm=TRUE)
      }
      c(pic,wleft,wright,length(which(x>=wleft & x<=wright)))
    }) ))
    names(temp) <- c("peak","left born","right born","size")
    L$"peaks x size" <- temp


    if(deltaApproach==2 &!is.null( minDeltaBetweenPeaks ) & L$'nb big peaks'>1){
      Lini <- L
      names(Lini) <- paste("initial",names(L))

      while(any( diff(L$"x abciss big peaks") < minDeltaBetweenPeaks) & L$"nb big peaks" >1){

        w <- which.min(abs(diff(L$"x abciss big peaks")) )


        lb <- min(L$"peaks x size"[w:(w+1),"left born"])
        rb <- max(L$"peaks x size"[w:(w+1),"right born"])
        si <- sum(L$"peaks x size"[w:(w+1),"size"])
        pe <- median(x[which(x>=lb & x<=rb)])

        L$"peaks x size"[w,] <- c(pe,lb,rb,si)
        L$"peaks x size" <- L$"peaks x size"[-(w+1),]

        L$"nb big peaks" <- L$"nb big peaks" - 1
        L$"x abciss big peaks"[w] <- pe
        L$"x abciss big peaks" <- L$"x abciss big peaks"[-(w+1)]
        L$"x abciss inter-peaks" <- L$"x abciss inter-peaks"[-(w+1)]
      }

      L <- c(Lini,L)
    }

    if(doplot){
      plot(v,...)
      abline(v=L$"x abciss big peaks",col="red")
      if(length(L$"x abciss inter-peaks")>0)abline(v=L$"x abciss inter-peaks",col="green",lty=3)
    }

    L
  }
