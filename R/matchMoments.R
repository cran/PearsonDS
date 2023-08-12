matchMoments <- function(mean,variance,skewness=NA,kurtosis=NA,type,moments,
                         skewness.sign=c("+","-"),return.distribution=FALSE) {
  sf <- if(match.arg(skewness.sign)=="+") 1 else -1
  suffix <- c("0","I","II","III","IV","V","VI","VII")
  if (is.numeric(type)) type <- suffix[type+1]
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  if (is.na(kkk)) {
    kkk <- switch(type,
      "0"   = 3,
      "III" = 3 + 1.5*sss^2,
      "V"   = if (sss^2<32) -3*(16+13*sss^2+2*sqrt(64+48*sss^2+12*sss^4+sss^6))/
                            (sss^2-32) else 
                stop("moments and distribution type do not fit"),
      stop("kurtosis must be provided for this type")
    )
  } 
  if (is.na(sss)) {
    sss <- switch(type,
      "0"   = 0,
      "II"  = 0,
      "VII" = 0,
      "III" = sf*sqrt((kkk-3)/1.5),
      "V"   = sf*1/72*(sss^4+78*sss^2-(sss^2+3)^1.5*sqrt(sss^2+147)-63),
      stop("skewness must be provided for this type")
    )
  }
  res <- pearsonFitM(mmm,vvv,sss,kkk)
  if (suffix[res$type+1]!=type) stop("moments and distribution type do not fit")
  if (return.distribution) res else 
    c(mean=mmm,variance=vvv,skewness=sss,kurtosis=kkk)
} 
