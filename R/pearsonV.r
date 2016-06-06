dpearsonV <- function(x,shape,location,scale,params,log=FALSE) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- 1/abs(scale); ssgn <- sign(scale)                                   # 2016-06-06: replaced abs(scale) with 1/abs(scale), bug reported by Andrew Beck
  if (log) {
    tmp <- dgamma(1/pmax(ssgn*(x-location),0),shape=shape,scale=gscale,
                  log=TRUE)-2*log(pmax(ssgn*(x-location),0))                           
    tmp[is.na(tmp)] <- -Inf
    tmp         
  } else {
    ifelse(x==location,0,dgamma(1/pmax(ssgn*(x-location),0),shape=shape,
                                scale=gscale)/((x-location)^2))
  }  
}

ppearsonV <- function(q,shape,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- 1/abs(scale); ssgn <- sign(scale)                                   # 2016-06-06: replaced abs(scale) with 1/abs(scale), bug reported by Andrew Beck
  ifelse(q==location,as.numeric(!(scale>0)),
         pgamma(1/pmax(ssgn*(q-location),0),shape=shape,scale=gscale,
           lower.tail=xor(scale>0,lower.tail),log.p=log.p))
}

qpearsonV <- function(p,shape,location,scale,params,lower.tail=TRUE,log.p=FALSE) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- 1/abs(scale); ssgn <- sign(scale)                                   # 2016-06-06: replaced abs(scale) with 1/abs(scale), bug reported by Andrew Beck
  ssgn/qgamma(p,shape=shape,scale=gscale,lower.tail=xor(scale>0,lower.tail),
              log.p=log.p)+location
}

rpearsonV <- function(n,shape,location,scale,params) {
  if (!missing(params)) { shape <- params[[1]]; 
                          location <- params[[2]]; scale <- params[[3]] }
#  stopifnot((shape>0)&&(scale!=0))
  gscale <- 1/abs(scale); ssgn <- sign(scale)                                   # 2016-06-06: replaced abs(scale) with 1/abs(scale), bug reported by Andrew Beck
  ssgn/rgamma(n,shape=shape,scale=gscale)+location
}
