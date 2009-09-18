.onLoad <- function(lib, pkg) {
	invisible(.C("FPUcheck",PACKAGE="PearsonDS"))
}

.onAttach <- function(lib, pkg) {
 	if ("gsl" %in% installed.packages()[,1]) {
    assign(".hasGSL",TRUE,pos=paste("package:",pkg,sep=""))
  } else {
    assign(".hasGSL",FALSE,pos=paste("package:",pkg,sep=""))
  }   
}
