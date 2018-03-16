Layer <- setClass("Layer",
                  slots=c(
                    m="complex",
                    r="numeric"
                    ))

Scatterer <- setClass("Scatterer",
               slots=c(
                 na="numeric",
                 lambda="numeric",
                 nt="numeric",
                 layers="list"),
               prototype=list(
                 na = 1.0,
                 lambda = 500,
                 nt=0
               )
               )

setMethod("show", signature(object="Layer"),function(object){
  cat("Layer Data \n")
  cat("Diameter:", 2*object@r," nm\n")
  cat("Complex Index:", object@m,"\n")
})

setMethod("show", signature(object="Scatterer"),function(object){
  cat("Scatterer Data \n")
  cat("Lambda:", object@lambda," nm\n")
  cat("Ambient Index", object@na, "\n")
})

setGeneric("na<-",function(x,value) standardGeneric("na<-"))
setGeneric("lambda<-",function(x,value) standardGeneric("lambda<-"))
setGeneric("nt<-",function(x,value) standardGeneric("nt<-"))
setGeneric("r<-",function(x,value) standardGeneric("r<-"))
setGeneric("m<-",function(x,value) standardGeneric("m<-"))


setReplaceMethod("na","Scatterer", function(x,value) {x@na <- value; validObject(x); x})
setReplaceMethod("lambda","Scatterer", function(x,value) {x@lambda <- value; validObject(x); x})
setReplaceMethod("nt","Scatterer", function(x,value) {x@nt <- value; validObject(x); x})
setReplaceMethod("r","Layer", function(x,value) {x@r <- value; validObject(x); x})
setReplaceMethod("m","Layer", function(x,value) {x@m <- value; validObject(x); x})

setGeneric("scattnlay", function(object) {standardGeneric("scattnlay")})

setMethod("scattnlay",signature(object="Scatterer"),function(object){
  Rpp<-S4_SCATTNLAY(object)
  return(Rpp)
})

setMethod("+", signature(e1="Scatterer",e2="Layer"), function(e1,e2){
  e1@layers <- c(e1@layers,e2)
  structure(e1,class="Scatterer")
})

