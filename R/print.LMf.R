#' Title output result of function "LMf"
#'
#' @param object an object inheriting from class LMf
#' @export
print.LMf<-function(object,...){
  print(object$Model,...)
}
