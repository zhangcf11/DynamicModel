#' Title output result of function "VS_LM"
#'
#' @param object an object inheriting from class VS_LM
#' @export
print.VS_LM<-function(object,...){
  print(object$Model,...)
}
