#' Title Individual dynamic prediction function.
#' @usage predict_LM(object,newdata)
#' @param object an objects inheriting from class "LMf" or "Vs_LM"
#' @param newdata data.frame
#' @return
#' \item{out} Predicted values.
#' @export
#' @examples
#' library(DynamicModel)
#' data(renal)
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=TRUE)
#' newdata<-renal[which(renal$id==5619),]
#' aa<-predict_LM(fit,newdata)
predict_LM<-function(object,newdata){
  if (!require("dynpred")) {
    install.packages("dynpred")
    library(dynpred)
  }
  fm<-object$Model
  time<-object$time
  status<-object$status
  data<-object$data$data
  sl<-object$tw$sl
  w<-object$tw$w
  rtime<-object$rtime
  fixed<-object$covs$fixed
  vary<-object$covs$vary
  func_covars<-object$func_covars
  func_lms<-object$func_lms
  tt<-seq(min(sl),max(sl),length=1001)
  inter<-findInterval(tt,newdata[[rtime]])
  out<-data.frame(time=tt,csurv=NA)
  for(i in 1:length(tt)){
    data_pre<-newdata[inter[i],]
    data_pre$LM<-tt[i]
    data_pred<-NULL
    data_pred$data<-data_pre
    data_pred$lm_col<-"LM"
    dt<-add_interactions(data_pred,c(fixed,vary),func_covars=func_covars,func_lms =func_lms,sl=sl)
    sfit<-summary(survfit(fm,newdata=dt$data))
    tmp<-dynpred::evalstep(sfit$time,sfit$surv,c(tt[i],tt[i]+w),subst=1)
    #tmp_lower<-evalstep(sfit$time,sfit$lower,c(tt[i],tt[i]+w),subst=1)
    #tmp_upper<-evalstep(sfit$time,sfit$upper,c(tt[i],tt[i]+w),subst=1)
    out[i,2]<-tmp[2]/tmp[1]
    #out[i,3]<-tmp_lower
    }
  out
}
