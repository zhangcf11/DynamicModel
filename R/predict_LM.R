#' Predict Survival Probabilities from Landmarking Models
#'
#' This function generates survival probability predictions from a fitted landmarking model
#' at specified landmark time points for new data.
#' @usage predict_LM(object,newdata)
#' @param object an objects inheriting from class "LMf" or "Vs_LM"
#' @param newdata A data.frame containing the new data for which predictions are made
#' @return A data.frame with two columns:
#' \itemize{
#'   \item{time: The landmark time points}
#'   \item{csurv: The conditional survival probabilities from each landmark time point
#'                to the landmark time plus the prediction window}
#' }
#' @importFrom dynpred evalstep
#' @importFrom survival survfit
#' @export
#' @examples
#' library(DynamicModel)
#' data(renal)
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=TRUE)
#' a1<-VS_LM(fit,cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")))
#' newdata<-renal[which(renal$id==5619),]
#' aa<-predict_LM(a1,newdata)
predict_LM<-function(object,newdata){
  if (!require("dynpred")) {
    install.packages("dynpred")
    library(dynpred)
  }
  # Input validation
if (!inherits(object, "LMf") && !inherits(object, "VS_LM")) {
  stop("object must be of class 'LMf' or 'VS_LM'")
}

  if (missing(newdata) || !is.data.frame(newdata)) {
    stop("newdata must be a data.frame")
  }
  # Extract model components
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

  #tt<-seq(min(sl),max(sl),length=1001)
  tt<-sl

  # For models with time-varying covariates, find appropriate intervals
  if (!is.null(rtime)) {
    inter <- findInterval(tt, newdata[[rtime]])
  }
  out<-data.frame(time=tt,csurv=NA)
  for(i in seq_along(tt)){
    # Prepare data for current landmark time point
    if (is.null(rtime)) {
      # For fixed covariates only, use the entire newdata
      data_pre <- newdata
    } else {
      # For time-varying covariates, select the appropriate row
      if (inter[i] > 0 && inter[i] <= nrow(newdata)) {
        data_pre <- newdata[inter[i], ]
      } else {
        warning(paste("No appropriate data found for landmark time", tt[i]))
        next
      }
    }
    data_pre$LM<-tt[i]
    data_pred<-list()
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

