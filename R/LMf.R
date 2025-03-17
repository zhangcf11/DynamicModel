#' Title The Construction of dynamic prediction model based on landmark
#' @usage LMf (data,tw=list(sl=sl,w=w),cov=list(fixed=fixed,vary=vary),id=id,rtime=rtime,time=time,status=status,inter=T,func_covars=c("linear","quadratic"),func_lms=c("linear","quadratic"))
#' @param data data.frames.
#' @param tw Time point options, including landmark time points (sl) and prediction windows (w).
#' @param cov Covariate selection, including baseline covariates: fixed and longitudinal covariates: vary.
#' @param id Individual ID
#' @param rtime Observation time column name
#' @param time Individual's survival time
#' @param status Individual's survival status
#' @param inter Indicator variables: whether covariates in the model need to include dynamic effects.
#' @param func_covars Functional forms for the dynamic effects of covariates. It mainly includes linear functions, quadratic functions, exponential functions, and logarithmic functions, with the default being linear and quadratic functions.
#' @param func_lms Functional forms for the dynamic effects of baseline hazard rates.It mainly includes linear functions, quadratic functions, exponential functions, and logarithmic functions, with the default being linear and quadratic functions.
#'
#' @return A list of class LMf with components:
#'  \describe{
#'    \item{Model}{The result of parameter estimation.}
#'    \item{data}{The Landmark dataset.}
#'    \item{func_covars}{Functional forms for the dynamic effects of covariates.}
#'    \item{func_lms}{Functional forms for the dynamic effects of baseline hazard rates.}
#'    \item{cov}{A list of covariates type. The components of this list are:
#'      \itemize{
#'        \item{fixed}{Baseline covariates.}
#'        \item{vary}{Longitudinal covariates.}
#'      }
#'    }
#'    \item{tw}{A list of time points. The components of this list are:
#'      \itemize{
#'        \item{sl}{Landmark time points.}
#'        \item{w}{Prediction window.}
#'      }
#'    }
#'    \item{time}{Individual's survival time.}
#'    \item{status}{Individual's survival status.}
#'    \item{id}{Individual ID.}
#'    \item{rtime}{Observation time column name.}
#'  }
#' @export
#' @examples
#' library(DynamicModel)
#' data(renal)
#' head(renal)
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=T)
#' print(fit)
LMf<-function(data,tw=list(sl=sl,w=w),cov=list(fixed=fixed,vary=vary),id=id,rtime=rtime,time=time,status=status,
              inter=T,func_covars=c("linear","quadratic"),func_lms=c("linear","quadratic")){
  LMdata<-NULL
  for (j in seq(along=tw$sl)){
    LM<-dynpred::cutLM(data,outcome=list(time=time,status=status),
              LM=tw$sl[j],horizon=tw$sl[j]+tw$w,covs=list(fixed=cov$fixed,varying=cov$vary),format="long",id=id,rtime=rtime,right=F)
    #LM$pseudo=pseudomean(time=LM[[time]]-tw$sl[j],event=LM[[status]],tmax=tw$w)
    LMdata<-rbind(LMdata,LM)
  }
  LMdata<-LMdata[order(LMdata[[id]]),]
  LMdata1<-NULL
  LMdata1$data<-LMdata
  LMdata1$lm_col="LM"
  if (inter==T){
    LMdata_i<-add_interactions(LMdata1,c(cov$fixed,cov$vary),func_covars=func_covars,func_lms =func_lms,sl=tw$sl)
    found_cols<-c()
    for (col in c(cov$fixed,cov$vary)) {
      matched_cols <- grep(paste0("^", col, "_"), colnames(LMdata_i$data), value = TRUE)
      found_cols <- c(found_cols, matched_cols)
    }
    LM_cols<-grep(paste0("^","LM","_"),colnames(LMdata_i$data),value=TRUE)
    f<-paste(c(cov$fixed,cov$vary,found_cols,LM_cols),collapse='+')
  }
  else{
    LMdata_i<-NULL
    LMdata_i$data<-LMdata1$data
    f<-paste(c(cov$fixed,cov$vary),collapse='+')
  }

  f1<-as.formula(paste("Surv(LM",",",time,",",status,")~",f,"+","cluster(",id,")"))
  Model<-coxph(f1,data=LMdata_i$data,method="breslow")

  object<-NULL
  object$Model<-Model
  object$data<-LMdata_i
  if (inter==T)
  {object$func_covars<-func_covars
  object$func_lms<-func_lms}
  object$covs$fixed<-cov$fixed
  object$covs$vary<-cov$vary
  #object$tw<-NULL
  object$tw$sl<-tw$sl
  object$tw$w<-tw$w
  object$time<-time
  object$status<-status
  object$id<-id
  object$rtime<-rtime
  object
}
add_interactions <- function(lmdata, lm_covs, func_covars, func_lms, sl,lm_col,
                             keep = T){
  if (missing(lm_col)){
    lm_col <- lmdata$lm_col
  }
  if (lm_col %in% func_covars){
    stop(paste0("arg lm_col (given as/inferred as ",lm_col,
                ") should not be in arg func_covars."))
  }
  data <- lmdata$data

  f1 <- function(t) t
  f2 <- function(t) t^2
  f3 <- function(t) log(1 + t)
  f4 <- function(t) exp(t)

  if (missing(func_covars)) func_covars <- list(f1, f2)
  if (missing(func_lms)) func_lms <- list(f1, f2)
  if (inherits(func_covars, "character")){
    funcs <- list()
    if ("linear" %in% func_covars) funcs <- c(funcs, list(f1))
    if ("quadratic" %in% func_covars) funcs <- c(funcs, list(f2))
    if ("log" %in% func_covars) funcs <- c(funcs, list(f3))
    if ("exp" %in% func_covars) funcs <- c(funcs, list(f4))
    func_covars <- funcs
  }
  if (inherits(func_lms, "character")){
    funcs <- list()
    if ("linear" %in% func_lms) funcs <- c(funcs, list(f1))
    if ("quadratic" %in% func_lms) funcs <- c(funcs, list(f2))
    if ("log" %in% func_lms) funcs <- c(funcs, list(f3))
    if ("exp" %in% func_lms) funcs <- c(funcs, list(f4))
    func_lms <- funcs
  }

  all_covs <- c(lm_covs)
  data_LM <- data[[lm_col]]
  # Add func_covarss: covariate LM interactions
  for(i in 1:length(lm_covs)){
    for (j in 1:length(func_covars)){
      f <- func_covars[[j]]
      name <- paste0(lm_covs[i],"_",j)
      data[[name]]  <- data[[lm_covs[i]]]*f(data_LM/(max(sl)-min(sl)))
      all_covs <- c(all_covs, name)
    }
  }
  # Add func_lms: LM interactions
  for (k in 1:length(func_lms)){
    g <- func_lms[[k]]
    name <- paste0("LM_",k)
    data[[name]]  <- g(data_LM/(max(sl)-min(sl)))
    all_covs <- c(all_covs, name)
  }

  if(!keep){
    remaining = colnames(data)[! colnames(data)  %in% lm_covs]
    data <- data[remaining]
  }
  lmdata$data <- data

  lmdata$func_covars <- func_covars
  lmdata$func_lms <- func_lms
  lmdata$lm_covs <- lm_covs
  lmdata$all_covs <- unique(all_covs)
  lmdata$lm_col <- lm_col

  return(lmdata)
}






