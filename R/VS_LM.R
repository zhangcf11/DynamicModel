#' Title Function for variable selection in landmark-based dynamic prediction models.
#' @usage VS_LM(object,cov=list(fixed=fixed,vary=vary))
#' @param object an objects inheriting from class "LMf"
#' @param cov A list of covariates selection options, including baseline covariates: fixed and longitudinal covariates: vary.
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
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=TRUE)
#' a1<-VS_LM(fit,cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")))
#' print(a1)
VS_LM<-function(object,cov=list(fixed=fixed,vary=vary)){
  if (!require("MASS")) {
    install.packages("MASS")
    library(MASS)
  }
  cl <- match.call()
  fm<-object$Model
  LMdata_i<-object$data
  f1<-paste(c(cov$fixed,cov$vary),collapse='+')
  id<-object$id
  sl<-object$tw$sl
  w<-object$tw$w
  time<-object$time
  status<-object$status
  if (is.null(cov$vary)){
    NULL
  }
  else{
    rtime<-object$rtime
  }
  func_covars<-object$func_covars
  func_lms<-object$func_lms
  f<-as.formula(paste("Surv(LM,",time,",",status,")","~",f1))
  scope<-list(lower=f)
  fit<-MASS::stepAIC(fm,scope=scope, direction="backward")
  object<-NULL
  object$Model<-fit
  object$data<-LMdata_i
  object$func_covars<-func_covars
  object$func_lms<-func_lms
  object$covs$fixed<-cov$fixed
  object$covs$vary<-cov$vary
  object$tw<-NULL
  object$tw$sl<-sl
  object$tw$w<-w
  object$time<-time
  object$status<-status
  object$id<-id
  if (is.null(cov$vary)){
    object$rtime<-NULL
  }
  else{
    object$rtime<-rtime
  }
  object$call <- cl
  class(object) <- "VS_LM"
  object
}
