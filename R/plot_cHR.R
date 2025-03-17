#' Title Function for dynamic effect hazard ratio (HR) plots.
#' @usage plot_cHR(object,covars,conf_int=TRUE,IDlegend=NULL,legend=NULL,xlab="landmark time points",ylab="Dynamic HR",ylim,...)
#' @param object an objects inheriting from class "LMf" or "Vs_LM"
#' @param covars Covariate names.
#' @param conf_int Indicator variable. It indicates whether confidence intervals are included.
#' @param IDlegend Annotations or descriptions of the plot.
#' @param legend Label of the plot.
#' @param xlab Label of the horizontal axis
#' @param ylab Label of the vertical axis.
#' @param ylim Scale of the vertical axis
#' @param ... extra aguments; currently none is used.
#' @export
#' @examples
#' library(DynamicModel)
#' data(renal)
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=TRUE)
#' plot_cHR(fit,c("age"),conf_int=TRUE,IDlegend="a.",legend="(Per 10 years)",xlab="landmark time points(s)",ylab="Dynamic HR (w=5)")
plot_cHR<-function(object,covars,conf_int=TRUE,IDlegend=NULL,legend=NULL,xlab="landmark time points",ylab="Dynamic HR",ylim,...){
  fm<-object$Model
  bet<-fm$coefficients
  func_covars <- object$data$func_covars
  if (conf_int) sig <- stats::vcov(fm)
  sl<-object$tw$sl
  w<-object$tw$w
  t <- seq(min(sl),max(sl) , by = 0.1)
  for (i in 1:length(covars)){
    idx <- startsWith(names(bet), covars[i])
    bet_var <- bet[idx]
    coef <- sapply(t/max(sl)-min(sl), function(x) { # eval HR over times x in t
      sum(sapply(1:length(bet_var), function(j){
        var = bet_var[j]
        name = names(bet_var)[j]
        if (name == covars[i]) { return(var) }
        else {
          idx <- as.numeric(sub(".*\\D+", "\\1", name))
          return(func_covars[[idx]](x) * var)
        }
      })) # bet0 + bet1*x + bet2*x^2 + ...
    })
    find_se_log <- function(t, coefs, covar, func_covars){
      if (!requireNamespace("msm", quietly = TRUE)) {
        stop("Package \"msm\" must be installed to use function find_se_log", call. = FALSE)}
      form <- "x1"
      if (length(coefs) > 1){
        for (i in 2:length(coefs)){
          form <- paste(form, sprintf("%s * %f", paste("x",i,sep=""), func_covars[[i-1]](t)),sep=" + ")
        }
      }
      form <- paste("~", form)
      se <- msm::deltamethod(stats::as.formula(form), coefs, covar)
      return(se)
    }
    se <- sapply(t/max(sl)-min(sl), find_se_log, bet_var, sig[idx, idx], func_covars)
    if (missing(ylim)) set_ylim <- TRUE
    else set_ylim <- FALSE
    vari <- exp(coef)
    se<-se
    lower<-exp(coef-1.96*se)
    upper<-exp(coef+1.96*se)
    if(set_ylim) ylim<-c(min(lower),max(upper))
    plot(t,vari,type="l",lwd=6,xlim=c(min(sl),max(sl)),
         ylim=ylim,bty="l",xlab="",ylab="",cex.axis=1.5,...)
    #title(list(covars[i],font=7,cex=2.5))
    mtext(IDlegend, side = 3, line = 1, adj = 0, cex = 1.5,font=7)
    title(main=list(covars[i],font=7,cex=1.5),line=1)
    title(xlab=xlab,font.lab=7,cex.lab=1.5,line=2.6)
    title(ylab=ylab,font.lab=7,cex.lab=1.5,line=2.6)
    #title(list(seq,font=7,cex=2),adj=0,line=1)
    lines(t,lower,type="l",lty=2,lwd=6)
    lines(t,upper,type="l",lty=2,lwd=6)
    if (is.null(legend)==FALSE)
      legend("topright", legend = legend,lty=1,lwd=6, bty = "n")
    abline(h=1,lwd=1,lty=3,col=2)

  }
}
