#' Title Dynamic Effects Plot of Covariates in a Joint Model-Based Dynamic Prediction Model
#' @usage plot_dynamicHR(object,variable,form_splines,IDnames=NULL,xlab = "Follow-up Time (years)",ylab = "Dynamic Hazard Ratio",ylim =NULL, col = c("black", "black", "black"), lty = c(1, 2, 2),lwd=6,...)
#' @param object an objects inheriting from class "jm"
#' @param variable Name of longitudinal covariates
#' @param form_splines Form of the dynamic effects function
#' @param tt Follow-up time with the maximum of interest
#' @param IDnames Label of the plot.
#' @param xlab Label of the horizontal axis
#' @param ylab Label of the vertical axis.
#' @param ylim Scale of the vertical axis
#' @param col Line Color Options
#' @param lty Line Type Options
#' @param lwd Line Thickness Options
#' @param ... extra aguments; currently none is used.
#' @export
#' @examples
#' library(JMbayes2)
#' library(DynamicModel)
#' data(renal)
#' renal.id<-renal[!duplicated(renal$id),]
#' fit1<-lme(GFR~age+weight+ns(yearse,3), random =~yearse | id,data = renal,control = list(opt = "optim"))
#' fit2<-lme(proteinuria~age+ns(yearse,3), random = ~yearse | id,  data =renal,control = list(opt = "optim"))
#' cox<-coxph(Surv(time,status)~age+sex+weight,data=renal.id,x=TRUE)
#' form_splines <- list("GFR"=~ value(GFR)* (yearse +I(yearse^2)),
#' "proteinuria"=~ value(proteinuria) * (yearse+I(yearse^2)))
#' jointFit1 <- jm(cox, list(fit1,fit2), time_var = "yearse",seed=1,n_chains = 1L, functional_forms = form_splines)
#' plot_dynamiceffect(jointFit1,"GFR",form_splines="quadratic",tt=15,IDnames="b.",ylim=c(0,0.5),xlim=c(0,15),cex.axis=1.5,cex.lab=1.5)
plot_dynamicHR<-function(object,variable,form_splines,tt,IDnames=NULL,xlab = "Follow-up Time (years)",ylab = "Dynamic Hazard Ratio",ylim =NULL, col = c("black", "black", "black"), lty = c(1, 2, 2),lwd=6,...){
  x_times <- seq(0.001, tt, length = 501)
  if (form_splines=="quadratic")
    x<-cbind(1,x_times,x_times^2)
  if(form_splines=="spline")
    x<-cbind(1,ns(x_times,3))
  x1<-as.matrix(object$mcmc$alphas[1])
  names<-colnames(x1)
  matches <- grep(variable, names, value = TRUE)
  mcmc_alpha<-x1[,matches]
  log_hr<-x %*% t(mcmc_alpha)
  log_hr_mean <- rowMeans(log_hr)
  log_hr_low <- apply(log_hr, 1, quantile, probs = 0.025)
  log_hr_upp <- apply(log_hr, 1, quantile, probs = 0.975)
  if(is.null(ylim)) ylim<-c(min(exp(log_hr_low)),max(exp(log_hr_upp)))
  matplot(x_times, cbind(exp(log_hr_mean), exp(log_hr_low), exp(log_hr_upp)),
          type = "l", col = col, lty =lty, lwd =lwd,xlab="",ylab="",
          ylim =ylim,main="",xaxt="n",yaxt="n",bty="n",...)
  axis(1,las=1,pos=0,tcl=-0.4,padj=-0.3,lwd=1.6,...)
  axis(2,las=1,pos=ylim[1],tcl=-0.4,hadj=0.9,lwd=1.6,...)
  title(xlab=xlab,font.lab=7,line=2.6,...)
  title(ylab=ylab,font.lab=7,line=2.6,...)
  abline(h = 1, lty = 2)
  title(main=list(variable,font=7,cex=1.5),line=1)
  mtext(IDnames, side = 3, line =1, adj =0, cex =1.5,font=7)
}
