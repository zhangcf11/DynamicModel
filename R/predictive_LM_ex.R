#' Title External Validation Function for Landmarking-Based Dynamic Prediction Models
#'
#' @param object an objects inheriting from class "LMf" or "Vs_LM".
#' @param ex_data An external validation dataset with variable names identical to those in the model-building dataset.
#'
#' @return A list of class Predictive_LM_ex with components:
#' \describe{
#' \item {time} landmark time points
#' \item {cindex_corr} C-index
#' \item {BS_corr} Brier score
#' \item {AUC_corr} AUC
#' }
#' @export
Predictive_LM_ex<-function(object,ex_data){
  if (!require("survivalROC")) {
    install.packages("survivalROC")
    library(survivalROC)
  }
  if (!require("dynpred")) {
    install.packages("dynpred")
    library(dynpred)
  }
  TSet<-object$data
  Model<-object$Model
  sl<-object$tw$sl
  nsl<-length(sl)
  id<-object$id
  w<-object$tw$w
  time<-object$time
  status<-object$status
  rtime<-object$rtime
  cov<-object$cov
  Vset<-NULL
  for (j in seq(along=sl)){
    LM<-dynpred::cutLM(ex_data,outcome=list(time=time,status=status),
                       LM=sl[j],horizon=sl[j]+w,covs=list(fixed=cov$fixed,varying=cov$vary),format="long",id=id,rtime=rtime,right=F)
    Vset<-rbind(Vset,LM)
  }
  Vset<-Vset[order(Vset[[id]]),]
  cindex<-rep(NA,nsl)
  score<-rep(NA,nsl)
  auc<-rep(NA,nsl)
  for (i in 1:nsl){
    Vdata<-VSet[VSet$LM==sl[i],]
    cindex[i]<-cal_cindex(model=Model,data=Vdata,time,status)
    score[i]<-cal_brierscore(model=Model,TSet,Vdata,width=w, tout=sl[i],time,status)
    auc[i]<-cal_auc(model=Model,data=Vdata,pred.t=sl[i]+w,time,status)$auc
    #all3<-cbind(sl[i],cal_auc(model=Model3,data=Vdata,pred.t=sl[i]+w)$all_max)
    #all2<-rbind(all2,all3)
  }
  return(list(time=sl,cindex=cindex,brierscore=score,auc=auc))
}
