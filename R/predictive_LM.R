#' Title Internal cross-validation function for landmark-based dynamic prediction models.
#' @usage predictive_LM(object,method,n,K,seed,r)
#' @param object an objects inheriting from class "LMf" or "Vs_LM"
#' @param method Types of internal cross validation.There are three types of internal validation methods: random splitting:sample, K-fold cross-validation:CV, and bootstrap:boot.
#' @param n Number of repetitions for internal cross-validation
#' @param K The specific number of folds in the K-fold cross-validation method.
#' @param seed The seed number for sampling.
#' @param r The proportion for random splitting.
#'
#' @return A list of class predictive_LM with components:
#' \describe{
#' \item {time} landmark time points
#' \item {cindex_corr} C-index
#' \item {BS_corr} Brier score
#' \item {AUC_corr} AUC
#' }
#' @export
#' @examples
#' library(DynamicModel)
#' data(renal)
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=TRUE)
#' LM_pre<-predictive_LM(fit,method="CV",n=1,seed=1,K=5)
predictive_LM<-function(object,method,n,K,seed,r){
  if (!require("caret")) {
    install.packages("caret")
    library(caret)
  }
  if (!require("survivalROC")) {
    install.packages("survivalROC")
    library(survivalROC)
  }
  fm<-object$Model
  id<-object$id
  data<-object$data$data
  sl<-object$tw$sl
  w<-object$tw$w
  time<-object$time
  status<-object$status
  if (method=="CV"){
    cindex_sm<-data.frame(array(NA,dim=c(n,length(sl))))
    BS_sm<-data.frame(array(NA,dim=c(n,length(sl))))
    AUC_sm<-data.frame(array(NA,dim=c(n,length(sl))))
    nrun<-0
    for(i in 1:(5*n)){
      tryCatch({
        set.seed(seed+10*i )
        folds<-caret::createFolds(unique(data[[id]]),K)
        cindex<-data.frame(array(NA,dim=c(K,length(sl))))
        BS<-data.frame(array(NA,dim=c(K,length(sl))))
        AUC<-data.frame(array(NA,dim=c(K,length(sl))))
        for(j in 1:K) {
          index<-unique(data[[id]])[folds[[j]]]
          train_data<-data[-which(data[[id]] %in% index),]
          test_data<-data[which(data[[id]] %in% index),]
          fit<-coxph(fm$formula,data=train_data)
          for(k in 1:length(sl)){
            test_data1<-test_data[test_data$LM==sl[k],]
            cindex[j,k]<-cal_cindex(model=fit,data=test_data1,time,status)
            BS[j,k]<-cal_brierscore(model=fit,Tdata=train_data,Vdata=test_data1,width=w,tout=sl[k],time,status)
            AUC[j,k]<-cal_auc(model=fit,data=test_data1,pred.t=sl[k]+w,time,status)
          }
        }
        cindex_sm[i,]<-colMeans(cindex,na.rm=TRUE)
        BS_sm[i,]<-colMeans(BS,na.rm=TRUE)
        AUC_sm[i,]<-colMeans(AUC,na.rm=TRUE)
        nrun<-nrun+1
      }, error = function(e) {
        # 捕获到错误后不输出任何信息，直接跳到下一个循环
        return(NULL)
      })
      if (nrun== n) {
        message("Reached maximum number of runs: ", nrun)
        break  # 退出循环
      }
    }
    cindex_corr<-colMeans(cindex_sm,na.rm=TRUE)
    BS_corr<-colMeans(BS_sm,na.rm=TRUE)
    AUC_corr<-colMeans(AUC_sm,na.rm=TRUE)
  }
  if (method=="boot"){
    cindex_app<-BS_app<-AUC_app<-c()
    for(k in 1:length(sl)){
      data1<-data[data$LM==sl[k],]
      cindex_app[k]<-cal_cindex(model=fm,data=data1,time,status)
      BS_app[k]<-cal_brierscore(model=fm,Tdata=data1,Vdata=data1,width=w,tout=sl[k],time,status)
      AUC_app[k]<-cal_auc(model=fm,data=data1,pred.t=sl[k]+w,time,status)
    }
    cindex_opt<-BS_opt<-AUC_opt<-data.frame(array(NA,dim=c(n,length(sl))))
    nrun<-0
    for (i in 1:(5*n)){
      tryCatch({
        set.seed(seed+10*i)
        train_index<-sample(unique(data[[id]]),length(unique(data[[id]])),replace=TRUE)
        train_data<-NULL
        for(j in seq_along(train_index)){
          t_data<-data[which(data[[id]]==train_index[j]),]
          train_data<-rbind(t_data,train_data)
        }
        train_data<-train_data[order(train_data[[id]]),]
        fit<-coxph(fm$formula,data=train_data)
        for(k in 1:length(sl)){
          train_data1<-train_data[train_data$LM==sl[k],]
          cindex_boot<-cal_cindex(model=fit,data=train_data1,time,status)
          BS_boot<-cal_brierscore(model=fit,Tdata=train_data1,Vdata=train_data1,width=w,tout=sl[k],time,status)
          AUC_boot<-cal_auc(model=fit,data=train_data1,pred.t=sl[k]+w,time,status)

          data1<-data[data$LM==sl[k],]
          cindex_orig<-cal_cindex(model=fit,data=data1,time,status)
          BS_orig<-cal_brierscore(model=fit,Tdata=train_data1,Vdata=data1,width=w,tout=sl[k],time,status)
          AUC_orig<-cal_auc(model=fit,data=data1,pred.t=sl[k]+w,time,status)

          cindex_opt[i,k] <- cindex_boot - cindex_orig
          BS_opt[i,k] <- BS_boot - BS_orig
          AUC_opt[i,k] <- AUC_boot - AUC_orig
          nrun<-nrun+1
        }
      }, error = function(e) {
        # 捕获到错误后不输出任何信息，直接跳到下一个循环
        return(NULL)
      })
      if (nrun == n) {
        message("Reached maximum number of runs: ", nrun)
        break  # 退出循环
      }

    }
    cindex_corr<-cindex_app-colMeans(cindex_opt,na.rm=TRUE)
    BS_corr<-BS_app+colMeans(BS_opt,na.rm=TRUE)
    AUC_corr<-AUC_app-colMeans(AUC_opt,na.rm=TRUE)
  }
  if (method=="sample"){
    cindex<-BS<-AUC<-data.frame(array(NA,dim=c(n,length(sl))))
    nrun<-0
    for(i in 1:(5*n)){
      tryCatch({
        set.seed(seed+10*i)
        index<-sample(unique(data[[id]]),length(unique(data[[id]]))*r)
        train_data<-data[which(data[[id]]%in% index),]
        test_data<-data[-which(data[[id]]%in% index),]
        fit<-coxph(fm$formula,data=train_data)
        for(k in 1:length(sl)){
          test_data1<-test_data[which(test_data$LM==sl[k]),]
          cindex[i,k]<-cal_cindex(model=fit,data=test_data1,time,status)
          BS[i,k]<-cal_brierscore(model=fit,Tdata=train_data,Vdata=test_data1,width=w,tout=sl[k],time,status)
          AUC[i,k]<-cal_auc(model=fit,data=test_data1,pred.t=sl[k]+w,time,status)
        }
        nrun<-nrun+1
      }, error = function(e) {
        # 捕获到错误后不输出任何信息，直接跳到下一个循环
        return(NULL)
      })
      if (nrun == n) {
        message("Reached maximum number of runs: ", nrun)
        break  # 退出循环
      }
    }
    cindex_corr<-colMeans(cindex,na.rm=TRUE)
    BS_corr<-colMeans(BS,na.rm=TRUE)
    AUC_corr<-colMeans(AUC,na.rm=TRUE)
  }

  return(list(time=sl,cindex_corr=cindex_corr,BS_corr=BS_corr,AUC_corr=AUC_corr))

}

cal_cindex<-function(model,data,time,status){
  nt<-length(data[[time]])
  ord<-order(data[[time]],-data[[status]])
  time<-data[[time]][ord]
  status<-data[[status]][ord]
  risk<-predict(model,newdata=data,type="lp")[ord]
  wh<-which(status==1)
  wh<-wh[which(wh<=nt-1)]
  total<-con<-0
  for (m in wh){
    for (n in ((m + 1):nt)) {
      if (time[n] > time[m]) {
        total <- total + 1
        if (risk[n] < risk[m])
          con <- con + 1
        if (risk[n] == risk[m])
          con <- con + 0.5
      }
    }
  }
  return(con/total)
}

cal_brierscore<-function(model,Tdata,Vdata,width,tout,Time,Status)
{
  ord<-order(Vdata[[Time]],-Vdata[[Status]])
  time<-Vdata[[Time]][ord]
  status<-Vdata[[Status]][ord]
  risk<-predict(model,newdata=Vdata,type="lp")[ord]
  cox1<-coxph(Surv(Vdata$LM[ord],Vdata[[Time]][ord],
                   Vdata[[Status]][ord])~risk)
  if (sum(risk^2)==0)
  {sf<-survfit(cox1,newdata=data.frame(risk=risk),
               type="kalbfl")
  }else sf<-survfit(cox1,newdata=data.frame(risk=risk))
  tt<-sf$time
  survmat<-sf$surv
  if (tt[1]>0){
    tsurv<-c(0,tt)
    survmat<-rbind(rep(1,nrow(Vdata)),survmat)
  } else tsurv<-tt
  coxcens<-coxph(Surv(Tdata$LM,Tdata[[Time]],1-Tdata[[Status]])~1)
  xcens<-predict(coxcens,newdata=Vdata,type="lp")[ord]
  coxcens<-coxph(Surv(Vdata$LM,Vdata[[Time]],1-Vdata[[Status]])~xcens)
  if (sum(xcens^2)==0)
  {sfcens<-survfit(coxcens,newdata=data.frame(xcens=xcens),
                   type="kalbfl")
  }else {sfcens<-survfit(coxcens,newdata=data.frame(xcens=xcens))}
  tcens<-sfcens$time
  censmat<-sfcens$surv
  if (tcens[1]>0) {
    tcens<-c(0,tcens)
    censmat<-rbind(rep(1,nrow(Vdata)),censmat)
  }
  res<-pew(time,status,tsurv,survmat,tcens,censmat,width,"Brier",tout)
  return(res[1,2])
}

cal_auc<-function(model,data,pred.t,time,status){
  risk<-predict(model,newdata=data)
  auc<-survivalROC::survivalROC(Stime=data[[time]],status=data[[status]],marker=risk,
                   predict.time = pred.t,
                   method="KM")

  auc0<-auc$AUC
  return(auc0)
}

