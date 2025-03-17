#' Title
#' @usage predictive_JM(jointFit,data,method,formulas,random_effect_formulas,form_splines,surv_formula,seq_len,w,n,K,seed,r)
#' @param jointFit an objects inheriting from class "jm"
#' @param data Dataset for fitting the joint model (JM).
#' @param method Types of internal cross validation. There are three types of internal validation methods: random splitting:sample, K-fold cross-validation:CV, and bootstrap:boot.
#' @param formulas A list of functional forms for fixed effects.
#' @param random_effect_formulas A list of functional forms for random effects.
#' @param form_splines A list of functional forms for the dynamic effects of longitudinal covariates.
#' @param surv_formula Functional form of the survival submodel.
#' @param seq_len a numeric vector of future times to calculate predictions
#' @param w prediction window
#' @param n umber of repetitions for internal cross-validation
#' @param K The specific number of folds in the K-fold cross-validation method.
#' @param seed The seed number for sampling.
#' @param r The proportion for random splitting.
#' @export
#' @return A list of class LMf with components:
#' \describe{
#' \item {sl} landmark time points
#' \item {cindex_corr} C-index
#' \item {BS_corr} Brier score
#' \item {AUC_corr} AUC
#' }
#'
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
#' aa<-predictive_JM(jointFit1,renal,method="CV",formulas=list(GFR~weight+age+ns(yearse,3),proteinuria~age+ns(yearse,3)),random_effect_formulas=list(~yearse|id,~yearse|id),surv_formula=Surv(time,status)~age+sex+weight,form_splines=form_splines,seq_len=seq(0,10,by=1),w=5,n=1,K=5,seed=1)
predictive_JM<-function(jointFit,data,method,formulas,random_effect_formulas,form_splines,surv_formula,seq_len,w,n,K,seed,r){

  id<-jointFit$model_info$var_names$idVar
  data.id<-data[!duplicated(data[[id]]),]
  if (method=="sample"){
    auc<-brier<-cindex<-data.frame(array(NA,dim=c(n,length(seq_len))))
    nrun<-0
    for (i in 1:(5*n)){
      tryCatch({
        set.seed(seed+10*i)
        ID<-sample(data.id[[id]],nrow(data.id)*r)
        data_train<-data[which(data[[id]] %in% ID),]
        data_train.id<-data.id[which(data.id[[id]] %in% ID),]
        data_test<-data[-which(data[[id]] %in% ID),]
        data_test.id<-data.id[-which(data.id[[id]] %in% ID),]
        #environment(fit_multiple_lme)<- environment()
        fit<-fit_multiple_lme(data_train,formulas,random_effect_formulas)
        surv_formula<-as.formula(surv_formula)
        cox<-coxph(surv_formula,data=data_train.id)
        jointFit1 <- jm(cox, fit, time_var =jointFit$model_info$var_names$time_var,functional_forms = form_splines,seed=1,n_chains = 1L)
        for (j in seq_len){
          auc[i,j+1]<-tvAUC(jointFit1,data_test, Tstart = j, Dt =w)$auc
          cindex[i,j+1]<-tvC_index(jointFit1,data_test, Tstart = j, Dt =w)
          brier[i,j+1]<-tvBrier(jointFit1,data_test, Tstart = j, Dt =w)$Brier
        }
        nrun<-nrun+1
      }, error = function(e) {
        message("Error in iteration ", i, ": ", e$message)
        # 捕获到错误后不输出任何信息，直接跳到下一个循环
        return(NULL)
      })
      if (nrun == n) {
        message("Reached maximum number of runs: ", n)
        break  # 退出循环
      }
    }
    cindex_corr<-colMeans(cindex,na.rm=TRUE)
    AUC_corr<-colMeans(auc,na.rm=TRUE)
    BS_corr<-colMeans(brier,na.rm=TRUE)
  }

  if (method=="CV"){
    cindex_sm<-data.frame(array(NA,dim=c(n,length(seq_len))))
    BS_sm<-data.frame(array(NA,dim=c(n,length(seq_len))))
    AUC_sm<-data.frame(array(NA,dim=c(n,length(seq_len))))
    nrun<-0
    for(i in 1:(5*n)){
      tryCatch({
        set.seed(seed+10*i)
        folds<-caret::createFolds(unique(data[[id]]),K)
        cindex<-data.frame(array(NA,dim=c(K,length(seq_len))))
        BS<-data.frame(array(NA,dim=c(K,length(seq_len))))
        AUC<-data.frame(array(NA,dim=c(K,length(seq_len))))
        for(j in 1:K) {
          index<-unique(data[[id]])[folds[[j]]]
          train_data<-data[-which(data[[id]] %in% index),]
          train_data.id<-data.id[-which(data.id[[id]] %in% index),]
          test_data<-data[which(data[[id]] %in% index),]
          test.id_data<-data.id[which(data.id[[id]] %in% index),]
          fit<-fit_multiple_lme(train_data,formulas,random_effect_formulas)
          surv_formula<-as.formula(surv_formula)
          cox<-coxph(surv_formula,data=train_data.id,x=TRUE)
          jointFit1 <- jm(cox, fit, time_var =jointFit$model_info$var_names$time_var,seed=1,functional_forms = form_splines, n_chains = 1L)
          for(k in 1:length(seq_len)){
            cindex[j,k]<-tvC_index(jointFit1,test_data,Tstart =seq_len[k], Dt =w)
            BS[j,k]<-tvBrier(jointFit1,test_data, Tstart = seq_len[k], Dt=w)$Brier
            AUC[j,k]<-tvAUC(jointFit1,test_data, Tstart = seq_len[k], Dt =w)$auc
          }
        }
        cindex_sm[i,]<-colMeans(cindex,na.rm=TRUE)
        BS_sm[i,]<-colMeans(BS,na.rm=TRUE)
        AUC_sm[i,]<-colMeans(AUC,na.rm=TRUE)

        nrun<-nrun+1
      }, error = function(e) {
        message("Error in iteration ", i, ": ", e$message)
        # 捕获到错误后不输出任何信息，直接跳到下一个循环
        return(NULL)
      })
      if (nrun == n) {
        message("Reached maximum number of runs: ", n)
        break  # 退出循环
      }
    }
    cindex_corr<-colMeans(cindex_sm,na.rm=TRUE)
    BS_corr<-colMeans(BS_sm,na.rm=TRUE)
    AUC_corr<-colMeans(AUC_sm,na.rm=TRUE)
  }

  if (method=="boot"){
    cindex_app<-BS_app<-AUC_app<-c()
    for(k in seq_len){
      cindex_app[k+1]<-tvC_index(jointFit,data, Tstart = k, Dt =w)
      BS_app[k+1]<-tvBrier(jointFit,data, Tstart = k, Dt =w)$Brier
      AUC_app[k+1]<-tvAUC(jointFit1,data, Tstart = k, Dt =w)$auc
    }
    cindex_opt<-BS_opt<-AUC_opt<-data.frame(array(NA,dim=c(n,length(seq_len))))
    nrun<-0
    for (i in 1:(5*n)){
      tryCatch({
        set.seed(seed+10*i)
        train_index<-sample(unique(data[[id]]),length(unique(data[[id]])),replace=TRUE)
        train_data<-NULL
        train_data.id<-NULL
        for(j in seq_along(train_index)){
          t_data<-data[which(data[[id]]==train_index[j]),]
          t_data[[id]]<-rep(j,nrow(t_data))
          t_data.id<-data.id[which(data.id[[id]]==train_index[j]),]
          t_data.id[[id]]<-j
          train_data<-rbind(t_data,train_data)
          train_data.id<-rbind(t_data.id,train_data.id)
        }
        train_data<-train_data[order(train_data[[id]]),]
        train_data.id<-train_data.id[order(train_data.id[[id]]),]

        fit<-fit_multiple_lme( train_data,formulas,random_effect_formulas)
        surv_formula<-as.formula(surv_formula)
        cox<-coxph(surv_formula,data= train_data.id,x=TRUE)
        jointFit1 <- jm(cox, fit, time_var = jointFit$model_info$var_names$time_var,seed=1,functional_forms = form_splines, n_chains = 1L)

        for(k in seq_len){
          cindex_boot<-tvC_index(jointFit1,train_data, Tstart = k, Dt =w)
          BS_boot<-tvBrier(jointFit1,train_data, Tstart = k, Dt =w)$Brier
          AUC_boot<-tvAUC(jointFit1,train_data, Tstart = k, Dt =w)$auc

          cindex_orig<-tvC_index(jointFit1,data, Tstart = k, Dt =w)
          BS_orig<-tvBrier(jointFit1,data, Tstart = k, Dt =w)$Brier
          AUC_orig<-tvAUC(jointFit1,data, Tstart = k, Dt =w)$auc

          cindex_opt[i,k+1] <- cindex_boot - cindex_orig
          BS_opt[i,k+1] <- BS_boot - BS_orig
          AUC_opt[i,k+1] <- AUC_boot - AUC_orig
        }
        nrun<-nrun+1
      }, error = function(e) {
        message("Error in iteration ", i, ": ", e$message)
        # 捕获到错误后不输出任何信息，直接跳到下一个循环
        return(NULL)
      })
      if (nrun == n) {
        message("Reached maximum number of runs: ", n)
        break  # 退出循环
      }

    }
    cindex_corr<-cindex_app-colMeans(cindex_opt,na.rm=TRUE)
    BS_corr<-BS_app+colMeans(BS_opt,na.rm=TRUE)
    AUC_corr<-AUC_app-colMeans(AUC_opt,na.rm=TRUE)

  }
  return(list(time=seq_len,cindex_corr=cindex_corr,BS_corr=BS_corr,AUC_corr=AUC_corr))
}
tvC_index<-function(object,newdata,Tstart,Thoriz = NULL, Dt = NULL){
  Tstart <- Tstart + 1e-06
  Thoriz <- Tstart + Dt
  id_var <- object$model_info$var_names$idVar
  time_var <- object$model_info$var_names$time_var
  Time_var <- object$model_info$var_names$Time_var
  event_var <- object$model_info$var_names$event_var
  tt<-newdata[[Time_var]]
  newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]
  id <- newdata[[id_var]]
  id <- match(id, unique(id))
  tt <- ave(tt, id, FUN = function (t) rep(tail(t, 1L) > Tstart, length(t)))
  newdata <- newdata[as.logical(tt), ]
  newdata <- newdata[newdata[[time_var]] <= Tstart, ]
  newdata2 <- newdata
  newdata2[[Time_var]] <- Tstart
  newdata2[[event_var]] <- 0
  predd<-predict(object,newdata=newdata2, process = "event",times=Thoriz)
  risk<-predd$pred[seq(2, length(predd$pred), by = 2)]
  names(risk)<-predd$id[seq(2, length(predd$id), by = 2)]

  id <- newdata[[id_var]]
  Time <- newdata[[Time_var]]
  event <- newdata[[event_var]]
  f <- factor(id, levels = unique(id))
  Time <- tapply(Time, f, tail, 1L)
  event <- tapply(event, f, tail, 1L)
  names(Time) <- names(event) <- as.character(unique(id))


  nt<-length(Time)
  ord<-order(Time,-event)
  time<-Time[ord]
  status<-event[ord]
  risk<-risk[ord]
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
fit_multiple_lme <-function(data,formulas,random_effect_formulas){
  models<-list()
  # 检查公式和随机效应公式的长度是否一致
  if (length(formulas) != length(random_effect_formulas)) {
    stop("The number of formulas must be the same as the number of random effect formulas.")
  }
  for (i in 1:length(formulas)){
    formula <- formulas[[i]]  # 获取当前公式
    random_effect_formula <- random_effect_formulas[[i]]  # 获取对应的随机效应公式
    model_name <- paste("model", i, sep = "_")  # 创建模型名称

    # 拟合lme模型
    models[[model_name]] <- lme(formula, random = random_effect_formula, data = data,control = list(opt = "optim"))
  }
  return(models)
}
