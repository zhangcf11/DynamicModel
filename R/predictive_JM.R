#' Predictive Performance Assessment for Joint Models
#'
#' This function evaluates the predictive performance of joint models using
#' different validation methods including data splitting, K-fold cross-validation,
#' and bootstrap validation.
#' @usage predictive_JM(jointFit,data,method,formulas,random_effect_formulas,form_splines,surv_formula,seq_len,w,n,K,seed,r)
#' @param jointFit A fitted joint model object from the `jm` function
#' @param data A data.frame containing the longitudinal and survival data
#' @param method Character string specifying validation method: "split", "K-fold", or "boot"
#' @param formulas A list of functional forms for fixed effects of longitudinal submodel.
#' @param random_effect_formulas A list of functional forms for random effects of longitudinal submodel.
#' @param form_splines A list of functional forms for the dynamic effects of longitudinal covariates，with default NULL.
#' @param surv_formula Functional form of the survival submodel.
#' @param seq_len Numeric vector of time points for prediction evaluation
#' @param w prediction window
#' @param n Number of repetitions for internal cross-validation. For random splitting, it refers to the number of random splitting iterations; for K-fold cross-validation, it refers to the number of outer loop iterations; for bootstrap, it refers to the number of resampling iterations.
#' @param K Integer, number of folds for K-fold cross-validation (only used when method="K-fold")
#' @param seed The seed number for sampling.
#' @param r Numeric, proportion of data to use for training (only used when method="split")
#' @param n_group number of groups (or points) for the calibration curve.Default is 10.
#' @importFrom caret createFolds
#' @importFrom survival coxph
#' @importFrom JMbayes2 jm tvAUC tvC_index tvBrier
#' @export
#' @return A list containing:
#' \itemize{
#'   \item{time: The sequence of time points used for evaluation}
#'   \item{cindex_corr: Corrected concordance index with confidence intervals}
#'   \item{BS_corr: Corrected Brier score with confidence intervals}
#'   \item{AUC_corr: Corrected AUC with confidence intervals}
#'   \item{calibration_slpoe_with_ci: Calibration slope with confidence intervals}
#'   \item{all_cal_autual: Complete calibration results}
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
#' JM_pre<-predictive_JM(jointFit1,renal,method="K-fold",formulas=list(GFR~weight+age+ns(yearse,3),proteinuria~age+ns(yearse,3)),random_effect_formulas=list(~yearse|id,~yearse|id),surv_formula=Surv(time,status)~age+sex+weight,form_splines=form_splines,seq_len=seq(0,10,by=1),w=5,n=1,K=5,seed=1)
predictive_JM <- function(jointFit, data, method, formulas, random_effect_formulas,
                          form_splines = NULL, surv_formula, seq_len, w, n, K, seed, r,
                          n_group = 10) {

  # Input validation
  if (!method %in% c("split", "K-fold", "boot")) {
    stop("method must be one of: 'split', 'K-fold', 'boot'")
  }

  if (!require("caret", quietly = TRUE)) {
    install.packages("caret")
    library(caret)
  }

  # Extract key variables from joint model
  id <- jointFit$model_info$var_names$idVar
  data.id <- data[!duplicated(data[[id]]), ]
  Time_var <- jointFit$model_info$var_names$Time_var
  event_var <- jointFit$model_info$var_names$event_var

  # Initialize result containers
  cindex_corr <- BS_corr <- AUC_corr <- calibration_slpoe_with_ci <- all_cal_autual <- NULL

  if (method == "split") {
    # Data splitting method implementation
    auc <- brier <- cindex <- data.frame(array(NA, dim = c((5 * n), length(seq_len))))
    nrun <- 0
    cal_pred <- data.frame()

    for (i in 1:(5 * n)) {
      tryCatch({
        set.seed(seed + 10 * i)
        ID <- sample(data.id[[id]], nrow(data.id) * r)
        data_train <- data[which(data[[id]] %in% ID), ]
        data_train.id <- data.id[which(data.id[[id]] %in% ID), ]
        data_test <- data[-which(data[[id]] %in% ID), ]
        data_test.id <- data.id[-which(data.id[[id]] %in% ID), ]

        # Fit longitudinal models
        fit <- fit_multiple_lme(data_train, formulas, random_effect_formulas)
        surv_formula <- as.formula(surv_formula)
        cox <- coxph(surv_formula, data = data_train.id)

        # Fit joint model
        jointFit1 <- jm(cox, fit, time_var = jointFit$model_info$var_names$time_var,
                        functional_forms = form_splines, seed = 1, n_chains = 1L)

        all_pred <- data.frame()
        for (j in seq_along(seq_len)) {
          time_point <- seq_len[j]
          auc[i, j] <- tvAUC(jointFit1, data_test, Tstart = time_point, Dt = w)$auc
          cindex[i, j] <- tvC_index(jointFit1, data_test, Tstart = time_point, Dt = w)
          brier[i, j] <- tvBrier(jointFit1, data_test, Tstart = time_point, Dt = w)$Brier

          all_pred1 <- individual_JM_predict(jointFit1, data_test, Tstart = time_point, Dt = w)
          all_surv <- cbind(all_pred1, rep(i, nrow(all_pred1)), rep(time_point, nrow(all_pred1)))
          colnames(all_surv)[c(4, 5)] <- c("num", "time_points")
          all_pred <- rbind(all_pred, all_surv)
        }
        nrun <- nrun + 1
        cal_pred <- rbind(cal_pred, all_pred)
      }, error = function(e) {
        message("Error in iteration ", i, ": ", e$message)
        return(NULL)
      })

      if (nrun == n) {
        message("Reached maximum number of successful runs: ", n)
        break
      }
    }

    # Calculate calibration metrics
    all_cal_autual <- data.frame()
    for (timepoints in unique(cal_pred$time_points)) {
      data_time <- cal_pred[which(cal_pred$time_points == timepoints), ]
      cal_re <- calibration_re_JM(Tstart = timepoints, Dt = w, all_pred = data_time,
                                  n_groups = n_group)$cal_results
      all_cal_autual <- rbind(all_cal_autual, cal_re)
    }

    calibration_slpoe_with_ci <- calculate_calibration_slope(all_cal_autual)
    cindex_corr <- calculate_metric_ci(metric_matrix = cindex, time_points = seq_len,
                                       method = "quantile")
    BS_corr <- calculate_metric_ci(metric_matrix = brier, time_points = seq_len,
                                   method = "quantile")
    AUC_corr <- calculate_metric_ci(metric_matrix = auc, time_points = seq_len,
                                    method = "quantile")

  } else if (method == "K-fold") {
    # K-fold cross-validation implementation
    cindex_sm <- BS_sm <- AUC_sm <- list()
    cal_pred <- list()
    nrun <- 0

    for (i in 1:(5 * n)) {
      tryCatch({
        set.seed(seed + 10 * i)
        folds <- caret::createFolds(unique(data[[id]]), K)
        cindex <- BS <- AUC <- data.frame(array(NA, dim = c(K, length(seq_len))))
        all_pred <- data.frame()

        for (j in 1:K) {
          index <- unique(data[[id]])[folds[[j]]]
          train_data <- data[-which(data[[id]] %in% index), ]
          train_data.id <- data.id[-which(data.id[[id]] %in% index), ]
          test_data <- data[which(data[[id]] %in% index), ]

          # Fit models
          fit <- fit_multiple_lme(train_data, formulas, random_effect_formulas)
          surv_formula <- as.formula(surv_formula)
          cox <- coxph(surv_formula, data = train_data.id, x = TRUE)
          jointFit1 <- jm(cox, fit, time_var = jointFit$model_info$var_names$time_var,
                          seed = 1, functional_forms = form_splines, n_chains = 1L)

          for (k in seq_along(seq_len)) {
            time_point <- seq_len[k]
            cindex[j, k] <- tvC_index(jointFit1, test_data, Tstart = time_point, Dt = w)
            BS[j, k] <- tvBrier(jointFit1, test_data, Tstart = time_point, Dt = w)$Brier
            AUC[j, k] <- tvAUC(jointFit1, test_data, Tstart = time_point, Dt = w)$auc

            all_pred1 <- individual_JM_predict(jointFit1, test_data, Tstart = time_point, Dt = w)
            idd<-rownames(all_pred1)
            all_surv <- cbind(all_pred1, rep(i, nrow(all_pred1)), rep(time_point, nrow(all_pred1)),idd)
            colnames(all_surv)[c(4, 5,6)] <- c("num", "time_points","id")
            all_pred <- rbind(all_pred, all_surv)
          }
        }

        cindex_sm[[i]] <- cindex
        BS_sm[[i]] <- BS
        AUC_sm[[i]] <- AUC
        cal_pred[[i]] <- all_pred
        nrun <- nrun + 1

      }, error = function(e) {
        message("Error in iteration ", i, ": ", e$message)
        return(NULL)
      })

      if (nrun == n) {
        message("Reached maximum number of successful runs: ", n)
        break
      }
    }

    # Process calibration results
    cal_array <- lapply(cal_pred, function(mat) {
      order_idx <- order(mat[, "id"], mat[, "time_points"])
      mat_sorted <- mat[order_idx, ]
      return(mat_sorted)
    })

    cal_pred_array <- matrix(NA, nrow = nrow(cal_array[[1]]), ncol = length(cal_array))
    for (i in seq_along(cal_array)) {
      cal_pred_array[, i] <- cal_array[[i]][, "preds"]
    }

    cal_pred_mean <- rowMeans(cal_pred_array)
    cal_pred_corr <- cbind(cal_array[[1]][, "Time"], cal_array[[1]][, "event"],
                           cal_array[[1]][, "time_points"], cal_pred_mean)
    cal_pred_corr <- as.data.frame(cal_pred_corr)
    names(cal_pred_corr) <- c("auctual_time", "auctual_status", "time_points", "preds")

    all_cal_autual <- data.frame()
    for (timepoints in unique(cal_pred_corr$time_points)) {
      data_time <- cal_pred_corr[which(cal_pred_corr$time_points == timepoints), ]
      data_time$auctual_status <- ifelse(data_time$auctual_time < data_time$time_points + w &
                                           data_time$auctual_status, 1, 0)
      cal_re <- calibration_re_JM(Tstart = timepoints, Dt = w, all_pred = data_time,
                                  n_groups = n_group)$cal_results
      all_cal_autual <- rbind(all_cal_autual, cal_re)
    }

    calibration_slpoe_with_ci <- calculate_calibration_slope(all_cal_autual)

    # Calculate performance metrics with confidence intervals
    cindex_array <- mean_ignore_na(cindex_sm)
    BS_array <- mean_ignore_na(BS_sm)
    AUC_array <- mean_ignore_na(AUC_sm)

    cindex_corr <- do.call(rbind, apply(cindex_array, 2, calculate_cv_ci_t))
    BS_corr <- do.call(rbind, apply(BS_array, 2, calculate_cv_ci_t))
    AUC_corr <- do.call(rbind, apply(AUC_array, 2, calculate_cv_ci_t))

  } else if (method == "boot") {
    # Bootstrap validation implementation
    cindex_app <- BS_app <- AUC_app <- numeric(length(seq_len))
    all_calibration <- data.frame()
    cuts <- list()

    # Apparent performance
    for (k in seq_along(seq_len)) {
      time_point <- seq_len[k]
      cindex_app[k] <- tvC_index(jointFit, data, Tstart = time_point, Dt = w)
      BS_app[k] <- tvBrier(jointFit, data, Tstart = time_point, Dt = w)$Brier
      AUC_app[k] <- tvAUC(jointFit, data, Tstart = time_point, Dt = w)$auc

      cal_re <- calibration_re_JM(object = jointFit, data = data, Tstart = time_point,
                                  Dt = w, n_groups = n_group)
      all_calibration <- rbind(all_calibration, cal_re$cal_results)
      cuts[[k]] <- cal_re$cuts
    }

    # Bootstrap correction
    cindex_opt <- BS_opt <- AUC_opt <- cindex_cor <- BS_cor <- AUC_cor <-
      data.frame(array(NA, dim = c((5 * n), length(seq_len))))
    cal_opt <- data.frame()
    nrun <- 0

    for (i in 1:(5 * n)) {
      tryCatch({
        set.seed(seed + 10 * i)
        train_index <- sample(unique(data[[id]]), length(unique(data[[id]])), replace = TRUE)
        train_data <- train_data.id <- NULL

        for (j in seq_along(train_index)) {
          t_data <- data[which(data[[id]] == train_index[j]), ]
          t_data[[id]] <- rep(j, nrow(t_data))
          t_data.id <- data.id[which(data.id[[id]] == train_index[j]), ]
          t_data.id[[id]] <- j
          train_data <- rbind(train_data, t_data)
          train_data.id <- rbind(train_data.id, t_data.id)
        }

        train_data <- train_data[order(train_data[[id]]), ]
        train_data.id <- train_data.id[order(train_data.id[[id]]), ]

        # Fit bootstrap model
        fit <- fit_multiple_lme(train_data, formulas, random_effect_formulas)
        surv_formula <- as.formula(surv_formula)
        cox <- coxph(surv_formula, data = train_data.id, x = TRUE)
        jointFit1 <- jm(cox, fit, time_var = jointFit$model_info$var_names$time_var,
                        seed = 1, functional_forms = form_splines, n_chains = 1L)

        cal_re_opt <- data.frame()
        for (k in seq_along(seq_len)) {
          time_point <- seq_len[k]

          # Bootstrap performance
          cindex_boot <- tvC_index(jointFit1, train_data, Tstart = time_point, Dt = w)
          BS_boot <- tvBrier(jointFit1, train_data, Tstart = time_point, Dt = w)$Brier
          AUC_boot <- tvAUC(jointFit1, train_data, Tstart = time_point, Dt = w)$auc

          # Original data performance
          cindex_orig <- tvC_index(jointFit1, data, Tstart = time_point, Dt = w)
          BS_orig <- tvBrier(jointFit1, data, Tstart = time_point, Dt = w)$Brier
          AUC_orig <- tvAUC(jointFit1, data, Tstart = time_point, Dt = w)$auc

          # Calculate optimism
          cindex_opt[i, k] <- cindex_boot - cindex_orig
          cindex_cor[i, k] <- cindex_app[k] - cindex_opt[i, k]

          BS_opt[i, k] <- BS_boot - BS_orig
          BS_cor[i, k] <- BS_app[k] - BS_opt[i, k]

          AUC_opt[i, k] <- AUC_boot - AUC_orig
          AUC_cor[i, k] <- AUC_app[k] - AUC_opt[i, k]

          # Calibration optimism
          cal_re_boot <- calibration_re_JM(object = jointFit1, data = train_data,
                                           Tstart = time_point, Dt = w, n_groups = n_group,
                                           probs = cuts[[k]])$cal_results
          cal_re_orig <- calibration_re_JM(object = jointFit1, data = data,
                                           Tstart = time_point, Dt = w, n_groups = n_group,
                                           probs = cuts[[k]])$cal_results

          ideal_groups <- levels(cut(0, breaks = cuts[[k]]))
          cal_opt1 <- data.frame(
            num = i,
            time_point = time_point,
            pred_group = ideal_groups,
            stringsAsFactors = FALSE
          )

          get_value <- function(data, group, column) {
            idx <- which(as.character(data$pred_group) == as.character(group))
            if (length(idx) > 0) data[idx, column] else 0
          }

          cal_re_boot1 <- sapply(ideal_groups, function(g) get_value(cal_re_boot, g, "actual_survival"))
          cal_re_orig1 <- sapply(ideal_groups, function(g) get_value(cal_re_orig, g, "actual_survival"))
          cal_opt1$cal_opt <- cal_re_boot1 - cal_re_orig1
          cal_opt1$cal_obs_corr <- all_calibration[which(all_calibration$time_point == time_point), ]$actual_survival - cal_opt1$cal_opt
          cal_re_opt <- rbind(cal_re_opt, cal_opt1)
        }

        nrun <- nrun + 1
        cal_opt <- rbind(cal_opt, cal_re_opt)

      }, error = function(e) {
        message("Error in iteration ", i, ": ", e$message)
        return(NULL)
      })

      if (nrun == n) {
        message("Reached maximum number of successful runs: ", n)
        break
      }
    }

    # Calculate final calibrated metrics
    cal_means <- aggregate(cal_obs_corr ~ time_point + pred_group, data = cal_opt,
                           FUN = mean, na.rm = TRUE)
    cal_means <- cal_means[order(cal_means$time_point, cal_means$pred_group), ]

    cal_lower <- aggregate(cal_obs_corr ~ time_point + pred_group, data = cal_opt,
                           FUN = function(x) quantile(x, probs = 0.025, na.rm = TRUE))

    cal_upper <- aggregate(cal_obs_corr ~ time_point + pred_group, data = cal_opt,
                           FUN = function(x) quantile(x, probs = 0.975, na.rm = TRUE))

    cal_stats <- Reduce(function(x, y) merge(x, y, by = c("time_point", "pred_group")),
                        list(cal_means, cal_lower, cal_upper))

    names(cal_stats) <- c("time_point", "pred_group", "actual_survival", "lower_survival", "upper_survival")
    all_cal_autual <- cbind(cal_stats, all_calibration[, c("pred_survival", "n_events")])

    calibration_slpoe_with_ci <- calculate_calibration_slope(all_cal_autual)
    cindex_corr <- calculate_metric_ci(cindex_cor, time_points = seq_len)
    BS_corr <- calculate_metric_ci(BS_cor, time_points = seq_len)
    AUC_corr <- calculate_metric_ci(AUC_cor, time_points = seq_len)
  }

  # Return results
  return(list(
    time = seq_len,
    cindex_corr = cindex_corr,
    BS_corr = BS_corr,
    AUC_corr = AUC_corr,
    calibration_slpoe_with_ci = calibration_slpoe_with_ci,
    all_cal_autual = all_cal_autual
  ))
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

  if (length(formulas) != length(random_effect_formulas)) {
    stop("The number of formulas must be the same as the number of random effect formulas.")
  }
  for (i in 1:length(formulas)){
    formula <- formulas[[i]]
    random_effect_formula <- random_effect_formulas[[i]]
    model_name <- paste("model", i, sep = "_")


    models[[model_name]] <- lme(formula, random = random_effect_formula, data = data,control = list(opt = "optim"))
  }
  return(models)
}

individual_JM_predict<-function(object,data,Tstart,Dt=NULL,Thoriz=NULL){
  if (is.null(Thoriz)) Thoriz <- Tstart + Dt
  Tstart <- Tstart + 1e-06
  Thoriz <- Thoriz + 1e-06
  id_var <- object$model_info$var_names$idVar
  time_var <- object$model_info$var_names$time_var
  Time_var <- object$model_info$var_names$Time_var
  event_var <- object$model_info$var_names$event_var
  newdata <-data[data[[Time_var]] > Tstart, ]
  newdata <- newdata[ newdata[[time_var]] <= Tstart, ]
  newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]
  test1 <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
  if (!any(test1))
    stop("it seems that there are no events in the interval [Tstart, Thoriz).")
  newdata2 <- newdata
  newdata2[[Time_var]] <- Tstart
  newdata2[[event_var]] <- 0
  preds <- predict(object, newdata = newdata2, process = "event",
                   times = Thoriz)
  pi_u_t <- preds$pred
  names(pi_u_t) <- preds$id
  pi_u_t <- pi_u_t[preds$times > Tstart]

  id <- newdata[[id_var]]
  Time <- newdata[[Time_var]]
  event <- newdata[[event_var]]
  f <- factor(id, levels = unique(id))
  Time <- tapply(Time, f, tail, 1L)
  event <- tapply(event, f, tail, 1L)
  names(Time) <- names(event) <- as.character(unique(id))
  cal_DF <- data.frame(Time = Time, event = event, preds = 1-pi_u_t[names(Time)])
  return(cal_DF)
}

calibration_re_JM<-function(object=NULL,data=NULL,Tstart,Dt=NULL,Thoriz=NULL,all_pred=NULL,n_groups,probs=NULL){
  if (is.null(Thoriz)) Thoriz <- Tstart + Dt
  if (is.null(all_pred)){
    all_pred<-individual_JM_predict(object,data,Tstart,Dt=NULL,Thoriz)
  }
  if(is.null(probs)){
    cuts=quantile(all_pred$preds,probs=seq(0,1,length=n_groups+1))
    all_pred$pred_group<-cut(all_pred$preds,breaks=cuts)
  }
  else{
    cuts <- probs
    all_pred$pred_group<-cut(all_pred$preds,breaks=cuts)
  }

  cal_results <- data.frame()
  for(group in levels(all_pred$pred_group)){
    data_group<-all_pred[which(all_pred$pred_group==group),]
    if(nrow(data_group) < 3) {
      warning
      (paste("original", group, "has too few observations:", nrow(data_group)))
      next
    }
    formula<-as.formula(paste("Surv(",names(data_group)[1],",",names(data_group)[2],")~1"))
    km_fit <- survfit(formula, data = data_group)
    km_summary <- summary(km_fit, times = Thoriz,extend=TRUE)
    actual_survival <- km_summary$surv[1]
    se_survival <- km_summary$std.err[1]
    lower_survival <- max(0, actual_survival - 1.96* se_survival)
    upper_survival <- min(1, actual_survival + 1.96 * se_survival)
    cal_results<-rbind(cal_results,data.frame(time_point=Tstart,pred_group = group,n_patients = nrow(data_group),mean_predicted_survival = mean(data_group$preds),actual_survival=actual_survival,lower_survival=lower_survival,upper_survival=upper_survival))
  }
  return(list(cal_results=cal_results,cuts=cuts))
}

