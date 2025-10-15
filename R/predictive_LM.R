#' Predictive Performance Assessment for Landmarking Models
#'
#' This function evaluates the predictive performance of landmarking models using
#' different validation methods including K-fold cross-validation, bootstrap validation,
#' and data splitting.
#' @usage predictive_LM(object,method,n,K,seed,r)
#' @param object an objects inheriting from class "LMf" or "Vs_LM"
#' @param method Character string specifying validation method: "K-fold", "boot", or "split"
#' @param n Number of repetitions for internal cross-validation. For random splitting, it refers to the number of random splitting iterations; for K-fold cross-validation, it refers to the number of outer loop iterations; for bootstrap, it refers to the number of resampling iterations.
#' @param K Integer, number of folds for K-fold cross-validation
#' @param seed Integer, random seed for reproducibility
#' @param r Numeric, proportion of data to use for training (only used when method="split")
#' @param n_group Integer, number of groups for calibration (default=10)
#'
#' @return A list containing:
#' \itemize{
#'   \item{time: The sequence of landmark times used for evaluation}
#'   \item{cindex_corr: Corrected concordance index with confidence intervals}
#'   \item{BS_corr: Corrected Brier score with confidence intervals}
#'   \item{AUC_corr: Corrected AUC with confidence intervals}
#'   \item{calibration_slpoe_with_ci: Calibration slope with confidence intervals}
#'   \item{all_cal_autual: Complete calibration results}
#' }
#' @export
#' @examples
#' library(DynamicModel)
#' data(renal)
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=TRUE)
#' a1<-VS_LM(fit,cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")))
#' LM_pre<-predictive_LM(a1,method="K-fold",n=1,seed=1,K=5)
predictive_LM <- function(object, method, n, K, seed, r, n_group = 10) {

  # Input validation
  if (!method %in% c("K-fold", "boot", "split")) {
    stop("method must be one of: 'K-fold', 'boot', 'split'")
  }
  if (!require("caret")) {
    install.packages("caret")
    library(caret)
  }
  if (!require("survivalROC")) {
    install.packages("survivalROC")
    library(survivalROC)
  }
  if (!require("Hmisc", quietly = TRUE)) {
    install.packages("Hmisc")
    library(Hmisc)
  }
  # Extract model components
  fm <- object$Model
  id <- object$id
  data <- object$data$data
  sl <- object$tw$sl
  w <- object$tw$w
  time <- object$time
  status <- object$status

  # Initialize result containers
  cindex_corr <- BS_corr <- AUC_corr <- calibration_slpoe_with_ci <- all_cal_autual <- NULL

  if (method == "K-fold") {
    # K-fold cross-validation implementation
    cindex_sm <- BS_sm <- AUC_sm <- list()
    cal_pred <- list()
    nrun <- 0

    for (i in 1:(5 * n)) {
      tryCatch({
        set.seed(seed + 10 * i)
        folds <- caret::createFolds(unique(data[[id]]), K)
        cindex <- BS <- AUC <- matrix(NA, nrow = K, ncol = length(sl))
        all_pred <- matrix(NA, ncol = 6, nrow = 0)
        colnames(all_pred) <- c("id", "surv", "time_points", "auctual_time", "auctual_status", "temp")

        for (j in 1:K) {
          index <- unique(data[[id]])[folds[[j]]]
          train_idx <- which(!data[[id]] %in% index)
          test_idx <- which(data[[id]] %in% index)

          # Create formula with cluster term
          formula_text <- paste(deparse(fm$formula), collapse = "")
          formula_text <- gsub("\\s+", " ", formula_text)
          new_formula_text <- paste0(formula_text, " + cluster(", id, ")")
          new_formula <- as.formula(new_formula_text)

          # Fit Cox model
          fit <- coxph(new_formula, data = data[train_idx, ])

          for (k in seq_along(sl)) {
            # Get test data at current landmark time
            test_idx1 <- test_idx[data$LM[test_idx] == sl[k]]
            if (length(test_idx1) == 0) next

            # Calculate performance metrics
            cindex[j, k] <- cal_cindex(model = fit, data = data[test_idx1, ], time, status)
            BS[j, k] <- cal_brierscore(model = fit, Tdata = data[train_idx, ],
                                       Vdata = data[test_idx1, ], width = w,
                                       tout = sl[k], time, status)
            AUC[j, k] <- cal_auc(model = fit, data = data[test_idx1, ],
                                 pred.t = sl[k] + w, time, status)

            # Get individual predictions
            all_pred1 <- individual_predict(model = fit, data = data[test_idx1, ], sl[k], w)
            all_surv <- cbind(
              data[test_idx1, ][[id]],
              all_pred1,
              rep(sl[k], length(all_pred1)),
              data[test_idx1, ][[time]] - sl[k],
              data[test_idx1, ][[status]],
              rep(i, length(all_pred1))  # iteration identifier
            )
            colnames(all_surv) <- c("id", "surv", "time_points", "auctual_time", "auctual_status", "num")
            all_pred <- rbind(all_pred, all_surv)
          }
        }

        # Store results for this iteration
        cindex_sm[[i]] <- cindex
        BS_sm[[i]] <- BS
        AUC_sm[[i]] <- AUC
        cal_pred[[i]] <- all_pred
        nrun <- nrun + 1

      }, error = function(e) {
        message("Error in iteration ", i, ": ", conditionMessage(e))
        return(NULL)
      })

      if (nrun == n) {
        message("Reached maximum number of successful runs: ", nrun)
        break
      }
    }

    # Calculate mean performance metrics
    cindex_array <- mean_ignore_na(cindex_sm)
    BS_array <- mean_ignore_na(BS_sm)
    AUC_array <- mean_ignore_na(AUC_sm)

    # Calculate confidence intervals
    cindex_corr <- do.call(rbind, apply(cindex_array, 2, calculate_cv_ci_t))
    BS_corr <- do.call(rbind, apply(BS_array, 2, calculate_cv_ci_t))
    AUC_corr <- do.call(rbind, apply(AUC_array, 2, calculate_cv_ci_t))

    # Process calibration predictions
    cal_array <- lapply(cal_pred, function(mat) {
      order_idx <- order(mat[, "id"], mat[, "time_points"])
      mat_sorted <- mat[order_idx, ]
      return(mat_sorted)
    })

    # Calculate mean predictions across iterations
    cal_pred_array <- matrix(NA, nrow = nrow(cal_array[[1]]), ncol = length(cal_array))
    for (i in seq_along(cal_array)) {
      cal_pred_array[, i] <- cal_array[[i]][, "surv"]
    }

    cal_pred_mean <- rowMeans(cal_pred_array)
    cal_pred_corr <- cbind(
      cal_array[[1]][, "id"],
      cal_pred_mean,
      cal_array[[1]][, "time_points"],
      cal_array[[1]][, "auctual_time"],
      cal_array[[1]][, "auctual_status"]
    )
    cal_pred_corr <- as.data.frame(cal_pred_corr)
    names(cal_pred_corr) <- c("id", "surv", "time_points", "auctual_time", "auctual_status")

    # Calculate occurrence status
    cal_pred_corr$occur_status <- ifelse(cal_pred_corr$auctual_time < w &
                                           cal_pred_corr$auctual_status == 1, 1, 0)

    # Calculate calibration
    all_cal_autual <- data.frame()
    for (i in unique(cal_pred_corr$time_points)) {
      cal_pred_sub <- cal_pred_corr[which(cal_pred_corr$time_points == i), ]
      cal_pred_sub$pred_group <- Hmisc::cut2(cal_pred_sub$surv, g = n_group)

      # Check data adequacy for grouping
      if (nrow(cal_pred_sub) < n_group) {
        warning(paste("Insufficient data at time", i, "for binning:",
                      nrow(cal_pred_sub), "obs vs required", n_group))
        next
      }

      # Check survival probability variability
      if (length(unique(cal_pred_sub$surv)) < 2) {
        warning(paste("Insufficient surv value variation at time", i, "- skipping"))
        next
      }

      bin_results <- data.frame()
      for (bin in levels(cal_pred_sub$pred_group)) {
        bin_data <- cal_pred_sub[which(cal_pred_sub$pred_group == bin), ]
        if (nrow(bin_data) < 3) {
          warning(paste("Bin", bin, "has too few observations:", nrow(bin_data)))
          next
        }

        km_fit <- survfit(Surv(auctual_time, auctual_status) ~ 1, data = bin_data)
        km_summary <- summary(km_fit, times = w)
        actual_survival <- ifelse(length(km_summary$surv) > 0, km_summary$surv[1], NA)
        se_survival <- ifelse(length(km_summary$std.err) > 0, km_summary$std.err[1], NA)

        if (!is.na(actual_survival) && !is.na(se_survival)) {
          lower_survival <- max(0, actual_survival - 1.96 * se_survival)
          upper_survival <- min(1, actual_survival + 1.96 * se_survival)
        } else {
          lower_survival <- upper_survival <- NA
        }

        bin_results <- rbind(bin_results, data.frame(
          time_point = i,
          pred_group = bin,
          n_patients = nrow(bin_data),
          mean_predicted_survival = mean(bin_data$surv),
          actual_survival = actual_survival,
          lower_survival = lower_survival,
          upper_survival = upper_survival
        ))
      }
      all_cal_autual <- rbind(all_cal_autual, bin_results)
    }

    # Calculate calibration slope
    calibration_slpoe_with_ci <- calculate_calibration_slope(all_cal_autual)

  } else if (method == "boot") {
    # Bootstrap validation implementation
    cindex_app <- BS_app <- AUC_app <- numeric(length(sl))
    all_calibration <- data.frame()
    cuts <- list()

    # Calculate apparent performance
    for (k in seq_along(sl)) {
      data1 <- data[data$LM == sl[k], ]
      cindex_app[k] <- cal_cindex(model = fm, data = data1, time, status)
      BS_app[k] <- cal_brierscore(model = fm, Tdata = data1, Vdata = data1,
                                  width = w, tout = sl[k], time, status)
      AUC_app[k] <- cal_auc(model = fm, data = data1, pred.t = sl[k] + w, time, status)

      cal_re <- calibration_re(model = fm, data = data1, time, status, sl[k], w, n_group)
      all_calibration <- rbind(all_calibration, cal_re$cal_results)
      cuts[[k]] <- cal_re$cuts
    }

    # Bootstrap correction
    cindex_opt <- BS_opt <- AUC_opt <- cindex_cor <- BS_cor <- AUC_cor <-
      data.frame(array(NA, dim = c((5 * n), length(sl))))
    cal_opt <- data.frame()
    nrun <- 0

    for (i in 1:(5 * n)) {
      tryCatch({
        set.seed(seed + 10 * i)
        train_index <- sample(unique(data[[id]]), length(unique(data[[id]])), replace = TRUE)
        train_data <- NULL

        # Create bootstrap sample
        for (j in seq_along(train_index)) {
          t_data <- data[which(data[[id]] == train_index[j]), ]
          train_data <- rbind(train_data, t_data)
        }
        train_data <- train_data[order(train_data[[id]]), ]

        # Fit model on bootstrap sample
        formula_text <- paste(deparse(fm$formula), collapse = "")
        formula_text <- gsub("\\s+", " ", formula_text)
        new_formula_text <- paste0(formula_text, " + cluster(", id, ")")
        new_formula <- as.formula(new_formula_text)
        fit <- coxph(new_formula, data = train_data)

        cal_re_opt <- data.frame()
        for (k in seq_along(sl)) {
          train_data1 <- train_data[train_data$LM == sl[k], ]
          data1 <- data[data$LM == sl[k], ]

          # Calculate bootstrap performance
          cindex_boot <- cal_cindex(model = fit, data = train_data1, time, status)
          BS_boot <- cal_brierscore(model = fit, Tdata = train_data1, Vdata = train_data1,
                                    width = w, tout = sl[k], time, status)
          AUC_boot <- cal_auc(model = fit, data = train_data1, pred.t = sl[k] + w, time, status)

          # Calculate performance on original data
          cindex_orig <- cal_cindex(model = fit, data = data1, time, status)
          BS_orig <- cal_brierscore(model = fit, Tdata = train_data1, Vdata = data1,
                                    width = w, tout = sl[k], time, status)
          AUC_orig <- cal_auc(model = fit, data = data1, pred.t = sl[k] + w, time, status)

          # Calculate calibration
          cal_re_boot <- calibration_re(model = fit, data = train_data1, time, status,
                                        sl[k], w, n_group, cuts[[k]])$cal_results
          cal_re_orig <- calibration_re(model = fit, data = data1, time, status,
                                        sl[k], w, n_group, cuts[[k]])$cal_results

          # Calculate optimism
          cindex_opt[i, k] <- cindex_boot - cindex_orig
          cindex_cor[i, k] <- cindex_app[k] - cindex_opt[i, k]

          BS_opt[i, k] <- BS_boot - BS_orig
          BS_cor[i, k] <- BS_app[k] - BS_opt[i, k]

          AUC_opt[i, k] <- AUC_boot - AUC_orig
          AUC_cor[i, k] <- AUC_app[k] - AUC_opt[i, k]

          # Calculate calibration optimism
          ideal_groups <- levels(cut(0, breaks = cuts[[k]]))
          cal_opt1 <- data.frame(
            num = i,
            time_point = sl[k],
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
          cal_opt1$cal_obs_corr <- all_calibration[which(all_calibration$time_point == sl[k]), ]$actual_survival - cal_opt1$cal_opt
          cal_re_opt <- rbind(cal_re_opt, cal_opt1)
        }

        nrun <- nrun + 1
        cal_opt <- rbind(cal_opt, cal_re_opt)

      }, error = function(e) {
        message("Error in iteration ", i, ": ", conditionMessage(e))
        return(NULL)
      })

      if (nrun == n) {
        message("Reached maximum number of successful runs: ", nrun)
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
    all_cal_autual <- cbind(cal_stats, all_calibration[, c("mean_predicted_survival", "n_patients")])

    calibration_slpoe_with_ci <- calculate_calibration_slope(all_cal_autual)
    cindex_corr <- calculate_metric_ci(cindex_cor, time_points = sl)
    BS_corr <- calculate_metric_ci(BS_cor, time_points = sl)
    AUC_corr <- calculate_metric_ci(AUC_cor, time_points = sl)

  } else if (method == "split") {
    # Data splitting implementation
    cindex <- BS <- AUC <- data.frame(array(NA, dim = c((5 * n), length(sl))))
    nrun <- 0
    cal_pred <- data.frame()

    for (i in 1:(5 * n)) {
      tryCatch({
        set.seed(seed + 10 * i)
        index <- sample(unique(data[[id]]), length(unique(data[[id]])) * r)
        train_data <- data[which(data[[id]] %in% index), ]
        test_data <- data[-which(data[[id]] %in% index), ]

        # Fit model
        formula_text <- paste(deparse(fm$formula), collapse = "")
        formula_text <- gsub("\\s+", " ", formula_text)
        new_formula_text <- paste0(formula_text, " + cluster(", id, ")")
        new_formula <- as.formula(new_formula_text)
        fit <- coxph(new_formula, data = train_data)

        all_pred <- data.frame()
        for (k in seq_along(sl)) {
          test_data1 <- test_data[which(test_data$LM == sl[k]), ]
          if (nrow(test_data1) == 0) next

          cindex[i, k] <- cal_cindex(model = fit, data = test_data1, time, status)
          BS[i, k] <- cal_brierscore(model = fit, Tdata = train_data, Vdata = test_data1,
                                     width = w, tout = sl[k], time, status)
          AUC[i, k] <- cal_auc(model = fit, data = test_data1, pred.t = sl[k] + w, time, status)

          all_pred1 <- individual_predict(model = fit, data = test_data1, sl[k], w)
          all_surv <- cbind(
            rep(i, nrow(test_data1)),
            test_data1[[id]],
            all_pred1,
            rep(sl[k], nrow(test_data1)),
            test_data1[[time]] - sl[k],
            test_data1[[status]]
          )
          colnames(all_surv) <- c("num", "id", "surv", "time_points", "auctual_time", "auctual_status")
          all_pred <- rbind(all_pred, all_surv)
        }
        nrun <- nrun + 1
        cal_pred <- rbind(cal_pred, all_pred)

      }, error = function(e) {
        message("Error in iteration ", i, ": ", conditionMessage(e))
        return(NULL)
      })

      if (nrun == n) {
        message("Reached maximum number of successful runs: ", nrun)
        break
      }
    }

    # Calculate calibration
    cal_pred$occur_status <- ifelse(cal_pred$auctual_time < w & cal_pred$auctual_status == 1, 1, 0)

    all_cal_autual <- data.frame()
    for (i in unique(cal_pred$time_points)) {
      cal_pred_sub <- cal_pred[which(cal_pred$time_points == i), ]
      cal_pred_sub$pred_group <- Hmisc::cut2(cal_pred_sub$surv, g = n_group)

      # Check data adequacy
      if (nrow(cal_pred_sub) < n_group) {
        warning(paste("Insufficient data at time", i, "for binning:",
                      nrow(cal_pred_sub), "obs vs required", n_group))
        next
      }

      if (length(unique(cal_pred_sub$surv)) < 2) {
        warning(paste("Insufficient surv value variation at time", i, "- skipping"))
        next
      }

      bin_results <- data.frame()
      for (bin in levels(cal_pred_sub$pred_group)) {
        bin_data <- cal_pred_sub[which(cal_pred_sub$pred_group == bin), ]
        if (nrow(bin_data) < 3) {
          warning(paste("Bin", bin, "has too few observations:", nrow(bin_data)))
          next
        }

        km_fit <- survfit(Surv(auctual_time, auctual_status) ~ 1, data = bin_data)
        km_summary <- summary(km_fit, times = w, extend = TRUE)
        actual_survival <- ifelse(length(km_summary$surv) > 0, km_summary$surv[1], NA)
        se_survival <- ifelse(length(km_summary$std.err) > 0, km_summary$std.err[1], NA)

        if (!is.na(actual_survival) && !is.na(se_survival)) {
          lower_survival <- max(0, actual_survival - 1.96 * se_survival)
          upper_survival <- min(1, actual_survival + 1.96 * se_survival)
        } else {
          lower_survival <- upper_survival <- NA
        }

        bin_results <- rbind(bin_results, data.frame(
          time_point = i,
          pred_group = bin,
          n_patients = nrow(bin_data),
          mean_predicted_survival = mean(bin_data$surv),
          actual_survival = actual_survival,
          lower_survival = lower_survival,
          upper_survival = upper_survival
        ))
      }
      all_cal_autual <- rbind(all_cal_autual, bin_results)
    }

    calibration_slpoe_with_ci <- calculate_calibration_slope(all_cal_autual)
    cindex_corr <- calculate_metric_ci(metric_matrix = cindex, time_points = sl, method = "quantile")
    BS_corr <- calculate_metric_ci(metric_matrix = BS, time_points = sl, method = "quantile")
    AUC_corr <- calculate_metric_ci(metric_matrix = AUC, time_points = sl, method = "quantile")
  }

  # Return results
  return(list(
    time = sl,
    cindex_corr = cindex_corr,
    BS_corr = BS_corr,
    AUC_corr = AUC_corr,
    calibration_slpoe_with_ci = calibration_slpoe_with_ci,
    all_cal_autual = all_cal_autual
  ))
}
calculate_cv_ci_t <- function(cv_results, alpha = 0.05) {

  n <- length(cv_results)
  mean_value <- mean(cv_results, na.rm = TRUE)
  se_value <- sd(cv_results, na.rm = TRUE) / sqrt(n)

  t_critical <- qt(1 - alpha/2, df = n - 1)

  ci_lower <- mean_value - t_critical * se_value
  ci_upper <- mean_value + t_critical * se_value

  return(
    data.frame(estimate = mean_value,
               ci_lower = ci_lower,
               ci_upper = ci_upper)
  )
}

mean_ignore_na <- function(matrices) {
  result <- matrix(NA, nrow = nrow(matrices[[1]]), ncol = ncol(matrices[[1]]))
  for(i in 1:nrow(result)) {
    for(j in 1:ncol(result)) {
      values <- sapply(matrices, function(x) x[i, j])
      result[i, j] <- mean(values, na.rm = TRUE)
    }
  }

  return(result)
}
calculate_calibration_slope<-function(calibration_data){
  calibration_metrics <- data.frame()
  for(i in unique(calibration_data$time_point)){
    time_dt<-calibration_data[which(calibration_data$time_point==i),]
    fit<-lm(actual_survival~mean_predicted_survival,data= time_dt,weights = n_patients)
    slope<-coef(fit)[2]
    vcov_matrix <- vcov(fit)
    slope_se <- sqrt(vcov_matrix[2, 2])
    slope_lower <- slope - 1.96 * slope_se
    slope_upper <- slope + 1.96 * slope_se
    calibration_metrics  <- rbind(calibration_metrics, data.frame( time_point = i, calibration_slope = slope, slope_lower_ci = slope_lower, slope_upper_ci  = slope_upper))
  }
  return(calibration_metrics)
}
calibration_re<-function(model,data,time,status,s,w,n_groups,probs=NULL){
  data$pred<-individual_predict(model=model,data=data,s,w)
  if(is.null(probs)){
    cuts=quantile(data$pred,probs=seq(0,1,length=n_groups+1))
    data$pred_group<-cut(data$pred,breaks=cuts)
  }
  else{
    cuts <- probs
    data$pred_group<-cut(data$pred,breaks=cuts)
  }

  cal_results <- data.frame()
  for(group in levels(data$pred_group)){
    data_group<-data[which(data$pred_group==group),]
    if(nrow(data_group) < 3) {
      warning
      (paste("original", group, "has too few observations:", nrow(data_group)))
      next
    }
    formula<-as.formula(paste("Surv(",time,",",status,")~1"))
    km_fit <- survfit(formula, data = data_group)
    km_summary <- summary(km_fit, times = s+w)
    actual_survival <- km_summary$surv[1]
    #se_survival <- km_summary$std.err[1]
    #lower_survival <- max(0, actual_survival - 1.96* se_survival)
    #upper_survival <- min(1, actual_survival + 1.96 * se_survival)
    cal_results<-rbind(cal_results,data.frame(time=s,pred_group = group,n_patients = nrow(data_group),mean_predicted_survival = mean(data_group$pred),actual_survival=actual_survival))
  }
  return(list(cal_results=cal_results,cuts=cuts))
}
individual_predict<-function(model,data,s,w){
  n<-nrow(data)
  predicted_conditional_survival<-c()
  for (i in 1:n) {
    individual_data <- data[i, ]
    sfit <- survfit(model, newdata = individual_data)
    sfit_summary <- summary(sfit)
    tmp <- dynpred::evalstep(sfit_summary$time, sfit_summary$surv, c(s, s + w), subst = 1)
    if (length(tmp) == 2 && !any(is.na(tmp))) {
      predicted_conditional_survival[i] <- tmp[2] / tmp[1]
    } else {
      predicted_conditional_survival[i] <- NA
    }
  }
  return(predicted_conditional_survival)
}
calculate_metric_ci <- function(metric_matrix, time_points,
                                alpha = 0.05, method = "quantile", n_boot = 1000) {

  result_df <- data.frame()


  for (k in 1:length(time_points)) {

    values <- metric_matrix[, k]


    values <- na.omit(values)


    if (length(values) < 5) {
      warning(paste("Insufficient data for time point", time_points[k],
                    "only", length(values), "observations available"))
      next
    }


    mean_val <- mean(values)
    se_val <- sd(values) / sqrt(length(values))


    if (method == "quantile") {

      lower <- quantile(values, alpha/2, na.rm = TRUE)
      upper <- quantile(values, 1 - alpha/2, na.rm = TRUE)
    } else if (method == "bootstrap") {

      tryCatch({
        boot_results <- boot::boot(
          data = values,
          statistic = function(x, i) mean(x[i]),
          R = n_boot
        )

        boot_ci <- boot::boot.ci(boot_results, conf = 1 - alpha, type = "perc")
        lower <- boot_ci$percent[4]
        upper <- boot_ci$percent[5]
      }, error = function(e) {
        message("Bootstrap failed for time point ", time_points[k],
                " in metric ", metric_name, ": ", conditionMessage(e))

        z <- qnorm(1 - alpha/2)
        lower <- mean_val - z * se_val
        upper <- mean_val + z * se_val
      })
    } else {

      z <- qnorm(1 - alpha/2)
      lower <- mean_val - z * se_val
      upper <- mean_val + z * se_val
    }


    result_df <- rbind(result_df, data.frame(
      estimate = mean_val,
      ci_lower = lower,
      ci_upper = upper
    ))
  }

  rownames(result_df) <- NULL
  return(as.data.frame(result_df))
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

