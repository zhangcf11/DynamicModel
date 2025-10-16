#' External Validation for Landmarking Models
#'
#' This function performs external validation of landmarking models on an external dataset.
#' It evaluates discrimination (C-index, AUC), calibration, and overall performance (Brier score)
#' at multiple landmark times.
#'
#' @param object an objects inheriting from class "LMf" or "Vs_LM".
#' @param ex_data An external validation dataset with variable names identical to those in the model-building dataset.
#' @param n_group Integer, number of groups for calibration (default=10)
#'
#' @return A list containing:
#' \itemize{
#'   \item{time: The sequence of landmark times used for evaluation}
#'   \item{cindex: Concordance index values for each landmark time}
#'   \item{BS: Brier score values for each landmark time}
#'   \item{AUC: AUC values for each landmark time}
#'   \item{calibration_slope_with_ci: Calibration slope with confidence intervals}
#'   \item{all_cal_autual: Complete calibration results}
#' }
#' @importFrom dynpred cutLM
#' @importFrom survival survfit Surv
#' @importFrom Hmisc cut2
#' @export
Predictive_LM_ex <- function(object, ex_data, n_group = 10) {

  # Load required packages
  if (!require("survivalROC", quietly = TRUE)) {
    install.packages("survivalROC")
    library(survivalROC)
  }

  if (!require("dynpred", quietly = TRUE)) {
    install.packages("dynpred")
    library(dynpred)
  }

  if (!require("Hmisc", quietly = TRUE)) {
    install.packages("Hmisc")
    library(Hmisc)
  }

  # Extract model components
  TSet <- object$data
  Model <- object$Model
  sl <- object$tw$sl
  nsl <- length(sl)
  id <- object$id
  w <- object$tw$w
  time <- object$time
  status <- object$status
  rtime <- object$rtime
  cov <- object$cov
  func_covars <- object$func_covars
  func_lms <- object$func_lms

  # Initialize Vset
  Vset <- NULL

  # Prepare landmark data
  if (is.null(cov$vary)) {
    # Case with only fixed covariates
    fixed1 <- c(id, cov$fixed)
    for (j in seq_along(sl)) {
      LM <- dynpred::cutLM(
        data = ex_data,
        outcome = list(time = time, status = status),
        LM = sl[j],
        horizon = sl[j] + w,
        covs = list(fixed = fixed1, varying = cov$vary)
      )
      Vset <- rbind(Vset, LM)
    }
  } else {
    # Case with time-varying covariates
    for (j in seq_along(sl)) {
      LM <- dynpred::cutLM(
        data = ex_data,
        outcome = list(time = time, status = status),
        LM = sl[j],
        horizon = sl[j] + w,
        covs = list(fixed = cov$fixed, varying = cov$vary),
        format = "long",
        id = id,
        rtime = rtime,
        right = FALSE
      )
      Vset <- rbind(Vset, LM)
    }
  }

  # Order data by ID
  Vset <- Vset[order(Vset[[id]]), ]

  # Prepare data structure for adding interactions
  Vset1 <- list()
  Vset1$data <- Vset
  Vset1$lm_col <- "LM"

  # Add interaction terms
  Vset_i <- add_interactions(
    Vset1,
    c(cov$fixed, cov$vary),
    func_covars = func_covars,
    func_lms = func_lms,
    sl = sl
  )
  Vset_2 <- Vset_i$data

  # Initialize performance metrics
  cindex <- score <- auc <- rep(NA, nsl)
  cal_pred <- data.frame()

  # Evaluate performance at each landmark time
  for (i in seq_len(nsl)) {
    Vdata <- Vset_2[Vset_2$LM == sl[i], ]

    # Skip if no data at this landmark time
    if (nrow(Vdata) == 0) {
      warning(paste("No data available at landmark time", sl[i]))
      next
    }

    # Calculate performance metrics
    cindex[i] <- cal_cindex(model = Model, data = Vdata, time, status)
    score[i] <- cal_brierscore(
      model = Model,
      Tdata = TSet$data,
      Vdata = Vdata,
      width = w,
      tout = sl[i],
      time, status
    )
    auc[i] <- cal_auc(
      model = Model,
      data = Vdata,
      pred.t = sl[i] + w,
      time, status
    )

    # Get individual predictions for calibration
    all_pred1 <- individual_predict(model = Model, data = Vdata, sl[i], w)
    all_surv <- cbind(
      Vdata[[id]],
      all_pred1,
      rep(sl[i], nrow(Vdata)),
      Vdata[[time]] - sl[i],
      Vdata[[status]]
    )
    colnames(all_surv) <- c("id", "surv", "time_points", "auctual_time", "auctual_status")
    cal_pred <- rbind(cal_pred, all_surv)
  }

  # Calculate calibration
  cal_pred$occur_status <- ifelse(
    cal_pred$auctual_time < w & cal_pred$auctual_status == 1, 1, 0
  )

  all_cal_autual <- data.frame()
  for (i in unique(cal_pred$time_points)) {
    cal_pred_sub <- cal_pred[which(cal_pred$time_points == i), ]

    # Skip if insufficient data
    if (nrow(cal_pred_sub) < n_group) {
      warning(paste("Insufficient data at time", i, "for binning:",
                    nrow(cal_pred_sub), "obs vs required", n_group))
      next
    }

    if (length(unique(cal_pred_sub$surv)) < 2) {
      warning(paste("Insufficient surv value variation at time", i, "- skipping"))
      next
    }

    # Create prediction groups
    cal_pred_sub$pred_group <- Hmisc::cut2(cal_pred_sub$surv, g = n_group)

    bin_results <- data.frame()
    for (bin in levels(cal_pred_sub$pred_group)) {
      bin_data <- cal_pred_sub[which(cal_pred_sub$pred_group == bin), ]

      if (nrow(bin_data) < 3) {
        warning(paste("Bin", bin, "has too few observations:", nrow(bin_data)))
        next
      }

      # Calculate Kaplan-Meier estimates
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

  # Calculate calibration slope
  calibration_slope_with_ci <- calculate_calibration_slope(all_cal_autual)

  # Return results
  return(list(
    time = sl,
    cindex = cindex,
    BS = score,
    AUC = auc,
    calibration_slope_with_ci = calibration_slope_with_ci,
    all_cal_autual = all_cal_autual
  ))
}


