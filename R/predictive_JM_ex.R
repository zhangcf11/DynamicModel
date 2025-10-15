#' External Validation for Joint Models
#'
#' This function performs external validation of joint models on an external dataset.
#' It evaluates discrimination (C-index, AUC), calibration, and overall performance (Brier score)
#' at multiple time points.
#'
#' @param object A fitted joint model object from the `jm` function
#' @param newdata A data.frame containing the external validation data
#' @param seq_len Numeric vector of time points for prediction evaluation
#' @param w Numeric, prediction window length (Dt parameter)
#' @param n_group Integer, number of groups for calibration (default=10)
#'
#' @return A list containing:
#' \itemize{
#'   \item{time: The sequence of time points used for evaluation}
#'   \item{cindex: Concordance index values for each time point}
#'   \item{BS: Brier score values for each time point}
#'   \item{AUC: AUC values for each time point}
#'   \item{calibration_slope_with_ci: Calibration slope with confidence intervals}
#'   \item{all_cal_actual: Complete calibration results}
#' }
#' @export
predictive_JM_ex<-function(object,newdata,seq_len,w, n_group = 10){
  c_index<-Brier_score<-AUC<-c()
  all_cal_autual<-data.frame()
  for (i in seq_along(seq_len) ){
    timepoint<-seq_len[i]
    c_index[i]<-tvC_index(object,newdata, Tstart =timepoint, Dt =w)
    Brier_score[i]<-tvBrier(object,newdata, Tstart =timepoint, Dt =w)$Brier
    AUC[i]<-tvAUC(object,newdata, Tstart =timepoint, Dt =w)$auc
    cal_autual_all<-calibration_re_JM(object = object, data = newdata,
                                      Tstart =timepoint, Dt = w, n_groups = n_group)$cal_results
    all_cal_autual<-rbind(all_cal_autual,cal_autual_all)
  }
  calibration_slpoe_with_ci <- calculate_calibration_slope(all_cal_autual)
  # Return results
  return(list(
    time = seq_len,
    cindex = c_index,
    BS =Brier_score,
    AUC = AUC,
    calibration_slope_with_ci = calibration_slope_with_ci,
    all_cal_autual = all_cal_autual
  ))
}
