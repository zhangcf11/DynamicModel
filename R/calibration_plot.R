#' Calibration Plot for Survival or Event Rate Predictions
#'
#' This function creates calibration plots for survival analysis models,
#' displaying the relationship between predicted and observed rates.
#' It can plot either survival rates or event rates (1 - survival rate)
#' with confidence intervals and an ideal calibration line.
#'
#' @param calbrationdata A data frame containing calibration data. Must include
#'   columns: `time_point`, `mean_predicted_survival`, `actual_survival`,
#'   `lower_survival`, and `upper_survival`.
#' @param time_points Numeric value specifying the time point for which to
#'   create the calibration plot.
#' @param xlim Numeric vector of length 2 specifying the x-axis limits.
#'   Default is c(0, 1).
#' @param ylim Numeric vector of length 2 specifying the y-axis limits.
#'   Default is c(0, 1).
#' @param xlab Character string for the x-axis label. If not specified and
#'   `inverse = TRUE`, defaults to "Predicted event rate".
#' @param ylab Character string for the y-axis label. If not specified and
#'   `inverse = TRUE`, defaults to "Observed event rate".
#' @param add_ideal Logical value indicating whether to add the ideal
#'   calibration line (diagonal). Default is TRUE.
#' @param pch Point character type. Default is 19 (solid circle).
#' @param col Point color. Default is "blue".
#' @param inverse Logical value indicating whether to plot survival rates
#'   (FALSE) or event rates (TRUE). Event rates are calculated as 1 - survival rate.
#'   Default is FALSE.
#' @param font_lab Font style for axis labels (1=normal, 2=bold, 3=italic, 4=bold italic).
#'   Default is 2 (bold).
#' @param font_main Font style for main title. Default is 2 (bold).
#' @param font_axis Font style for axis tick labels. Default is 1 (normal).
#' @param cex_lab Font size multiplier for axis labels. Default is 1.
#' @param cex_main Font size multiplier for main title. Default is 1.2.
#' @param cex_axis Font size multiplier for axis tick labels. Default is 1.
#' @param cex_point Point size multiplier. Default is 1.2.
#' @param line_lab Distance of axis labels from the plot. Default is 2.6.
#' @param line_main Distance of main title from the plot. Default is 1.
#' @param arrow_col Color for error bar arrows. Default is semi-transparent blue.
#' @param arrow_lwd Line width for arrows. Default is 1.5.
#' @param arrow_length Length of arrow heads. Default is 0.05.
#' @param ideal_col Color for ideal calibration line. Default is "red".
#' @param ideal_lwd Line width for ideal calibration line. Default is 2.
#' @param ideal_lty Line type for ideal calibration line. Default is 2 (dashed).
#' @param main_title Character string for the main title. Default is empty string.
#'
#' @return No return value. The function is called for its side effect of
#'   generating a calibration plot.
#' @export
#'
#' @examples
#' library(DynamicModel)
#' data(renal)
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=TRUE)
#' a1<-VS_LM(fit,cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")))
#' LM_pre<-predictive_LM(a1,method="CV",n=1,seed=1,K=5)
#' calibration_plot(LM_pre$all_cal_autual,5)
calibration_plot <- function(calbrationdata, time_points,
                             xlim = c(0, 1), ylim = c(0, 1),
                             xlab = "Predicted survival rate",
                             ylab = "Observed survival rate",
                             add_ideal = TRUE,
                             pch = 19,
                             col = "blue",
                             # New: whether to calculate event rate (1 - survival rate)
                             inverse = FALSE,      # FALSE: survival rate, TRUE: event rate
                             # Font related parameters
                             font_lab = 2,
                             font_main = 2,
                             font_axis = 1,
                             cex_lab = 1,
                             cex_main = 1.2,
                             cex_axis = 1,
                             cex_point = 1.2,
                             # Line related parameters
                             line_lab = 2.6,
                             line_main = 1,
                             # Arrow related parameters
                             arrow_col = rgb(0, 0, 1, 0.4),
                             arrow_lwd = 1.5,
                             arrow_length = 0.05,
                             # Ideal line related parameters
                             ideal_col = "red",
                             ideal_lwd = 2,
                             ideal_lty = 2,
                             # Title
                             main_title = "") {

  # Filter data for specific time point
  data <- calbrationdata[which(calbrationdata$time_point == time_points), ]

  # Determine whether to use survival rate or event rate based on inverse parameter
  if (inverse) {
    # Calculate event rate (1 - survival rate)
    x_values <- 1 - data$mean_predicted_survival
    y_values <- 1 - data$actual_survival

    # For confidence intervals, note that upper and lower bounds need to be swapped
    # because upper event rate = 1 - lower survival rate
    # and lower event rate = 1 - upper survival rate
    lower_values <- 1 - data$upper_survival
    upper_values <- 1 - data$lower_survival

    # Update axis labels if using default values
    if (missing(xlab)) xlab <- "Predicted event rate"
    if (missing(ylab)) ylab <- "Observed event rate"
  } else {
    # Use survival rate (default)
    x_values <- data$mean_predicted_survival
    y_values <- data$actual_survival
    lower_values <- data$lower_survival
    upper_values <- data$upper_survival
  }

  # Create base plot
  plot(x_values, y_values,
       pch = pch, col = col, cex = cex_point,
       xlab = "", ylab = "",
       ylim = ylim, xlim = xlim, main = main_title,
       bty = "l", cex.main = cex_main, font.main = font_main)

  # Add axis labels with specified formatting
  title(xlab = xlab, font.lab = font_lab, line = line_lab, cex.lab = cex_lab)
  title(ylab = ylab, font.lab = font_lab, line = line_lab, cex.lab = cex_lab)

  # Add ideal calibration line (diagonal)
  if (add_ideal) {
    abline(0, 1, col = ideal_col, lwd = ideal_lwd, lty = ideal_lty)
  }

  # Draw arrows only for non-zero length to avoid warnings
  non_zero_arrows <- which(lower_values != upper_values)
  if (length(non_zero_arrows) > 0) {
    arrows(x0 = x_values[non_zero_arrows],
           y0 = lower_values[non_zero_arrows],
           x1 = x_values[non_zero_arrows],
           y1 = upper_values[non_zero_arrows],
           angle = 90, code = 3, length = arrow_length,
           col = arrow_col, lwd = arrow_lwd)
  }
}
