#' Plot Landmarking Model Predictions
#'
#' This function creates customized plots for landmarking model predictions,
#' allowing visualization of survival probabilities over landmark time points.
#' @usage plot_predict(object,newdata,ylim,main="Patient",xlab="landmark time points (s)",ylab="5-year survival",legend=NULL,...)
#' @param object An objects inheriting from class "LMf" or "Vs_LM".
#' @param newdata A data.frame containing the data for which predictions are made
#' @param ylim Numeric vector of length 2 specifying y-axis limits (optional)
#' @param main Character string for the main title (default="Patient")
#' @param xlab Character string for x-axis label (default="landmark time points (s)")
#' @param ylab Character string for y-axis label (default="5-year survival")
#' @param legend Character vector for legend labels (optional, default=NULL)
#' @param line_type Line type for the prediction curve (default=1)
#' @param line_width Line width for the prediction curve (default=6)
#' @param cex_main Character expansion for main title (default=1.5)
#' @param cex_lab Character expansion for axis labels (default=1.5)
#' @param cex_axis Character expansion for axis text (default=1.5)
#' @param axis_lwd Line width for axes (default=1.6)
#' @param tcl Tick length (default=-0.4)
#' @param padj_x Vertical adjustment for x-axis labels (default=-0.3)
#' @param hadj_y Horizontal adjustment for y-axis labels (default=0.9)
#' @param font_main Font for main title (default=7 for bold)
#' @param font_lab Font for axis labels (default=7 for bold)
#' @param line_main Line position for main title (default=1)
#' @param line_xlab Line position for x-axis label (default=2.6)
#' @param line_ylab Line position for y-axis label (default=2.6)
#' @param legend_pos Position of legend (default="topright")
#' @param legend_lty Line type in legend (default=1)
#' @param legend_lwd Line width in legend (default=6)
#' @param legend_bty Box type for legend (default="n" for no box)
#' @param ... Additional graphical parameters passed to plot, axis, and title functions
#' @param ... extra aguments; currently none is used.
#' @importFrom graphics plot axis title legend par
#' @export
#' @examples
#' library(DynamicModel)
#' data(renal)
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=TRUE)
#' a1<-VS_LM(fit,cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")))
#' newdata<-renal[which(renal$id==5619),]
#' plot_predict(object=a1,newdata=newdata,main="")
plot_predict <- function(object, newdata, ylim, main = "Patient",
                         xlab = "landmark time points (s)",
                         ylab = "5-year survival", legend = NULL,
                         line_type = 1, line_width = 6,
                         cex_main = 1.5, cex_lab = 1.5, cex_axis = 1.5,
                         axis_lwd = 1.6, tcl = -0.4,
                         padj_x = -0.3, hadj_y = 0.9,
                         font_main = 7, font_lab = 7,
                         line_main = 1, line_xlab = 2.6, line_ylab = 2.6,
                         legend_pos = "topright", legend_lty = 1,
                         legend_lwd = 6, legend_bty = "n", ...) {

  # Generate predictions
  pp <- predict_LM(object, newdata)

  # Validate prediction structure
  if (ncol(pp) < 2) {
    stop("Prediction matrix must have at least 2 columns (time and survival probability)")
  }

  # Set y-axis limits if not provided
  if (missing(ylim)) {
    ylim <- c(0, 1)
  }

  # Set graphical parameters
  old_par <- par(no.readonly = TRUE)  # Save current parameters
  on.exit(par(old_par))  # Restore parameters on exit

  par(xaxs = "i", yaxs = "i", mar = c(5, 5, 3, 2))

  # Create main plot
  plot(pp[, 1], pp[, 2],
       type = "l",
       lwd = line_width,
       lty = line_type,
       xlim = c(min(pp[, 1]), max(pp[, 1])),
       ylim = ylim,
       bty = "l",  # L-shaped box
       xaxt = "n", yaxt = "n",  # Suppress default axes
       xlab = "", ylab = "",  # Empty labels for custom ones
       ...)  # Pass additional graphical parameters

  # Add custom x-axis
  axis(1, las = 1, pos = 0,
       cex.axis = cex_axis,
       tcl = tcl,
       padj = padj_x,
       lwd = axis_lwd, ...)

  # Add custom y-axis
  axis(2, las = 1, pos = 0,
       cex.axis = cex_axis,
       tcl = tcl,
       hadj = hadj_y,
       lwd = axis_lwd, ...)

  # Add titles
  title(main = list(main, font = font_main, cex = cex_main),
        line = line_main)
  title(xlab = xlab, font.lab = font_lab, cex.lab = cex_lab,
        line = line_xlab, ...)
  title(ylab = ylab, font.lab = font_lab, cex.lab = cex_lab,
        line = line_ylab, ...)

  # Add legend if requested - FIXED: condition was inverted
  if (!is.null(legend)) {
    legend(legend_pos,
           legend = legend,
           lty = legend_lty,
           lwd = legend_lwd,
           bty = legend_bty,
           ...)
  }

}
