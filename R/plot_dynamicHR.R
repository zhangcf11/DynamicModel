#'#' This function creates plots of dynamic hazard ratios over time from joint models,
#' showing the time-varying effect of covariates on the hazard function.
#' @usage plot_dynamicHR(object,variable,form_splines)
#' @param object An objects inheriting from class "jm".
#' @param variable Character string specifying the covariate for which to plot the hazard ratio
#' @param form_splines Character string specifying the functional form: "quadratic" or "spline"
#' @param tt Numeric, the maximum follow-up time for the plot
#' @param IDnames Character string for patient identifier to display (optional, default=NULL)
#' @param legend Character vector for legend labels (optional, default=NULL)
#' @param xlab Character string for x-axis label (default="Follow-up Time (years)")
#' @param ylab Character string for y-axis label (default="Dynamic Hazard Ratio")
#' @param ylim Numeric vector of length 2 specifying y-axis limits (optional)
#' @param col Vector of colors for the mean, lower and upper bounds (default=c("black", "black", "black"))
#' @param lty Vector of line types for the mean, lower and upper bounds (default=c(1, 2, 2))
#' @param lwd Line width for the plot (default=6)
#' @param axis_lwd Line width for axes (default=1.6)
#' @param axis_tcl Tick length for axes (default=-0.4)
#' @param axis_padj_x Vertical adjustment for x-axis labels (default=-0.3)
#' @param axis_hadj_y Horizontal adjustment for y-axis labels (default=0.9)
#' @param axis_las Orientation of axis labels (0=parallel, 1=horizontal, 2=perpendicular, 3=vertical) (default=1)
#' @param title_font_lab Font for axis labels (default=7 for bold)
#' @param title_line_xlab Line position for x-axis label (default=2.6)
#' @param title_line_ylab Line position for y-axis label (default=2.6)
#' @param main_title_font Font for main title (default=7 for bold)
#' @param main_title_cex Character expansion for main title (default=1.5)
#' @param main_title_line Line position for main title (default=1)
#' @param abline_lty Line type for reference line at HR=1 (default=2 for dashed)
#' @param idnames_side Side to display ID names (1=bottom, 2=left, 3=top, 4=right) (default=3)
#' @param idnames_line Line position for ID names (default=1)
#' @param idnames_adj Adjustment for ID names (0=left, 1=right) (default=0)
#' @param idnames_cex Character expansion for ID names (default=1.5)
#' @param idnames_font Font for ID names (default=7 for bold)
#' @param legend_pos Position of legend (default="topright")
#' @param legend_lty Line type in legend (default=1)
#' @param legend_lwd Line width in legend (default=6)
#' @param legend_bty Box type for legend (default="n" for no box)
#' @param mar Numeric vector of length 4 specifying margins (bottom, left, top, right) (default=c(5, 5, 4, 2))
#' @param ... extra aguments; currently none is used.
#' @importFrom splines ns
#' @importFrom graphics matplot axis title abline mtext legend par
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
#' plot_dynamicHR(jointFit1,"GFR",form_splines="quadratic",tt=15,IDnames="b.",ylim=c(0,0.5),xlim=c(0,15),cex.axis=1.5,cex.lab=1.5)
plot_dynamicHR <- function(object, variable, form_splines, tt, IDnames = NULL,
                           legend = NULL, xlab = "Follow-up Time (years)",
                           ylab = "Dynamic Hazard Ratio", ylim = NULL,
                           col = c("black", "black", "black"),
                           lty = c(1, 2, 2), lwd = 6,
                           axis_lwd = 1.6,
                           axis_tcl = -0.4,
                           axis_padj_x = -0.3,
                           axis_hadj_y = 0.9,
                           axis_las = 1,
                           title_font_lab = 7,
                           title_line_xlab = 2.6,
                           title_line_ylab = 2.6,
                           main_title_font = 7,
                           main_title_cex = 1.5,
                           main_title_line = 1,
                           abline_lty = 2,
                           idnames_side = 3,
                           idnames_line = 1,
                           idnames_adj = 0,
                           idnames_cex = 1.5,
                           idnames_font = 7,
                           legend_pos = "topright",
                           legend_lty = 1,
                           legend_lwd = 6,
                           legend_bty = "n",
                           mar = c(5, 5, 4, 2),
                           ...) {

  # Load required packages
  if (!require("splines", quietly = TRUE)) {
    install.packages("splines")
    library(splines)
  }

  # Input validation
  if (missing(object) || is.null(object$mcmc$alphas)) {
    stop("object must be a joint model with MCMC samples in object$mcmc$alphas")
  }

  if (missing(variable) || !is.character(variable)) {
    stop("variable must be a character string specifying the covariate")
  }

  if (missing(form_splines) || !form_splines %in% c("quadratic", "spline")) {
    stop("form_splines must be either 'quadratic' or 'spline'")
  }

  if (missing(tt) || tt <= 0) {
    stop("tt must be a positive number")
  }

  # Generate time points
  x_times <- seq(0.001, tt, length = 501)

  # Create design matrix based on functional form
  if (form_splines == "quadratic") {
    x <- cbind(1, x_times, x_times^2)
  } else if (form_splines == "spline") {
    x <- cbind(1, splines::ns(x_times, 3))
  }

  # Extract MCMC samples for the specified variable
  x1 <- as.matrix(object$mcmc$alphas[1])
  names <- colnames(x1)
  matches <- grep(variable, names, value = TRUE)

  if (length(matches) == 0) {
    stop("Variable '", variable, "' not found in MCMC samples")
  }

  mcmc_alpha <- x1[, matches, drop = FALSE]

  # Calculate log hazard ratios
  log_hr <- x %*% t(mcmc_alpha)
  log_hr_mean <- rowMeans(log_hr)
  log_hr_low <- apply(log_hr, 1, quantile, probs = 0.025, na.rm = TRUE)
  log_hr_upp <- apply(log_hr, 1, quantile, probs = 0.975, na.rm = TRUE)

  # Set y-axis limits if not provided
  if (is.null(ylim)) {
    ylim <- c(min(exp(log_hr_low), na.rm = TRUE), max(exp(log_hr_upp), na.rm = TRUE))
  }

  # Save current graphical parameters
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Set graphical parameters
  par(mar = mar)

  # Create main plot
  matplot(x_times, cbind(exp(log_hr_mean), exp(log_hr_low), exp(log_hr_upp)),
          type = "l", col = col, lty = lty, lwd = lwd,
          xlab = "", ylab = "", ylim = ylim,
          main = "", xaxt = "n", yaxt = "n", bty = "n", ...)

  # Add custom axes with parameterized settings
  axis(1, las = axis_las, pos = 0, tcl = axis_tcl, padj = axis_padj_x, lwd = axis_lwd, ...)
  axis(2, las = axis_las, pos = ylim[1], tcl = axis_tcl, hadj = axis_hadj_y, lwd = axis_lwd, ...)

  # Add labels and titles with parameterized settings
  title(xlab = xlab, font.lab = title_font_lab, line = title_line_xlab, ...)
  title(ylab = ylab, font.lab = title_font_lab, line = title_line_ylab, ...)

  # Add reference line at HR = 1 with parameterized settings
  abline(h = 1, lty = abline_lty)

  # Add main title with parameterized settings
  title(main = list(variable, font = main_title_font, cex = main_title_cex),
        line = main_title_line)

  # Add ID names if provided with parameterized settings
  if (!is.null(IDnames)) {
    mtext(IDnames, side = idnames_side, line = idnames_line, adj = idnames_adj,
          cex = idnames_cex, font = idnames_font)
  }

  # Add legend if requested with parameterized settings
  if (!is.null(legend)) {
    legend(legend_pos,
           legend = legend,
           lty = legend_lty,
           lwd = legend_lwd,
           bty = legend_bty,
           ...)
  }

  # Return results invisibly
  invisible(list(
    times = x_times,
    hr_mean = exp(log_hr_mean),
    hr_lower = exp(log_hr_low),
    hr_upper = exp(log_hr_upp)
  ))
}
