#' Plot Conditional Hazard Ratios from Landmarking Models
#'
#' This function plots conditional hazard ratios over landmark time points from
#' landmarking models, with optional confidence intervals.
#' @usage plot_cHR(object,covars,conf_int=TRUE,IDlegend=NULL,legend=NULL,xlab="landmark time points",ylab="cHR",ylim,...)
#' @param object An objects inheriting from class "LMf" or "VS_LM"
#' @param covars Character vector of covariate names for which to plot hazard ratios
#' @param conf_int Logical indicating whether to plot confidence intervals (default=TRUE)
#' @param IDlegend Character string for identifier legend (optional, default=NULL)
#' @param legend Character vector for legend labels (optional, default=NULL)
#' @param xlab Character string for x-axis label (default="landmark time points")
#' @param ylab Character string for y-axis label (default="cHR")
#' @param ylim Numeric vector of length 2 specifying y-axis limits (optional)
#' @param IDlegend_side Side to display ID legend (1=bottom, 2=left, 3=top, 4=right) (default=3)
#' @param IDlegend_line Line position for ID legend (default=1)
#' @param IDlegend_adj Adjustment for ID legend (0=left, 1=right) (default=0)
#' @param IDlegend_cex Character expansion for ID legend (default=1.5)
#' @param IDlegend_font Font for ID legend (default=7 for bold)
#' @param main_title_font Font for main title (default=7 for bold)
#' @param main_title_cex Character expansion for main title (default=1.5)
#' @param main_title_line Line position for main title (default=1)
#' @param xlab_font Font for x-axis label (default=7 for bold)
#' @param xlab_cex Character expansion for x-axis label (default=1.5)
#' @param xlab_line Line position for x-axis label (default=2.6)
#' @param ylab_font Font for y-axis label (default=7 for bold)
#' @param ylab_cex Character expansion for y-axis label (default=1.5)
#' @param ylab_line Line position for y-axis label (default=2.6)
#' @param conf_int_lty Line type for confidence intervals (default=2 for dashed)
#' @param conf_int_lwd Line width for confidence intervals (default=6)
#' @param legend_pos Position of legend (default="topright")
#' @param legend_lty Line type in legend (default=1)
#' @param legend_lwd Line width in legend (default=6)
#' @param legend_bty Box type for legend (default="n" for no box)
#' @param abline_lwd Line width for reference line (default=1)
#' @param abline_lty Line type for reference line (default=3 for dotted)
#' @param abline_col Color for reference line (default=2 for red)
#' @param ... extra aguments; currently none is used.
#' @importFrom msm deltamethod
#' @importFrom stats vcov as.formula
#' @importFrom graphics plot lines mtext title legend abline
#' @export
#' @examples
#' library(DynamicModel)
#' data(renal)
#' fit<-LMf(renal,tw=list(sl=seq(0, 10, by=0.1),w=5),cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")),id="id",rtime="yearse",time="time",status="status",inter=TRUE)
#' a1<-VS_LM(fit,cov=list(fixed=c("age", "weight", "sex"),vary=c("GFR","proteinuria")))
#' plot_cHR(a1,c("age"),conf_int=TRUE,IDlegend="a.",legend="(Per 10 years)",xlab="landmark time points(years)",ylab="cHR")
plot_cHR <- function(object, covars, conf_int = TRUE, IDlegend = NULL,
                     legend = NULL, xlab = "landmark time points",
                     ylab = "cHR", ylim,
                     IDlegend_side = 3,
                     IDlegend_line = 1,
                     IDlegend_adj = 0,
                     IDlegend_cex = 1.5,
                     IDlegend_font = 7,
                     main_title_font = 7,
                     main_title_cex = 1.5,
                     main_title_line = 1,
                     xlab_font = 7,
                     xlab_cex = 1.5,
                     xlab_line = 2.6,
                     ylab_font = 7,
                     ylab_cex = 1.5,
                     ylab_line = 2.6,
                     conf_int_lty = 2,
                     conf_int_lwd = 6,
                     legend_pos = "topright",
                     legend_lty = 1,
                     legend_lwd = 6,
                     legend_bty = "n",
                     abline_lwd = 1,
                     abline_lty = 3,
                     abline_col = 2,
                     ...){
  # Load required packages
  if (!require("msm", quietly = TRUE)) {
    install.packages("msm")
    library(msm)
  }

  # Input validation
if (!inherits(object, "LMf") && !inherits(object, "VS_LM")) {
  stop("object must be of class 'LMf' or 'VS_LM'")
}

  if (missing(covars) || length(covars) == 0) {
    stop("covars must be specified with at least one covariate name")
  }

  # Extract model components
  fm <- object$Model
  bet <- fm$coefficients
  func_covars <- object$data$func_covars

  if (conf_int) {
    sig <- stats::vcov(fm)
  }

  sl <- object$tw$sl
  w <- object$tw$w

  # Create time sequence for plotting
  t <- seq(min(sl), max(sl), by = 0.1)

  # Initialize result list
  results <- list()

  # Save current graphical parameters
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Set graphical parameters for multiple plots if needed
  if (length(covars) > 1) {
    n_plots <- length(covars)
    n_col <- min(2, n_plots)
    n_row <- ceiling(n_plots / n_col)
    par(mfrow = c(n_row, n_col))
  }

  # Extract model components
  fm <- object$Model
  bet <- fm$coefficients
  func_covars <- object$data$func_covars

  if (conf_int) {
    sig <- stats::vcov(fm)
  }

  sl <- object$tw$sl
  w <- object$tw$w

  # Create time sequence for plotting
  t <- seq(min(sl), max(sl), by = 0.1)

  # Initialize result list
  results <- list()

  # Save current graphical parameters
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Set graphical parameters for multiple plots if needed
  if (length(covars) > 1) {
    n_plots <- length(covars)
    n_col <- min(2, n_plots)
    n_row <- ceiling(n_plots / n_col)
    par(mfrow = c(n_row, n_col))
  }

  # Process each covariate
  for (i in seq_along(covars)) {
    covar_name <- covars[i]

    # Find coefficients for this covariate
    idx <- startsWith(names(bet), covar_name)

    if (sum(idx) == 0) {
      warning(paste("Covariate", covar_name, "not found in model coefficients"))
      next
    }

    bet_var <- bet[idx]

    # Calculate hazard ratios over time
    coef <- sapply(t / max(sl) - min(sl), function(x) {
      sum(sapply(seq_along(bet_var), function(j) {
        var <- bet_var[j]
        name <- names(bet_var)[j]

        if (name == covar_name) {
          return(var)  # Main effect
        } else {
          # Extract index for functional form
          idx_num <- as.numeric(sub(".*\\D+", "\\1", name))
          if (!is.na(idx_num) && idx_num <= length(func_covars)) {
            return(func_covars[[idx_num]](x) * var)
          } else {
            return(0)
          }
        }
      }))
    })

    # Calculate standard errors if confidence intervals are requested
    se <- NULL
    if (conf_int) {
      find_se_log <- function(t, coefs, covar, func_covars) {
        if (!requireNamespace("msm", quietly = TRUE)) {
          stop("Package \"msm\" must be installed to use function find_se_log",
               call. = FALSE)
        }

        form <- "x1"
        if (length(coefs) > 1) {
          for (k in 2:length(coefs)) {
            form <- paste(form,
                          sprintf("%s * %f", paste("x", k, sep = ""),
                                  func_covars[[k - 1]](t)),
                          sep = " + ")
          }
        }
        form <- paste("~", form)
        se_val <- msm::deltamethod(stats::as.formula(form), coefs, covar)
        return(se_val)
      }

      se <- sapply(t / max(sl) - min(sl), find_se_log,
                   bet_var, sig[idx, idx], func_covars)
    }

    # Calculate hazard ratios and confidence bounds
    vari <- exp(coef)

    if (conf_int && !is.null(se)) {
      lower <- exp(coef - 1.96 * se)
      upper <- exp(coef + 1.96 * se)
    } else {
      lower <- upper <- NULL
    }

    # Set y-axis limits if not provided
    if (missing(ylim)) {
      if (conf_int && !is.null(lower) && !is.null(upper)) {
        ylim <- c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
      } else {
        ylim <- c(min(vari, na.rm = TRUE), max(vari, na.rm = TRUE))
      }
    }

    # Create plot
    plot(t, vari, type = "l", lwd = 6,
         xlim = c(min(sl), max(sl)), ylim = ylim,
         bty = "l", xlab = "", ylab = "", cex.axis = 1.5, ...)

    # Add ID legend if provided with parameterized settings
    if (!is.null(IDlegend)) {
      mtext(IDlegend,
            side = IDlegend_side,
            line = IDlegend_line,
            adj = IDlegend_adj,
            cex = IDlegend_cex,
            font = IDlegend_font)
    }

    # Add main title with parameterized settings
    title(main = list(covar_name,
                      font = main_title_font,
                      cex = main_title_cex),
          line = main_title_line)

    # Add axis labels with parameterized settings
    title(xlab = xlab,
          font.lab = xlab_font,
          cex.lab = xlab_cex,
          line = xlab_line)
    title(ylab = ylab,
          font.lab = ylab_font,
          cex.lab = ylab_cex,
          line = ylab_line)

    # Add confidence intervals if requested with parameterized settings
    if (conf_int && !is.null(lower) && !is.null(upper)) {
      lines(t, lower, type = "l", lty = conf_int_lty, lwd = conf_int_lwd)
      lines(t, upper, type = "l", lty = conf_int_lty, lwd = conf_int_lwd)
    }

    # Add legend if requested with parameterized settings
    if (!is.null(legend)) {
      legend(legend_pos,
             legend = legend,
             lty = legend_lty,
             lwd = legend_lwd,
             bty = legend_bty)
    }

    # Add reference line at HR = 1 with parameterized settings
    abline(h = 1, lwd = abline_lwd, lty = abline_lty, col = abline_col)

  }

  # Reset graphical parameters if we changed them
  if (length(covars) > 1) {
    par(mfrow = c(1, 1))
  }

}


