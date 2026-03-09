#' Generate Publication-Ready QARDL Tables
#'
#' Creates formatted tables suitable for academic publications from
#' QARDL estimation results.
#'
#' @param x An object of class \code{"qardl"}.
#' @param type Character. Type of table: \code{"latex"}, \code{"html"},
#'   or \code{"text"}. Default is \code{"text"}.
#' @param include Character vector. Which parameters to include:
#'   \code{"beta"}, \code{"phi"}, \code{"gamma"}, \code{"rho"}.
#'   Default is \code{c("beta", "gamma")}.
#' @param stars Logical. Include significance stars. Default is \code{TRUE}.
#' @param digits Integer. Number of decimal places. Default is 3.
#' @param caption Character. Table caption. Default is \code{NULL}.
#' @param label Character. LaTeX label. Default is \code{NULL}.
#'
#' @return Character string containing the formatted table.
#'
#' @examples
#' data(qardl_sim)
#' fit <- qardl(y ~ x1 + x2, data = qardl_sim, tau = c(0.25, 0.50, 0.75), p = 2, q = 2)
#' cat(qardl_table(fit, type = "text"))
#'
#' @export
qardl_table <- function(x, type = c("text", "latex", "html"),
                        include = c("beta", "gamma"),
                        stars = TRUE, digits = 3,
                        caption = NULL, label = NULL) {

  if (!inherits(x, "qardl")) {
    stop("'x' must be of class 'qardl'", call. = FALSE)
  }

  type <- match.arg(type)
  ntau <- length(x$tau)

  # Build table content
  lines <- character()

  if (type == "text") {
    lines <- c(lines, build_text_table(x, include, stars, digits))
  } else if (type == "latex") {
    lines <- c(lines, build_latex_table(x, include, stars, digits, caption, label))
  } else {
    lines <- c(lines, build_html_table(x, include, stars, digits, caption))
  }

  result <- paste(lines, collapse = "\n")
  return(result)
}


#' Build Text Table
#' @keywords internal
build_text_table <- function(x, include, stars, digits) {

  lines <- character()
  ntau <- length(x$tau)
  fmt <- paste0("%.", digits, "f")

  # Header
  header <- sprintf("%-15s", "")
  for (t_idx in seq_len(ntau)) {
    header <- paste0(header, sprintf("%15s", paste0("tau=", x$tau[t_idx])))
  }
  lines <- c(lines, header)
  lines <- c(lines, strrep("-", 15 + 15 * ntau))

  # Beta
  if ("beta" %in% include) {
    lines <- c(lines, "Long-Run (beta)")
    for (i in seq_len(x$k)) {
      row <- sprintf("  %-13s", x$indepvars[i])
      for (t_idx in seq_len(ntau)) {
        val <- x$beta[i, t_idx]
        se <- x$beta_se[i, t_idx]
        pval <- 2 * (1 - stats::pnorm(abs(val / se)))
        star_str <- ifelse(stars, get_stars(pval), "")
        row <- paste0(row, sprintf("%14s%s", sprintf(fmt, val), star_str))
      }
      lines <- c(lines, row)

      # SE row
      se_row <- sprintf("  %-13s", "")
      for (t_idx in seq_len(ntau)) {
        se_row <- paste0(se_row, sprintf("(%13s)", sprintf(fmt, x$beta_se[i, t_idx])))
      }
      lines <- c(lines, se_row)
    }
  }

  # Gamma
  if ("gamma" %in% include) {
    lines <- c(lines, "")
    lines <- c(lines, "Short-Run (gamma)")
    for (i in seq_len(x$k)) {
      row <- sprintf("  %-13s", x$indepvars[i])
      for (t_idx in seq_len(ntau)) {
        val <- x$gamma[i, t_idx]
        se <- x$gamma_se[i, t_idx]
        pval <- 2 * (1 - stats::pnorm(abs(val / se)))
        star_str <- ifelse(stars, get_stars(pval), "")
        row <- paste0(row, sprintf("%14s%s", sprintf(fmt, val), star_str))
      }
      lines <- c(lines, row)

      se_row <- sprintf("  %-13s", "")
      for (t_idx in seq_len(ntau)) {
        se_row <- paste0(se_row, sprintf("(%13s)", sprintf(fmt, x$gamma_se[i, t_idx])))
      }
      lines <- c(lines, se_row)
    }
  }

  # Phi
  if ("phi" %in% include) {
    lines <- c(lines, "")
    lines <- c(lines, "AR Coefficients (phi)")
    for (i in seq_len(x$p)) {
      row <- sprintf("  %-13s", paste0("Lag ", i))
      for (t_idx in seq_len(ntau)) {
        val <- x$phi[i, t_idx]
        se <- x$phi_se[i, t_idx]
        pval <- 2 * (1 - stats::pnorm(abs(val / se)))
        star_str <- ifelse(stars, get_stars(pval), "")
        row <- paste0(row, sprintf("%14s%s", sprintf(fmt, val), star_str))
      }
      lines <- c(lines, row)

      se_row <- sprintf("  %-13s", "")
      for (t_idx in seq_len(ntau)) {
        se_row <- paste0(se_row, sprintf("(%13s)", sprintf(fmt, x$phi_se[i, t_idx])))
      }
      lines <- c(lines, se_row)
    }
  }

  # Rho
  if ("rho" %in% include) {
    lines <- c(lines, "")
    lines <- c(lines, "ECM Coefficient (rho)")
    row <- sprintf("  %-13s", "rho")
    for (t_idx in seq_len(ntau)) {
      val <- x$rho[t_idx]
      se <- x$rho_se[t_idx]
      pval <- 2 * (1 - stats::pnorm(abs(val / se)))
      star_str <- ifelse(stars, get_stars(pval), "")
      row <- paste0(row, sprintf("%14s%s", sprintf(fmt, val), star_str))
    }
    lines <- c(lines, row)

    se_row <- sprintf("  %-13s", "")
    for (t_idx in seq_len(ntau)) {
      se_row <- paste0(se_row, sprintf("(%13s)", sprintf(fmt, x$rho_se[t_idx])))
    }
    lines <- c(lines, se_row)
  }

  lines <- c(lines, strrep("-", 15 + 15 * ntau))
  lines <- c(lines, sprintf("Observations: %d", x$nobs))
  if (stars) {
    lines <- c(lines, "*** p<0.01, ** p<0.05, * p<0.10")
  }

  return(lines)
}


#' Build LaTeX Table
#' @keywords internal
build_latex_table <- function(x, include, stars, digits, caption, label) {

  lines <- character()
  ntau <- length(x$tau)
  fmt <- paste0("%.", digits, "f")

  # Begin table
  col_spec <- paste0("l", paste(rep("c", ntau), collapse = ""))
  lines <- c(lines, "\\begin{table}[htbp]")
  lines <- c(lines, "\\centering")
  if (!is.null(caption)) {
    lines <- c(lines, sprintf("\\caption{%s}", caption))
  }
  if (!is.null(label)) {
    lines <- c(lines, sprintf("\\label{%s}", label))
  }
  lines <- c(lines, sprintf("\\begin{tabular}{%s}", col_spec))
  lines <- c(lines, "\\hline\\hline")

  # Header
  header <- ""
  for (t_idx in seq_len(ntau)) {
    header <- paste0(header, sprintf(" & $\\tau = %.2f$", x$tau[t_idx]))
  }
  header <- paste0(header, " \\\\")
  lines <- c(lines, header)
  lines <- c(lines, "\\hline")

  # Beta
  if ("beta" %in% include) {
    lines <- c(lines, "\\multicolumn{", ntau + 1, "}{l}{\\textit{Long-Run Parameters ($\\beta$)}} \\\\")
    for (i in seq_len(x$k)) {
      row <- x$indepvars[i]
      for (t_idx in seq_len(ntau)) {
        val <- x$beta[i, t_idx]
        se <- x$beta_se[i, t_idx]
        pval <- 2 * (1 - stats::pnorm(abs(val / se)))
        star_str <- ifelse(stars, get_stars_latex(pval), "")
        row <- paste0(row, sprintf(" & %s%s", sprintf(fmt, val), star_str))
      }
      row <- paste0(row, " \\\\")
      lines <- c(lines, row)

      # SE row
      se_row <- ""
      for (t_idx in seq_len(ntau)) {
        se_row <- paste0(se_row, sprintf(" & (%s)", sprintf(fmt, x$beta_se[i, t_idx])))
      }
      se_row <- paste0(se_row, " \\\\")
      lines <- c(lines, se_row)
    }
  }

  # Gamma
  if ("gamma" %in% include) {
    lines <- c(lines, "\\multicolumn{", ntau + 1, "}{l}{\\textit{Short-Run Parameters ($\\gamma$)}} \\\\")
    for (i in seq_len(x$k)) {
      row <- x$indepvars[i]
      for (t_idx in seq_len(ntau)) {
        val <- x$gamma[i, t_idx]
        se <- x$gamma_se[i, t_idx]
        pval <- 2 * (1 - stats::pnorm(abs(val / se)))
        star_str <- ifelse(stars, get_stars_latex(pval), "")
        row <- paste0(row, sprintf(" & %s%s", sprintf(fmt, val), star_str))
      }
      row <- paste0(row, " \\\\")
      lines <- c(lines, row)

      se_row <- ""
      for (t_idx in seq_len(ntau)) {
        se_row <- paste0(se_row, sprintf(" & (%s)", sprintf(fmt, x$gamma_se[i, t_idx])))
      }
      se_row <- paste0(se_row, " \\\\")
      lines <- c(lines, se_row)
    }
  }

  lines <- c(lines, "\\hline")
  lines <- c(lines, sprintf("Observations & \\multicolumn{%d}{c}{%d} \\\\", ntau, x$nobs))
  lines <- c(lines, "\\hline\\hline")
  lines <- c(lines, "\\end{tabular}")

  if (stars) {
    lines <- c(lines, "\\begin{tablenotes}")
    lines <- c(lines, "\\small")
    lines <- c(lines, "\\item Standard errors in parentheses. $^{***}$ p$<$0.01, $^{**}$ p$<$0.05, $^{*}$ p$<$0.10")
    lines <- c(lines, "\\end{tablenotes}")
  }

  lines <- c(lines, "\\end{table}")

  return(lines)
}


#' Build HTML Table
#' @keywords internal
build_html_table <- function(x, include, stars, digits, caption) {

  lines <- character()
  ntau <- length(x$tau)
  fmt <- paste0("%.", digits, "f")

  lines <- c(lines, "<table class='qardl-table'>")
  if (!is.null(caption)) {
    lines <- c(lines, sprintf("<caption>%s</caption>", caption))
  }

  # Header
  lines <- c(lines, "<thead><tr>")
  lines <- c(lines, "<th></th>")
  for (t_idx in seq_len(ntau)) {
    lines <- c(lines, sprintf("<th>&tau; = %.2f</th>", x$tau[t_idx]))
  }
  lines <- c(lines, "</tr></thead>")

  lines <- c(lines, "<tbody>")

  # Beta
  if ("beta" %in% include) {
    lines <- c(lines, sprintf("<tr><td colspan='%d'><em>Long-Run (&beta;)</em></td></tr>", ntau + 1))
    for (i in seq_len(x$k)) {
      lines <- c(lines, "<tr>")
      lines <- c(lines, sprintf("<td>%s</td>", x$indepvars[i]))
      for (t_idx in seq_len(ntau)) {
        val <- x$beta[i, t_idx]
        se <- x$beta_se[i, t_idx]
        pval <- 2 * (1 - stats::pnorm(abs(val / se)))
        star_str <- ifelse(stars, get_stars_html(pval), "")
        lines <- c(lines, sprintf("<td>%s%s<br/>(%s)</td>",
                                   sprintf(fmt, val), star_str,
                                   sprintf(fmt, se)))
      }
      lines <- c(lines, "</tr>")
    }
  }

  # Gamma
  if ("gamma" %in% include) {
    lines <- c(lines, sprintf("<tr><td colspan='%d'><em>Short-Run (&gamma;)</em></td></tr>", ntau + 1))
    for (i in seq_len(x$k)) {
      lines <- c(lines, "<tr>")
      lines <- c(lines, sprintf("<td>%s</td>", x$indepvars[i]))
      for (t_idx in seq_len(ntau)) {
        val <- x$gamma[i, t_idx]
        se <- x$gamma_se[i, t_idx]
        pval <- 2 * (1 - stats::pnorm(abs(val / se)))
        star_str <- ifelse(stars, get_stars_html(pval), "")
        lines <- c(lines, sprintf("<td>%s%s<br/>(%s)</td>",
                                   sprintf(fmt, val), star_str,
                                   sprintf(fmt, se)))
      }
      lines <- c(lines, "</tr>")
    }
  }

  lines <- c(lines, "</tbody>")
  lines <- c(lines, sprintf("<tfoot><tr><td colspan='%d'>Observations: %d</td></tr></tfoot>",
                             ntau + 1, x$nobs))
  lines <- c(lines, "</table>")

  return(lines)
}


#' Get LaTeX Stars
#' @keywords internal
get_stars_latex <- function(pval) {
  if (is.na(pval)) return("")
  if (pval < 0.01) return("$^{***}$")
  if (pval < 0.05) return("$^{**}$")
  if (pval < 0.10) return("$^{*}$")
  return("")
}


#' Get HTML Stars
#' @keywords internal
get_stars_html <- function(pval) {
  if (is.na(pval)) return("")
  if (pval < 0.01) return("<sup>***</sup>")
  if (pval < 0.05) return("<sup>**</sup>")
  if (pval < 0.10) return("<sup>*</sup>")
  return("")
}


#' Predict Method for QARDL
#'
#' Generate predictions from a fitted QARDL model.
#'
#' @param object An object of class \code{"qardl"}.
#' @param newdata Optional data frame for prediction. If \code{NULL},
#'   uses the original data.
#' @param tau Quantile(s) for prediction. Default uses all fitted quantiles.
#' @param ... Additional arguments (unused).
#'
#' @return Matrix of predicted quantiles.
#'
#' @export
predict.qardl <- function(object, newdata = NULL, tau = NULL, ...) {

  if (is.null(tau)) {
    tau <- object$tau
  }

  if (is.null(newdata)) {
    # Return fitted values from quantreg fits
    ntau <- length(object$tau)
    fitted_vals <- matrix(NA, nrow = object$nobs, ncol = ntau)
    colnames(fitted_vals) <- paste0("tau=", object$tau)

    for (t_idx in seq_along(object$tau)) {
      if (!is.null(object$qr_fits[[t_idx]])) {
        fitted_vals[, t_idx] <- fitted(object$qr_fits[[t_idx]])
      }
    }

    return(fitted_vals)
  }

  # For new data, construct design matrix and predict
  message("Prediction with new data requires reconstructing the model design. ",
          "Using fitted values for original data.")
  return(predict.qardl(object, newdata = NULL, tau = tau, ...))
}
