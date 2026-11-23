#' @export
print.drlate <- function(x, digits = 7, ...) {
  title <- if (x$estimand == "late") "Local average treatment effect"
           else "Local average treatment effect on the treated"
  meth <- c(ipwra = "IPWRA", ipw = "IPW", aipw = "AIPW", ra = "RA")[[x$method]]
  if (x$method %in% c("ipw", "aipw"))
    meth <- paste0(meth, if (x$normalized) " (normalized)" else " (unnormalized)")
  mname <- function(m) c(linear = "linear", logit = "logit", poisson = "Poisson")[[m]]
  om <- if (x$method == "ipw") "weighted mean" else mname(x$omodel)
  tm <- if (x$method == "ipw") "weighted mean" else mname(x$tmodel)
  im <- if (x$method == "ra") "none" else
    c(logit = "logit (MLE)", cbps = "logit (CBPS)", ipt = "logit (IPT)")[[x$ivmodel]]

  cat("\n", title, sep = "")
  cat(sprintf("%*s%d\n", max(2, 60 - nchar(title)), "Number of obs = ", x$N))
  if (!is.null(x$N_clust)) cat(sprintf("%60s%d\n", "Number of clusters = ", x$N_clust))
  cat(sprintf("%-17s: %s\n", "Estimator", meth))
  cat(sprintf("%-17s: %s\n", "Outcome model", om))
  cat(sprintf("%-17s: %s\n", "Treatment model", tm))
  cat(sprintf("%-17s: %s\n", "Instrument model", im))
  cat("\n")
  print(summary_table(x, digits = digits))
  se_type <- if (x$cluster) "Cluster-robust" else "Robust"
  cat(sprintf("\n%s standard errors.\n", se_type))
  invisible(x)
}

summary_table <- function(x, digits = 7) {
  b <- x$coefficients
  se <- sqrt(diag(x$vcov))
  z <- b / se
  p <- 2 * stats::pnorm(-abs(z))
  q <- stats::qnorm(0.975)
  tab <- data.frame(
    Coefficient = signif(b, digits), `Std. err.` = signif(se, digits),
    z = round(z, 2), `P>|z|` = round(p, 3),
    `[95% conf.` = signif(b - q * se, digits), `interval]` = signif(b + q * se, digits),
    check.names = FALSE
  )
  rownames(tab) <- names(b)
  tab
}

#' @export
coef.drlate <- function(object, ...) object$coefficients

#' @export
vcov.drlate <- function(object, ...) object$vcov
