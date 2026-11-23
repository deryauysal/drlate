#' Doubly robust estimation of the LATE and the LATT
#'
#' \code{drlate} estimates the local average treatment effect (LATE) and the
#' local average treatment effect on the treated (LATT) with a binary
#' instrument. The treatment can be binary, continuous, or a count, and so can
#' the outcome. The causal parameter is the ratio of two treatment-effect-type
#' parameters with the instrument as the "treatment": the effect of Z on Y
#' (numerator) and the effect of Z on D (denominator).
#'
#' All parameters (instrument propensity score, outcome and treatment models,
#' numerator, denominator, and the LATE or LATT) are estimated jointly as one
#' just-identified system of moment conditions, so the standard errors account
#' for the estimation of all nuisance models. This is the R version of the
#' Stata command \code{drlate} and uses the same moment conditions.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{"ipwra"} (default)}{Inverse probability weighted regression
#'   adjustment. Outcome and treatment models are fitted within each instrument
#'   group with inverse probability weights; the parameters are averages of the
#'   predictions. Doubly robust. Always normalized.}
#'   \item{\code{"ipw"}}{Inverse probability weighting. No outcome or treatment
#'   model: \code{outcome} and \code{treatment} must have no covariates
#'   (e.g. \code{y ~ 1}).}
#'   \item{\code{"aipw"}}{Augmented inverse probability weighting. Outcome and
#'   treatment models are fitted without weights, and a weighted residual
#'   correction is added. Doubly robust.}
#'   \item{\code{"ra"}}{Regression adjustment. No instrument propensity score:
#'   \code{instrument} must have no covariates (e.g. \code{z ~ 1}).}
#' }
#'
#' @section Normalization:
#' For \code{"ipw"} and \code{"aipw"}, \code{normalized = TRUE} (the default)
#' uses weights that sum to one within each instrument group. For
#' \code{"ipwra"} and \code{"ra"} the option has no effect. With
#' \code{ivmodel = "ipt"} the weights are normalized by construction, so the
#' normalized and unnormalized estimators coincide.
#'
#' @section Overlap:
#' \code{drlate} stops with an error if an estimated instrument propensity
#' score lies outside \code{[pstolerance, 1 - pstolerance]}. The error object
#' has a component \code{violators} with the row numbers (in \code{data}) of
#' the observations concerned.
#'
#' @param outcome formula for the outcome model, e.g. \code{y ~ x1 + x2}.
#'   Use \code{y ~ 1} for no covariates.
#' @param treatment formula for the treatment model, e.g. \code{d ~ x1 + x2}.
#' @param instrument formula for the instrument propensity score model, e.g.
#'   \code{z ~ x1 + x2}. The instrument must be 0/1.
#' @param data a data frame.
#' @param estimand \code{"late"} (default) or \code{"latt"}.
#' @param method \code{"ipwra"} (default), \code{"ipw"}, \code{"aipw"}, or
#'   \code{"ra"}.
#' @param omodel outcome model: \code{"linear"} (default), \code{"logit"}
#'   (outcome must be 0/1), or \code{"poisson"} (outcome must be
#'   non-negative).
#' @param tmodel treatment model: \code{"logit"} (default; treatment must be
#'   0/1), \code{"linear"}, or \code{"poisson"}.
#' @param ivmodel instrument propensity score model, a logit estimated by
#'   \code{"logit"} (maximum likelihood, default), \code{"cbps"} (covariate
#'   balancing, Imai and Ratkovic, 2014; not available for the LATT), or
#'   \code{"ipt"} (inverse probability tilting, Graham, Pinto, and Egel, 2012).
#' @param normalized logical; normalized (default) or unnormalized moment
#'   conditions for \code{"ipw"} and \code{"aipw"}.
#' @param weights optional sampling weights (treated like Stata's pweights).
#' @param cluster optional cluster identifier for cluster-robust standard
#'   errors; a vector of length \code{nrow(data)} or the name of a variable in
#'   \code{data}.
#' @param pstolerance tolerance for the overlap check (default \code{1e-5}).
#'
#' @return An object of class \code{drlate} with components
#' \item{coefficients}{the causal estimate, the numerator, and the denominator}
#' \item{vcov}{their variance-covariance matrix}
#' \item{N, N_clust}{number of observations and clusters}
#' \item{dmeanz1, dmeanz0}{mean of the treatment for Z = 1 and Z = 0}
#' \item{theta, V}{all parameters of the joint system and their covariance}
#' and the settings used.
#'
#' @references
#' Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2022). "Doubly Robust
#' Estimation of Local Average Treatment Effects Using Inverse Probability
#' Weighted Regression Adjustment." arXiv:2208.01300.
#'
#' Graham, B. S., C. C. de Xavier Pinto, and D. Egel (2012). "Inverse
#' Probability Tilting for Moment Condition Models with Missing Data."
#' Review of Economic Studies 79(3), 1053-1079.
#'
#' Imai, K., and M. Ratkovic (2014). "Covariate Balancing Propensity Score."
#' Journal of the Royal Statistical Society Series B 76(1), 243-263.
#'
#' @examples
#' # SIPP data, as in the Stata help file (requires the haven package)
#' if (requireNamespace("haven", quietly = TRUE)) {
#'   sipp <- as.data.frame(haven::read_dta(system.file("extdata", "sipp.dta",
#'                                                     package = "drlate")))
#'   sipp <- sipp[!(is.na(sipp$kwage) | is.na(sipp$educ) | sipp$rsncode == 999), ]
#'   sipp$lwage <- log(sipp$kwage)
#'
#'   # LATE, default settings (IPWRA, logit instrument and treatment models,
#'   # linear outcome model)
#'   drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = sipp)
#'
#'   # LATE with IPT instrument propensity score
#'   drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = sipp,
#'          ivmodel = "ipt")
#'
#'   # LATT with IPT instrument propensity score
#'   drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = sipp,
#'          ivmodel = "ipt", estimand = "latt")
#'
#'   # IPW
#'   drlate(lwage ~ 1, nvstat ~ 1, rsncode ~ age_5, data = sipp, method = "ipw")
#'
#'   # AIPW, unnormalized
#'   drlate(lwage ~ age_5, nvstat ~ age_5, rsncode ~ age_5, data = sipp,
#'          method = "aipw", normalized = FALSE)
#'
#'   # Regression adjustment, Poisson model for the wage in levels
#'   drlate(kwage ~ age_5, nvstat ~ age_5, rsncode ~ 1, data = sipp,
#'          method = "ra", omodel = "poisson")
#' }
#'
#' @importFrom stats model.frame model.matrix model.response plogis pnorm qnorm
#' @export
drlate <- function(outcome, treatment, instrument, data,
                   estimand = c("late", "latt"),
                   method = c("ipwra", "ipw", "aipw", "ra"),
                   omodel = c("linear", "logit", "poisson"),
                   tmodel = c("logit", "linear", "poisson"),
                   ivmodel = c("logit", "cbps", "ipt"),
                   normalized = TRUE, weights = NULL, cluster = NULL,
                   pstolerance = 1e-5) {
  cl_call <- match.call()
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  omodel <- match.arg(omodel)
  tmodel <- match.arg(tmodel)
  ivmodel <- match.arg(ivmodel)
  if (!is.logical(normalized) || length(normalized) != 1 || is.na(normalized))
    stop("'normalized' must be TRUE or FALSE")
  if (!is.numeric(pstolerance) || length(pstolerance) != 1 ||
      pstolerance < 0 || pstolerance >= 1)
    stop("'pstolerance' must be >= 0 and < 1")
  if (ivmodel == "cbps" && estimand == "latt") stop("CBPS not possible for the LATT")
  if (!is.data.frame(data)) stop("'data' must be a data frame")

  # ---- sample: rows with no missing values in any model variable ----
  fl <- list(outcome = outcome, treatment = treatment, instrument = instrument)
  for (nm in names(fl)) {
    if (!inherits(fl[[nm]], "formula") || length(fl[[nm]]) != 3)
      stop(sprintf("'%s' must be a two-sided formula, e.g. y ~ x1 + x2", nm))
  }
  allvars <- unique(unlist(lapply(fl, all.vars)))
  missing_vars <- setdiff(allvars, names(data))
  if (length(missing_vars)) stop("variables not found in data: ", paste(missing_vars, collapse = ", "))
  if (!is.null(cluster) && is.character(cluster) && length(cluster) == 1) {
    if (!cluster %in% names(data)) stop("cluster variable not found in data")
    cluster <- data[[cluster]]
  }
  if (!is.null(cluster) && length(cluster) != nrow(data))
    stop("'cluster' must have length nrow(data)")
  if (!is.null(weights) && length(weights) != nrow(data))
    stop("'weights' must have length nrow(data)")
  ok <- stats::complete.cases(data[, allvars, drop = FALSE])
  if (!is.null(weights)) ok <- ok & !is.na(weights)
  if (!is.null(cluster)) ok <- ok & !is.na(cluster)
  rows <- which(ok)
  d <- data[rows, , drop = FALSE]
  w <- if (is.null(weights)) rep(1, length(rows)) else as.numeric(weights[rows])
  if (any(w < 0)) stop("negative weights not allowed")
  cl <- if (is.null(cluster)) NULL else cluster[rows]

  get_resp <- function(f) {
    mf <- stats::model.frame(f, d, na.action = NULL)
    as.numeric(stats::model.response(mf))
  }
  get_X <- function(f) {
    X <- stats::model.matrix(stats::delete.response(stats::terms(f)), d)
    if (!"(Intercept)" %in% colnames(X)) stop("models must include an intercept")
    .drlate_std(X, w)
  }
  Y <- get_resp(outcome)
  D <- get_resp(treatment)
  Z <- get_resp(instrument)
  Xy <- get_X(outcome)
  Xd <- get_X(treatment)
  Xz <- get_X(instrument)
  N <- length(Y)

  # ---- input checks (as in the Stata version) ----
  if (!all(sort(unique(Z)) == c(0, 1))) stop("instrument must be 0/1 with both values present")
  if (method == "ipw" && (ncol(Xy) > 1 || ncol(Xd) > 1))
    stop("outcome and treatment models must have no covariates with method = 'ipw' (use y ~ 1, d ~ 1)")
  if (method == "ra" && ncol(Xz) > 1)
    stop("instrument model must have no covariates with method = 'ra' (use z ~ 1)")
  if (omodel == "logit" && !all(Y %in% c(0, 1))) stop("outcome must be 0/1 for omodel = 'logit'")
  if (omodel == "poisson" && any(Y < 0)) stop("outcome must be non-negative for omodel = 'poisson'")
  if (tmodel == "logit" && !all(D %in% c(0, 1))) stop("treatment must be 0/1 for tmodel = 'logit'")
  if (tmodel == "poisson" && any(D < 0)) stop("treatment must be non-negative for tmodel = 'poisson'")
  dmeanz1 <- mean(D[Z == 1])
  dmeanz0 <- mean(D[Z == 0])
  if (dmeanz1 == dmeanz0 && dmeanz1 %in% c(0, 1)) stop("treatment does not vary with the instrument")

  sys <- .drlate_system()

  # ---- instrument propensity score ----
  ps_names <- character(0)
  p1fun <- p0fun <- NULL
  if (method != "ra") {
    start <- .drlate_glm(Xz, Z, w, "logit", "logit model for the instrument")
    mom_logit <- function(g) (Z - plogis(as.vector(Xz %*% g))) * Xz
    mom_cbps <- function(g) { p <- plogis(as.vector(Xz %*% g)); (Z / p - (1 - Z) / (1 - p)) * Xz }
    mom_ipt1 <- function(g) { p <- plogis(as.vector(Xz %*% g)); (Z / p - 1) * Xz }
    mom_ipt0 <- function(g) { p <- plogis(as.vector(Xz %*% g)); ((1 - Z) / (1 - p) - 1) * Xz }
    solve_ps <- function(mom, what) {
      .drlate_solve(function(g) colSums(w * mom(g)) / sum(w), start, what)
    }
    if (ivmodel == "logit") {
      g <- start
      sys$add("gamma", g, function(par) mom_logit(par$gamma))
      p1fun <- p0fun <- function(par) plogis(as.vector(Xz %*% par$gamma))
    } else if (ivmodel == "cbps") {
      g <- solve_ps(mom_cbps, "CBPS estimation")
      sys$add("gamma", g, function(par) mom_cbps(par$gamma))
      p1fun <- p0fun <- function(par) plogis(as.vector(Xz %*% par$gamma))
    } else if (estimand == "late") {
      g1 <- solve_ps(mom_ipt1, "IPT estimation (Z=1)")
      g0 <- solve_ps(mom_ipt0, "IPT estimation (Z=0)")
      sys$add("gamma1", g1, function(par) mom_ipt1(par$gamma1))
      sys$add("gamma0", g0, function(par) mom_ipt0(par$gamma0))
      p1fun <- function(par) plogis(as.vector(Xz %*% par$gamma1))
      p0fun <- function(par) plogis(as.vector(Xz %*% par$gamma0))
    } else {
      # LATT with IPT: odds-weighted controls balance the treated (= IPT for Z=0)
      g <- solve_ps(mom_ipt0, "IPT estimation")
      sys$add("gamma", g, function(par) mom_ipt0(par$gamma))
      p1fun <- p0fun <- function(par) plogis(as.vector(Xz %*% par$gamma))
    }
    # overlap check
    par0 <- sys$values
    ips <- ifelse(Z == 1, p1fun(par0), p0fun(par0))
    viol <- ips < pstolerance | ips > 1 - pstolerance
    if (any(viol)) {
      msg <- sprintf(paste0("%d observation%s violate%s the overlap assumption ",
                            "(instrument propensity score outside [%g, 1-%g]); ",
                            "see the 'violators' component of this error for their row numbers"),
                     sum(viol), if (sum(viol) > 1) "s" else "", if (sum(viol) > 1) "" else "s",
                     pstolerance, pstolerance)
      cond <- structure(class = c("drlate_overlap", "error", "condition"),
                        list(message = msg, call = cl_call,
                             violators = as.integer(rownames(data)[rows[viol]])))
      stop(cond)
    }
  }
  par0 <- sys$values

  # ---- one potential-outcome-type mean for variable v ----
  # Adds the blocks for E[v(k)] (LATE) or E[v(k)|Z=1] (LATT) and returns a
  # function of the parameter list giving its value (or a constant).
  add_mean <- function(v, X, model, k, label, const = NULL) {
    nm <- paste0("mu_", label, k)
    if (!is.null(const)) return(function(par) const)
    zk <- if (k == 1) Z else 1 - Z
    if (estimand == "late") {
      om <- if (k == 1) function(par) Z / p1fun(par) else function(par) (1 - Z) / (1 - p0fun(par))
      wsum <- sum(w)
      if (method == "ra" || method == "ipwra") {
        bw <- if (method == "ra") zk else om(par0)
        b <- .drlate_glm(X, v, w * bw, model, sprintf("%s model (Z=%d)", label, k))
        bn <- paste0("b_", label, k)
        if (method == "ra") {
          sys$add(bn, b, function(par) zk * (v - .drlate_mean(X, par[[bn]], model)) * X)
        } else {
          sys$add(bn, b, function(par) om(par) * (v - .drlate_mean(X, par[[bn]], model)) * X)
        }
        mu <- sum(w * .drlate_mean(X, b, model)) / wsum
        sys$add(nm, mu, function(par) matrix(.drlate_mean(X, par[[bn]], model) - par[[nm]]))
      } else if (method == "ipw") {
        o <- om(par0)
        if (normalized) {
          mu <- sum(w * o * v) / sum(w * o)
          sys$add(nm, mu, function(par) matrix(om(par) * (v - par[[nm]])))
        } else {
          mu <- sum(w * o * v) / wsum
          sys$add(nm, mu, function(par) matrix(om(par) * v - par[[nm]]))
        }
      } else { # aipw
        b <- .drlate_glm(X, v, w * zk, model, sprintf("%s model (Z=%d)", label, k))
        bn <- paste0("b_", label, k)
        sys$add(bn, b, function(par) zk * (v - .drlate_mean(X, par[[bn]], model)) * X)
        m <- .drlate_mean(X, b, model)
        o <- om(par0)
        rn <- paste0("r_", label, k)
        if (normalized) {
          r <- sum(w * o * (v - m)) / sum(w * o)
          sys$add(rn, r, function(par) matrix(om(par) * (v - .drlate_mean(X, par[[bn]], model) - par[[rn]])))
        } else {
          r <- sum(w * o * (v - m)) / wsum
          sys$add(rn, r, function(par) matrix(om(par) * (v - .drlate_mean(X, par[[bn]], model)) - par[[rn]]))
        }
        mu <- sum(w * m) / wsum + r
        sys$add(nm, mu, function(par) matrix(.drlate_mean(X, par[[bn]], model) + par[[rn]] - par[[nm]]))
      }
    } else { # latt
      n1 <- sum(w * Z)
      if (k == 1) {
        mu <- sum(w * Z * v) / n1
        sys$add(nm, mu, function(par) matrix(Z * (v - par[[nm]])))
        return(function(par) par[[nm]])
      }
      odds <- function(par) { p <- p0fun(par); p / (1 - p) }
      if (method == "ra" || method == "ipwra") {
        bw <- if (method == "ra") (1 - Z) else (1 - Z) * odds(par0)
        b <- .drlate_glm(X, v, w * bw, model, sprintf("%s model (Z=0)", label))
        bn <- paste0("b_", label, 0)
        if (method == "ra") {
          sys$add(bn, b, function(par) (1 - Z) * (v - .drlate_mean(X, par[[bn]], model)) * X)
        } else {
          sys$add(bn, b, function(par) (1 - Z) * odds(par) * (v - .drlate_mean(X, par[[bn]], model)) * X)
        }
        mu <- sum(w * Z * .drlate_mean(X, b, model)) / n1
        sys$add(nm, mu, function(par) matrix(Z * (.drlate_mean(X, par[[bn]], model) - par[[nm]])))
      } else if (method == "ipw") {
        o <- (1 - Z) * odds(par0)
        if (normalized) {
          mu <- sum(w * o * v) / sum(w * o)
          sys$add(nm, mu, function(par) matrix((1 - Z) * odds(par) * (v - par[[nm]])))
        } else {
          mu <- sum(w * o * v) / n1
          sys$add(nm, mu, function(par) matrix((1 - Z) * odds(par) * v - Z * par[[nm]]))
        }
      } else { # aipw
        b <- .drlate_glm(X, v, w * (1 - Z), model, sprintf("%s model (Z=0)", label))
        bn <- paste0("b_", label, 0)
        sys$add(bn, b, function(par) (1 - Z) * (v - .drlate_mean(X, par[[bn]], model)) * X)
        m <- .drlate_mean(X, b, model)
        o <- (1 - Z) * odds(par0)
        rn <- paste0("r_", label, 0)
        if (normalized) {
          r <- sum(w * o * (v - m)) / sum(w * o)
          sys$add(rn, r, function(par) matrix((1 - Z) * odds(par) * (v - .drlate_mean(X, par[[bn]], model) - par[[rn]])))
        } else {
          r <- sum(w * o * (v - m)) / n1
          sys$add(rn, r, function(par) matrix((1 - Z) * odds(par) * (v - .drlate_mean(X, par[[bn]], model)) - Z * par[[rn]]))
        }
        mu <- sum(w * Z * m) / n1 + r
        sys$add(nm, mu, function(par) matrix(Z * (.drlate_mean(X, par[[bn]], model) + par[[rn]] - par[[nm]])))
      }
    }
    function(par) par[[nm]]
  }

  muY1 <- add_mean(Y, Xy, omodel, 1, "y")
  muY0 <- add_mean(Y, Xy, omodel, 0, "y")
  muD1 <- add_mean(D, Xd, tmodel, 1, "d", const = if (dmeanz1 %in% c(0, 1)) dmeanz1 else NULL)
  muD0 <- add_mean(D, Xd, tmodel, 0, "d", const = if (dmeanz0 %in% c(0, 1)) dmeanz0 else NULL)

  pv <- sys$values
  num <- muY1(pv) - muY0(pv)
  den <- muD1(pv) - muD0(pv)
  if (abs(den) < 1e-12) stop("zero denominator: the LATE is not defined")
  sys$add("num", num, function(par) matrix(muY1(par) - muY0(par) - par$num, N, 1))
  sys$add("denom", den, function(par) matrix(muD1(par) - muD0(par) - par$denom, N, 1))
  sys$add("late", num / den, function(par) matrix(par$late - par$num / par$denom, N, 1))

  # ---- joint system and sandwich variance ----
  skel <- sys$values
  theta <- unlist(skel)
  fns <- sys$fns
  psi <- function(th) {
    par <- utils::relist(th, skel)
    do.call(cbind, lapply(fns, function(f) f(par)))
  }
  chk <- colSums(w * psi(theta)) / sum(w)
  if (max(abs(chk)) > 1e-6) stop("internal error: moment conditions not satisfied at the estimates")
  V <- .drlate_vcov(psi, theta, w, cl)

  k <- length(theta)
  idx <- c(k, k - 2, k - 1)
  b <- theta[idx]
  vc <- V[idx, idx]
  lab <- if (estimand == "late") c("LATE: D on Y", "ATE: Z on Y", "ATE: Z on D")
         else c("LATT: D on Y", "ATT: Z on Y", "ATT: Z on D")
  names(b) <- lab
  dimnames(vc) <- list(lab, lab)

  structure(list(
    coefficients = b, vcov = vc, N = N,
    N_clust = if (is.null(cl)) NULL else length(unique(cl)),
    dmeanz1 = dmeanz1, dmeanz0 = dmeanz0,
    estimand = estimand, method = method, omodel = omodel, tmodel = tmodel,
    ivmodel = ivmodel, normalized = normalized || ivmodel == "ipt",
    weights = !is.null(weights), cluster = !is.null(cl),
    theta = theta, V = V, call = cl_call
  ), class = "drlate")
}
