# Internal helpers for drlate (not exported).

# Canonical-link GLM family for a model name.
.drlate_family <- function(model) {
  switch(model,
    linear  = stats::gaussian(),
    logit   = stats::quasibinomial(),
    poisson = stats::quasipoisson(),
    stop(sprintf("model '%s' not recognized; use 'linear', 'logit', or 'poisson'", model))
  )
}

# Fit a canonical-link GLM with prior weights; return the coefficient vector.
# For canonical links the estimating equations are sum_i wt_i (y_i - mu_i) x_i = 0,
# which is what Stata's regress / logit / poisson solve with [pw] weights.
.drlate_glm <- function(X, y, wt, model, what) {
  keep <- wt > 0
  fit <- suppressWarnings(stats::glm.fit(X[keep, , drop = FALSE], y[keep],
                                         weights = wt[keep], family = .drlate_family(model)))
  b <- fit$coefficients
  if (any(is.na(b))) {
    stop(sprintf("%s: collinear covariates (coefficient%s %s not identified)", what,
                 if (sum(is.na(b)) > 1) "s" else "", paste(names(b)[is.na(b)], collapse = ", ")))
  }
  if (!isTRUE(fit$converged)) stop(sprintf("%s: convergence not achieved", what))
  b
}

# Mean function of the canonical-link GLM.
.drlate_mean <- function(X, b, model) {
  eta <- as.vector(X %*% b)
  switch(model, linear = eta, logit = stats::plogis(eta), poisson = exp(eta))
}

# Newton's method with step halving for a just-identified system f(x) = 0.
.drlate_solve <- function(f, start, what, tol = 1e-10, maxit = 200) {
  x <- start
  fx <- f(x)
  for (it in seq_len(maxit)) {
    if (all(is.finite(fx)) && max(abs(fx)) < tol) return(x)
    J <- numDeriv::jacobian(f, x)
    step <- tryCatch(solve(J, fx), error = function(e) NULL)
    if (is.null(step)) break
    t <- 1
    repeat {
      xn <- x - t * step
      fn <- f(xn)
      if (all(is.finite(fn)) && sum(fn^2) < sum(fx^2)) break
      t <- t / 2
      if (t < 1e-10) break
    }
    if (t < 1e-10) break
    x <- xn
    fx <- fn
  }
  if (all(is.finite(fx)) && max(abs(fx)) < 1e-7) return(x)
  stop(sprintf("convergence not achieved for %s", what))
}

# Standardize non-constant columns with more than two distinct values
# (weighted mean and sd). The column space is unchanged because every model
# has an intercept, so the causal estimates are not affected.
.drlate_std <- function(X, w) {
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    if (length(unique(x)) > 2) {
      m <- sum(w * x) / sum(w)
      s <- sqrt(sum(w * (x - m)^2) / (sum(w) - 1))
      X[, j] <- (x - m) / s
    }
  }
  X
}

# Collects parameter blocks and their per-observation moment functions.
# Each block's moment function takes the full named parameter list and
# returns an N x length(value) matrix.
.drlate_system <- function() {
  env <- new.env()
  env$values <- list()
  env$fns <- list()
  env$add <- function(name, value, fn) {
    env$values[[name]] <- value
    env$fns[[name]] <- fn
    invisible(NULL)
  }
  env
}

# Sandwich variance for a just-identified M-estimator.
# psi(theta) returns the N x q matrix of unweighted moments; w are sampling
# weights (pweights), cl an optional cluster identifier.
.drlate_vcov <- function(psi, theta, w, cl = NULL) {
  gbar <- function(th) colSums(w * psi(th))
  A <- numDeriv::jacobian(gbar, theta)
  S <- w * psi(theta)
  if (is.null(cl)) {
    B <- crossprod(S)
  } else {
    Sg <- rowsum(S, cl)
    B <- crossprod(Sg)
  }
  Ainv <- solve(A)
  V <- Ainv %*% B %*% t(Ainv)
  dimnames(V) <- list(names(theta), names(theta))
  V
}
