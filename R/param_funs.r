par_bounds <- function(object, ...) {
  UseMethod("par_bounds")
}

par_bounds.egpd <- function(object, ...) {
  cbind(
    lower = c(kappa = 0, sigma = 0, xi = -Inf),
    upper = c(kappa = Inf, sigma = Inf, xi = Inf)
  )
}

par_bounds.copula <- function(object, ...) {
  cbind(
    lower = object@param.lowbnd,
    upper = object@param.upbnd
  )
}

par_bounds.copula_markov <- function(object, ...) {
  margin_bounds <- par_bounds(object$margin)
  copula_bounds <- par_bounds(object$copula)

  rbind(margin_bounds, copula_bounds)
}

transf_par <- function(object, par, ...) {
  UseMethod("transf_par")
}

untransf_par <- function(object, par, ...) {
  UseMethod("untransf_par")
}

transf_par.default <- function(object, par, ...) {
  bounds <- par_bounds(object)

  vapply(seq_along(par), function(i) {
    lower <- bounds[i, "lower"]
    upper <- bounds[i, "upper"]
    par_i <- par[i]

    if (is.infinite(lower) && is.infinite(upper)) {
      return(par_i)
    }

    if (is.finite(lower) && is.infinite(upper)) {
      return(log(par_i - lower))
    }

    if (is.infinite(lower) && is.finite(upper)) {
      return(log(upper - par_i))
    }

    qlogis((par_i - lower) / (upper - lower))
  }, numeric(1))
}

untransf_par.default <- function(object, par, ...) {
  bounds <- par_bounds(object)

  vapply(seq_along(par), function(i) {
    lower <- bounds[i, "lower"]
    upper <- bounds[i, "upper"]
    par_i <- par[i]

    if (is.infinite(lower) && is.infinite(upper)) {
      return(par_i)
    }

    if (is.finite(lower) && is.infinite(upper)) {
      return(lower + exp(par_i))
    }

    if (is.infinite(lower) && is.finite(upper)) {
      return(upper - exp(par_i))
    }

    lower + (upper - lower) * plogis(par_i)
  }, numeric(1))
}
