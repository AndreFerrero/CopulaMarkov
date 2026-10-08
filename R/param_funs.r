#' Parameter Bounds
#'
#' Retrieve the lower and upper bounds for the parameters of an object.
#'
#' `par_bounds()` is an S3 generic used to obtain the parameter domain
#' associated with an object. Methods return a two-column matrix with
#' columns `lower` and `upper`, and one row for each parameter.
#'
#' The bounds are used by the parameter transformation functions
#' [transf_par()] and [untransf_par()] to map parameters between their
#' natural and unconstrained representations.
#'
#' @param object An object for which parameter bounds are required.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric matrix with columns `lower` and `upper`. Row names
#'   correspond to parameter names.
#'
#' @export
par_bounds <- function(object, ...) {
  UseMethod("par_bounds")
}


#' Parameter Bounds for an EGPD Margin
#'
#' Return the parameter bounds for an EGPD margin.
#'
#' The EGPD shape parameter `kappa` and scale parameter `sigma` are
#' constrained to be positive, while the shape parameter `xi` is
#' unrestricted.
#'
#' @param object An `egpd` margin object.
#' @param ... Additional arguments passed to the method.
#'
#' @return A numeric matrix with rows `kappa`, `sigma`, and `xi`, and
#'   columns `lower` and `upper`.
#'
#' @export
par_bounds.egpd <- function(object, ...) {
  cbind(
    lower = c(kappa = 0, sigma = 0, xi = -Inf),
    upper = c(kappa = Inf, sigma = Inf, xi = Inf)
  )
}


#' Parameter Bounds for a Copula
#'
#' Return the parameter bounds defined by a `copula` package copula object.
#'
#' The lower and upper parameter bounds are obtained from the
#' `param.lowbnd` and `param.upbnd` slots of the S4 copula object.
#'
#' @param object A copula object from the `copula` package.
#' @param ... Additional arguments passed to the method.
#'
#' @return A numeric matrix with columns `lower` and `upper`, and one
#'   row for each copula parameter.
#'
#' @export
par_bounds.copula <- function(object, ...) {
  cbind(
    lower = object@param.lowbnd,
    upper = object@param.upbnd
  )
}


#' Parameter Bounds for a Copula Markov Model
#'
#' Return the combined parameter bounds for the margin and copula
#' components of a copula Markov model.
#'
#' The parameter bounds of the margin are followed by the parameter
#' bounds of the copula. The ordering therefore corresponds to the
#' ordering of parameters expected by the model likelihood.
#'
#' @param object A `copula_markov` model object.
#' @param ... Additional arguments passed to the method.
#'
#' @return A numeric matrix with columns `lower` and `upper`. Rows
#'   corresponding to margin parameters are followed by rows
#'   corresponding to copula parameters.
#'
#' @export
par_bounds.copula_markov <- function(object, ...) {
  margin_bounds <- par_bounds(object$margin)
  copula_bounds <- par_bounds(object$copula)

  rbind(margin_bounds, copula_bounds)
}


#' Transform Parameters to an Unconstrained Scale
#'
#' Transform parameters from their natural, constrained scale to an
#' unconstrained scale.
#'
#' `transf_par()` uses the parameter bounds returned by [par_bounds()]
#' to determine an appropriate transformation for each parameter.
#' The resulting parameters can therefore be passed to an unconstrained
#' optimization routine such as [stats::optim()] with method `"BFGS"`.
#'
#' The transformations are determined by the parameter bounds:
#'
#' \describe{
#'   \item{`(-Inf, Inf)`}{Identity transformation.}
#'   \item{`[a, Inf)`}{`log(par - a)`.}
#'   \item{`(-Inf, b]`}{`log(b - par)`.}
#'   \item{`[a, b]`}{Logit transformation of the rescaled parameter.}
#' }
#'
#' For finite lower and upper bounds, the transformation maps the
#' interior of the interval `(a, b)` to the entire real line.
#'
#' @param object An object for which parameter bounds are defined.
#' @param par Numeric vector of parameters on their natural scale.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric vector of parameters on the unconstrained scale.
#'
#' @seealso [untransf_par()], [par_bounds()]
#'
#' @export
transf_par <- function(object, par, ...) {
  UseMethod("transf_par")
}


#' Transform Parameters to Their Natural Scale
#'
#' Transform parameters from an unconstrained scale to their natural,
#' constrained scale.
#'
#' `untransf_par()` uses the parameter bounds returned by
#' [par_bounds()] to determine the inverse transformation for each
#' parameter. It is the inverse of [transf_par()].
#'
#' The inverse transformations are determined by the parameter bounds:
#'
#' \describe{
#'   \item{`(-Inf, Inf)`}{Identity transformation.}
#'   \item{`[a, Inf)`}{`a + exp(par)`.}
#'   \item{`(-Inf, b]`}{`b - exp(par)`.}
#'   \item{`[a, b]`}{`a + (b - a) * plogis(par)`.}
#' }
#'
#' @param object An object for which parameter bounds are defined.
#' @param par Numeric vector of parameters on the unconstrained scale.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric vector of parameters on their natural scale.
#'
#' @seealso [transf_par()], [par_bounds()]
#'
#' @export
untransf_par <- function(object, par, ...) {
  UseMethod("untransf_par")
}


#' Transform Parameters to an Unconstrained Scale
#'
#' Default method for [transf_par()].
#'
#' The transformation applied to each parameter is determined by the
#' corresponding lower and upper bounds returned by [par_bounds()].
#'
#' @param object An object for which parameter bounds are defined.
#' @param par Numeric vector of parameters on their natural scale.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric vector of parameters on the unconstrained scale.
#'
#' @export
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


#' Transform Parameters to Their Natural Scale
#'
#' Default method for [untransf_par()].
#'
#' The inverse transformation applied to each parameter is determined
#' by the corresponding lower and upper bounds returned by
#' [par_bounds()].
#'
#' @param object An object for which parameter bounds are defined.
#' @param par Numeric vector of parameters on the unconstrained scale.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric vector of parameters on their natural scale.
#'
#' @export
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