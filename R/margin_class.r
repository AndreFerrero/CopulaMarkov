#' Cumulative distribution function of a margin
#' @param margin A `margin` object.
#' @param x Numeric vector of quantiles.
#' @param ... Passed to methods.
#' @return Numeric vector of probabilities in (0, 1).
#' @export
p_margin <- function(margin, x, ...) UseMethod("p_margin")

#' Log-density of a margin
#' @param margin A `margin` object.
#' @param x Numeric vector of quantiles.
#' @param ... Passed to methods.
#' @return Numeric vector of log-density values.
#' @export
d_margin <- function(margin, x, log = FALSE, ...) UseMethod("d_margin")

#' Quantile function of a margin
#' @param margin A `margin` object.
#' @param p Numeric vector of probabilities in (0, 1).
#' @param ... Passed to methods.
#' @return Numeric vector of quantiles.
#' @export
q_margin <- function(margin, p, ...) UseMethod("q_margin")

#' Random generation from a margin
#' @param margin A `margin` object.
#' @param n Number of observations to draw.
#' @param ... Passed to methods.
#' @return Numeric vector of length `n`.
#' @export
r_margin <- function(margin, n, ...) UseMethod("r_margin")

#' Low-level constructor for margin objects
#'
#' Internal helper used by family constructors such as [egpd_margin()].
#' It performs no validation; family constructors are responsible for that.
#'
#' @param param Named numeric vector of margin parameters.
#' @param ... Additional fields stored in the object.
#' @param class Character vector of subclass names, prepended to `"margin"`.
#' @return An object of class `c(class, "margin")`.
#' @keywords internal
new_margin <- function(param, ..., class) {
  structure(list(param = param, ...), class = c(class, "margin"))
}

#' Utility to update margin parameters
#' @export
param_update <- function(margin, new_param, ...) {
  stopifnot(length(new_param) == length(margin$param))

  param_names <- names(margin$param)
  names(new_param) <- param_names

  margin$param <- new_param
  margin
}

#' Print a margin
#' @param x A `margin` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.margin <- function(x, ...) {
  print(class(x))
  print(x$param)
}
