#' Extended Generalized Pareto margin
#'
#' Constructs an EGPD margin with parameters `kappa`, `sigma` and `xi`,
#' backed by the \pkg{egpd} package.
#'
#' @param par Named parameter vector (kappa, sigma, xi)
#' @return An object of class `c("egpd", "margin")`.
#' @export
egpd_margin <- function(par = NULL) {
  if(is.null(par)) return(new_margin(par, class = "egpd"))

  kappa <- par[1]
  sigma <- par[2]
  xi <- par[3]

  stopifnot(kappa > 0, sigma > 0)
  new_margin(c(kappa = kappa, sigma = sigma, xi = xi), class = "egpd")
}

#' @describeIn p_margin EGPD margin, via `egpd::pegpd()`.
#' @export
p_margin.egpd <- function(margin, x, ...) {
  param <- margin$param
  prob <- egpd::pegpd(
    x,
    kappa = param["kappa"],
    sigma = param["sigma"],
    xi = param["xi"]
  )

  unname(prob)
}

#' @describeIn q_margin EGPD margin, via `egpd::qegpd()`.
#' @export
q_margin.egpd <- function(margin, p, ...) {
  param <- margin$param
  quant <- egpd::qegpd(
    p,
    kappa = param["kappa"],
    sigma = param["sigma"],
    xi = param["xi"]
  )

  unname(quant)
}

#' @describeIn d_margin EGPD margin, via `egpd::degpd_density()`.
#' @export
d_margin.egpd <- function(margin, x, log = FALSE, ...) {
  param <- margin$param

  dens <- egpd::degpd_density(
    x,
    kappa = param["kappa"],
    sigma = param["sigma"],
    xi = param["xi"],
    log = log
  )

  unname(dens)
}

#' @describeIn r_margin EGPD margin, via `egpd::regpd()`.
#' @export
r_margin.egpd <- function(margin, n, ...) {
  param <- margin$param
  egpd::regpd(n, kappa = param["kappa"], sigma = param["sigma"], xi = param["xi"])
}
