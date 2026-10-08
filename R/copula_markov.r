#' Copula Markov Model
#'
#' Construct a copula Markov model from its margin and copula
#'
#' @param margin Package-defined margin class.
#' @param copula Copula object from the `copula` package.
#' @return A `copula_markov` object.
#' @export
copula_markov <- function(margin, copula) {
  stopifnot(inherits(margin, "margin"))

  structure(
    list(
      margin = margin,
      copula = copula
    ),
    class = "copula_markov"
  )
}


#' Simulate Copula Markov Model using Inverse Rosenblatt Transform
#'
#' @param object Copula Markov model.
#' @param nsim Sample size.
#' @param seed Optional random seed.
#' @param ... Passed to methods.
#' @return A list containing data on the margin scale (`X`) and copula
#'   scale (`U`).
#' @export
simulate.copula_markov <- function(object, nsim, seed = NULL, ...) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  V <- runif(nsim, 1e-5, 1 - 1e-5)

  U <- numeric(nsim)
  U[1] <- V[1]

  for (t in 2:nsim) {
    U[t] <- copula::cCopula(
      cbind(U[t - 1], V[t]),
      copula = object$copula,
      inverse = TRUE
    )[, 2]
  }

  X <- q_margin(object$margin, U)

  list(
    X = X,
    U = U
  )
}


#' Log-likelihood of a model
#'
#' @param object A model object.
#' @param par Parameter vector on the natural parameter scale.
#' @param x Observed data.
#' @param ... Passed to methods.
#' @return Numeric log-likelihood.
#' @export
log_lik <- function(object, par, x, ...) {
  UseMethod("log_lik")
}


#' Log-likelihood of a Copula Markov Model
#'
#' @param object Copula Markov Model object.
#' @param par Parameter vector on the natural parameter scale.
#' @param x Data vector.
#' @param ... Passed to methods.
#' @return Numeric log-likelihood.
#' @export
log_lik.copula_markov <- function(object, par, x, ...) {

  # Number of margin parameters
  n_margin <- length(object$margin$param)

  # Split parameter vector
  margin_par <- par[seq_len(n_margin)]
  copula_par <- par[-seq_len(n_margin)]

  # Update margin
  margin <- param_update(
    object$margin,
    margin_par
  )

  # Update copula
  cop <- copula::setTheta(
    object$copula,
    copula_par
  )

  # Transform observations to copula scale
  nsim <- length(x)

  log_f <- d_margin(
    margin,
    x,
    log = TRUE
  )

  u <- p_margin(
    margin,
    x
  )

  # Copula contribution for consecutive observations
  log_c <- copula::dCopula(
    cbind(u[-nsim], u[-1]),
    cop,
    log = TRUE
  )

  ll <- sum(log_f) + sum(log_c)

  if (!is.finite(ll)) {
    return(-Inf)
  }

  ll
}


#' Fit a Copula Markov Model
#'
#' @param object A Copula Markov Model object.
#' @param init Initial parameter vector on the natural parameter scale.
#' @param x Observed data.
#' @param ... Passed to methods.
#' @return Fitted parameter vector
#' @export
fit <- function(object, init, x, ...) {
  UseMethod("fit")
}


#' @export
fit.copula_markov <- function(object, init, x, ...) {

  # Check that the initial parameter vector has the
  # expected length
  bounds <- par_bounds(object)

  if (length(init) != nrow(bounds)) {
    stop(
      "`init` must contain ",
      nrow(bounds),
      " parameters."
    )
  }

  # Convert natural parameters to unconstrained parameters
  init_unconstrained <- transf_par(
    object,
    init
  )

  # Objective function operates on the unconstrained scale
  objective <- function(par_unconstrained) {

    # Convert back to natural parameter scale
    par <- untransf_par(
      object,
      par_unconstrained
    )

    # Negative log-likelihood for minimization
    -log_lik(
      object,
      par,
      x
    )
  }

  # Unconstrained optimization
  opt <- optim(
    par = init_unconstrained,
    fn = objective,
    method = "BFGS",
    ...
  )

  par_hat <- untransf_par(object, opt$par)

  par_hat
}


#' @export
print.copula_markov <- function(x, ...) {

  cat("Copula-Markov model\n")
  cat(" margin :", class(x$margin)[1], "\n")
  cat(" copula :", class(x$copula)[1], "\n")

  invisible(x)
}