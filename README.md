# CopulaMarkov

**Copula-based Markov models in R.**

`CopulaMarkov` is an R package for constructing, simulating, and fitting copula-based Markov models. It builds on the [`copula`](https://CRAN.R-project.org/package=copula) package and supports flexible combinations of copulas and marginal distributions.

The package separates the dependence structure, represented by a copula, from the marginal distribution, allowing users to specify and study a range of univariate Markov models.

Currently, the package supports the Extended Generalized Pareto Distribution (EGPD) as a marginal distribution, using functionality from the [`egpd`](https://github.com/sdwfrost/egpd) package.

## Features

- **Model construction:** Specify copula-based Markov models using a copula and a marginal distribution.
- **Simulation:** Generate random samples from fitted or specified models.
- **Parameter estimation:** Fit model parameters to observed data.
- **Flexible margins:** Work with marginal distributions through a common interface, with EGPD support currently implemented.
- **Copula-based dependence:** Use the copula framework to represent dependence between successive observations.

## Installation

```r
# install.packages("remotes")
remotes::install_github("AndreFerrero/CopulaMarkov")
```

## Basic usage

The following example illustrates how to construct a copula-based Markov model, simulate a time series from it, and fit the model to the simulated observations.

```r
library(CopulaMarkov)

# Define the marginal distribution with its parameters
m <- egpd_margin(c(2, 1, 0.2))

# Define the dependence structure
cop <- copula::gumbelCopula(2)

# Construct the copula Markov model
mod <- copula_markov(m, cop)

# Simulate a time series
# (warning: Numerical inversion failures for some copulas might require multiple runs)
sim <- simulate(mod, nsim = 2000)

# Fit the model to the simulated observations
init <- c(2.5, 3, 0.6, 5)
fit_mod <- fit(mod, init, sim$X)

# Inspect the fitted model parameters
fit_mod
```

In this example, an EGPD marginal distribution is combined with a Gumbel copula to define a Markov model. A time series of 2,000 observations is simulated, and the model is then fitted to those observations using an initial parameter vector.

## Dependencies

`CopulaMarkov` builds on the following R packages:

- [`copula`](https://CRAN.R-project.org/package=copula) for copula modelling.
- [`egpd`](https://CRAN.R-project.org/package=egpd) for the Extended Generalized Pareto Distribution.

## Development status

The package is under active development. Its functionality and API may change as additional marginal distributions, model-fitting procedures, and features are implemented.

## License

See the `LICENSE` file for licensing information.
