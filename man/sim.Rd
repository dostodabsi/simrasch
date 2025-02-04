% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/main.R
\name{sim}
\alias{sim}
\title{Runs the main simulation}
\usage{
sim(n, items, times, model_sim, ..., Sigma = NULL, weights = NULL,
  parallel = TRUE, cores = detectCores())
}
\arguments{
\item{n,}{items: numeric vectors}

\item{times:}{numeric}

\item{model_sim:}{function to simulate Rasch data}

\item{Sigma:}{factor correlation matrix for multidimensional data}

\item{weights:}{item loadings of the factor structure}

\item{parallel:}{if simulation should be run in parallel}

\item{cores:}{number of cores to run in parallel}

\item{...:}{additional arguments to model_sim}
}
\value{
returns a list that stores objects of class 'sim_res'
each array stores #times simulations of each possible n x items combination
}
\description{
Runs the main simulation
}
\examples{
sim(c(100, 400), c(10, 50), 1000, sim.2pl, .5, parallel = TRUE, cores = 3)
}

