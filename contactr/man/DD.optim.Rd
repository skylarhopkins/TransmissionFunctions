% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DD.optim.R
\name{DD.optim}
\alias{DD.optim}
\title{Minimize the negative log likelihood to find the best-fitting density-dependent transmission function}
\usage{
DD.optim(
  FoI,
  N,
  gamma,
  datasets,
  initial.inf,
  initial.sus,
  samp.sizes,
  pops,
  time.outs,
  time.samps
)
}
\arguments{
\item{FoI}{The true force of infection, which is used with nref to find a truebeta as a mean for picking random starting values}

\item{N}{The size of the mid-sized or reference population (Nref), which is used to scale the FOI to a truebeta for the DD function}

\item{gamma}{The true value of gamma, which is used as a mean for picking random starting values}

\item{datasets}{Dataframe with dataset for each population as a column}

\item{initial.inf}{vector of initial number of infected individuals in each population}

\item{initial.sus}{vector of initial number of infected individuals in each population, which should be pops - initial.S}

\item{samp.sizes}{vector of sample sizes, in number of individuals, that were sampled to estimate prevalence at each time point}

\item{pops}{vector of total population sizes in each population, which should be initial.S + initial.I}

\item{time.outs}{a single value indicating how many time units (e.g., 21 days days) the ODE should be simulated for - should be as long as largest time.samp}

\item{time.samps}{vector of the times when the populations were sampled, such as 0, 7, 14, and 21 days}
}
\description{
'DD.optim' uses optim() to minimize the nll and fit the density-dependent transmission function
given random starting parameters centered around the true/known values for Beta and gamma.
}
\details{
You can use this function on it's own, but it is meant to go inside the optimization loop for all models
}
