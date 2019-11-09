#'Minimize the negative log likelihood to find the best-fitting density-dependent transmission function
#'
#''DD.optim' uses optim() to minimize the nll and fit the density-dependent transmission function
#'given random starting parameters centered around the true/known values for Beta and gamma.
#'
#'You can use this function on it's own, but it is meant to go inside the optimization loop for all models
#'
#'@param FoI The true force of infection, which is used with nref to find a truebeta as a mean for picking random starting values
#'@param N The size of the mid-sized or reference population (Nref), which is used to scale the FOI to a truebeta for the DD function
#'@param gamma The true value of gamma, which is used as a mean for picking random starting values
#'@param datasets Dataframe with dataset for each population as a column
#'@param initial.inf vector of initial number of infected individuals in each population
#'@param initial.sus vector of initial number of infected individuals in each population, which should be pops - initial.S
#'@param samp.sizes vector of sample sizes, in number of individuals, that were sampled to estimate prevalence at each time point
#'@param pops vector of total population sizes in each population, which should be initial.S + initial.I
#'@param time.outs a single value indicating how many time units (e.g., 21 days days) the ODE should be simulated for - should be as long as largest time.samp
#'@param time.samps vector of the times when the populations were sampled, such as 0, 7, 14, and 21 days
#'
#' @export
DD.optim <- function(FoI,
                     N,
                     gamma,
                     datasets,
                     initial.inf,
                     initial.sus,
                     samp.sizes, # sample sizes
                     pops,
                     time.outs,
                     time.samps) # the times they were observed at
{
  nll.DD.fn <- function(pars) # log(B) & log(gamma) in a vector, must be named
  {
    pars <- c(exp(pars[1]), exp(pars[2]), K=1) #K=1 for DD
    nlls<-rep(NA, length(pops))
    for(i in 1:length(pops)) {
      #Create epidemic time series from the parameters
      ts.sir.temp <- data.frame(deSolve::ode(y = c(S=initial.sus[i], I=initial.inf[i]), times = time.outs,
                                    func = sir, parms = c(pars, N=pops[i])))
      ts.sir.temp$P <- ts.sir.temp$I / pops[i]
      #Pull out the prevalences at the random time points
      prev.samp.temp <- ts.sir.temp$P[ts.sir.temp$time %in% time.samps]
      prev.samp.temp <- pmin(pmax(prev.samp.temp, 0.0001), 0.9999)
      #calculate the negative log likelihood for these params in this population
      nlls[i]<- -sum(stats::dbinom(datasets[,i], samp.sizes, prev.samp.temp, log = TRUE))
    }
    nll<-sum(nlls)
    nll
  }
  startpar= c(beta = log(stats::rlnorm(1, mean=log((FoI*N)/(N^1)), sdlog=1)), gamma = log(stats::rlnorm(1, mean=log(gamma), sdlog=1)))
  outDD<-stats::optim(startpar, nll.DD.fn, control = list(trace = 0, maxit = 1000), method = "Nelder-Mead")
  #print(outNL)
}
