#'Negative log likelihood when fitting the nonlinear transmission function
#'
#''nll.NL.fn' returns the negative log likelihood (nll) from fitting a nonlinear transmission function
#'with given parameters to a prevalence dataset. It calculates the nll for one population/epidemic
#'at a time, and later the loop adds and minimizes the total nll.
#'
#'You can use this function on it's own, but it is meant to have optim pass it
#'parameters during the optimization procedure.
#'
#'@param pars named numeric vector. # log(B), log(gamma), and log(K) in a vector, must be named - usually passed from optim()
#'@param initial.inf vector of initial number of infected individuals in each population
#'@param initial.sus vector of initial number of infected individuals in each population, which should be pops - initial.S
#'@param n vector of sample sizes, in number of individuals, that were sampled to estimate prevalence at each time point
#'@param popsizes vector of total population sizes in each population, which should be initial.S + initial.I
#'@param time.outs a single value indicating how many time units (e.g., 21 days days) the ODE should be simulated for - should be as long as largest time.samp
#'@param time.samps vector of the times when the populations were sampled, such as 0, 7, 14, and 21 days
#'
#'@examples
#'#startpar = c(beta = log(0.05), gamma = log(0.05), K= log(0.5))
#'#initial.I<-c(1,2,5,10,15,20)
#'#pops<-c(100, 200, 500, 1000, 1500, 2000)
#'#initial.S<-pops-initial.I
#'#time.out <- seq(0,150,by = 1) ##for simulating data from ODEs - output same with 0.1 step
#'#time.samp <- seq(0,133, by = 7) ##will sample 20 time points
#'#samp.sizes <- rep(100, length(time.samp)) ##sample 100 indiv each time
#'#nll.NL.fn(pars=startpar, initial.inf = initial.I,initial.sus = initial.S,n = samp.sizes,datasets = data,popsizes = pops,time.outs = time.out,time.samps = time.samp)


nll.NL.fn <- function(pars,                # log(B), log(gamma), and log(K) in a vector, must be named
                      datasets=data,
                      initial.inf=initial.I,
                      initial.sus=initial.S,
                      n = samp.sizes, # sample sizes
                      popsizes = pops,
                      time.outs = time.out,
                      time.samps = time.samp) # the times they were observed at
{
  pars <- c(exp(pars[1]), exp(pars[2]), exp(pars[3]))
  nlls<-rep(NA, length(popsizes))
  for(i in 1:length(popsizes)) {
    #Create epidemic time series from the parameters
    ts.sir.temp <- data.frame(ode(y = c(S=initial.sus[i], I=initial.inf[i]), times = time.outs,
                                  func = sir, parms = c(pars, N=popsizes[i])))
    ts.sir.temp$P <- ts.sir.temp$I / popsizes[i]
    #Pull out the prevalences at the random time points
    prev.samp.temp <- ts.sir.temp$P[ts.sir.temp$time %in% time.samps]
    prev.samp.temp <- pmin(pmax(prev.samp.temp, 0.0001), 0.9999)
    #calculate the negative log likelihood for these params in this population
    nlls[i]<- -sum(dbinom(datasets[,i], n, prev.samp.temp, log = TRUE))
  }
  nll<-sum(nlls)
  nll
}
