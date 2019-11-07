#'Minimize the negative log likelihood to find the best-fitting frequency-dependent transmission function
#'
#''FD.optim' uses optim() to minimize the nll and fit the frequency-dependent transmission function
#'given random starting parameters centered around the true/known value for gamma. Given convergence issues,
#'the range of starting values for Beta is set between -2 and 0 (on the logscale).
#'
#'You can use this function on it's own, but it is meant to go inside the optimization loop for all models
#'
#'@param randombeta A random starting value for Beta on the logscale - works best between -2 and 0, but can be changed
#'@param gamma The true value of gamma, which is used as a mean for picking random starting values

FD.optim <- function(randombeta=runif(1, -2, 0), gamma=truegamma)
{
  startpar= c(beta = randombeta, gamma = log(rlnorm(1, mean=log(gamma), sdlog=1)))
  outDD<-optim(startpar, nll.FD.fn, control = list(trace = 0, maxit = 1000), method = "Nelder-Mead")
  #print(outNL)
}
