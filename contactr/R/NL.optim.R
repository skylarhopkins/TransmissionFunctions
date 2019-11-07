#'Minimize the negative log likelihood to find the best-fitting nonlinear transmission function
#'
#''NL.optim' uses optim() to minimize the nll and fit the nonlinear transmission function
#'given random starting parameters centered around the true/known values for Beta and gamma. The random starting
#'value for K is somewhere between ~0 and ~1 (on the logscale).
#'
#'You can use this function on it's own, but it is meant to go inside the optimization loop for all models
#'
#'@param beta The true value of Beta, which is used as a mean for picking random starting values
#'@param gamma The true value of gamma, which is used as a mean for picking random starting values

NL.optim <- function(beta=truebeta, gamma=truegamma)
{
  startpar = c(beta = log(rlnorm(1, mean=log(beta), sdlog=1)), gamma = log(rlnorm(1, mean=log(gamma), sdlog=1)), K=runif(1, -2.3, -0.1))
  outNL<-optim(startpar, nll.NL.fn, control = list(trace = 0, maxit = 1000), method = "Nelder-Mead")
  #print(outNL)
}
