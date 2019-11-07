#'Minimize the negative log likelihood to find the best-fitting density-dependent transmission function
#'
#''DD.optim' uses optim() to minimize the nll and fit the density-dependent transmission function
#'given random starting parameters centered around the true/known values for Beta and gamma.
#'
#'You can use this function on it's own, but it is meant to go inside the optimization loop for all models
#'
#'@param FOI The true force of infection, which is used with nref to find a truebeta as a mean for picking random starting values
#'@param N The size of the mid-sized or reference population (Nref), which is used to scale the FOI to a truebeta for the DD function
#'@param gamma The true value of gamma, which is used as a mean for picking random starting values

DD.optim <- function(FoI=FOI, N=Nref, gamma=truegamma)
{
  startpar= c(beta = log(rlnorm(1, mean=log((FOI*Nref)/(Nref^1)), sdlog=1)), gamma = log(rlnorm(1, mean=log(gamma), sdlog=1)))
  outDD<-optim(startpar, nll.DD.fn, control = list(trace = 0, maxit = 1000), method = "Nelder-Mead")
  #print(outNL)
}
