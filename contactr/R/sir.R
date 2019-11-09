#'Simulate an SIR model using ordinary differential equations a constant host population size
#'
#''sir' returns the instantaneous rates of change in the number of susceptible and infected
#'individuals, dSdt and dIdt, respectively
#'
#''sir' is used with the ode function to simulate an epidemic through time assuming constant area and
#'given known starting population sizes of susceptible and infected individuals, beta: a transmission rate,
#'gamma: a recovery rate, and K: a unitless density dependence parameter that allows the
#'transmission function to take on any form. When K=0, the transmission function is frequency-
#'dependent. When K=1, the transmission function is density-dependent. And when 0<K<1, the
#'transmission function is an increasing, saturating function of host density.
#'
#'@param t numeric vector. The time at which the instantaneous rates are evalated.
#'@param y named numeric vector. Initial population sizes of Susceptible and Infected individuals, labeled with S and I
#'@param parms  named numeric vector. List the parameter values being used for beta, gamma, K, and N, the total population size.
#'@return Returns the instantaneous rates of change in the number of susceptible and infected
#'individuals, dSdt and dIdt, respectively
#'@examples
#'#initial.S <- 99
#'#initial.I <- 1
#'#y<-c(S = initial.S, I = initial.I)
#'#parameters<-c(beta=0.0005, gamma=0.05, K=1, N=100)
#'#sir(t=1, y=y, parms=parameters)
#'
#' @export
sir <- function(t,y,parms){
  with(c(as.list(y),parms),{
    dSdt <- -beta*(N^K)*S*I/N
    dIdt <- beta*(N^K)*S*I/N - gamma*I
    #Note: Population size is constant, so don't need to specify dRdt
    list(c(dSdt,dIdt))
  })
}
