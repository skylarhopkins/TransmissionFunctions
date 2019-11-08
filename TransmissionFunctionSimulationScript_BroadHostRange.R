#This script defines all the functions in the contactr package, a tool to simulate epidemics
#in populations with different sizes (and thus densities, assuming constant area) using SIR
#epidemiological models with variable amounts of density-dependence in the transmission function, ranging from
#the density-dependent (DD) transmission function (K=1) to the frequency-dependent (FD) transmission function (K=0).
#After simulating, the package fits three models to the simulated datasets to compare the accuracy of their
#parameter estimates: the DD function, the FD function, and a flexible nonlinear function.

library(deSolve)

#################################################################################
####################System of Diffeqs Function#####################################
#################################################################################
sir <- function(t,y,parms){
  with(c(as.list(y),parms),{
    dSdt <- -beta*(N^K)*S*I/N
    dIdt <- beta*(N^K)*S*I/N - gamma*I
    #Note: Population size is constant, so don't need to specify dRdt
    list(c(dSdt,dIdt))
  })
}

################################################################################
##########NL Neg Log Likelihood and optimization functions#######################
#################################################################################
NL.optim <- function(beta,
                     gamma,
                     datasets,
                     initial.inf,
                     initial.sus,
                     samp.sizes, # sample sizes
                     pops,
                     time.outs,
                     time.samps) # the times they were observed at
{
  nll.NL.fn <- function(pars) # log(B) & log(gamma) in a vector, must be named
  {
    pars <- c(exp(pars[1]), exp(pars[2]), exp(pars[3]))
    nlls<-rep(NA, length(pops))
    for(i in 1:length(pops)) {
      #Create epidemic time series from the parameters
      ts.sir.temp <- data.frame(ode(y = c(S=initial.sus[i], I=initial.inf[i]), times = time.outs,
                                    func = sir, parms = c(pars, N=pops[i])))
      ts.sir.temp$P <- ts.sir.temp$I / pops[i]
      #Pull out the prevalences at the random time points
      prev.samp.temp <- ts.sir.temp$P[ts.sir.temp$time %in% time.samps]
      prev.samp.temp <- pmin(pmax(prev.samp.temp, 0.0001), 0.9999)
      #calculate the negative log likelihood for these params in this population
      nlls[i]<- -sum(dbinom(datasets[,i], samp.sizes, prev.samp.temp, log = TRUE))
    }
    nll<-sum(nlls)
    nll
  }
  startpar = c(beta = log(rlnorm(1, mean=log(beta), sdlog=1)), gamma = log(rlnorm(1, mean=log(gamma), sdlog=1)), K=log(rlnorm(1, mean=log(0.1), sdlog=1)))
  outNL<-optim(startpar, nll.NL.fn, control = list(trace = 0, maxit = 1000), method = "Nelder-Mead")
  #print(outNL)
}
ks.<-0.1
truebeta.<-FOI.*(Nref.)/(Nref.^ks.) #use known FOI to define truebeta for simulation and random parameter picking
data<-as.data.frame(just.sim(FOI=FOI., truegamma=truegamma., pops=pops., Nref=Nref., ks=ks., initial.I=initial.I., initial.S=initial.S., time.out=time.out., time.samp=time.samp., samp.sizes=samp.sizes., nrestarts=nrestarts., ndatasets=ndatasets.))
data.<-data
outNL<-NL.optim(beta=truebeta., gamma = truegamma., datasets = data., initial.inf = initial.I., initial.sus = initial.S., samp.sizes = samp.sizes., pops = pops., time.outs = time.out., time.samps = time.samp.)
exp(as.numeric(outNL$par))

##############################################################################
##########DD Neg Log Likelihood and optimization functions#######################
###############################################################################
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
      ts.sir.temp <- data.frame(ode(y = c(S=initial.sus[i], I=initial.inf[i]), times = time.outs,
                                    func = sir, parms = c(pars, N=pops[i])))
      ts.sir.temp$P <- ts.sir.temp$I / pops[i]
      #Pull out the prevalences at the random time points
      prev.samp.temp <- ts.sir.temp$P[ts.sir.temp$time %in% time.samps]
      prev.samp.temp <- pmin(pmax(prev.samp.temp, 0.0001), 0.9999)
      #calculate the negative log likelihood for these params in this population
      nlls[i]<- -sum(dbinom(datasets[,i], samp.sizes, prev.samp.temp, log = TRUE))
    }
    nll<-sum(nlls)
    nll
  }
  startpar= c(beta = log(rlnorm(1, mean=log((FoI*N)/(N^1)), sdlog=1)), gamma = log(rlnorm(1, mean=log(gamma), sdlog=1)))
  outDD<-optim(startpar, nll.DD.fn, control = list(trace = 0, maxit = 1000), method = "Nelder-Mead")
  #print(outNL)
}

data.<-as.data.frame(matrix(seq(1,20), nrow = 20, ncol=6))
outDD<-DD.optim(FoI = FOI., N = Nref., gamma = truegamma., datasets = data., initial.inf = initial.I., initial.sus = initial.S., samp.sizes = samp.sizes., pops = pops., time.outs = time.out., time.samps = time.samp.)

#################################################################################
##########FD Neg Log Likelihood and optimization functions#######################
################################################################################
FD.optim <- function(gamma,
                     datasets,
                     initial.inf,
                     initial.sus,
                     samp.sizes, # sample sizes
                     pops,
                     time.outs,
                     time.samps) # the times they were observed at
{
  nll.FD.fn <- function(pars) # log(B) & log(gamma) in a vector, must be named
  {
    pars <- c(exp(pars[1]), exp(pars[2]), K=0) #K=1 for FD
    nlls<-rep(NA, length(pops))
    for(i in 1:length(pops)) {
      #Create epidemic time series from the parameters
      ts.sir.temp <- data.frame(ode(y = c(S=initial.sus[i], I=initial.inf[i]), times = time.outs,
                                    func = sir, parms = c(pars, N=pops[i])))
      ts.sir.temp$P <- ts.sir.temp$I / pops[i]
      #Pull out the prevalences at the random time points
      prev.samp.temp <- ts.sir.temp$P[ts.sir.temp$time %in% time.samps]
      prev.samp.temp <- pmin(pmax(prev.samp.temp, 0.0001), 0.9999)
      #calculate the negative log likelihood for these params in this population
      nlls[i]<- -sum(dbinom(datasets[,i], samp.sizes, prev.samp.temp, log = TRUE))
    }
    nll<-sum(nlls)
    nll
  }
  startpar=c(beta=runif(1, -2, 0), gamma = log(rlnorm(1, mean=log(gamma), sdlog=1)))
  outFD<-optim(startpar, nll.FD.fn, control = list(trace = 0, maxit = 1000), method = "Nelder-Mead")
  #print(outNL)
}
data.<-as.data.frame(matrix(seq(1,20), nrow = 20, ncol=6))
outFD<-FD.optim(gamma = truegamma., datasets = data., initial.inf = initial.I., initial.sus = initial.S., samp.sizes = samp.sizes., pops = pops., time.outs = time.out., time.samps = time.samp.)

###############################################################################
###############Simulation and Optimization Function############################
###############################################################################

just.sim<-function(FOI, truegamma, pops, Nref, ks, initial.I, initial.S, time.out, time.samp, samp.sizes, nrestarts, ndatasets) {
  for (k in 1:length(ks)) {
    ##set up dataframe for output of the fitting process - one per K value
    factors<-expand.grid(fittingattempt=seq(1, nrestarts,1), dataset=seq(1, ndatasets, 1), K=ks[k])
    L<-length(factors$K)
    compareests<-cbind(factors, data.frame(betaNL=rep(NA,L),gammaNL=rep(NA,L), KNL=rep(NA,L), betaDD=rep(NA,L), gammaDD=rep(NA,L), betaFD=rep(NA,L), gammaFD=rep(NA,L), lik=rep(NA,L), conv=rep(NA,L), likDD=rep(NA,L), convDD=rep(NA,L), likFD=rep(NA,L), convFD=rep(NA,L), truebeta=rep(NA, L)))
    #create each simulated dataset and try to fit to it nrestart times
    for (j in 1:ndatasets) {
      truebeta<-FOI*(Nref)/(Nref^ks[k]) #use known FOI to define truebeta for simulation and random parameter picking
      #loop to get each sample dataset as a column in a dataframe
      data<-data.frame(matrix(vector(), nrow=length(time.samp), ncol=length(pops))) #empty dataframe for data
      for(e in 1:length(pops)) {
        ts.sir <- data.frame(ode(
          y = c(S=initial.S[e], I=initial.I[e]),               # Initial conditions for population
          times = time.out,             # Timepoints for evaluation
          func = sir,                   # Function to evaluate
          parms = c(beta=truebeta, gamma=truegamma, N=pops[e], K=ks[k]),
          method="lsoda"))               # Vector of parameters
        ts.sir$P <- ts.sir$I / pops[e]
        prev.samp <- ts.sir$P[ts.sir$time %in% time.samp]
        data[,e] <- rbinom(length(time.samp), size = samp.sizes, prob = prev.samp)
      }
      #print(data)
    }
  }
  data
}

just.sim(FOI=FOI., truegamma=truegamma., pops=pops., Nref=Nref., ks=ks., initial.I=initial.I., initial.S=initial.S., time.out=time.out., time.samp=time.samp., samp.sizes=samp.sizes., nrestarts=nrestarts., ndatasets=ndatasets.)

outputlocation.="~/Documents/Transmission Function Literature Review"
sim.and.opt<-function(FOI, truegamma, pops, Nref, ks, initial.I, initial.S, time.out, time.samp, samp.sizes, nrestarts, ndatasets, outputlocation) {
  for (k in 1:length(ks)) {
    ##set up dataframe for output of the fitting process - one per K value
    factors<-expand.grid(fittingattempt=seq(1, nrestarts,1), dataset=seq(1, ndatasets, 1), K=ks[k])
    L<-length(factors$K)
    compareests<-cbind(factors, data.frame(betaNL=rep(NA,L),gammaNL=rep(NA,L), KNL=rep(NA,L), betaDD=rep(NA,L), gammaDD=rep(NA,L), betaFD=rep(NA,L), gammaFD=rep(NA,L), lik=rep(NA,L), conv=rep(NA,L), likDD=rep(NA,L), convDD=rep(NA,L), likFD=rep(NA,L), convFD=rep(NA,L), truebeta=rep(NA, L)))
    #create each simulated dataset and try to fit to it nrestart times
    for (j in 1:ndatasets) {
      truebeta<-FOI*(Nref)/(Nref^ks[k]) #use known FOI to define truebeta for simulation and random parameter picking
      #loop to get each sample dataset as a column in a dataframe
      data<-data.frame(matrix(vector(), nrow=length(time.samp), ncol=length(pops))) #empty dataframe for data
      for(e in 1:length(pops)) {
        ts.sir <- data.frame(ode(
          y = c(S=initial.S[e], I=initial.I[e]),               # Initial conditions for population
          times = time.out,             # Timepoints for evaluation
          func = sir,                   # Function to evaluate
          parms = c(beta=truebeta, gamma=truegamma, N=pops[e], K=ks[k]),
          method="lsoda"))               # Vector of parameters
        ts.sir$P <- ts.sir$I / pops[e]
        prev.samp <- ts.sir$P[ts.sir$time %in% time.samp]
        data[,e] <- rbinom(length(time.samp), size = samp.sizes, prob = prev.samp)
      }
      for (i in 1:nrestarts) {
        print(c("NEW PARAMS", "K=", ks[k], "Dataset #", j, "Restart #", i)) #print to see loop progress
        ##optimization procedures - set random restarts and re-try once if they immediately produce errors
        tryNL <- try(outNL<-NL.optim(beta = truebeta, gamma=truebeta, datasets = data, initial.inf = initial.I, initial.sus = initial.S, time.outs = time.out, samp.sizes = samp.sizes, pops = pops, time.samps = time.samp))
        if (class(tryNL) == "try-error") {
          outNL<-NL.optim(beta = truebeta, gamma=truebeta, datasets = data, initial.inf = initial.I, initial.sus = initial.S, time.outs = time.out, samp.sizes = samp.sizes, pops = pops, time.samps = time.samp)
        }
        tryDD <- try(outDD<-DD.optim(FoI = FOI, N=Nref, gamma=truebeta, datasets = data, initial.inf = initial.I, initial.sus = initial.S, time.outs = time.out, samp.sizes = samp.sizes, pops = pops, time.samps = time.samp))
        if (class(tryDD) == "try-error") {
          outDD<-DD.optim(FoI = FOI, N=Nref, gamma=truebeta, datasets = data, initial.inf = initial.I, initial.sus = initial.S, time.outs = time.out, samp.sizes = samp.sizes, pops = pops, time.samps = time.samp)
        }
        tryFD <- try(outFD<-FD.optim(gamma=truebeta, datasets = data, initial.inf = initial.I, initial.sus = initial.S, time.outs = time.out, samp.sizes = samp.sizes, pops = pops, time.samps = time.samp))
        if (class(tryFD) == "try-error") {
          outFD<-FD.optim(gamma=truebeta, datasets = data, initial.inf = initial.I, initial.sus = initial.S, time.outs = time.out, samp.sizes = samp.sizes, pops = pops, time.samps = time.samp)
        }
        ##save output
        AssignmentRows<-which(compareests$K==ks[k] & compareests$fittingattempt==i & compareests$dataset==j)
        compareests[AssignmentRows, c(4:10)]<-exp(c(as.numeric(outNL$par), as.numeric(outDD$par), as.numeric(outFD$par)))
        compareests$lik[AssignmentRows]<-outNL$value  #negative log likelihood value
        compareests$conv[AssignmentRows]<-outNL$convergence  #if not 0, didn't converge
        compareests$likDD[AssignmentRows]<-outDD$value  #negative log likelihood value
        compareests$convDD[AssignmentRows]<-outDD$convergence  #if not 0, didn't converge
        compareests$likFD[AssignmentRows]<-outFD$value  #negative log likelihood value
        compareests$convFD[AssignmentRows]<-outFD$convergence  #if not 0, didn't converge
        compareests$truebeta[AssignmentRows]<-truebeta
      }
    }
    #output one dataframe per K value to a CSV to save
    write.csv(compareests, paste(outputlocation,"/compareests","K",ks[k],"FOI",FOI,"gamma",truegamma,".csv",sep=""), row.names=FALSE)
  }
}

#################################################################################
#########################Set Global Variables#####################################
#################################################################################
#If you assign these outside of the function arguments, you cannot use the same names as the
#argument names, or else you'll get a recursive argument error. Sticking a period after each
#name solves this problem

#number of fits of each model to each dataset with different random starting parameters
nrestarts. = 1
#number of sample datasets for each value of K that you're simulating over
ndatasets. = 1
#The values of K (unitless density-dependence parameter) that you're simulating over
ks.<-seq(0.0, 1.0, 0.1)
#The Force of Infection value that you are simulating. If you want to do more than one at a time, you'll need to write a loop
FOI.<-0.0001
#The gamma value that you are simulating. If you want to do more than one at a time, you'll need to write a loop
truegamma.<-0.1
#A list of constant population sizes in different populations that will experience simultaneous epidemics; you can
#include as few/many populations as desired
pops.<-c(100, 200, 500, 1000, 1500, 2000)
#To help with comparing across populations and picking relevant starting values for optimization, specify
#a reference population size around the middle of your population sizes
Nref.<-1000
#Initial number of infected individuals in each population. To keep the initial prevalence roughly constant, you
#can scale these relative to the population sizes
initial.I.<-c(1,2,5,10,15,20)
#Initial number of infected individuals in each population. S+I should equal N
initial.S.<-pops.-initial.I.
#When simulating data from ODEs, what time steps do you want simulation output for? EX: each day for 150 days
time.out. <- seq(0,150,by = 1)
#When do you want to sample each population during the epidemic? EX: One per week for 133 days
time.samp. <- seq(0,133, by = 7)
#How many individuals will you sample in each population during each sampling event?
#Must be the same number of individuals in each population each time, but it can change across sampling events
samp.sizes. <- rep(100, length(time.samp.))
#Where should the CSV output file for each K be saved?
outputlocation.<-getwd()
outputlocation.<-"~/Documents/Transmission Function Literature Review"

#################################################################################
########################Run the tool##############################################
################################################################################
#WARNING: this can take a very long time to run depending on ndatasets, nrestarts, and length of ks,
#so you might want to estimate run times (end_time - start_time) on a smaller subset first

start_time <- Sys.time()

sim.and.opt(FOI=FOI., truegamma=truegamma., pops=pops., Nref=Nref., ks=ks., initial.I=initial.I., initial.S=initial.S., time.out=time.out., time.samp=time.samp., samp.sizes=samp.sizes., nrestarts=nrestarts., ndatasets=ndatasets., outputlocation=outputlocation.)

end_time <- Sys.time()
end_time - start_time