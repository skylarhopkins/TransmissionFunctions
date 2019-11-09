#'Simulate epidemics with different transmission functions and fit multiple transmission functions to the simulated data
#'
#''sim.and.opt' simulates epidemics in populations with different sizes (and thus densities, assuming constant area) using SIR
#'epidemiological models with variable amounts of density-dependence in the transmission function, ranging from
#'the density-dependent (DD) transmission function (K=1) to the frequency-dependent (FD) transmission function (K=0).
#'After simulating, the package fits three models to the simulated datasets to compare the accuracy of their
#'parameter estimates: the DD function, the FD function, and a flexible nonlinear function.
#'
#'This function is the wrapper for the entire simulation and optimization process, using all other functions in the
#'contactr package. It runs simulations for a single FOI, a single gamma, and a range of Ks. You could write a loop
#'to run more than one FOI or gamma at a time, but given the long code run times, this is not advisable. The function
#'will save one output CSV file containing the parameter estimates for each dataset and each random restart for each
#'values of K.
#'
#'@param FOI The Force of Infection value that you are simulating. If you want to do more than one at a time, you'll need to write a loop
#'@param truegamma The gamma value that you are simulating. If you want to do more than one at a time, you'll need to write a loop
#'@param pops A vector of constant population sizes in different populations that will experience simultaneous epidemics; you can include as few/many populations as desired
#'@param Nref Specify a reference population size around the middle of your population sizes for scaling across transmission functions
#'@param ks A vector containing the values of K (unitless density-dependence parameter) that you're simulating over
#'@param initial.I Vector of initial numbers of infected individuals in each pop. To keep the initial prevalence constant, you can scale these relative to pop sizes
#'@param initial.S Vector of initial numbers of infected individuals in each population. S+I should equal N.
#'@param time.out a single value indicating how many time units (e.g., 21 days days) the ODE should be simulated for - should be as long as largest time.samp
#'@param time.samp vector of the times when the populations were sampled, such as 0, 7, 14, and 21 days
#'@param samp.sizes vector of how many individuals will you sample in each population during each sampling event; same length as time.samp
#'@param nrestarts number of fits of each model to each dataset with different random starting parameters
#'@param ndatasets number of sample datasets for each value of K that you're simulating over
#'@param outputlocation The root directory where the output files should be stored
#'
#'@examples
#'###Set Global Variables
#'##If you assign these variables outside of the function arguments, you cannot use the
#'##same names as the argument names, or else you'll get a recursive argument error.
#'##Sticking a period after each name solves this problem
#'#FOI.<-0.0001
#'#truegamma.<-0.1
#'#pops.<-c(100, 200, 500, 1000, 1500, 2000)
#'#Nref.<-1000
#'#ks.<-seq(0.0, 1.0, 0.1)
#'#initial.I.<-c(1,2,5,10,15,20)
#'#initial.S.<-pops.-initial.I.
#'#time.out. <- seq(0,150,by = 1)
#'#time.samp. <- seq(0,133, by = 7)
#'#samp.sizes. <- rep(100, length(time.samp.))
#'#nrestarts. = 1
#'#ndatasets. = 1
#'#outputlocation.<-getwd()
#'###Run the tool
#'##WARNING: this can take a very long time to run depending on ndatasets, nrestarts,
#'##and length of ks,so you might want to estimate run times (end_time - start_time)
#'##on a smaller subset first
#'#start_time <- Sys.time()
#'#sim.and.opt(FOI=FOI., truegamma=truegamma., pops=pops., Nref=Nref., ks=ks., initial.I=initial.I.,
#'#initial.S=initial.S., time.out=time.out., time.samp=time.samp., samp.sizes=samp.sizes.,
#'#nrestarts=nrestarts., ndatasets=ndatasets., outputlocation=outputlocation.)
#'#end_time <- Sys.time()
#'#end_time - start_time
#'
#' @export
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
        ts.sir <- data.frame(deSolve::ode(
          y = c(S=initial.S[e], I=initial.I[e]),               # Initial conditions for population
          times = time.out,             # Timepoints for evaluation
          func = sir,                   # Function to evaluate
          parms = c(beta=truebeta, gamma=truegamma, N=pops[e], K=ks[k]),
          method="lsoda"))               # Vector of parameters
        ts.sir$P <- ts.sir$I / pops[e]
        prev.samp <- ts.sir$P[ts.sir$time %in% time.samp]
        data[,e] <- stats::rbinom(length(time.samp), size = samp.sizes, prob = prev.samp)
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
    utils::write.csv(compareests, paste(outputlocation,"/compareests","K",ks[k],"FOI",FOI,"gamma",truegamma,".csv",sep=""), row.names=FALSE)
  }
}
