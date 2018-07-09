##5 restarts, 100 datasets, k=0-1, N=100, 200, 400, 800, 1600, 3200
#We added a wider density range to address a reviewer's comment
##5 July 2018

library(deSolve) 

#################################################################################
###########################Global Variables#####################################
nrestarts = 5 ##number of fits of each model to each dataset with different starting parameters
ndatasets = 100 ##number of sample datasets for each K
FOI<-0.0005 #0.001, 0.0005, 0.0001
gamma<-0.05 #0.02, 0.05, 0.1
ks<-seq(0, 1, 0.1) 

#Epidemics will be generated in 6 populations with densities that are constant in time 
#(ie no demography) and then combined to create a dataset where density varies in space
N0<-100
N02<-200
N03<-500
N04<-1000
N05<-1500
N06<-2000

#initial population sizes - introduce X-X infected individuals to each population
#such that each population starts with 1% infection prevalence
initial.SI1 <- c(S = (N0-1), I = 1)
initial.SI2 <- c(S = (N02-2), I = 2)
initial.SI3 <- c(S = (N03-5), I = 5)
initial.SI4 <- c(S = (N04-10), I = 10)
initial.SI5 <- c(S = (N05-15), I = 15)
initial.SI6 <- c(S = (N06-20), I = 20)
time.out <- seq(0,150,by = 1) ##for simulating data from ODEs - output same with 0.1 step
time.samp <- seq(0,133, by = 7) ##will sample 20 time points
samp.sizes <- rep(100, length(time.samp)) ##sample 100 indiv each time

#################################################################################
##########################System of Diffeqs#####################################
#################################################################################
sir <- function(t,y,parms){
  with(c(as.list(y),parms),{
    dSdt <- -beta*(N^K)*S*I/N
    dIdt <- beta*(N^K)*S*I/N - gamma*I
    #Note: Population size is constant, so don't need to specify dRdt
    list(c(dSdt,dIdt))
  })
}

#################################################################################
#############NONLINEAR Neg Log Likelihood Function################################
#################################################################################
nll.fn <- function(pars,                # log(B) & log(gamma) in a vector, must be named
                   data1 = num.inf.samp1,             # number positive in sample
                   data2 = num.inf.samp2,
                   data3 = num.inf.samp3,
                   data4 = num.inf.samp4,
                   data5 = num.inf.samp5,
                   data6 = num.inf.samp6,
                   n = samp.sizes,                        # sample sizes
                   pop.size1 = N0,
                   pop.size2 = N02,
                   pop.size3 = N03,
                   pop.size4 = N04,
                   pop.size5 = N05,
                   pop.size6 = N06,
                   times = time.samp) # the times they were observed at
{
  pars <- c(exp(pars[1]), exp(pars[2]), exp(pars[3])) ##exponentiate parameters
  
  ##Epidemic 1
  ts.sir.temp1 <- data.frame(ode(y = initial.SI1, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size1)))
  ts.sir.temp1$P <- ts.sir.temp1$I / pop.size1
  prev.samp.temp1 <- ts.sir.temp1$P[ts.sir.temp1$time %in% times]
  prev.samp.temp1[prev.samp.temp1<0.0001]<-0.0001
  prev.samp.temp1[prev.samp.temp1>0.9999]<-0.9999
  
  ##Epidemic 2
  ts.sir.temp2 <- data.frame(ode(y = initial.SI2, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size2)))
  ts.sir.temp2$P <- ts.sir.temp2$I / pop.size2
  prev.samp.temp2 <- ts.sir.temp2$P[ts.sir.temp2$time %in% times]
  prev.samp.temp2[prev.samp.temp2<0.0001]<-0.0001
  prev.samp.temp2[prev.samp.temp2>0.9999]<-0.9999
  
  ##Epidemic 3
  ts.sir.temp3 <- data.frame(ode(y = initial.SI3, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size3)))
  ts.sir.temp3$P <- ts.sir.temp3$I / pop.size3
  prev.samp.temp3 <- ts.sir.temp3$P[ts.sir.temp3$time %in% times]
  prev.samp.temp3[prev.samp.temp3<0.0001]<-0.0001
  prev.samp.temp3[prev.samp.temp3>0.9999]<-0.9999
  
  ##Epidemic 4
  ts.sir.temp4 <- data.frame(ode(y = initial.SI4, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size4)))
  ts.sir.temp4$P <- ts.sir.temp4$I / pop.size4
  prev.samp.temp4 <- ts.sir.temp4$P[ts.sir.temp4$time %in% times]
  prev.samp.temp4[prev.samp.temp4<0.0001]<-0.0001
  prev.samp.temp4[prev.samp.temp4>0.9999]<-0.9999
  
  ##Epidemic 5
  ts.sir.temp5 <- data.frame(ode(y = initial.SI5, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size5)))
  ts.sir.temp5$P <- ts.sir.temp5$I / pop.size5
  prev.samp.temp5 <- ts.sir.temp5$P[ts.sir.temp5$time %in% times]
  prev.samp.temp5[prev.samp.temp5<0.0001]<-0.0001
  prev.samp.temp5[prev.samp.temp5>0.9999]<-0.9999
  
  ##Epidemic 6
  ts.sir.temp6 <- data.frame(ode(y = initial.SI6, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size6)))
  ts.sir.temp6$P <- ts.sir.temp6$I / pop.size6
  prev.samp.temp6 <- ts.sir.temp6$P[ts.sir.temp6$time %in% times]
  prev.samp.temp6[prev.samp.temp6<0.0001]<-0.0001
  prev.samp.temp6[prev.samp.temp6>0.9999]<-0.9999
  
  nll<-sum((-sum(dbinom(data1, n, prev.samp.temp1, log = TRUE))) +
             (-sum(dbinom(data2, n, prev.samp.temp2, log = TRUE))) +
             (-sum(dbinom(data3, n, prev.samp.temp3, log = TRUE))) +
             (-sum(dbinom(data4, n, prev.samp.temp4, log = TRUE))) +
             (-sum(dbinom(data5, n, prev.samp.temp5, log = TRUE))) +
             (-sum(dbinom(data6, n, prev.samp.temp6, log = TRUE)))
  )
  nll
}

##############################################################################
##########################DD NLL Function######################################
###############################################################################
nllDD.fn <- function(pars,                # log(B) & log(gamma) in a vector, must be named
                     data1 = num.inf.samp1,             # number positive in sample
                     data2 = num.inf.samp2,
                     data3 = num.inf.samp3,
                     data4 = num.inf.samp4,
                     data5 = num.inf.samp5,
                     data6 = num.inf.samp6,
                     n = samp.sizes,                        # sample sizes
                     pop.size1 = N0,
                     pop.size2 = N02,
                     pop.size3 = N03,
                     pop.size4 = N04,
                     pop.size5 = N05,
                     pop.size6 = N06,
                     times = time.samp) # the times they were observed at
{
  pars <- c(exp(pars[1]), exp(pars[2]), K=1) ##this is the only difference from the NL function (K=1)
  ##Epidemic 1
  ts.sir.temp1 <- data.frame(ode(y = initial.SI1, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size1)))
  ts.sir.temp1$P <- ts.sir.temp1$I / pop.size1
  prev.samp.temp1 <- ts.sir.temp1$P[ts.sir.temp1$time %in% times]
  prev.samp.temp1[prev.samp.temp1<0.0001]<-0.0001
  prev.samp.temp1[prev.samp.temp1>0.9999]<-0.9999
  ##Epidemic 2
  ts.sir.temp2 <- data.frame(ode(y = initial.SI2, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size2)))
  ts.sir.temp2$P <- ts.sir.temp2$I / pop.size2
  prev.samp.temp2 <- ts.sir.temp2$P[ts.sir.temp2$time %in% times]
  prev.samp.temp2[prev.samp.temp2<0.0001]<-0.0001
  prev.samp.temp2[prev.samp.temp2>0.9999]<-0.9999
  ##Epidemic 3
  ts.sir.temp3 <- data.frame(ode(y = initial.SI3, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size3)))
  ts.sir.temp3$P <- ts.sir.temp3$I / pop.size3
  prev.samp.temp3 <- ts.sir.temp3$P[ts.sir.temp3$time %in% times]
  prev.samp.temp3[prev.samp.temp3<0.0001]<-0.0001
  prev.samp.temp3[prev.samp.temp3>0.9999]<-0.9999
  ##Epidemic 4
  ts.sir.temp4 <- data.frame(ode(y = initial.SI4, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size4)))
  ts.sir.temp4$P <- ts.sir.temp4$I / pop.size4
  prev.samp.temp4 <- ts.sir.temp4$P[ts.sir.temp4$time %in% times]
  prev.samp.temp4[prev.samp.temp4<0.0001]<-0.0001
  prev.samp.temp4[prev.samp.temp4>0.9999]<-0.9999
  ##Epidemic 5
  ts.sir.temp5 <- data.frame(ode(y = initial.SI5, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size5)))
  ts.sir.temp5$P <- ts.sir.temp5$I / pop.size5
  prev.samp.temp5 <- ts.sir.temp5$P[ts.sir.temp5$time %in% times]
  prev.samp.temp5[prev.samp.temp5<0.0001]<-0.0001
  prev.samp.temp5[prev.samp.temp5>0.9999]<-0.9999
  ##Epidemic 6
  ts.sir.temp6 <- data.frame(ode(y = initial.SI6, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size6)))
  ts.sir.temp6$P <- ts.sir.temp6$I / pop.size6
  prev.samp.temp6 <- ts.sir.temp6$P[ts.sir.temp6$time %in% times]
  prev.samp.temp6[prev.samp.temp6<0.0001]<-0.0001
  prev.samp.temp6[prev.samp.temp6>0.9999]<-0.9999
  
  nll<-sum((-sum(dbinom(data1, n, prev.samp.temp1, log = TRUE))) +
             (-sum(dbinom(data2, n, prev.samp.temp2, log = TRUE))) +
             (-sum(dbinom(data3, n, prev.samp.temp3, log = TRUE))) +
             (-sum(dbinom(data4, n, prev.samp.temp4, log = TRUE))) +
             (-sum(dbinom(data5, n, prev.samp.temp5, log = TRUE))) +
             (-sum(dbinom(data6, n, prev.samp.temp6, log = TRUE)))
  )
  
  nll
}

#######################################################################
##########################FD FIT################################
#######################################################################
nllFD.fn <- function(pars,                # log(B) & log(gamma) in a vector, must be named
                     data1 = num.inf.samp1,             # number positive in sample
                     data2 = num.inf.samp2,
                     data3 = num.inf.samp3,
                     data4 = num.inf.samp4,
                     data5 = num.inf.samp5,
                     data6 = num.inf.samp6,
                     n = samp.sizes,                        # sample sizes
                     pop.size1 = N0,
                     pop.size2 = N02,
                     pop.size3 = N03,
                     pop.size4 = N04,
                     pop.size5 = N05,
                     pop.size6 = N06,
                     times = time.samp) # the times they were observed at
{
  pars <- c(exp(pars[1]), exp(pars[2]), K=0) ##this is the only difference from the NL function (K=0)
  ##Epidemic 1
  ts.sir.temp1 <- data.frame(ode(y = initial.SI1, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size1)))
  ts.sir.temp1$P <- ts.sir.temp1$I / pop.size1
  prev.samp.temp1 <- ts.sir.temp1$P[ts.sir.temp1$time %in% times]
  prev.samp.temp1[prev.samp.temp1<0.0001]<-0.0001
  prev.samp.temp1[prev.samp.temp1>0.9999]<-0.9999
  ##Epidemic 2
  ts.sir.temp2 <- data.frame(ode(y = initial.SI2, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size2)))
  ts.sir.temp2$P <- ts.sir.temp2$I / pop.size2
  prev.samp.temp2 <- ts.sir.temp2$P[ts.sir.temp2$time %in% times]
  prev.samp.temp2[prev.samp.temp2<0.0001]<-0.0001
  prev.samp.temp2[prev.samp.temp2>0.9999]<-0.9999
  ##Epidemic 3
  ts.sir.temp3 <- data.frame(ode(y = initial.SI3, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size3)))
  ts.sir.temp3$P <- ts.sir.temp3$I / pop.size3
  prev.samp.temp3 <- ts.sir.temp3$P[ts.sir.temp3$time %in% times]
  prev.samp.temp3[prev.samp.temp3<0.0001]<-0.0001
  prev.samp.temp3[prev.samp.temp3>0.9999]<-0.9999
  ##Epidemic 4
  ts.sir.temp4 <- data.frame(ode(y = initial.SI4, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size4)))
  ts.sir.temp4$P <- ts.sir.temp4$I / pop.size4
  prev.samp.temp4 <- ts.sir.temp4$P[ts.sir.temp4$time %in% times]
  prev.samp.temp4[prev.samp.temp4<0.0001]<-0.0001
  prev.samp.temp4[prev.samp.temp4>0.9999]<-0.9999
  ##Epidemic 5
  ts.sir.temp5 <- data.frame(ode(y = initial.SI5, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size5)))
  ts.sir.temp5$P <- ts.sir.temp5$I / pop.size5
  prev.samp.temp5 <- ts.sir.temp5$P[ts.sir.temp5$time %in% times]
  prev.samp.temp5[prev.samp.temp5<0.0001]<-0.0001
  prev.samp.temp5[prev.samp.temp5>0.9999]<-0.9999
  ##Epidemic 6
  ts.sir.temp6 <- data.frame(ode(y = initial.SI6, times = time.out,
                                 func = sir, parms = c(pars, N=pop.size6)))
  ts.sir.temp6$P <- ts.sir.temp6$I / pop.size6
  prev.samp.temp6 <- ts.sir.temp6$P[ts.sir.temp6$time %in% times]
  prev.samp.temp6[prev.samp.temp6<0.0001]<-0.0001
  prev.samp.temp6[prev.samp.temp6>0.9999]<-0.9999
  
  nll<-sum((-sum(dbinom(data1, n, prev.samp.temp1, log = TRUE))) +
             (-sum(dbinom(data2, n, prev.samp.temp2, log = TRUE))) +
             (-sum(dbinom(data3, n, prev.samp.temp3, log = TRUE))) +
             (-sum(dbinom(data4, n, prev.samp.temp4, log = TRUE))) +
             (-sum(dbinom(data5, n, prev.samp.temp5, log = TRUE))) +
             (-sum(dbinom(data6, n, prev.samp.temp6, log = TRUE)))
  )
  
  nll
}

###############################################################################
#######################Optimization Procedures###################################
###############################################################################
start_time <- Sys.time()
for (k in 1:length(ks)) {
  ##set up dataframe for output - one per K value
  factors<-expand.grid(fittingattempt=seq(1, nrestarts,1), dataset=seq(1, ndatasets, 1), K=ks[k])
  L<-length(factors$K)
  compareests<-cbind(factors, data.frame(initbetaNL=rep(NA,L), initgammaNL=rep(NA,L), initKNL=rep(NA,L), initbetaDD=rep(NA,L), initgammaDD=rep(NA,L), initbetaFD=rep(NA,L), initgammaFD=rep(NA,L),betaNL=rep(NA,L),gammaNL=rep(NA,L), KNL=rep(NA,L), betaDD=rep(NA,L), gammaDD=rep(NA,L), betaFD=rep(NA,L), gammaFD=rep(NA,L), lik=rep(NA,L), conv=rep(NA,L), likDD=rep(NA,L), convDD=rep(NA,L), likFD=rep(NA,L), convFD=rep(NA,L), NLLtrue=rep(NA, L), truebeta=rep(NA, L)))
  for (j in 1:ndatasets) {
    truebeta<-(FOI*N0)/((N0^ks[k])*1)
    ts.sir <- data.frame(ode(
      y = initial.SI1,               # Initial conditions for population
      times = time.out,             # Timepoints for evaluation
      func = sir,                   # Function to evaluate
      parms = c(beta=truebeta, gamma=gamma, N=N0, K=ks[k]),
      method="lsoda"))               # Vector of parameters
    ts.sir$P <- ts.sir$I / N0
    prev.samp1 <- ts.sir$P[ts.sir$time %in% time.samp]
    num.inf.samp1 <- rbinom(length(time.samp), size = samp.sizes, prob = prev.samp1)
    nll.true1 <- - sum(dbinom(num.inf.samp1, samp.sizes, prev.samp1, log = TRUE))
    ts.sir2 <- data.frame(ode(
      y = initial.SI2,               # Initial conditions for population
      times = time.out,             # Timepoints for evaluation
      func = sir,                   # Function to evaluate
      parms = c(beta=truebeta, gamma=gamma, N=N02, K=ks[k]),
      method="lsoda")) 
    ts.sir2$P <- ts.sir2$I / N02
    prev.samp2 <- ts.sir2$P[ts.sir2$time %in% time.samp]
    num.inf.samp2 <- rbinom(length(time.samp), size = samp.sizes, prob = prev.samp2)
    nll.true2 <- - sum(dbinom(num.inf.samp2, samp.sizes, prev.samp2, log = TRUE))
    ts.sir3 <- data.frame(ode(
      y = initial.SI3,               # Initial conditions for population
      times = time.out,             # Timepoints for evaluation
      func = sir,                   # Function to evaluate
      parms = c(beta=truebeta, gamma=gamma, N=N03, K=ks[k]),
      method="lsoda"))
    ts.sir3$P <- ts.sir3$I / N03  
    prev.samp3 <- ts.sir3$P[ts.sir3$time %in% time.samp]
    num.inf.samp3 <- rbinom(length(time.samp), size = samp.sizes, prob = prev.samp3)
    nll.true3 <- - sum(dbinom(num.inf.samp3, samp.sizes, prev.samp3, log = TRUE))
    ts.sir4 <- data.frame(ode(
      y = initial.SI4,               # Initial conditions for population
      times = time.out,             # Timepoints for evaluation
      func = sir,                   # Function to evaluate
      parms = c(beta=truebeta, gamma=gamma, N=N04, K=ks[k]),
      method="lsoda"))
    ts.sir4$P <- ts.sir4$I / N04  
    prev.samp4 <- ts.sir4$P[ts.sir4$time %in% time.samp]
    num.inf.samp4 <- rbinom(length(time.samp), size = samp.sizes, prob = prev.samp4)
    nll.true4 <- - sum(dbinom(num.inf.samp4, samp.sizes, prev.samp4, log = TRUE))
    ts.sir5 <- data.frame(ode(
      y = initial.SI5,               # Initial conditions for population
      times = time.out,             # Timepoints for evaluation
      func = sir,                   # Function to evaluate
      parms = c(beta=truebeta, gamma=gamma, N=N05, K=ks[k]),
      method="lsoda"))
    ts.sir5$P <- ts.sir5$I / N05  
    prev.samp5 <- ts.sir5$P[ts.sir5$time %in% time.samp]
    num.inf.samp5 <- rbinom(length(time.samp), size = samp.sizes, prob = prev.samp5)
    nll.true5 <- - sum(dbinom(num.inf.samp5, samp.sizes, prev.samp5, log = TRUE))
    ts.sir6 <- data.frame(ode(
      y = initial.SI6,               # Initial conditions for population
      times = time.out,             # Timepoints for evaluation
      func = sir,                   # Function to evaluate
      parms = c(beta=truebeta, gamma=gamma, N=N06, K=ks[k]),
      method="lsoda"))
    ts.sir6$P <- ts.sir6$I / N06  
    prev.samp6 <- ts.sir6$P[ts.sir6$time %in% time.samp]
    num.inf.samp6 <- rbinom(length(time.samp), size = samp.sizes, prob = prev.samp6)
    nll.true6 <- - sum(dbinom(num.inf.samp6, samp.sizes, prev.samp6, log = TRUE))
    for (i in 1:nrestarts) {
      print(c("NEW PARAMS", "K=", ks[k], "Dataset #", j, "Restart #", i))
      ##Random starts
      beta1<-log(rlnorm(1, mean=log(truebeta), sdlog=1))
      gamma1<-log(rlnorm(1, mean=log(gamma), sdlog=1))
      K1<-log(rlnorm(1, mean=log(0.1), sdlog=1));
      beta2<-log(rlnorm(1, mean=log((FOI*N0)/(N0^1)), sdlog=1))
      #beta2<-log(rlnorm(1, mean=log(truebeta), sdlog=1))
      gamma2<-log(rlnorm(1, mean=log(gamma), sdlog=1));
      beta3<-log(rlnorm(1, mean=log((FOI*N0)/(N0^0)), sdlog=1))
      gamma3<-log(rlnorm(1, mean=log(gamma), sdlog=1));
      startpar<-c(beta=beta1, gamma=gamma1, K=K1)
      startpar2<-c(beta=beta2, gamma=gamma2)
      startpar3<-c(beta=beta3, gamma=gamma3)
      ##optimization procedures
      outNL<-optim(startpar, nll.fn, control = list(trace = 0, maxit = 1000), method = "Nelder-Mead")
      outDD<-optim(startpar2, nllDD.fn, control = list(trace = 0, maxit = 1000), method= "Nelder-Mead")
      outFD<-optim(startpar3, nllFD.fn, control = list(trace = 0, maxit = 1000), method = "Nelder-Mead")
      ##save output
      AssignmentRows<-which(compareests$K==ks[k] & compareests$fittingattempt==i & compareests$dataset==j)
      compareests[AssignmentRows, c(4:10)]<-c(startpar, startpar2, startpar3) 
      compareests[AssignmentRows, c(11:17)]<-exp(c(as.numeric(outNL$par), as.numeric(outDD$par), as.numeric(outFD$par)))
      compareests$lik[AssignmentRows]<-outNL$value  #negative log likelihood value
      compareests$conv[AssignmentRows]<-outNL$convergence  #if not 0, didn't converge
      compareests$likDD[AssignmentRows]<-outDD$value  #negative log likelihood value
      compareests$convDD[AssignmentRows]<-outDD$convergence  #if not 0, didn't converge
      compareests$likFD[AssignmentRows]<-outFD$value  #negative log likelihood value
      compareests$convFD[AssignmentRows]<-outFD$convergence  #if not 0, didn't converge
      compareests$NLLtrue[AssignmentRows]<-sum(nll.true1, nll.true2, nll.true3, nll.true4, nll.true5, nll.true6)
      compareests$truebeta[AssignmentRows]<-truebeta
    }
  }
  #output one dataframe per K value to a CSV to save
  write.csv(compareests, paste("/Users/hopkins/Documents/Transmission Function Literature Review/TransmissionFunctions/compareests","K",ks[k],"FOI",FOI,"gamma",gamma,".csv",sep=""), row.names=F)
  #gc()
}
end_time <- Sys.time()
end_time - start_time 
