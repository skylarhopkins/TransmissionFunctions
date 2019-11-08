# Transmission functions in epidemiological models - when do they create parameter estimate biases?
This repo contains a model simulation tool that can illustrate bias in parameter estimates when linear transmission functions are used to approximate nonlinear transmission dynamics or when the linear, density-dependent transmission function is used with truly frequency-dependent transmission and vice versa. This simulation tool was created for an analysis in a manuscript in revision at Methods in Ecology and Evolution, and I later turned it into a package so that more flexible simulations could be run by people interested in evaluating bias in their own epidemiological models.

The older, primary script file, TransmissionFunctionSimulationScript_BroadHostRange.R, generates epidemics in six popululations with six different host densities, assuming that the underlying transmission-density relationship is frequency-dependent (K=0), density-dependent (K=1), or nonlinear (0<K<1). The script then fits the three models (FD, DD, NL) to all generated data and records the resulting parameter estimates (Beta, gamma, K) and model fits (AICs). These estimates can be compared to the true values to calculate bias. All relevant functions are defined inside this script file.

The subdirectory, contactr, contains a package with five functions that accomplish the same goal as TransmissionFunctionSimulationScript_BroadHostRange.R, but allow for more flexibility in choosing population sizes. The sim.and.opt() function is the wrapper that runs everything. Below is an example of the work flow, which is also included in the documentation for sim.and.opt():

## Load the package
install.packages("devtools")
library("devtools")

devtools::install_github("skylarhopkins/TransmissionFunctions", subdir="contactr") #does it work? I hope so.
library("contactr")

#If you assign these outside of the function arguments, you cannot use the same names as the
#argument names, or else you'll get a recursive argument error. Sticking a period after each
#name solves this problem

#number of fits of each model to each dataset with different random starting parameters
nrestarts. = 7
#number of sample datasets for each value of K that you're simulating over
ndatasets. = 100
#The values of K (unitless density-dependence parameter) that you're simulating over
ks.<-seq(0.0, 1.0, 0.1)
#The Force of Infection value that you are simulating. If you want to do more than one at a time, you'll need to write a loop
FOI.<-0.001
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

## Run the tool
#WARNING: this can take a very long time to run depending on ndatasets, nrestarts, and length of ks,
#so you might want to estimate run times (end_time - start_time) on a smaller subset first

start_time <- Sys.time()

sim.and.opt(FOI=FOI., truegamma=truegamma., pops=pops., Nref=Nref., ks=ks., initial.I=initial.I., initial.S=initial.S., time.out=time.out., time.samp=time.samp., samp.sizes=samp.sizes., nrestarts=nrestarts., ndatasets=ndatasets., outputlocation=outputlocation.)

end_time <- Sys.time()
end_time - start_time
