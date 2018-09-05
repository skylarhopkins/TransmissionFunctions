# Transmission functions in epidemiological models - when do they create parameter estimate biases?
This repo contains a model simulation tool that can illustrate bias in parameter estimates when linear transmission functions are used to approximate nonlinear transmission dynamics or when the linear, density-dependent transmission function is used with truly frequency-dependent transmission and vice versa.

The primary script file, TransmissionFunctionSimulationScript_BroadHostRange.R, generates epidemics in six popululations with six different host densities, assuming that the underlying transmission-density relationship is frequency-dependent (K=0), density-dependent (K=1), or nonlinear (0<K<1). The script then fits the three models (FD, DD, NL) to all generated data and records the resulting parameter estimates (Beta, gamma, K) and model fits (AICs). These estimates can be compared to the true values to calculate bias.

The many .txt files contain simulation runs over 9 combinations of FOI and gamma and 11 values of K (0 to 1 by intervals of 0.1). For each of these 99 parameter combinations (contained in 99 .txt files), we created simulated 100 datasets, for a total of 9900 'datasets'. The output in the 99 .txt files can be used to plot a figure illustrating model fits and parameter estimate biases. We will push the script that creates this figure to Github after publication.    

Finally, we also include a short script file to demonstrate how K, the density dependence parameter, influences the shape of the nonlinear transmission function that we use: B*N^K.
