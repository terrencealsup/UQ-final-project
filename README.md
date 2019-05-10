# Stochastic Modeling and UQ Final Project
Final project on MFCE for rare-event estimation for UQ class.
"Multifidelity Cross-Entropy Estimation of Rare Events in Steady Heat Conduction"

Authors: Frederick Law and Terrence Alsup

## Description of the code.
Most Matlab functions and scripts are contained in the directory titled "code". The main Matlab files are:

1. CrossEntropy.m
2. MFCE.m
3. CE_LogNormal.m
4. MFCE_LogNormal.m

Files 1. and 2. implement CE and MFCE with a family of Gaussian biasing densities.
Files 3. and 4. implement CE and MFCE with a family of Log-normal biasing densities.
To see how to run the code check the documentation of these files.

The "ellip" directory contains the necessary functions and files for partitioning the domain and solving the PDE using finite elements.  Data from our test runs are contained in the directory "sim data".


Several additional Matlab scripts are included for some test cases.  Here are some quick descriptions of the functions and scripts used to generate data and plots.

1. CE_biasing_density.m -- plots the biasing density found from cross-entropy as well as the optimal biasing density.
2. sqcov_plot.m -- reads in simulated data to make loglog plot of MFCE and CE SQCoV vs. runtime
3. heat_test.m -- builds the contour plots
4. heat.m -- builds the contour plots
5. CoeffVar_test.m -- estimates the SQCoV against a CE reference
6. CE_LN_PDE_test.m -- compare CE against the Monte Carlo references
7. build_reference.m -- build the Monte Carlo references
8. build_CE_reference.m -- build the CE references


