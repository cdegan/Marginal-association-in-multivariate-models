# Marginal association in multivariate models 

Code used for the manuscript "*Estimating the evolution of the marginal association between a continuous time-varying covariate and a longitudinal non-Gaussian outcome in a joint modeling framework*". 

## Overview

This repository contains the code used for the joint scaled models (JSM) analyses and the marginal association derivation described in the associated manuscript. The project includes simulation scripts, model estimation functions, and saved model objects used to reproduce the reported results.

## Repository structure

### R
In this folder we included all the codings files used for the simulations (and can also be used to reproduce the analysis in real data):
-   *Data_parameters.R* file containing all the fixed parameters used for the simulation (R file);
-   *Data_creation.R* file containing the functions used to create the datasets for the simulation (R file);
-   *models_preparation.R* file containing the functions used to prepare estimation of JMM and JSM (R file);
-   *JMM.R* file containing the specific functions used to estimate the JMM (R file);
-   *JSM.R* file containing the specific functions used to estimate the JSM (R file);
-   *Ass_coeff_derivation.R* file containing the functions used to estimate the association coefficient in both JMM and JSM (R file);
-   *Dependencies.R* file containing the list of packages useful to run the code.

### Validation marginal association

It contains all the files that are needed to compile the simulation that evaluate the validity of the proposed methods for estimating the marginal
association. This simulation is shown in the manuscript (Section 2.1 of the Supplementary Material):
-   *Data_parameters_gaussian.R* file containing all the fixed parameters used for this simulation (R file);
-   *Data_creation_gaussian.R* file containing the functions used to create the datasets for this simulation (R file);
-   *Simulated models* folder containing the files with the saved simulated estimated models (*models_JSM_gaussian_200.RData*).
-   *Simulation.Rmd* file containing the code used to implement the simulation and look at the result (R-markdown file).

### Simulation 65 vs 200

It contains all the files that are needed to compile the simulation that compare the models estimation when assuming small or big sample size. This simulation is shown in the manuscript (Section 2.2.2 of the Supplementary Material).
  -   *Simulated models* folder containing the files with the saved simulated estimated JSM models, distinguished between 65 and 200 patients (*models_JSM_beta_65.RData*, etc.).
  -   *Simulation_65vs200.Rmd* file containing the code used to implement the simulation and look at the result (R-markdown file).

### Simulation Default vs Informative

It contains all the files that are needed to compile the simulation that compare the models estimation with or without an informative hyperprior for the variance parameters. This simulation is also described in the manuscript (Section 2.2.3 of the Supplementary Material).
  -   *Simulated models* folder containing the files with the saved simulated estimated JSM models, distinguished between 65 and 200 patients using the informative hyperprior setting.
  -   *Simulation_DefaultvsInformative.Rmd* file containing the code used to implement the simulation and look at the result (R-markdown file).

## Reproducibility

To reproduce the simulations and analyses:

1. Download all the R scripts in the `R` folder and upload them all in the same folder in which you want to do the analysis;
2. Download and open the relevant simulation rmarkdown(s) in the simulation folders;
3. Run the code to reproduce the figures, model comparisons, and saved outputs.

## Requirements

This project is implemented in R, R-version 4.2.2. Before running the scripts, ensure that the required R packages (listed in */R/Dependences.R*) are installed. The list of the versions of the packages is listed in the file *main/sessionInfo.txt*.

## Citation

If you use this code in your work, please cite the corresponding manuscript.
  
