Sonde Oxygen
========

Analyses associated with the manuscript "Time-varying responses of lake metabolism to light and temperature."
-------

Joseph Phillips
Integrative Biology, UW-Madison


## Description

This repository contains analyses for inferring temporal variation in photosynthesis-irradiance (P-I) curves in the context of ecosystem metabolism for Lake Mývatn in northern Iceland. The analyses employ an extension of previous methods for inferring ecosystem metabolism from dissovled oxgyen data, with the main innovation being the explicit modeling of temporal variation in the parameters of the P-I curve (e.g. the maximum rate of GPP at high light). This takes advantage of the fact that the physical and biological processes governing ecosystem metabolism and other aspects of DO dynamics are correlated through time, which means that this shared information can be used to inform the parameter estimates across all time points. 

The model is fit in a Bayesian framework using Stan, run from R using the 'rstan' package. This repository contains all of the raw data and code for reproducing the analyses. 

## Key Contents

* `data`: raw data files and code for processing data 
  -`clean_data.R`: code for combining, cleaning, and filtering sensor data
  -`sonde_final.csv`: the final output from the data processing
  
* `analyses`: primary analyses
  -`model_fit`: primary model fitting
    -`input`: data files formatted for model fitting, with associated code
    -`fit_model`: code for fitting model using Stan (model file itself is contained in `model`)
    -`output`: main output of model fit