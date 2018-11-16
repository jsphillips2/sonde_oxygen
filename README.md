Sonde Oxygen
========

Analyses associated with the manuscript "Time-varying responses of lake metabolism to light and temperature."
-------

Joseph Phillips
Integrative Biology, UW-Madison


## Description

This repository contains analyses for inferring temporal variation in photosynthesis-irradiance (P-I) curves in the context of ecosystem metabolism for Lake Mývatn in northern Iceland. The analyses employ an extension of previous methods for inferring ecosystem metabolism from dissovled oxgyen data, with the main innovation being the explicit modeling of temporal variation in the parameters of the P-I curve (e.g. the maximum rate of GPP at high light). This takes advantage of the fact that the physical and biological processes governing ecosystem metabolism and other aspects of DO dynamics are correlated through time, which means that this shared information can be used to inform the parameter estimates across all time points. 

The model is fit in a Bayesian framework using Stan, run from R using the 'rstan' package. This repository contains all of the raw data and code for reproducing the analyses. 

## Contents

* `data`: Raw data files and code for processing data.

* `analyses`: Primary analyses, including code for fitting the model using Stan.

* `model`: Files for specifying the model in Stan. The version `o2_model.stan` takes the obervation error standard deviation as an input (set to a low value), while `o2_model_sig_obs.stan` estimates the standard deviation from the data.

* `simulations`: Code for generating simulated data and fitting the model to those data.

* `supplement`: Rmd file for generaitng supplementary materials associated with manuscript.   