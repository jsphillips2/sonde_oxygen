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

* `analyses`: Primary analyses, including code for figures, variance partitioning, and additional calculaitons.

* `model`: Files for specifying and fitting the model in Stan. Different input and output folders (`alt_k`,`fixed`,`main`,`sig_obs`, and `surface_par`) correspond to various model fits, differing either due to the input data add/or the Stan file used to specify the model.  All of the model fits can be run from the file `fit_model.R`. The folder `stan` contains files specifying the structure of different versions of the model.  The file `o2_model.stan` specifies the model used for the main analysis, including temporal variation in the parameters of the photosynthesis-irradiance curve and takes the obervation error standard deviation as an input (set to a low value). The file  `o2_model_sig_obs.stan` is identical to `o2_model.stan`, except that it includes the observation error standard deviation as a parameter to estimate from the data. The file  `o2_model_fixed.stan` fits the model with the parameters of the photosynthesis-irriadiance curve fixed through time.

* `supplement`: Rmd file for generaitng supplementary materials associated with manuscript.   
