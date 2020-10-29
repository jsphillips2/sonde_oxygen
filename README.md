Sonde Oxygen
========

Analyses associated with the manuscript "Time-varying responses of lake metabolism to light and temperature."
-------

Joseph Phillips, Integrative Biology, UW-Madison

## Update: 29 October 2020

A sign error in the code associated with the calculation of the temperature dependence 
of the Schmidt number has been corrected (see lines 364-368 of `data/clean_data.R`). 
The code now runs for the corrected version, 
but with the original preserved for future reference.
All downstream files associated with the update are marked as such,
and otherwise correspond to the originals.
The error had only a very minor effect on the model estimates 
and no effect on the biologcial inferences,
due to the (a) inclusion of process error and (b) the model's treatment
of biological processes as temporally autocorrelated. 
Deviations due to the erroneous temperature dependence of the Schmidt number were 
therefore not attributed to the biological processes of primary interest.
I have contacted the journal and a corrigendum is currently under consideration.
The error was discovered by a colleague, to whom I am greatly indebted,
who sought to reproduce a portion of the results.

## Description

This repository contains analyses for inferring temporal variation in photosynthesis-irradiance (P-I) curves in the context of ecosystem metabolism for Lake Mývatn in northern Iceland. The analyses employ an extension of previous methods for inferring ecosystem metabolism from dissolved oxygen (DO) data, with the main innovation being the explicit modeling of temporal variation in the parameters of the P-I curve (e.g. the maximum rate of GPP at high light). This takes advantage of the fact that the physical and biological processes governing ecosystem metabolism and other aspects of DO dynamics are correlated through time, which means that this shared information can be used to inform the parameter estimates across all time points. 

The model is fit in a Bayesian framework using Stan, run from R using the 'rstan' package. This repository contains all of the raw data and code for reproducing the analyses. 

## Contents

* `data`: Raw data files and code for processing data.

* `analyses`: Analyses of the model fits, including code for figures, variance partitioning, and additional calculations.

* `model`: Files for specifying and fitting the model in Stan. Different input and output folders (`alt_k`,`fixed`,`main`,`sig_obs`, and `surface_par`) correspond to various model fits, differing either due to the input data add/or the Stan file used to specify the model.  All of the model fits can be run from the file `fit_model.R`. The folder `stan` contains files specifying the structure of different versions of the model.  The file `o2_model.stan` specifies the model used for the main analysis, including temporal variation in the parameters of the photosynthesis-irradiance curve and takes the observation error standard deviation as an input (set to a low value). The file  `o2_model_sig_obs.stan` is identical to `o2_model.stan`, except that it includes the observation error standard deviation as a parameter to estimate from the data. The file  `o2_model_fixed.stan` fits the model with the parameters of the photosynthesis-irradiance curve fixed through time.

* `supplement`: Rmd file for generating supplementary materials associated with manuscript.   
