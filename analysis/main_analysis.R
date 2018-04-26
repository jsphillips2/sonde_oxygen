#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)

# import data and model fit
sonde_data = read_csv("data/sonde_prep.csv")
params_full = read_csv("model_output/o2_model/fixed_pars_full.csv")
post_pred = read_csv("model_output/o2_model/post_pred_full.csv")
model_fit = read_csv("model_output/o2_model/summary_clean.csv")
