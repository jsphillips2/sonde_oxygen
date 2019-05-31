#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)
library(nlme)





#==========
#========== Import data
#==========

# set paths for main and simulation analyses
main_path <- "model/output/"

# data
sonde_data <- read_csv(paste0("model/input/main/sonde_prep.csv"))

# main analysis
model_fit <- read_csv(paste0("model/output/main/summary_clean.csv"))

# midges
midges <- read_csv("data/midges.csv")





#==========
#========== Alpha and light compensation
#==========

# exctract alpha and mean daily light
alpha_light <- model_fit %>%
  filter(name %in% c("alpha")) %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
  mutate(middle = middle,
         lower16 = lower16,
         upper84 = upper84) %>% 
  left_join(sonde_data %>% 
              group_by(year, yday) %>%
              summarize(par = mean(par),
                        par_int = mean(par_int)))

# examine & export
cor.test(alpha_light$middle, alpha_light$par_int)






#==========
#========== Drivers of max GPP 
#==========

#### all phycocyanin observtions (no midges)

# prepare data 
beta0_phyc <- model_fit %>%
  filter(name %in% c("beta0")) %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              filter(pcyv < 0.3) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day),
                        pcyv = mean(pcyv))) %>%
  select(year, yday, middle, pcyv) %>%
  rename(beta0 = middle) 

# fit model 
m_phyc <- gls(log(beta0) ~ pcyv, correlation = corCAR1(form = ~ yday|year), data = beta0_phyc)
summary(mm_phyc1)
anova(m_phyc)



#### midges and phycocyanin (weekly measurements)
 
# prepare midge data
midges_summary = midges %>%
  filter(sta %in% c(3, 33)) %>%
  mutate(year = year(sampledate),
         yday = yday(sampledate)) %>%
  group_by(year, yday, coreid) %>%
  summarize(tanyt = sum(tanyt/fract_count),
            chiro = sum(chiro/fract_count),
            midges = tanyt + chiro) %>%
  group_by(year, yday) %>%
  summarize(tanyt = mean(tanyt, na.rm=T),
            chiro = mean(chiro, na.rm=T),
            midges = mean(midges, na.rm=T))

# combine 
beta0_phyc_midge = beta0_phyc %>%
  left_join(midges_summary) %>%
  na.omit()

# fit model 
m_phyc_midge = gls(log(beta0) ~ pcyv + midges, 
                   correlation = corCAR1(form = ~ yday|year), data = beta0_phyc_midge)
summary(m_phyc_midge)
anova(m_phyc_midge)





#==========
#========== Correlation between GPP and ER
#==========

# overall rates
model_fit %>%
  filter(name %in% c("GPP","ER")) %>%
  select(name, middle, day) %>%
  spread(name, middle) %>%
  {cor.test(.$GPP, .$ER)}

# max GPP and baseline ER (all years)
model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  select(name, middle, day) %>%
  spread(name, middle) %>%
  {cor.test(.$beta0, .$rho)}

# max GPP and baseline ER (by year)
model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  select(name, middle, day) %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
  spread(name, middle) %>%
  split(.$year) %>%
  lapply(function(x){cor.test(x$beta0, x$rho)})





#==========
#========== DO change explained by NEP
#==========

# prepare data
flux_d <- model_fit %>%
  filter(name %in% c("nep","air")) %>%
  select(name, day, index, middle) %>%
  spread(name, middle) %>%
  arrange(day, index) %>%
  select(-day, -index) %>%
  bind_cols(sonde_data %>% 
              na.omit() %>%
              arrange(year, yday, hour)) %>%
  select(year, yday, hour, par_int, do, nep, air) %>%
  arrange(year, yday, hour) %>%
  group_by(year, yday) %>%
  mutate(air = air/3300,
         nep = nep/3300,
         flux = c(diff(do), NA),
         error = flux - (nep + air)) %>%
  na.omit()

# calculate proportional contribution of nep
nep_var <- var(flux_d$nep) + cov(flux_d$nep, flux_d$air) + cov(flux_d$nep, flux_d$error) 
nep_var/var(flux_d$flux, na.rm=T)




#==========
#========== Variation in P-I parameters
#==========

model_fit %>% 
  filter(name %in% c("beta0","alpha","rho")) %>%
  group_by(name) %>%
  summarize(fold = max(middle)/min(middle))




#==========
#========== "Smoothed" days
#==========

fit_days <- model_fit %>% 
  filter(name == "beta0") %>% nrow()

expand_nep <- model_fit %>% 
  filter(name == "NEP") %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
  full_join(sonde_data %>% expand(year,yday,name=c("NEP"))) 

day_bounds <- expand_nep %>%
  na.omit() %>%
  group_by(year) %>%
  summarize(min = min(yday), max = max(yday))

tot_days <- expand_nep %>%
  group_by(year) %>%
  left_join(day_bounds) %>%
  filter(!(yday < min), !(yday > max)) %>% nrow()

tot_days-fit_days





#==========
#========== Model comparison
#==========

looic_main <- readRDS("model/output/main/loo.rds")
looic_fixed <- readRDS("model/output/fixed/loo.rds")





#==========
#========== Overall mean metabolism parameters
#==========

model_fit %>% 
  filter(name %in% c("beta0","alpha","rho")) %>%
  group_by(name) %>%
  summarize(middle = median(middle),
            lower16 = median(lower16),
            upper84 = median(upper84))





