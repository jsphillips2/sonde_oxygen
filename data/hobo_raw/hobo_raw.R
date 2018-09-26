#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)





#==========
#========== Clean data
#==========

# specify year
year = 2018

# import all files
data = list.files(paste0("data/hobo_raw/raw_data/",year), full.names = T) %>%
  lapply(read_csv, skip = 1)
data[[1]]

# rename and select relevant columns
data_slim = data %>%
  lapply(function(x){
    x %>%
      rename(date_time = grep("Date",names(data[[1]])),
             lux = grep("Intensity",names(data[[1]]))) %>%
      select(date_time, 
             lux) %>%
      return()
  }) %>%
  bind_rows()

# examine
data_slim

# plot lux
data_slim %>%
  ggplot(aes(date_time, lux))+
  geom_line()+
  theme_bw()

# export
write_csv(data_slim, paste0("data/hobo_raw/clean_data/hobo_clean_",year,".csv"))
