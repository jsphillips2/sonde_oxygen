#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)





#==========
#========== Clean data
#==========

# specify year
year = 2018

# import all files
data = list.files(paste0("data/sonde_raw/raw_data/",year), full.names = T) %>%
  lapply(read_lines)

# examine header to determine lines of meta data
head(data[[1]], 20)

# meta data lines
# meta_lines = c(-2,-3) #2012
# meta_lines = c(-2) # 2013
# meta_lines = c(-1:-12, -14:-15) # 2015-2017

# remove metadata
data_slim = data %>%
  lapply(function(x){
    x[meta_lines]
  })
head(data_slim[[1]],5)

# clean
data_csv = data_slim %>%
  lapply(function(x){
    # total number of columns (including blanks)
    ncols = length(str_split(x[1],",")[[1]])
    
    # logical vector for non-empty column positions
    inds = str_split(x[1],",")[[1]] %>%
      `!=`("\"\"")
    
    # clean 
    x_clean = x %>% 
      str_split(",") %>%
      # omit lines with fewer than the appropriate number of columns  
      discard(~length(.x) < ncols) %>%
      # select non-empty columns
      map(~.x[inds]) %>%
      map_chr(~str_c(.x, collapse = ","))
    
    # return as csv
    return(str_c(x_clean, collapse = "\n") %>%
             read_csv()
    )
  }) %>%
  bind_rows() %>%
  mutate(Date_Time = ymd_hms(paste(mdy(Date), Time))) %>%
  select(Date_Time, Temp, SpCond, TurbSC, `LDO%`, LDO, PCYV)

# examine
data_csv

# plot
data_csv %>%
  mutate(index = 1:nrow(data_csv)) %>%
  ggplot(aes(index, LDO))+
  geom_line()+
  theme_bw()

# export
# write_csv(data_csv, paste0("data/sonde_raw/clean_data/sonde_clean_",year,".csv"))



