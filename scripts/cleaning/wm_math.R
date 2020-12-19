# Working memory

# Scores for working memory are not homogenous, because Spanish speakers always obtain lower scores 
# (My hypothesis is that the lower scores are due to the need to remember gender)

# So, since anticipation is closely linked to processing speed, which is another component of working memory,
# we are going to use other OSpan measures to ensure that populations are homogeneous.
# Namely, we are going to use processing speed in match while ignoring storage (which is the score we tested first).

library(data.table)
library(tidyverse)

# Load all data and bind them into one dataframe (model: https://iamkbpark.com)
csv_files <- list.files (path       = "./data/wm", 
                         pattern    = "*.csv", 
                         full.names = T)

wm_ENES <- as_tibble (rbindlist (lapply (csv_files, fread)))

# Mean of RT for each participant
wm_RT <- wm_ENES %>%
  filter(., subj_form_resp == correct_resp) %>%
  group_by(., subject_id) %>%
  summarize(., mean_rt = mean(rt_formula)) %>%
  write.csv(., "./data/clean/wm_EN_ES.csv")


  


  
  