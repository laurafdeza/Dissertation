# Working memory

# Scores for working memory are not homogeneous—the highest mean score is for the EN, then MA, and then ES
# (My hypothesis is that the lower scores are due to the need to remember gender and concepts rather than just the word)

# So, since anticipation is closely linked to processing speed, which is another component of working memory,
# we are going to use other OSpan measures to ensure that populations are homogeneous.
# Namely, we are going to use processing speed in match while ignoring storage (which is the score we tested first).

# Load packages
library(data.table)
library(tidyverse)
library(lmerTest)

# Load visuospatial WM data and bind them into one dataframe
csv_files <- list.files (path       = "./data/corsi", 
                         pattern    = "*.csv", 
                         full.names = T)

corsi_speed <- as_tibble (rbindlist (lapply (csv_files, fread)))

# Select correct trials, rename participants' IDs and participant column
corsi_speed <- corsi_speed %>%
  filter(., correct == 1) %>%
  rename(., participant = subject_id)

corsi_speed$participant <- str_replace(corsi_speed$participant, "ae", "aes")
corsi_speed$participant <- str_replace(corsi_speed$participant, "ie", "ies")
corsi_speed$participant <- str_replace(corsi_speed$participant, "am", "ams")
corsi_speed$participant <- str_replace(corsi_speed$participant, "im", "ims")
corsi_speed$participant <- str_replace(corsi_speed$participant, "mo", "mon")

# Find random effects for each participant
corsi_glm <- lmer(time ~ level + (1 | participant),
                   data = corsi_speed)

corsi_ranef <- ranef(corsi_glm) %>% as_tibble()       

corsi_sel <- corsi_ranef %>%
  select(., grp, condval, condsd) %>%
  rename(., participant = grp,
         corsi_rt = condval,
         corsi_sd = condsd)



### Repeat process for verbal WM  

# Load data  
csv_files <- list.files (path       = "./data/ospan", 
                         pattern    = "*.csv", 
                         full.names = T)

ospan_speed <- as_tibble (rbindlist (lapply (csv_files, fread)))

# Select correct trials, rename participants' IDs 
ospan_speed <- ospan_speed %>%
  filter(., correct_resp == 1) 

ospan_speed$subject_id <- str_replace(ospan_speed$subject_id, "ae", "aes")
ospan_speed$subject_id <- str_replace(ospan_speed$subject_id, "ie", "ies")
ospan_speed$subject_id <- str_replace(ospan_speed$subject_id, "am", "ams")
ospan_speed$subject_id <- str_replace(ospan_speed$subject_id, "im", "ims")
ospan_speed$subject_id <- str_replace(ospan_speed$subject_id, "mo", "mon")

# Calculate random effects for each participant
ospan_glm <- lmer(rt_formula ~ seq_length + (1 | subject_id),
                  data = ospan_speed)

ospan_ranef <- ranef(ospan_glm) %>% as_tibble()       

ospan_sel <- ospan_ranef %>%
  select(., grp, condval, condsd) %>%
  rename(., participant = grp,
         ospan_rt = condval,
         ospan_sd = condsd)

### merge dataframes and save resulting one as a .csv for analysis
wm_speed <- merge(corsi_sel, ospan_sel, by="participant")

write.csv(wm_speed,'./data/clean/wm_processing_speed.csv', row.names = F)









# # Load all data and bind them into one dataframe (model: https://iamkbpark.com)
# csv_files <- list.files (path       = "./data/wm", 
#                          pattern    = "*.csv", 
#                          full.names = T)
# 
# wm_ENES <- as_tibble (rbindlist (lapply (csv_files, fread)))
# 
# # Mean of RT for each participant
# wm_RT <- wm_ENES %>%
#   filter(., subj_form_resp == correct_resp) %>%
#   group_by(., subject_id) %>%
#   summarize(., mean_rt = mean(rt_formula)) %>%
#   write.csv(., "./data/clean/wm_EN_ES.csv")


  


  
  