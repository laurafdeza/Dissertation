source(here::here("scripts", "00_load_libs.R"))

pitch <- read_csv(here("data", 'clean', "pitch_long.csv"))
pitch <- pitch %>%
  select(., -X1) %>%
  rename(., pitch_rt = rt)


# rt ~ cond + (1 | id)

pitch_glm <- lmer(pitch_rt ~ base_note + direction + (1 | participant),
                       data = pitch)

pitch_ranef <- ranef(pitch_glm) %>% as_tibble()                   


rhythm <- read_csv(here("data", 'clean', "rhythm_long.csv"))
rhythm <- rhythm %>%
  select(., -X1) %>%
  rename(., rhythm_time_dev = mean_dev,
         rhythm_cond = condition)

rhythm_glm <- lmer(rhythm_time_dev ~ rhythm_cond + (1 | participant),
                  data = rhythm)

rhythm_ranef <- ranef(rhythm_glm) %>% as_tibble() 

pitch_sel <- pitch_ranef %>%
  select(., grp, condval, condsd) %>%
  rename(., participant = grp,
         pitch_dev = condval,
         pitch_sd = condsd)

rhythm_sel <- rhythm_ranef %>%
  select(., grp, condval, condsd) %>%
  rename(., participant = grp,
         rhythm_dev = condval,
         rhythm_sd = condsd)

auditory <- merge(pitch_sel, rhythm_sel, by="participant")

write.csv(auditory,'./data/clean/auditory_scores.csv')




car <- read_csv(here("data", 'clean', 'car_sec_long.csv'))

car <-  car %>%
  select(., -X1)

car_glm <-  lmer(mean_dev_sc ~ speed + direction + (1 | subject_id),
                 data = car)

car_ranef <- ranef(car_glm) %>% as_tibble()  

car_sel <- car_ranef %>%
  select(., grp, condval, condsd) %>%
  rename(., participant = grp,
         car_dev = condval,
         car_sd = condsd)

car_sel$participant <- str_replace(car_sel$participant, "ae", "aes")
car_sel$participant <- str_replace(car_sel$participant, "ie", "ies")
car_sel$participant <- str_replace(car_sel$participant, "am", "ams")
car_sel$participant <- str_replace(car_sel$participant, "im", "ims")
car_sel$participant <- str_replace(car_sel$participant, "mo", "mon")

write.csv(car_sel,'./data/clean/vision_scores.csv')
