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

write.csv(auditory,'./data/clean/auditory_scores.csv', row.names = F)




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

write.csv(car_sel,'./data/clean/vision_scores.csv', row.names = F)


verbal_wm <- read_csv(here("data", 'pupurri_analysis.csv'))

agg <- separate(data = verbal_wm,
                col = group,
                into = c("prof", "l1"),
                sep = 1,
                remove = FALSE) 

ospan_stats_group <- agg %>%
  group_by(l1) %>%
  summarize(ospan_mean = mean(WM_set),
            ospan_sd = sd(WM_set))

print(ospan_stats_group)
# l1    ospan_mean  ospan_sd

# 1 es        8.89      2.11 EN 
# 2 ms        7.78      2.14 MA
# 3 on        6.2       2.72 ES


agg <- agg %>% mutate(ospan_mean = case_when(
  l1 == 'es' ~ 8.89,
  l1 == 'ms' ~ 7.78, 
  l1 == 'on' ~ 6.2),
  ospan_sd = case_when(
  l1 == 'es' ~ 2.11,
  l1 == 'ms' ~ 2.14, 
  l1 == 'on' ~ 2.72  
  )) %>%
  mutate(ospan = (WM_set - ospan_mean)/ospan_sd)

ver_wm <- agg %>% select(., participant, ospan) 

ver_wm$participant <- tolower(ver_wm$participant)

write.csv(ver_wm,'./data/clean/ospan_set_z_scores.csv', row.names = F)





vis_wm <- read_csv(here("data", 'clean', 'corsi.csv'))

vis_wm$subject_id <- str_replace(vis_wm$subject_id, "ae", "aes")
vis_wm$subject_id <- str_replace(vis_wm$subject_id, "ie", "ies")
vis_wm$subject_id <- str_replace(vis_wm$subject_id, "am", "ams")
vis_wm$subject_id <- str_replace(vis_wm$subject_id, "im", "ims")
vis_wm$subject_id <- str_replace(vis_wm$subject_id, "mo", "mon")

agg_vis <- vis_wm %>%
  select(., -X1) %>%
  rename(., participant = subject_id) %>%
  separate(., col = participant,
           into = c('group', 'id'),
           sep = 3,
           remove = FALSE) %>%
  separate(., col = group,
           into = c("prof", "l1"),
           sep = 1,
           remove = FALSE)

corsi_stats_group <- agg_vis %>%
  group_by(l1) %>%
  summarize(corsi_mean = mean(corsi_pt),
            corsi_sd = sd(corsi_pt))

print(corsi_stats_group)
# l1    ospan_mean  ospan_sd

# 1 es        4.77      1.16 EN 
# 2 ms        4.5       1.02 MA
# 3 on        4.68      1.04 ES


agg_vis <- agg_vis %>% mutate(corsi_mean = case_when(
  l1 == 'es' ~ 4.77,
  l1 == 'ms' ~ 4.5, 
  l1 == 'on' ~ 4.68),
  corsi_sd = case_when(
    l1 == 'es' ~ 1.16,
    l1 == 'ms' ~ 1.02, 
    l1 == 'on' ~ 1.04  
  )) %>%
  mutate(corsi = (corsi_pt - corsi_mean)/corsi_sd)

vis_wm <- agg_vis %>% select(., participant, corsi) 

write.csv(vis_wm,'./data/clean/corsi_z_scores.csv', row.names = F)
