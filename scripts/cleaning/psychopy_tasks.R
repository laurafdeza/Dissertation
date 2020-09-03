library(dplyr)
library(tidyr)
library(tidyverse)

setwd("./data/pitch/")

myfiles = list.files(pattern="*.csv", full.names=TRUE)
pitch_df <- plyr::ldply(myfiles, read_csv)
view(pitch_df)

setwd("..")
setwd("..")

lapply(pitch_df, class)

pitch_df$subject_id <- as.factor(pitch_df$subject_id)
pitch_df$trial_num <- as.factor(pitch_df$trial_num)
pitch_df$condition <- as.factor(pitch_df$condition)

pitch_df <- pitch_df %>%
  select(., -date) %>%
  filter(., correct == TRUE)

pitch_means <- pitch_df %>%
  group_by(subject_id, condition) %>%
  summarise(mean_rt = mean(rt)) %>%
#            sd_rt = sd(rt))
  spread(., condition, mean_rt)

write.csv(pitch_means,'./data/clean/pitch.csv')



setwd("./data/rhythm/")

myfiles = list.files(pattern="*.csv", full.names=TRUE)
rhythm_df <- plyr::ldply(myfiles, read_csv)
view(rhythm_df)

setwd("..")
setwd("..")

lapply(rhythm_df, class)

rhythm_df$subject_id <- as.factor(rhythm_df$subject_id)
rhythm_df$trial_num <- as.factor(rhythm_df$trial_num)
rhythm_df$condition <- as.factor(rhythm_df$condition)

rhythm_means <- rhythm_df %>%
  select(., -date) %>%
  group_by(subject_id, condition) %>%
  summarise(mean_dev = mean(deviation)) %>%
  spread(., condition, mean_dev)

write.csv(rhythm_means,'./data/clean/rhythm.csv')



setwd("./data/car/")

myfiles = list.files(pattern="*.csv", full.names=TRUE)
car_df <- plyr::ldply(myfiles, read_csv)
view(car_df)

setwd("..")
setwd("..")

lapply(car_df, class)

car_df$subject_id <- as.factor(car_df$subject_id)
car_df$trial <- as.factor(as.character(car_df$trial))
car_df$direction <- as.factor(car_df$direction)
car_df$speed <- as.factor(car_df$speed)

car_sec <- car_df %>%
  select(., -date) %>%
  unite(., condition, direction, speed, sep = "_", remove = TRUE) %>%
  group_by(subject_id, condition) %>%
  summarise(mean_dev_sc = mean(deviation_sec)) %>%
  spread(., condition, mean_dev_sc) %>%
  view()

write.csv(car_sec,'./data/clean/car_sec.csv')




setwd("./data/corsi/")

myfiles = list.files(pattern="*.csv", full.names=TRUE)
corsi_df <- plyr::ldply(myfiles, read_csv)
view(corsi_df)

setwd("..")
setwd("..")

lapply(corsi_df, class)

corsi_df$subject_id <- as.factor(corsi_df$subject_id)
corsi_df$level <- as.factor(as.character(corsi_df$level))


## mean time for all accurate responses
# corsi_time <- corsi_df %>%
#   filter(., correct == 1) %>%
#   select(., -date, -X1, -lang, -correct, -sequence, -response) %>%
#   group_by(subject_id) %>%
#   summarise(mean_time_corr = mean(time)) %>%
#   view()

# mean time for max level 3 completed
corsi_score <- corsi_df %>%
  filter(., correct == 1) %>%
  select(., -date, -X1, -lang, -correct, -sequence, -response) %>%
  group_by(subject_id, level) %>%
  summarise(n_corr = n(),
            mean_time_corr = mean(time)) %>%
  filter(., n_corr == 3)

corsi_score$level <- as.numeric(as.character(corsi_score$level))

corsi_score <- corsi_score %>%
  filter(., level == max(level)) %>%
  rename(., corsi_pt = level)

write.csv(corsi_score,'./data/clean/corsi.csv')
