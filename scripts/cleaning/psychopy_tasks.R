library(dplyr)
library(tidyr)
library(tidyverse)

setwd("./data/pitch/")

myfiles = list.files(pattern="*.csv", full.names=TRUE)
pitch_df <- plyr::ldply(myfiles, read_csv)
view(pitch_df)
# data is wide

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

pitch_means <- rename(pitch_means, participant = subject_id)

pitch_means$participant <- str_replace(pitch_means$participant, "ae", "aes")
pitch_means$participant <- str_replace(pitch_means$participant, "ie", "ies")
pitch_means$participant <- str_replace(pitch_means$participant, "am", "ams")
pitch_means$participant <- str_replace(pitch_means$participant, "im", "ims")
pitch_means$participant <- str_replace(pitch_means$participant, "mo", "mon")

write.csv(pitch_means,'./data/clean/pitch_wide.csv')

pitch_long <- pitch_means %>%
  gather(., condition, rt, DO_down:SOL_up, factor_key = T) %>%
  separate(., condition, c("base_note", "direction"), sep = "_")

write.csv(pitch_long,'./data/clean/pitch_long.csv')




setwd("./data/rhythm/")

myfiles = list.files(pattern="*.csv", full.names=TRUE)
rhythm_df <- plyr::ldply(myfiles, read_csv)
view(rhythm_df)
# data is long

setwd("..")
setwd("..")

lapply(rhythm_df, class)

rhythm_df$subject_id <- as.factor(rhythm_df$subject_id)
rhythm_df$trial_num <- as.factor(as.character(rhythm_df$trial_num))
rhythm_df$condition <- as.factor(rhythm_df$condition)

rhythm_means <- rhythm_df %>%
  select(., -date) %>%
  
 # filter responses after 150 ms because they're considered 
 # reaction times rather than anticipation (Pagliarini, 2016)
  filter(., deviation < .150) %>%
  rename(., participant = subject_id) %>%
  group_by(participant, condition) %>%
  summarise(mean_dev = mean(deviation)) 

rhythm_means$participant <- str_replace(rhythm_means$participant, "ae", "aes")
rhythm_means$participant <- str_replace(rhythm_means$participant, "ie", "ies")
rhythm_means$participant <- str_replace(rhythm_means$participant, "am", "ams")
rhythm_means$participant <- str_replace(rhythm_means$participant, "im", "ims")
rhythm_means$participant <- str_replace(rhythm_means$participant, "mo", "mon")

write.csv(rhythm_means,'./data/clean/rhythm_long.csv')

rhythm_wide <- spread(rhythm_means, condition, mean_dev)

write.csv(rhythm_wide,'./data/clean/rhythm_wide.csv')


setwd("./data/car/")

myfiles = list.files(pattern="*.csv", full.names=TRUE)
car_df <- plyr::ldply(myfiles, read_csv)
view(car_df)
# long

setwd("..")
setwd("..")

lapply(car_df, class)

car_df$subject_id <- as.factor(car_df$subject_id)
car_df$trial <- as.factor(as.character(car_df$trial))
car_df$direction <- as.factor(car_df$direction)
car_df$speed <- as.factor(car_df$speed)

car_sec <- car_df %>%
  select(., -date) %>%
  group_by(subject_id, direction, speed) %>%
  summarise(mean_dev_sc = mean(deviation_sec)) %>%
  view()

write.csv(car_sec,'./data/clean/car_sec_long.csv')

car_wide <- car_sec %>%
  unite(., condition, direction, speed, sep = "_", remove = TRUE) %>%
  spread(., condition, mean_dev_sc) %>%
  view()

write.csv(car_wide,'./data/clean/car_sec_wide.csv')




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
