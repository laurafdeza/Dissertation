# Morphosyntactic predictability: tidy stress data ----------------------------
#
# This script will load and tidy the raw eye tracking data
# with 50 ms bins and save the output to data/clean
#
# Last update: 06/11/2019 [working fine]
# Original script by Joseph Casillas
# Adapted to this project by Cristina Lozano-Argüelles
# -----------------------------------------------------------------------------


# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------
# load eye-tracking data
stress50 <- read_csv(here("data", "clean", "stress_50ms_final.csv"))
stress50 <- stress50 %>% 
  rename(., linx_stress = cond)


# Load pitch and rhythm data
pitch <- read_csv(here("data", 'clean', "pitch_long.csv"))
pitch <- pitch %>%
  select(., -X1) %>%
  rename(., pitch_rt = rt)


rhythm <- read_csv(here("data", 'clean', "rhythm_long.csv"))
rhythm <- rhythm %>%
  select(., -X1) %>%
  rename(., rhythm_time_dev = mean_dev,
         rhythm_cond = condition)
  

# Add pitch and rhythm score to eyetracking data frame
stress50_pi <- merge(x = stress50, y = pitch, by = "participant", all.x=TRUE)
stress50_mu <- merge(x = stress50_pi, y = rhythm, by = "participant", all.x=TRUE)

stress50_mu <- na.omit(stress50_mu)

write_csv(stress50_mu, here::here("data", 'clean', "music_50.csv"))
