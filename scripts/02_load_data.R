# Load data -------------------------------------------------------------------

source(here::here("scripts", "01_helpers.R"))

# Stress
#stress10 <- readRDS(here("data", "clean", "stress_10ms_clean.rds")) #if problems bc too big
stress10 <- read_csv(here::here("data", "clean", "stress_10ms_final.csv"))
stress50 <- read_csv(here::here("data", "clean", "stress_50ms_final.csv"))


# -----------------------------------------------------------------------------
