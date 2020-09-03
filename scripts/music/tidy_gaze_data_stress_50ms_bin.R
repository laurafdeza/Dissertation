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
# Merge bilingual and monlingual datasets

data <- read_csv(here("data", "clean", "stress50clean.csv"))

# Tidy data -------------------------------------------------------------------

# Read data
stress50 <- data %>%

  # Remove unnecessary columns
  select(., -X1) %>%


  # Create eLog variable and respective wts
  mutate(.,eLog = log((target_count + 0.5) / (50 - target_count + 0.5)),
         wts = 1 / (target_count + 0.5) + 1 / (50 - target_count + 0.5)) %>%

  # Create 'group' column and new 'id'
  # in order to match participant ids with WM df
  # separate(., col = participant, into = c("group", "id"), sep = -2, remove = FALSE) %>%
  # 
  # 
  # # Filter out incorrect responses
  # filter(ACCURACY == 1) %>%

  # Select necessary columns
  # Gather data to prepare for bin adjustment
  # Get suffix onset label and center at 0 for each
  # participant for each item
  dplyr::select(participant, group, target, cond, target, bin_adj,
                target_count, target_prop, eLog, wts, onset_c3) %>%
  gather(., landmark, lm_bin, -c(participant:wts)) %>%
  mutate(., lm_bin = (lm_bin / 50) %>% ceiling(.),
         t_onset = if_else(bin_adj == lm_bin, TRUE, FALSE)) %>%
  # Must remove word3_c2 for coda words (its the start of the coda, unnecessary)
  # and remove word3_c3 for non coda words (they dont have c3, lm_bin = 0)
# filter(coda == 1 & landmark == "s_c2_sonset" | coda == 0 & lm_bin != 0) %>%
  

  group_by(., participant, target) %>%
  mutate(., time_zero = onset_pupil(bin_adj, t_onset, event = c("TRUE"))) %>%
  ungroup(.) %>%
  write_csv(here("data", "clean", "stress_50ms_final.csv"))

# -----------------------------------------------------------------------------
stress50$cond <- factor(stress50$cond, levels = c("1", "2"), 
                        labels = c("Present", "Past"))

# Test plot
stress50 %>%
  filter(time_zero > -10) %>%
  ggplot(., aes(x = time_zero, y = target_prop, color = group)) +
  facet_grid(. ~ cond) +
  geom_vline(xintercept = 4, lty = 3) +
  geom_hline(yintercept = 0.5, color = "white", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  ggtitle("Time course per verbal tense") +
  xlab("Time in 50 ms bins (0 = ???)") +
  ylab("Proportion of fixations on target") +
  scale_color_discrete(name="Group",
                     breaks = c("aes", 'ams', 'ies', 'ims', 'mon'),
                     labels = c("Adv EN", 'Adv MA', "Int EN", 'Int MA', 'ES Sp'))



