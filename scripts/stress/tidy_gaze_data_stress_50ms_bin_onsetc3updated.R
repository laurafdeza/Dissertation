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
# load data
stress_50 <- read.delim("./data/stress_50bin.txt")

# Check gaze fixation columns have different values
unique(stress_50$AVERAGE_IA_1_SAMPLE_COUNT)  # looking at target according to IA_#_ID
unique(stress_50$AVERAGE_IA_2_SAMPLE_COUNT)  # looking at distractor
unique(stress_50$AVERAGE_IA_0_SAMPLE_COUNT)  # elsewhere

# How much data we lose by selecting only accurate trials
sum( ( stress_50$ACCURACY == 0 ) / length( stress_50$ACCURACY ) )
# 0.004198754


stress_50$t04[stress_50$id == "S10_C2"] <- 853 
stress_50$t04[stress_50$id == "S15_C1"] <- 1144
stress_50$t04[stress_50$id == "S15_C2"] <- 939
stress_50$t04[stress_50$id == "S16_C1"] <- 879




# Tidy data -------------------------------------------------------------------

# Read data
stress_50 <- stress_50 %>%
  
  # create variable group
  separate(., col = RECORDING_SESSION_LABEL,
           into = c("group", "group_member"),
           sep = 3,
           remove = FALSE) %>%
  
  #select and rename variables of interest
  select(., RECORDING_SESSION_LABEL, TRIAL_INDEX, BIN_INDEX,
         AVERAGE_IA_0_SAMPLE_COUNT, AVERAGE_IA_0_SAMPLE_COUNT_.,
         AVERAGE_IA_1_SAMPLE_COUNT, AVERAGE_IA_1_SAMPLE_COUNT_.,
         AVERAGE_IA_2_SAMPLE_COUNT, AVERAGE_IA_2_SAMPLE_COUNT_.,
         ACCURACY, RT, block, cond, 
         id, lex_freq, phonot_freq, 
         t01, t02, t03, t04, t05, t06, t07, target, version, group) %>%
  dplyr::rename(., participant = RECORDING_SESSION_LABEL,
                trial = TRIAL_INDEX, 
                bin = BIN_INDEX,
                target_count = AVERAGE_IA_1_SAMPLE_COUNT, 
                target_prop = AVERAGE_IA_1_SAMPLE_COUNT_.,
                offset_prev_word = t01,
                onset_v1 = t02,
                onset_c2 = t03,
                onset_c3 = t04,
                onset_v2 = t05,
                offset_target = t06,
                endSentence = t07,
                sentence_id = id) %>%
  
  # remove incorrect
  filter(., ACCURACY == 1) %>%
  
  # drop unused levels of factors
  droplevels(.) %>%

  # Create eLog variable and respective wts
  mutate(.,eLog = log((target_count + 0.5) / (50 - target_count + 0.5)),
         wts = 1 / (target_count + 0.5) + 1 / (50 - target_count + 0.5)) %>%
    
  # Select necessary columns
  # Gather data to prepare for bin adjustment
  # Get suffix onset label and center at 0 for each
  # participant for each item
  dplyr::select(participant, group, target, cond, target, bin,
                target_count, target_prop, eLog, wts, onset_c3) %>%   
  # change onset_c3 in previous line depending on what trigger we want to observe
  gather(., landmark, lm_bin, -c(participant:wts)) %>%
  mutate(., lm_bin = (lm_bin / 50) %>% ceiling(.),
         t_onset = if_else(bin == lm_bin, TRUE, FALSE)) %>%
  
  group_by(., participant, target) %>%
  mutate(., time_zero = onset_pupil(bin, t_onset, event = c("TRUE"))) %>%
  ungroup(.)
  
# Load verbal WM
dem <- read_csv(here("data", "pupurri_analysis.csv"))
dem <- dem %>%
  select(., participant,DELE, percent_l2_week)
  
dem$participant <- tolower(dem$participant)
dem$DELE <- as.numeric(dem$DELE)

  
# Add verbal wm score to eyetracking data frame
stress_50 <- merge(x = stress_50, y = dem, by = "participant", all.x=TRUE)


# Create L1 column
stress_50 <- separate(stress_50,
                    col = group,
                    into = c("prof", "l1"),
                    sep = 1,
                    remove = FALSE)
  

stress_50$l1 <- str_replace(stress_50$l1, "es", "en")
stress_50$l1 <- str_replace(stress_50$l1, "ms", "ma")
stress_50$l1 <- str_replace(stress_50$l1, "on", "es")


stress_50$DELE[is.na(stress_50$DELE) & stress_50$l1 == 'es'] <- 56
stress_50$percent_l2_week[is.na(stress_50$percent_l2_week) & stress_50$l1 == 'es'] <- 0

stress50 <- stress_50 %>%
  # mutate(DELE = DELE + runif(n(), min = -0.15, max = 0.15) * (n() > 1)) %>%
  mutate(., #ospan = (WM_set - mean(WM_set))/sd(WM_set),
         use_z = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week),
    DELE_z = (DELE - mean(DELE))/sd(DELE)
  ) %>%
  
  # filter participants to make ES use and proficiency homogeneous
  
  filter(participant != 'ies04' & participant != 'ies17' & 
           participant != 'aes32' & participant != 'ies28')
 

# change name of .csv if trigger checked different
write_csv(stress50, here("data", "clean", "stress_50ms_final_onsetc3updated.csv"))

# -----------------------------------------------------------------------------

# # stress50 <- read_csv(here("data", "clean", "stress_50ms_final_onsetc3updated.csv"))

# stress50$cond <- factor(stress50$cond, levels = c("1", "2"), 
#                         labels = c("Paroxytone (present)", "Oxytone (preterite)"))

stress50 <- stress50 %>% mutate(l1 = fct_relevel(l1, 'es'))

# Test plot
timecourse <- stress50 %>%
  filter(time_zero > -9 & time_zero < 9) %>%
  ggplot(., aes(x = time_zero, y = target_prop, color = l1, fill = l1)) +
  geom_vline(xintercept = 4, lty = 3) +
  geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
                                stroke = 0.5, pch = 21, show.legend = FALSE) +
  scale_color_discrete(name="L1",
                     breaks = c('es', 'en', 'ma'),
                     labels = c('Spanish', "English", 'Chinese')) +
  scale_x_continuous(breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8),
                     labels = c("-400", "-300", "-200", "-100", "0", "100", "200", "300", "400")) +
  # ggtitle("Time course per verbal tense") +
  # xlab("Time in 50 ms bins (0 = marker time before accounting for 200 ms processing)") +
  labs(x = "Time relative to final syllable onset (ms)",
       y = "Proportion of fixations on target",
       caption = "Mean +/- 95% CI") +
  # annotate("text", x = 3.65, y = 0.53, label = '200ms',
  #                       angle = 90, size = 3, hjust = 0) +
  theme_grey(base_size = 10, base_family = "Times") +
  theme(legend.position = 'bottom')


ggsave("./figs/stress/gca/LL_changes/timecourse.png", timecourse, width = 180,
       height = 120, units = "mm", dpi = 600)

