# Morphosyntactic predictability: tidy stress data ----------------------------
#
# This script will load and tidy the raw eye tracking data
# with 10 ms bins and save the output to data/clean
#
# Last update: 05/16/2019
#
# -----------------------------------------------------------------------------




# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------

data <- read_csv(here("data", "clean", "stress10clean.csv"))

# Tidy data -------------------------------------------------------------------

# Read data
stress10 <- data %>%
  
  # Select necessary columns
 select(., -X1) %>%
  #        RECORDING_SESSION_LABEL, TRIAL_INDEX, BIN_INDEX,
  #        AVERAGE_IA_0_SAMPLE_COUNT, `AVERAGE_IA_0_SAMPLE_COUNT_%`,
  #        AVERAGE_IA_1_SAMPLE_COUNT, `AVERAGE_IA_1_SAMPLE_COUNT_%`,
  #        AVERAGE_IA_2_SAMPLE_COUNT, `AVERAGE_IA_2_SAMPLE_COUNT_%`,
  #        IA_1_ID, IA_2_ID, IA_0_ID, ACCURACY, RT, cond, 
  #        distractor, id, lex_freq, phonot_freq, sentence, 
  #        t01, t02, t03, t04, t05, t06, t07, target, version) %>%
  
  # Rename some columns
  # rename(., participant = RECORDING_SESSION_LABEL,
  #        trial = TRIAL_INDEX, 
  #        bin = BIN_INDEX,
  #        target_count = AVERAGE_IA_1_SAMPLE_COUNT, 
  #        target_prop = `AVERAGE_IA_1_SAMPLE_COUNT_%`,
  #        offset_prev_word = t01,
  #        onset_v1 = t02,
  #        onset_c2 = t03,
  #        onset_c3 = t04,
  #        onset_v2 = t05,
  #        offset_target = t06,
  #        endSentence = t07) %>%
  
  
  # Create eLog variable and respective wts
  mutate(.,eLog = log((target_count + 0.5) / (10 - target_count + 0.5)),
         wts = 1 / (target_count + 0.5) + 1 / (10 - target_count + 0.5)) %>%
  
  # Create 'group' column and new 'id'
  # in order to match participant ids with WM df
  # separate(., col = participant,
  #          into = c("group", "id"),
  #          sep = 3,
  #          remove = FALSE) %>%
  
  # Filter out incorrect responses
  # filter(ACCURACY == 1) %>%
  
  # Select necessary columns
  # Gather data to prepare for bin adjustment
  # Get suffix onset label and center at 0 for each
  # participant for each item
  dplyr::select(participant, group, target, cond, target, bin_adj,
                target_count, target_prop, eLog, wts, onset_c3) %>%
  gather(., landmark, lm_bin, -c(participant:wts)) %>%
  mutate(., lm_bin = (lm_bin / 10) %>% ceiling(.),
         t_onset = if_else(bin_adj == lm_bin, TRUE, FALSE)) %>%
  # Must remove word3_c2 for coda words (its the start of the coda, unnecessary)
  # and remove word3_c3 for non coda words (they dont have c3, lm_bin = 0)
  # filter(coda == 1 & landmark == "s_c2_sonset" | coda == 0 & lm_bin != 0) %>%
  
  
  group_by(., participant, target) %>%
  mutate(., time_zero = onset_pupil(bin_adj, t_onset, event = c("TRUE"))) %>%
  ungroup(.) %>%
  write_csv(here("data", "clean", "stress_10ms_final.csv"))


# -----------------------------------------------------------------------------

# Test plot for stress_unrelated
stress10$cond <- factor(stress10$cond, levels = c("1", "2"), 
                  labels = c("Present", "Past"))

stress10 %>%
  filter(time_zero > -50, time_zero < 100) %>% 
  ggplot(., aes(x = time_zero, y = target_prop, color = group)) +
  facet_grid(. ~ cond) +
  geom_vline(xintercept = 20, lty = 3) +
  geom_hline(yintercept = 0.5, color = "white", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  ggtitle("Time course per verbal tense") +
  xlab("Time in 10 ms bins (0 = ???)") +
  ylab("Proportion of fixations on target") +
  scale_color_discrete(name="Group",
                      breaks = c("aes", 'ams', 'ies', 'ims', 'mon'),
                      labels = c("Adv EN", 'Adv MA', "Int EN", 'Int MA', 'ES Sp'))


