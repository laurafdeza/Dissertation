source(here::here("scripts", "00_load_libs.R"))

# # Load wm data and combine with stress_subset_0 (proportion data)
# # in order to add working memory as a covariate
# wm_df <- read_csv(here("data", "pupurri_analysis.csv")) %>%
#   # filter(., #!(group %in% c("HS", "L")),
#   #        !(participant %in% c('', '', '', '', '',        # add participants to remove
#   #                             '', '', '', '', '',
#   #                             '', '', '', '', '',
#   #                             '', '', '', '', '',
#   #                             '', '', ''))) %>%
#   select(., participant, group, WM_set)
# 
# scale_this <- function(x) as.vector(scale(x))



# I don't know where df_short_temp comes from
# Original df did not have bin start and end times

df_short <- read.delim(here("data", "bin_times_stress_10.txt"))

df_short_temp <- df_short %>%

  # Filter out incorrect responses
  filter(ACCURACY == 1) %>%
  
  # Select necessary columns
  select(., BIN_END_TIME, BIN_START_TIME, startsentence,
         RECORDING_SESSION_LABEL, TRIAL_INDEX, BIN_INDEX,
         AVERAGE_IA_0_SAMPLE_COUNT, AVERAGE_IA_0_SAMPLE_COUNT_.,
         AVERAGE_IA_1_SAMPLE_COUNT, AVERAGE_IA_1_SAMPLE_COUNT_.,
         AVERAGE_IA_2_SAMPLE_COUNT, AVERAGE_IA_2_SAMPLE_COUNT_.,
         IA_1_ID, IA_2_ID, IA_0_ID, RT, cond,
         distractor, id, lex_freq, phonot_freq,
         t01, t02, t03, t04, t05, t06, t07, target, version) %>%
  
  # Rename some columns
  rename(., participant = RECORDING_SESSION_LABEL,
         trial = TRIAL_INDEX,
         bin = BIN_INDEX,
         targetCount = AVERAGE_IA_1_SAMPLE_COUNT,
         targetProp = AVERAGE_IA_1_SAMPLE_COUNT_.,
         offset_prev_word = t01,
         onset_v1 = t02,
         onset_c2 = t03,
         onset_c3 = t04,
         onset_v2 = t05,
         offset_target = t06,
         endSentence = t07) %>%

  # Create eLog variable and respective wts
  mutate(.,eLog = log((targetCount + 0.5) / (10 - targetCount + 0.5)),
        wts = 1 / (targetCount + 0.5) + 1 / (10 - targetCount + 0.5)) %>%

# Create 'group' column
  separate(., col = participant,
           into = c("group", "proficiency"),
           sep = 1,
           remove = FALSE) %>%
  
  # Select necessary columns
  # Gather data to prepare for bin adjustment
  # Get suffix onset label and center at 0 for each
  # participant for each item
  dplyr::select(participant, group, target, cond, target, bin,
                targetCount, targetProp, eLog, wts, startsentence,
                offset_prev_word, onset_v1, onset_c2, onset_c3,
                onset_v2, offset_target, endSentence,
                BIN_START_TIME, BIN_END_TIME) %>%
  gather(., landmark, lm_bin, -c(participant:wts)) %>%
  mutate(., lm_bin = (lm_bin / 10) %>% ceiling(.),
         t_onset = if_else(bin == lm_bin, TRUE, FALSE)) %>%
  # Must remove word3_c2 for coda words (its the start of the coda, unnecessary)
  # and remove word3_c3 for non coda words (they dont have c3, lm_bin = 0)
  # filter(coda == 1 & landmark == "s_c2_sonset" | coda == 0 & lm_bin != 0) %>%
  
  
  group_by(., participant, target) %>%
  mutate(., time_zero = onset_pupil(bin, t_onset, event = c("TRUE"))) %>%
  ungroup(.)





################################################


glimpse(df_stress)
glimpse(df_short_temp)

# calculate proportion of target fixations by sub, by target, by stress cond, by group
# for each landmark in the time course

df_timecourse <- df_short_temp %>%
  select(., participant, group, bin, BIN_START_TIME, BIN_END_TIME, target,
         startsentence:target, cond, targetCount:targetProp, eLog) %>%
  mutate(., bin = bin * 10)

range(df_timecourse$bin)

df_timecourse_startsentence <- df_timecourse %>%
  filter(., startsentence + 1550 >= BIN_START_TIME, startsentence + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, cond:eLog) %>%
  mutate(., landmark = 'start_sentence')

df_timecourse_word2_c1v1 <- df_timecourse %>%
  filter(., word2_c1v1 + 1550 >= BIN_START_TIME, word2_c1v1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word2_c1v1')


df_timecourse_word3_20msafterv1 <- df_timecourse %>%
  filter(., word3_20msafterv1 + 1550 >= BIN_START_TIME, word3_20msafterv1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_20msafterv1')


df_timecourse_word3_c1v1 <- df_timecourse %>%
  filter(., word3_c1v1 + 1550 >= BIN_START_TIME, word3_c1v1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_c1v1')


df_timecourse_word3_c2 <- df_timecourse %>%
  filter(., word3_c2 + 1550 >= BIN_START_TIME, word3_c2 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_c2')

df_timecourse_word3_c3 <- df_timecourse %>%
  filter(., word3_c3 + 1550 >= BIN_START_TIME, word3_c3 + 1550 <= BIN_END_TIME)  %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_c3')

df_timecourse_word3_suffix <- df_timecourse %>%
  filter(., word3_suffix + 1550 >= BIN_START_TIME, word3_suffix + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_suffix')

df_timecourse_word4_c1v1 <- df_timecourse %>%
  filter(., word4_c1v1 + 1550 >= BIN_START_TIME, word4_c1v1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word4_c1v1')

df_timecourse_word5 <- df_timecourse %>%
  filter(., word5 + 1550 >= BIN_START_TIME, word5 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word5')

df_timecourse_end_sentence <- df_timecourse %>%
  filter(., end_sentence + 1550 >= BIN_START_TIME, end_sentence + 1550 <= BIN_END_TIME)  %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'end_sentence')

df_landmarks <- do.call("rbind", list(df_timecourse_startsentence,
                                      df_timecourse_word2_c1v1,
                                      df_timecourse_word3_20msafterv1,
                                      df_timecourse_word3_c1v1,
                                      df_timecourse_word3_c2,
                                      df_timecourse_word3_c3,
                                      df_timecourse_word3_suffix,
                                      df_timecourse_word4_c1v1,
                                      df_timecourse_word5,
                                      df_timecourse_end_sentence))

glimpse(df_landmarks)

df_landmarks <- mutate(df_landmarks,
                       landmark = factor(landmark, levels = c('start_sentence',
                                                              'word2_c1v1',
                                                              'word3_c1v1',
                                                              'word3_20msafterv1',
                                                              'word3_c2',
                                                              'word3_c3',
                                                              'word3_suffix',
                                                              'word4_c1v1',
                                                              'word5',
                                                              'end_sentence')))




df_landmarks %>%
  na.omit(.) %>%
  filter(., coda == 0,
         landmark %in% c('word3_c1v1', 'word3_20msafterv1', 'word3_suffix', 'word4_c1v1')) %>%
  group_by(., participant, target, group, landmark) %>%
  summarize(., target_fix = mean(targetProp)) %>%
  ggplot(., aes(x = landmark, y = target_fix, shape = group, dodge = group)) +
  geom_hline(yintercept = 0.5, color = 'black') +
  ylim(0, 1) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', position = position_dodge(width = 0.3)) +
  theme_bw()

df_landmarks %>%
  na.omit(.) %>%
  filter(., coda == 1,
         landmark %in% c('word3_c1v1', 'word3_20msafterv1', 'word3_c2', 'word3_suffix', 'word4_c1v1')) %>%
  group_by(., participant, target, group, landmark) %>%
  summarize(., target_fix = mean(targetProp)) %>%
  ggplot(., aes(x = landmark, y = target_fix, shape = group, dodge = group)) +
  geom_hline(yintercept = 0.5, color = 'black') +
  ylim(0, 1) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', position = position_dodge(width = 0.3)) +
  theme_bw()

