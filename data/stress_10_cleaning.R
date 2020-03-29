# load packages
library(plyr)
library(tidyr)
library(tidyverse)

# load data
stress_10 <- read_tsv("data/stress_10.txt")

# Check gaze fixation columns have different values
unique(stress_10$AVERAGE_IA_1_SAMPLE_COUNT)  # looking at target according to IA_#_ID
unique(stress_10$AVERAGE_IA_2_SAMPLE_COUNT)  # looking at distractor
unique(stress_10$AVERAGE_IA_0_SAMPLE_COUNT)  # elsewhere

# create variable group
ss_10 <- stress_10 %>%
  separate(., col = RECORDING_SESSION_LABEL,
           into = c("group", "group_member"),
           sep = 3,
           remove = FALSE) %>%
  
  #select and rename variables of interest
  select(., RECORDING_SESSION_LABEL, TRIAL_INDEX, BIN_INDEX,
         AVERAGE_IA_0_SAMPLE_COUNT, `AVERAGE_IA_0_SAMPLE_COUNT_%`,
         AVERAGE_IA_1_SAMPLE_COUNT, `AVERAGE_IA_1_SAMPLE_COUNT_%`,
         AVERAGE_IA_2_SAMPLE_COUNT, `AVERAGE_IA_2_SAMPLE_COUNT_%`,
         IA_1_ID, IA_2_ID, IA_0_ID, ACCURACY, RT, cond, 
         distractor, id, lex_freq, phonot_freq, sentence, 
         t01, t02, t03, t04, t05, t06, t07, target, version, group) %>%
  dplyr::rename(., participant = RECORDING_SESSION_LABEL,
         trial = TRIAL_INDEX, 
         bin = BIN_INDEX,
         target_count = AVERAGE_IA_1_SAMPLE_COUNT, 
         target_prop = `AVERAGE_IA_1_SAMPLE_COUNT_%`,
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
  droplevels(.)

# Adjust processing time in bin columns
ss_10$bin_adj <- ss_10$bin - 20


# find data row for bin at C2 onset (verb structure = C1V1C2.C3V2)


################################################################################
s10 <- ss_10

## create function to find time marker
# to debug options(error=recover)
on_c2 <- ldply(s10, function(s10) {
  for (i in s10$sentence_id) {
    bin_mark_i <- round(unique(s10[s10[sentence_id == i], "onset_c3"]) / 10)
    
    bin_sum_i <- filter(s10, s10[sentence_id == i], bin_adj == as.numeric(bin_mark_i))
  }
}
)



function(x, y) {
  
  bin_mark <- round(unique(df[df$x, y]) / 10)
  
  bin_sum_x <- filter(df, x, bin_adj == bin_num)

  rbind(bin_sum_x)
}


################################################################################


bin_S01_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S01_C1", "onset_c2"])) / 10) 
bin_S01_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S01_C2", "onset_c2"])) / 10)

st_S01c2 <- filter(ss_10, sentence_id == "S01_C1", bin_adj == as.numeric(bin_S01_stc2))
un_S01c2 <- filter(ss_10, sentence_id == "S01_C2", bin_adj == as.numeric(bin_S01_unc2)) 



bin_S02_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S02_C1", "onset_c2"])) / 10) 
bin_S02_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S02_C2", "onset_c2"])) / 10)

st_S02c2 <- filter(ss_10, sentence_id == "S02_C1", bin_adj == as.numeric(bin_S02_stc2))
un_S02c2 <- filter(ss_10, sentence_id == "S02_C2", bin_adj == as.numeric(bin_S02_unc2))



bin_S03_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S03_C1", "onset_c2"])) / 10) 
bin_S03_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S03_C2", "onset_c2"])) / 10)

st_S03c2 <- filter(ss_10, sentence_id == "S03_C1", bin_adj == as.numeric(bin_S03_stc2))
un_S03c2 <- filter(ss_10, sentence_id == "S03_C2", bin_adj == as.numeric(bin_S03_unc2))



bin_S04_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S04_C1", "onset_c2"])) / 10) 
bin_S04_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S04_C2", "onset_c2"])) / 10)

st_S04c2 <- filter(ss_10, sentence_id == "S04_C1", bin_adj == as.numeric(bin_S04_stc2))
un_S04c2 <- filter(ss_10, sentence_id == "S04_C2", bin_adj == as.numeric(bin_S04_unc2))



bin_S05_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S05_C1", "onset_c2"])) / 10) 
bin_S05_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S05_C2", "onset_c2"])) / 10)

st_S05c2 <- filter(ss_10, sentence_id == "S05_C1", bin_adj == as.numeric(bin_S05_stc2))
un_S05c2 <- filter(ss_10, sentence_id == "S05_C2", bin_adj == as.numeric(bin_S05_unc2))



bin_S06_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S06_C1", "onset_c2"])) / 10) 
bin_S06_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S06_C2", "onset_c2"])) / 10)

st_S06c2 <- filter(ss_10, sentence_id == "S06_C1", bin_adj == as.numeric(bin_S06_stc2))
un_S06c2 <- filter(ss_10, sentence_id == "S06_C2", bin_adj == as.numeric(bin_S06_unc2))



bin_S07_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S07_C1", "onset_c2"])) / 10) 
bin_S07_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S07_C2", "onset_c2"])) / 10)

st_S07c2 <- filter(ss_10, sentence_id == "S07_C1", bin_adj == as.numeric(bin_S07_stc2))
un_S07c2 <- filter(ss_10, sentence_id == "S07_C2", bin_adj == as.numeric(bin_S07_unc2))



bin_S08_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S08_C1", "onset_c2"])) / 10) 
bin_S08_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S08_C2", "onset_c2"])) / 10)

st_S08c2 <- filter(ss_10, sentence_id == "S08_C1", bin_adj == as.numeric(bin_S08_stc2))
un_S08c2 <- filter(ss_10, sentence_id == "S08_C2", bin_adj == as.numeric(bin_S08_unc2))



bin_S09_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S09_C1", "onset_c2"])) / 10) 
bin_S09_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S09_C2", "onset_c2"])) / 10)

st_S09c2 <- filter(ss_10, sentence_id == "S09_C1", bin_adj == as.numeric(bin_S09_stc2))
un_S09c2 <- filter(ss_10, sentence_id == "S09_C2", bin_adj == as.numeric(bin_S09_unc2))



bin_S10_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S10_C1", "onset_c2"])) / 10) 
bin_S10_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S10_C2", "onset_c2"])) / 10)

st_S10c2 <- filter(ss_10, sentence_id == "S10_C1", bin_adj == as.numeric(bin_S10_stc2))
un_S10c2 <- filter(ss_10, sentence_id == "S10_C2", bin_adj == as.numeric(bin_S10_unc2))



bin_S11_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S11_C1", "onset_c2"])) / 10) 
bin_S11_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S11_C2", "onset_c2"])) / 10)

st_S11c2 <- filter(ss_10, sentence_id == "S11_C1", bin_adj == as.numeric(bin_S11_stc2))
un_S11c2 <- filter(ss_10, sentence_id == "S11_C2", bin_adj == as.numeric(bin_S11_unc2))



bin_S12_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S12_C1", "onset_c2"])) / 10) 
bin_S12_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S12_C2", "onset_c2"])) / 10)

st_S12c2 <- filter(ss_10, sentence_id == "S12_C1", bin_adj == as.numeric(bin_S12_stc2))
un_S12c2 <- filter(ss_10, sentence_id == "S12_C2", bin_adj == as.numeric(bin_S12_unc2))



bin_S13_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S13_C1", "onset_c2"])) / 10) 
bin_S13_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S13_C2", "onset_c2"])) / 10)

st_S13c2 <- filter(ss_10, sentence_id == "S13_C1", bin_adj == as.numeric(bin_S13_stc2))
un_S13c2 <- filter(ss_10, sentence_id == "S13_C2", bin_adj == as.numeric(bin_S13_unc2))



bin_S14_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S14_C1", "onset_c2"])) / 10) 
bin_S14_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S14_C2", "onset_c2"])) / 10)

st_S14c2 <- filter(ss_10, sentence_id == "S14_C1", bin_adj == as.numeric(bin_S14_stc2))
un_S14c2 <- filter(ss_10, sentence_id == "S14_C2", bin_adj == as.numeric(bin_S14_unc2))



bin_S15_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S15_C1", "onset_c2"])) / 10) 
bin_S15_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S15_C2", "onset_c2"])) / 10)

st_S15c2 <- filter(ss_10, sentence_id == "S15_C1", bin_adj == as.numeric(bin_S15_stc2))
un_S15c2 <- filter(ss_10, sentence_id == "S15_C2", bin_adj == as.numeric(bin_S15_unc2))



bin_S16_stc2 <- round((unique(ss_10[ss_10$sentence_id == "S16_C1", "onset_c2"])) / 10) 
bin_S16_unc2 <- round((unique(ss_10[ss_10$sentence_id == "S16_C2", "onset_c2"])) / 10)

st_S16c2 <- filter(ss_10, sentence_id == "S16_C1", bin_adj == as.numeric(bin_S16_stc2))
un_S16c2 <- filter(ss_10, sentence_id == "S16_C2", bin_adj == as.numeric(bin_S16_unc2))

# bind individual rows into one single dataframe
st10_onC2 <- rbind(st_S01c2, un_S01c2, st_S02c2, un_S02c2, st_S03c2, un_S03c2, st_S04c2, un_S04c2,
                   st_S05c2, un_S05c2, st_S06c2, un_S06c2, st_S07c2, un_S07c2, st_S08c2, un_S08c2,
                   st_S09c2, un_S09c2, st_S10c2, un_S10c2, st_S12c2, un_S12c2, st_S12c2, un_S12c2,
                   st_S13c2, un_S13c2, st_S14c2, un_S14c2, st_S15c2, un_S15c2, st_S16c2, un_S16c2)

# Bonferroni correction for alpha level
# Divided by 5 because using same data for multiple tests (5 different groups)
print(alphaAdj <- 0.05/10)
# 0.005

# One Sample t-test
# proportion of fixations on target by monolinguals when tense = present
mon_pres <- st10_onC2 %>% 
  filter(., group == "mon", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pres$mean_prop, alternative = "greater", mu = 0.5)
# t = -1.7905, df = 29,
# p-value = 0.9581

# proportion of fixations on target by monolinguals when tense = preterit
mon_pret <- st10_onC2 %>% 
  filter(., group == "mon", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pret$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.7348, df = 29, p-value = 0.2342


# proportion of fixations on target by intermediate EN when tense = present
ies_pres <- st10_onC2 %>% 
  filter(., group == "ies", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pres$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.3635, df = 13, p-value = 0.9828


# proportion of fixations on target by intermediate EN when tense = preterit
ies_pret <- st10_onC2 %>% 
  filter(., group == "ies", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pret$mean_prop, alternative = "greater", mu = 0.5)
# t = -0.67472, df = 13, p-value = 0.7442


# proportion of fixations on target by advanced EN when tense = present
aes_pres <- st10_onC2 %>% 
  filter(., group == "aes", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pres$mean_prop, alternative = "greater", mu = 0.5)
# t = -3.8389, df = 20, p-value = 0.9995


# proportion of fixations on target by advanced EN when tense = preterit
aes_pret <- st10_onC2 %>% 
  filter(., group == "aes", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pret$mean_prop, alternative = "greater", mu = 0.5)
# t = -1.5333, df = 20, p-value = 0.9296



# proportion of fixations on target by intermediate CH when tense = present
ims_pres <- st10_onC2 %>% 
  filter(., group == "ims", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pres$mean_prop, alternative = "greater", mu = 0.5)
# t = -3.5678, df = 13, p-value = 0.9983


# proportion of fixations on target by intermediate CH when tense = preterit
ims_pret <- st10_onC2 %>% 
  filter(., group == "ims", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pret$mean_prop, alternative = "greater", mu = 0.5)
# t = -1.5384, df = 13, p-value = 0.926


# proportion of fixations on target by advanced CH when tense = present
ams_pres <- st10_onC2 %>% 
  filter(., group == "ams", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pres$mean_prop, alternative = "greater", mu = 0.5)
# t = -3.1269, df = 18, p-value = 0.9971


# proportion of fixations on target by advanced CH when tense = preterit
ams_pret <- st10_onC2 %>% 
  filter(., group == "ams", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pret$mean_prop, alternative = "greater", mu = 0.5)
# t = -1.8669, df = 18, p-value = 0.9609




########################################

# Prediction at onset C3?

bin_S01_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S01_C1", "onset_c3"])) / 10) 
bin_S01_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S01_C2", "onset_c3"])) / 10)

st_S01c3 <- filter(ss_10, sentence_id == "S01_C1", bin_adj == as.numeric(bin_S01_stc3))
un_S01c3 <- filter(ss_10, sentence_id == "S01_C2", bin_adj == as.numeric(bin_S01_unc3)) 



bin_S02_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S02_C1", "onset_c3"])) / 10) 
bin_S02_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S02_C2", "onset_c3"])) / 10)

st_S02c3 <- filter(ss_10, sentence_id == "S02_C1", bin_adj == as.numeric(bin_S02_stc3))
un_S02c3 <- filter(ss_10, sentence_id == "S02_C2", bin_adj == as.numeric(bin_S02_unc3))



bin_S03_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S03_C1", "onset_c3"])) / 10) 
bin_S03_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S03_C2", "onset_c3"])) / 10)

st_S03c3 <- filter(ss_10, sentence_id == "S03_C1", bin_adj == as.numeric(bin_S03_stc3))
un_S03c3 <- filter(ss_10, sentence_id == "S03_C2", bin_adj == as.numeric(bin_S03_unc3))



bin_S04_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S04_C1", "onset_c3"])) / 10) 
bin_S04_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S04_C2", "onset_c3"])) / 10)

st_S04c3 <- filter(ss_10, sentence_id == "S04_C1", bin_adj == as.numeric(bin_S04_stc3))
un_S04c3 <- filter(ss_10, sentence_id == "S04_C2", bin_adj == as.numeric(bin_S04_unc3))



bin_S05_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S05_C1", "onset_c3"])) / 10) 
bin_S05_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S05_C2", "onset_c3"])) / 10)

st_S05c3 <- filter(ss_10, sentence_id == "S05_C1", bin_adj == as.numeric(bin_S05_stc3))
un_S05c3 <- filter(ss_10, sentence_id == "S05_C2", bin_adj == as.numeric(bin_S05_unc3))



bin_S06_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S06_C1", "onset_c3"])) / 10) 
bin_S06_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S06_C2", "onset_c3"])) / 10)

st_S06c3 <- filter(ss_10, sentence_id == "S06_C1", bin_adj == as.numeric(bin_S06_stc3))
un_S06c3 <- filter(ss_10, sentence_id == "S06_C2", bin_adj == as.numeric(bin_S06_unc3))



bin_S07_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S07_C1", "onset_c3"])) / 10) 
bin_S07_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S07_C2", "onset_c3"])) / 10)

st_S07c3 <- filter(ss_10, sentence_id == "S07_C1", bin_adj == as.numeric(bin_S07_stc3))
un_S07c3 <- filter(ss_10, sentence_id == "S07_C2", bin_adj == as.numeric(bin_S07_unc3))



bin_S08_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S08_C1", "onset_c3"])) / 10) 
bin_S08_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S08_C2", "onset_c3"])) / 10)

st_S08c3 <- filter(ss_10, sentence_id == "S08_C1", bin_adj == as.numeric(bin_S08_stc3))
un_S08c3 <- filter(ss_10, sentence_id == "S08_C2", bin_adj == as.numeric(bin_S08_unc3))



bin_S09_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S09_C1", "onset_c3"])) / 10) 
bin_S09_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S09_C2", "onset_c3"])) / 10)

st_S09c3 <- filter(ss_10, sentence_id == "S09_C1", bin_adj == as.numeric(bin_S09_stc3))
un_S09c3 <- filter(ss_10, sentence_id == "S09_C2", bin_adj == as.numeric(bin_S09_unc3))



bin_S10_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S10_C1", "onset_c3"])) / 10) 
bin_S10_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S10_C2", "onset_c3"])) / 10)

st_S10c3 <- filter(ss_10, sentence_id == "S10_C1", bin_adj == as.numeric(bin_S10_stc3))
un_S10c3 <- filter(ss_10, sentence_id == "S10_C2", bin_adj == as.numeric(bin_S10_unc3))



bin_S11_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S11_C1", "onset_c3"])) / 10) 
bin_S11_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S11_C2", "onset_c3"])) / 10)

st_S11c3 <- filter(ss_10, sentence_id == "S11_C1", bin_adj == as.numeric(bin_S11_stc3))
un_S11c3 <- filter(ss_10, sentence_id == "S11_C2", bin_adj == as.numeric(bin_S11_unc3))



bin_S12_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S12_C1", "onset_c3"])) / 10) 
bin_S12_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S12_C2", "onset_c3"])) / 10)

st_S12c3 <- filter(ss_10, sentence_id == "S12_C1", bin_adj == as.numeric(bin_S12_stc3))
un_S12c3 <- filter(ss_10, sentence_id == "S12_C2", bin_adj == as.numeric(bin_S12_unc3))



bin_S13_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S13_C1", "onset_c3"])) / 10) 
bin_S13_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S13_C2", "onset_c3"])) / 10)

st_S13c3 <- filter(ss_10, sentence_id == "S13_C1", bin_adj == as.numeric(bin_S13_stc3))
un_S13c3 <- filter(ss_10, sentence_id == "S13_C2", bin_adj == as.numeric(bin_S13_unc3))



bin_S14_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S14_C1", "onset_c3"])) / 10) 
bin_S14_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S14_C2", "onset_c3"])) / 10)

st_S14c3 <- filter(ss_10, sentence_id == "S14_C1", bin_adj == as.numeric(bin_S14_stc3))
un_S14c3 <- filter(ss_10, sentence_id == "S14_C2", bin_adj == as.numeric(bin_S14_unc3))



bin_S15_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S15_C1", "onset_c3"])) / 10) 
bin_S15_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S15_C2", "onset_c3"])) / 10)

st_S15c3 <- filter(ss_10, sentence_id == "S15_C1", bin_adj == as.numeric(bin_S15_stc3))
un_S15c3 <- filter(ss_10, sentence_id == "S15_C2", bin_adj == as.numeric(bin_S15_unc3))



bin_S16_stc3 <- round((unique(ss_10[ss_10$sentence_id == "S16_C1", "onset_c3"])) / 10) 
bin_S16_unc3 <- round((unique(ss_10[ss_10$sentence_id == "S16_C2", "onset_c3"])) / 10)

st_S16c3 <- filter(ss_10, sentence_id == "S16_C1", bin_adj == as.numeric(bin_S16_stc3))
un_S16c3 <- filter(ss_10, sentence_id == "S16_C2", bin_adj == as.numeric(bin_S16_unc3))

# bind individual rows into one single dataframe
st10_onC3 <- rbind(st_S01c3, un_S01c3, st_S02c3, un_S02c3, st_S03c3, un_S03c3, st_S04c3, un_S04c3,
                   st_S05c3, un_S05c3, st_S06c3, un_S06c3, st_S07c3, un_S07c3, st_S08c3, un_S08c3,
                   st_S09c3, un_S09c3, st_S10c3, un_S10c3, st_S12c3, un_S12c3, st_S12c3, un_S12c3,
                   st_S13c3, un_S13c3, st_S14c3, un_S14c3, st_S15c3, un_S15c3, st_S16c3, un_S16c3)


write.csv(st10_onC3, "data/clean/fixations_onC3_10.csv")

# Bonferroni correction for alpha level
# Divided by 5 because using same data for multiple tests (5 different groups)
print(alphaAdj <- 0.05/10)
# 0.005

# One Sample t-test
# proportion of fixations on target by monolinguals when tense = present
mon_pres_c3 <- st10_onC3 %>% 
  filter(., group == "mon", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 4.3768, df = 29, p-value = 7.133e-05

# proportion of fixations on target by monolinguals when tense = preterit
mon_pret_c3 <- st10_onC3 %>% 
  filter(., group == "mon", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 4.9592, df = 29, p-value = 1.42e-05


# proportion of fixations on target by intermediate EN when tense = present
ies_pres_c3 <- st10_onC3 %>% 
  filter(., group == "ies", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.8065, df = 13, p-value = 0.9926


# proportion of fixations on target by intermediate EN when tense = preterit
ies_pret_c3 <- st10_onC3 %>% 
  filter(., group == "ies", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 1.6859, df = 13, p-value = 0.05783


# proportion of fixations on target by advanced EN when tense = present
aes_pres_c3 <- st10_onC3 %>% 
  filter(., group == "aes", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.8624, df = 20, p-value = 0.9952


# proportion of fixations on target by advanced EN when tense = preterit
aes_pret_c3 <- st10_onC3 %>% 
  filter(., group == "aes", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.27751, df = 20, p-value = 0.3921



# proportion of fixations on target by intermediate CH when tense = present
ims_pres_c3 <- st10_onC3 %>% 
  filter(., group == "ims", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.4941, df = 13, p-value = 0.9866


# proportion of fixations on target by intermediate CH when tense = preterit
ims_pret_c3 <- st10_onC3 %>% 
  filter(., group == "ims", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.42373, df = 13, p-value = 0.3393


# proportion of fixations on target by advanced CH when tense = present
ams_pres_c3 <- st10_onC3 %>% 
  filter(., group == "ams", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = -1.5285, df = 18, p-value = 0.9281


# proportion of fixations on target by advanced CH when tense = preterit
ams_pret_c3 <- st10_onC3 %>% 
  filter(., group == "ams", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.32532, df = 18, p-value = 0.3743



#####
# Gaze fixations at segmental disambiguator onset
bin_S01_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S01_C1", "onset_v2"])) / 10) 
bin_S01_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S01_C2", "onset_v2"])) / 10)

st_S01v2 <- filter(ss_10, sentence_id == "S01_C1", bin_adj == as.numeric(bin_S01_stv2))
un_S01v2 <- filter(ss_10, sentence_id == "S01_C2", bin_adj == as.numeric(bin_S01_unv2)) 



bin_S02_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S02_C1", "onset_v2"])) / 10) 
bin_S02_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S02_C2", "onset_v2"])) / 10)

st_S02v2 <- filter(ss_10, sentence_id == "S02_C1", bin_adj == as.numeric(bin_S02_stv2))
un_S02v2 <- filter(ss_10, sentence_id == "S02_C2", bin_adj == as.numeric(bin_S02_unv2))



bin_S03_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S03_C1", "onset_v2"])) / 10) 
bin_S03_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S03_C2", "onset_v2"])) / 10)

st_S03v2 <- filter(ss_10, sentence_id == "S03_C1", bin_adj == as.numeric(bin_S03_stv2))
un_S03v2 <- filter(ss_10, sentence_id == "S03_C2", bin_adj == as.numeric(bin_S03_unv2))



bin_S04_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S04_C1", "onset_v2"])) / 10) 
bin_S04_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S04_C2", "onset_v2"])) / 10)

st_S04v2 <- filter(ss_10, sentence_id == "S04_C1", bin_adj == as.numeric(bin_S04_stv2))
un_S04v2 <- filter(ss_10, sentence_id == "S04_C2", bin_adj == as.numeric(bin_S04_unv2))



bin_S05_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S05_C1", "onset_v2"])) / 10) 
bin_S05_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S05_C2", "onset_v2"])) / 10)

st_S05v2 <- filter(ss_10, sentence_id == "S05_C1", bin_adj == as.numeric(bin_S05_stv2))
un_S05v2 <- filter(ss_10, sentence_id == "S05_C2", bin_adj == as.numeric(bin_S05_unv2))



bin_S06_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S06_C1", "onset_v2"])) / 10) 
bin_S06_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S06_C2", "onset_v2"])) / 10)

st_S06v2 <- filter(ss_10, sentence_id == "S06_C1", bin_adj == as.numeric(bin_S06_stv2))
un_S06v2 <- filter(ss_10, sentence_id == "S06_C2", bin_adj == as.numeric(bin_S06_unv2))



bin_S07_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S07_C1", "onset_v2"])) / 10) 
bin_S07_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S07_C2", "onset_v2"])) / 10)

st_S07v2 <- filter(ss_10, sentence_id == "S07_C1", bin_adj == as.numeric(bin_S07_stv2))
un_S07v2 <- filter(ss_10, sentence_id == "S07_C2", bin_adj == as.numeric(bin_S07_unv2))



bin_S08_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S08_C1", "onset_v2"])) / 10) 
bin_S08_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S08_C2", "onset_v2"])) / 10)

st_S08v2 <- filter(ss_10, sentence_id == "S08_C1", bin_adj == as.numeric(bin_S08_stv2))
un_S08v2 <- filter(ss_10, sentence_id == "S08_C2", bin_adj == as.numeric(bin_S08_unv2))



bin_S09_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S09_C1", "onset_v2"])) / 10) 
bin_S09_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S09_C2", "onset_v2"])) / 10)

st_S09v2 <- filter(ss_10, sentence_id == "S09_C1", bin_adj == as.numeric(bin_S09_stv2))
un_S09v2 <- filter(ss_10, sentence_id == "S09_C2", bin_adj == as.numeric(bin_S09_unv2))



bin_S10_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S10_C1", "onset_v2"])) / 10) 
bin_S10_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S10_C2", "onset_v2"])) / 10)

st_S10v2 <- filter(ss_10, sentence_id == "S10_C1", bin_adj == as.numeric(bin_S10_stv2))
un_S10v2 <- filter(ss_10, sentence_id == "S10_C2", bin_adj == as.numeric(bin_S10_unv2))



bin_S11_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S11_C1", "onset_v2"])) / 10) 
bin_S11_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S11_C2", "onset_v2"])) / 10)

st_S11v2 <- filter(ss_10, sentence_id == "S11_C1", bin_adj == as.numeric(bin_S11_stv2))
un_S11v2 <- filter(ss_10, sentence_id == "S11_C2", bin_adj == as.numeric(bin_S11_unv2))



bin_S12_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S12_C1", "onset_v2"])) / 10) 
bin_S12_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S12_C2", "onset_v2"])) / 10)

st_S12v2 <- filter(ss_10, sentence_id == "S12_C1", bin_adj == as.numeric(bin_S12_stv2))
un_S12v2 <- filter(ss_10, sentence_id == "S12_C2", bin_adj == as.numeric(bin_S12_unv2))



bin_S13_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S13_C1", "onset_v2"])) / 10) 
bin_S13_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S13_C2", "onset_v2"])) / 10)

st_S13v2 <- filter(ss_10, sentence_id == "S13_C1", bin_adj == as.numeric(bin_S13_stv2))
un_S13v2 <- filter(ss_10, sentence_id == "S13_C2", bin_adj == as.numeric(bin_S13_unv2))



bin_S14_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S14_C1", "onset_v2"])) / 10) 
bin_S14_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S14_C2", "onset_v2"])) / 10)

st_S14v2 <- filter(ss_10, sentence_id == "S14_C1", bin_adj == as.numeric(bin_S14_stv2))
un_S14v2 <- filter(ss_10, sentence_id == "S14_C2", bin_adj == as.numeric(bin_S14_unv2))



bin_S15_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S15_C1", "onset_v2"])) / 10) 
bin_S15_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S15_C2", "onset_v2"])) / 10)

st_S15v2 <- filter(ss_10, sentence_id == "S15_C1", bin_adj == as.numeric(bin_S15_stv2))
un_S15v2 <- filter(ss_10, sentence_id == "S15_C2", bin_adj == as.numeric(bin_S15_unv2))



bin_S16_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S16_C1", "onset_v2"])) / 10) 
bin_S16_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S16_C2", "onset_v2"])) / 10)

st_S16v2 <- filter(ss_10, sentence_id == "S16_C1", bin_adj == as.numeric(bin_S16_stv2))
un_S16v2 <- filter(ss_10, sentence_id == "S16_C2", bin_adj == as.numeric(bin_S16_unv2))

# bind individual rows into one single dataframe
st10_onv2 <- rbind(st_S01v2, un_S01v2, st_S02v2, un_S02v2, st_S03v2, un_S03v2, st_S04v2, un_S04v2,
                   st_S05v2, un_S05v2, st_S06v2, un_S06v2, st_S07v2, un_S07v2, st_S08v2, un_S08v2,
                   st_S09v2, un_S09v2, st_S10v2, un_S10v2, st_S12v2, un_S12v2, st_S12v2, un_S12v2,
                   st_S13v2, un_S13v2, st_S14v2, un_S14v2, st_S15v2, un_S15v2, st_S16v2, un_S16v2)


# Bonferroni correction for alpha level
# Divided by 5 because using same data for multiple tests (5 different groups)
print(alphaAdj <- 0.05/10)
# 0.005

# One Sample t-test
# proportion of fixations on target by monolinguals when tense = present
mon_pres_v2 <- st10_onv2 %>% 
  filter(., group == "mon", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 7.8217, df = 29, p-value = 6.303e-09

# proportion of fixations on target by monolinguals when tense = preterit
mon_pret_v2 <- st10_onv2 %>% 
  filter(., group == "mon", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 5.5193, df = 29, p-value = 2.998e-06


# proportion of fixations on target by intermediate EN when tense = present
ies_pres_v2 <- st10_onv2 %>% 
  filter(., group == "ies", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.6194, df = 13, p-value = 0.9894


# proportion of fixations on target by intermediate EN when tense = preterit
ies_pret_v2 <- st10_onv2 %>% 
  filter(., group == "ies", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.85843, df = 13, p-value = 0.2031


# proportion of fixations on target by advanced EN when tense = present
aes_pres_v2 <- st10_onv2 %>% 
  filter(., group == "aes", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.4616, df = 20, p-value = 0.9885


# proportion of fixations on target by advanced EN when tense = preterit
aes_pret_v2 <- st10_onv2 %>% 
  filter(., group == "aes", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.93483, df = 20, p-value = 0.1805



# proportion of fixations on target by intermediate CH when tense = present
ims_pres_v2 <- st10_onv2 %>% 
  filter(., group == "ims", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = -1.569, df = 13, p-value = 0.9297


# proportion of fixations on target by intermediate CH when tense = preterit
ims_pret_v2 <- st10_onv2 %>% 
  filter(., group == "ims", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 1.3879, df = 13, p-value = 0.09424


# proportion of fixations on target by advanced CH when tense = present
ams_pres_v2 <- st10_onv2 %>% 
  filter(., group == "ams", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = -0.90233, df = 18, p-value = 0.8106


# proportion of fixations on target by advanced CH when tense = preterit
ams_pret_v2 <- st10_onv2 %>% 
  filter(., group == "ams", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 1.1074, df = 18, p-value = 0.1414







