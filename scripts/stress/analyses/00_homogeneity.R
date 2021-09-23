library(tidyverse); library(TOSTER)

# Load data
#dem_all <- read_csv(here::here("data", "clean", "ospan_set_z_scores.csv"))
dem_all <- read_csv(here::here("data", "pupurri_analysis.csv"))


# Remove participants to make the groups homogeneous in L2 proficiency (DELE)

unique(dem_all$participant)

dem_all <- dem_all %>%
  # separate(., col = participant,
  #          into = c('group', 'id'),
  #          sep = 3,
  #          remove = FALSE) %>%
  separate(., col = group,
           into = c("proficiency", "l1"), # es = EN speaker, ms = MA speaker, on = ES speaker
           sep = 1,
           remove = FALSE)

dem_all$l1 <- str_replace(dem_all$l1, "es", "en")
dem_all$l1 <- str_replace(dem_all$l1, "ms", "ma")
dem_all$l1 <- str_replace(dem_all$l1, "on", "es")


dem_all$participant <- tolower(dem_all$participant)
dem_all$DELE <- as.numeric(as.character(dem_all$DELE))
dem_all$percent_l2_week <- as.numeric(as.character(dem_all$percent_l2_week))


dem_all %>%
  group_by(l1) %>%
  summarize(n = n()) # en = English speakers, ma = Mandarin Chinese speakers, es = Spanish speakers
#   l1        n
# 1 en       65
# 2 es       30
# 3 ma       64


dem_all %>%
  #filter(., group %in% c("hs", "l2")) %>%
  group_by(., l1) %>%
  summarise(mean_perc_week_Spa = mean(percent_l2_week),
            sd_perc_week_Spa = sd(percent_l2_week),
            mean_DELE = mean(DELE),
            sd_DELE = sd(DELE),
            ospan_mean = round(mean(WM_set),2),   # ospan instead of WM_set if using z scores
            ospan_sd = round(sd(WM_set),2),          # ospan instead of WM_set if using z scores
            n = length(unique(participant))) %>%
  knitr::kable()
# |l1 | mean_perc_week_Spa| sd_perc_week_Spa| mean_DELE|  sd_DELE| ospan_mean| ospan_sd|  n|
# |:--|------------------:|----------------:|---------:|--------:|----------:|--------:|--:|
# |en |           33.30769|         17.43911|  38.47692| 8.185558|       8.89|     2.11| 65|
# |es |           88.93333|          9.43191|        NA|       NA|       6.20|     2.72| 30|
# |ma |           41.64062|         21.65851|  39.17188| 7.564652|       7.78|     2.14| 64|


t.test(percent_l2_week ~ l1, data = dem_all %>% filter(l1 != 'es' & participant != 'ies04' & participant != 'ies17' & 
                                                         participant != 'aes32' & participant != 'ies28'), var.equal = TRUE)
# t = -1.9516, df = 123, p-value = 0.05326
bartlett.test(percent_l2_week ~ l1, data = dem_all %>% filter(l1 != 'es' & participant != 'ies04' & participant != 'ies17' & 
                                                                participant != 'aes32' & participant != 'ies28'))
# K-squared = 3.6896, df = 1, p-value = 0.05475
t.test(DELE ~ l1, data = dem_all %>% filter(l1 != 'es' & participant != 'ies04' & participant != 'ies17' & 
                                              participant != 'aes32' & participant != 'ies28'), var.equal = TRUE)
# t = -0.14658, df = 123, p-value = 0.8837
bartlett.test(DELE ~ l1, data = dem_all %>% filter(l1 != 'es' & participant != 'ies04' & participant != 'ies17' & 
                                                     participant != 'aes32' & participant != 'ies28'))
# K-squared = 0.23176, df = 1, p-value = 0.6302





bartlett.test(WM_set ~ l1, data = dem_all %>% filter(participant != 'ies04' & participant != 'ies17' & 
                                                       participant != 'aes32' & participant != 'ies28'))
# Bartlett's K-squared = 2.8135, df = 2, p-value = 0.2449


t.test(WM_set ~ l1, data = dem_all %>% filter(l1 != 'es' & participant != 'ies04' & participant != 'ies17' & 
                                                           participant != 'aes32' & participant != 'ies28'), var.equal = TRUE)
# t = 2.9647, df = 127, p-value = 0.002919

t.test(WM_set ~ l1, data = dem_all %>% filter(l1 != "en" & participant != 'ies04' & participant != 'ies17' & 
                                                participant != 'aes32' & participant != 'ies28'), var.equal = TRUE)
# t = -3.0541, df = 92, p-value = 0.002953

t.test(WM_set ~ l1, data = dem_all %>% filter(l1 != "ma" & participant != 'ies04' & participant != 'ies17' & 
                                                participant != 'aes32' & participant != 'ies28'), var.equal = TRUE)
# t = 5.2551, df = 93, p-value = 1.127e-06



dem_all %>%
  filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
           participant != 'aes32' & participant != 'ies28') %>%
  group_by(., l1) %>%
  summarise(mean_perc_week_Spa = mean(percent_l2_week),
            sd_perc_week_Spa = sd(percent_l2_week),
            mean_DELE = mean(DELE),
            sd_DELE = sd(DELE),
            n = length(unique(participant))) %>%
  knitr::kable()
# |l1 | mean_perc_week_Spa| sd_perc_week_Spa| mean_DELE|  sd_DELE|  n|
# |:--|------------------:|----------------:|---------:|--------:|--:|
# |en |           34.83607|         16.90580|  38.96721| 8.045635| 61|
# |ma |           41.64062|         21.65851|  39.17188| 7.564652| 64|

# Proficiency equivalence test
TOSTtwo(m1 = 38.97, sd1 = 8.05, n1 = 61, # EN
        m2 = 39.17, sd2 = 7.57, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
  # null hypothesis test non significant
# The equivalence test was non-significant, t(121.54) = 1.532, p = 0.064, given equivalence bounds of -2.344 and 2.344 (on a raw scale) and an alpha of 0.05.
# 
# The null hypothesis test was non-significant, t(121.54) = -0.143, p = 0.887, given an alpha of 0.05.


# L2 use equivalence test
TOSTtwo(m1 = 34.84, sd1 = 16.91, n1 = 61, # EN
        m2 = 41.64, sd2 = 21.66, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
  # null hypothesis test non significant
# The equivalence test was non-significant, t(118.47) = -0.280, p = 0.610, given equivalence bounds of -5.829 and 5.829 (on a raw scale) and an alpha of 0.05.
# 
# The null hypothesis test was non-significant, t(118.47) = -1.962, p = 0.0522, given an alpha of 0.05.

dem_all %>%
  filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
           participant != 'aes32' & participant != 'ies28') %>%
  summarise(mean_perc_week_Spa = mean(percent_l2_week),
            sd_perc_week_Spa = sd(percent_l2_week),
            mean_DELE = mean(DELE),
            sd_DELE = sd(DELE),
            n = length(unique(participant))) %>%
  knitr::kable()
# | mean_perc_week_Spa| sd_perc_week_Spa| mean_DELE|  sd_DELE|   n|
# |------------------:|----------------:|---------:|--------:|---:|
# |              38.32|         19.70492|    39.072| 7.772134| 125|


dem_all$age_fluent_L2 <- as.numeric(dem_all$age_fluent_L2)
dem_all %>%
  filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
           participant != 'aes32' & participant != 'ies28') %>% 
  group_by(l1) %>%
  summarise(mean_AoA = mean(AoA_L2),
            sd_AoA = sd(AoA_L2),
            mean_fluent = mean(age_fluent_L2),
            sd_fluent = sd(age_fluent_L2),
            mean_abroad = mean(mo_ES_country),
            sd_abroad = sd(mo_ES_country),
            n = length(unique(participant)))
#   l1    mean_AoA sd_AoA mean_fluent sd_fluent mean_abroad sd_abroad    n
# 1 en        16.0   5.41        22.0      4.22        38.1      34.0   61
# 2 ma        18.9   3.59        21.2      3.91        40.8      45.5   64



dem_all %>%
  filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
           participant != 'aes32' & participant != 'ies28') %>% 
  group_by(l1, l2_school) %>%
  summarise(educ_es = n())
# l1    l2_school  educ_es
# 1 en    no              24
# 2 en    U.              37
# 3 ma    H., U.           4
# 4 ma    M., H., U.       2
# 5 ma    no               2
# 6 ma    U.              56


# Two onesided test of equivalence for months spent in Spain (null hypothesis test non-significant)
TOSTtwo(m1 = 38.1, sd1 = 34.0, n1 = 61, # EN
        m2 = 40.8, sd2 = 45.5, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# The equivalence test was non-significant, t(116.47) = 1.305, p = 0.0972, given equivalence bounds of -12.049 and 12.049 (on a raw scale) and an alpha of 0.05.
# 
# The null hypothesis test was non-significant, t(116.47) = -0.377, p = 0.707, given an alpha of 0.05.




bartlett.test(mo_ES_country ~ l1, data = dem_all %>% filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
                                                       participant != 'aes32' & participant != 'ies28'))
# Bartlett's K-squared = 5.0992, df = 1, p-value = 0.02394

t.test(mo_ES_country ~ l1, data = dem_all %>% filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
                                                participant != 'aes32' & participant != 'ies28'), var.equal = FALSE)
# t = -0.38728, df = 116.44, p-value = 0.6993



