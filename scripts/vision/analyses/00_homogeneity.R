library(tidyverse); library(TOSTER)

# Load data

demographics <- read_csv(here::here("data", "pupurri_analysis.csv"))
wm <- read_csv('./data/clean/wm_processing_speed.csv')

demographics$participant <- tolower(demographics$participant)

dem_all <- merge(demographics, wm, by = 'participant')

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



dem_all$DELE <- as.numeric(as.character(dem_all$DELE))
dem_all$percent_l2_week <- as.numeric(as.character(dem_all$percent_l2_week))


dem_all %>%
  group_by(l1) %>%
  summarize(n = n()) # en = English speakers, ma = Mandarin Chinese speakers, es = Spanish speakers
#   l1        n
# 1 en       65
# 2 es       29
# 3 ma       64


dem_all %>%
  #filter(., group %in% c("hs", "l2")) %>%
  group_by(., l1) %>%
  summarise(mean_perc_week_Spa = mean(percent_l2_week),
            sd_perc_week_Spa = sd(percent_l2_week),
            mean_DELE = mean(DELE),
            sd_DELE = sd(DELE),
            ospan_mean = round(mean(ospan_rt),2),   
            ospan_sd = round(sd(ospan_sd),2),  
            corsi_mean = round(mean(corsi_rt),2),   
            corsi_sd = round(sd(corsi_sd),2), 
            n = length(unique(participant))) %>%
  knitr::kable()
# |l1 | mean_perc_week_Spa|  sd_perc_week_Spa| mean_DELE|  sd_DELE|  ospan_mean| ospan_sd| corsi_mean| corsi_sd|  n|
# |:--|------------------:|-----------------:|---------:|--------:|-----------:|--------:|----------:|--------:|--:|
# |en |           33.30769|         17.439111|  38.47692| 8.185558|       -0.08|        0|      -0.34|     0.03| 65|
# |es |           88.72414|          9.527766|        NA|       NA|        0.23|        0|       0.37|     0.07| 29|
# |ma |           41.64062|         21.658509|  39.17188| 7.564652|       -0.03|        0|       0.17|     0.03| 64|


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


bartlett.test(ospan_rt ~ l1, data = dem_all %>% filter(participant != 'ies04' & participant != 'ies17' & 
                                                       participant != 'aes32' & participant != 'ies28'))
# Bartlett's K-squared = 13.113, df = 2, p-value = 0.001421


t.test(ospan_rt ~ l1, data = dem_all %>% filter(l1 != 'es' & participant != 'ies04' & participant != 'ies17' & 
                                                           participant != 'aes32' & participant != 'ies28'), var.equal = TRUE)
# t = -0.81529, df = 123, p-value = 0.4165

t.test(corsi_rt ~ l1, data = dem_all %>% filter(l1 != "en" & participant != 'ies04' & participant != 'ies17' & 
                                                participant != 'aes32' & participant != 'ies28'), var.equal = TRUE)
# t = 0.99238, df = 91, p-value = 0.3236

bartlett.test(corsi_rt ~ l1, data = dem_all %>% filter(participant != 'ies04' & participant != 'ies17' & 
                                                         participant != 'aes32' & participant != 'ies28'))
# Bartlett's K-squared = 19.732, df = 2, p-value = 5.192e-05




dem_all %>%
  filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
           participant != 'aes32' & participant != 'ies28' & participant != 'ams09') %>%
  group_by(., l1) %>%
  summarise(mean_perc_week_Spa = mean(percent_l2_week),
            sd_perc_week_Spa = sd(percent_l2_week),
            mean_DELE = mean(DELE),
            sd_DELE = sd(DELE),
            n = length(unique(participant))) %>%
  knitr::kable()
# With ams09
# |l1 | mean_perc_week_Spa| sd_perc_week_Spa| mean_DELE|  sd_DELE|  n|
# |:--|------------------:|----------------:|---------:|--------:|--:|
# |en |           34.83607|         16.90580|  38.96721| 8.045635| 61|
# |ma |           41.64062|         21.65851|  39.17188| 7.564652| 64|


# Without ams09
# |l1 | mean_perc_week_Spa| sd_perc_week_Spa| mean_DELE|  sd_DELE|  n|
# |:--|------------------:|----------------:|---------:|--------:|--:|
# |en |           34.83607|         16.90580|  38.96721| 8.045635| 61|
# |ma |           41.50794|         21.80624|  39.15873| 7.624676| 63|

# Proficiency equivalence test (w/ ams09)
TOSTtwo(m1 = 38.97, sd1 = 8.05, n1 = 61, # EN
        m2 = 39.17, sd2 = 7.57, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
  # null hypothesis test non significant
# The equivalence test was non-significant, t(121.54) = 1.532, p = 0.064, given equivalence bounds of -2.344 and 2.344 (on a raw scale) and an alpha of 0.05.
# 
# The null hypothesis test was non-significant, t(121.54) = -0.143, p = 0.887, given an alpha of 0.05.


# Proficiency equivalence test (w/o ams09)
TOSTtwo(m1 = 38.97, sd1 = 8.05, n1 = 61, # EN
        m2 = 39.16, sd2 = 7.62, n2 = 63, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# The equivalence test was non-significant, t(121.08) = 1.534, p = 0.0638, given equivalence bounds of -2.351 and 2.351 (on a raw scale) and an alpha of 0.05.
# 
# The null hypothesis test was non-significant, t(121.08) = -0.135, p = 0.893, given an alpha of 0.05.


# L2 use equivalence test (w/ ams09)
TOSTtwo(m1 = 34.84, sd1 = 16.91, n1 = 61, # EN
        m2 = 41.64, sd2 = 21.66, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
  # null hypothesis test non significant
# The equivalence test was non-significant, t(118.47) = -0.280, p = 0.610, given equivalence bounds of -5.829 and 5.829 (on a raw scale) and an alpha of 0.05.
# 
# The null hypothesis test was non-significant, t(118.47) = -1.962, p = 0.0522, given an alpha of 0.05.

# w/o
TOSTtwo(m1 = 34.84, sd1 = 16.91, n1 = 61, # EN
        m2 = 41.51, sd2 = 21.81, n2 = 63, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# The equivalence test was non-significant, t(116.49) = -0.233, p = 0.592, given equivalence bounds of -5.854 and 5.854 (on a raw scale) and an alpha of 0.05.
# 
# The null hypothesis test was non-significant, t(116.49) = -1.907, p = 0.059, given an alpha of 0.05.



dem_all$age_fluent_L2 <- as.numeric(dem_all$age_fluent_L2)
dem_all %>%
  filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
           participant != 'aes32' & participant != 'ies28' & participant != 'ams09') %>% 
  group_by(l1) %>%
  summarise(mean_AoA = mean(AoA_L2),
            sd_AoA = sd(AoA_L2),
            mean_fluent = mean(age_fluent_L2),
            sd_fluent = sd(age_fluent_L2),
            mean_abroad = mean(mo_ES_country),
            sd_abroad = sd(mo_ES_country),
            n = length(unique(participant)))
# with ams09
#   l1    mean_AoA sd_AoA mean_fluent sd_fluent mean_abroad sd_abroad    n
# 1 en        16.0   5.41        22.0      4.22        38.1      34.0   61
# 2 ma        18.9   3.59        21.2      3.91        40.8      45.5   64

# without ams09
#   l1    mean_AoA sd_AoA mean_fluent sd_fluent mean_abroad sd_abroad     n
# 1 en        16.0   5.41        22.0      4.22        38.1      34.0    61
# 2 ma        18.8   3.61        21.2      3.94        40.7      45.8    63



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
# with ams09
TOSTtwo(m1 = 38.1, sd1 = 34.0, n1 = 61, # EN
        m2 = 40.8, sd2 = 45.5, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# The equivalence test was non-significant, t(116.47) = 1.305, p = 0.0972, given equivalence bounds of -12.049 and 12.049 (on a raw scale) and an alpha of 0.05.
# 
# The null hypothesis test was non-significant, t(116.47) = -0.377, p = 0.707, given an alpha of 0.05.

# without ams09
TOSTtwo(m1 = 38.1, sd1 = 34.0, n1 = 61, # EN
        m2 = 40.7, sd2 = 45.8, n2 = 63, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
#The equivalence test was non-significant, t(114.37) = 1.314, p = 0.0957, given equivalence bounds of -12.100 and 12.100 (on a raw scale) and an alpha of 0.05.


bartlett.test(mo_ES_country ~ l1, data = dem_all %>% filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
                                                       participant != 'aes32' & participant != 'ies28'))
# Bartlett's K-squared = 5.0992, df = 1, p-value = 0.02394

t.test(mo_ES_country ~ l1, data = dem_all %>% filter(l1 != "es" & participant != 'ies04' & participant != 'ies17' & 
                                                participant != 'aes32' & participant != 'ies28'), var.equal = FALSE)
# t = -0.38728, df = 116.44, p-value = 0.6993




dem_all %>%
  filter(participant != 'ies04' & participant != 'ies17' & 
           participant != 'aes32' & participant != 'ies28' & participant != 'ams09') %>%
  group_by(., l1) %>%
  summarise( ospan_mean = round(mean(ospan_rt),2),   
             ospan_sd = round(sd(ospan_sd),2),  
             corsi_mean = round(mean(corsi_rt),2),   
             corsi_sd = round(sd(corsi_sd),2),
            n = length(unique(participant))) %>%
  knitr::kable()
# with ams09
# |l1 | ospan_mean| ospan_sd| corsi_mean| corsi_sd|  n|
# |:--|----------:|--------:|----------:|--------:|--:|
# |en |      -0.09|        0|      -0.32|     0.03| 61|
# |es |       0.23|        0|       0.37|     0.07| 29|
# |ma |      -0.03|        0|       0.17|     0.03| 64|

# without ams09
# |l1 | ospan_mean| ospan_sd| corsi_mean| corsi_sd|  n|
# |:--|----------:|--------:|----------:|--------:|--:|
# |en |      -0.09|        0|      -0.32|     0.03| 61|
# |es |       0.23|        0|       0.37|     0.07| 29|
# |ma |      -0.03|        0|       0.18|     0.03| 63|


# corsi equivalence test
TOSTtwo(m1 = -0.32, sd1 = 0.03, n1 = 61, # EN
        m2 = 0.17, sd2 = 0.03, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# null hypothesis test non significant
# The equivalence test was non-significant, t(122.71) = -89.603, p = 1.000, given equivalence bounds of -0.009 and 0.009 (on a raw scale) and an alpha of 0.05.

# Null Hypothesis Test Result:
# The null hypothesis test was significant, t(122.71) = -91.280, p = 1.19e-114, given an alpha of 0.05.

TOSTtwo(m1 = -0.32, sd1 = 0.03, n1 = 61, # EN
        m2 = 0.37, sd2 = 0.07, n2 = 29, # ES
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# Equivalence Test Result:
#   The equivalence test was non-significant, t(32.99) = -49.714, p = 1.000, given equivalence bounds of -0.0162 and 0.0162 (on a raw scale) and an alpha of 0.05.
# 
# Null Hypothesis Test Result:
#   The null hypothesis test was significant, t(32.99) = -50.906, p = 0.00000000000000000000000000000000622, given an alpha of 0.05.

TOSTtwo(m1 = 0.17, sd1 = 0.03, n1 = 64, # MA
        m2 = 0.37, sd2 = 0.07, n2 = 29, # ES
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# Equivalence Test Result:
#   The equivalence test was non-significant, t(32.75) = -13.589, p = 1.000, given equivalence bounds of -0.0162 and 0.0162 (on a raw scale) and an alpha of 0.05.
# 
# Null Hypothesis Test Result:
#   The null hypothesis test was significant, t(32.75) = -14.783, p = 0.000000000000000475, given an alpha of 0.05.
