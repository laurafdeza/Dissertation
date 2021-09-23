library(tidyverse); library(TOSTER)

# Load data
#dem_all <- read_csv(here::here("data", "clean", "ospan_set_z_scores.csv"))
dem_all <- read_csv(here::here("data", "pupurri_analysis.csv"))

# Clean data

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

# Remove participants to make the groups homogeneous in L2 proficiency (DELE)

dem_all %>%
  filter(l1 != 'ma') %>%
  group_by(l1) %>%
  summarize(n = n()) # en = English speakers, ma = Mandarin Chinese speakers, es = Spanish speakers
#   l1        n
# 1 en       65
# 2 es       30



dem_all %>%
  filter(l1 != 'ma') %>%
  group_by(., l1) %>%
  summarise(
            age_mean = mean(age),
            age_sd = sd(age),
            ospan_mean = round(mean(WM_set),2),   # ospan instead of WM_set if using z scores
            ospan_sd = round(sd(WM_set),2),       # ospan instead of WM_set if using z scores
            n = length(unique(participant))) %>%
  knitr::kable()
# |l1 | age_mean|   age_sd| ospan_mean| ospan_sd|  n|
# |:--|--------:|--------:|----------:|--------:|--:|
# |en | 26.75385| 4.541581|       8.89|     2.11| 65|
# |es | 26.16667| 8.824274|       6.20|     2.72| 30|

  
  

dem_all$age_fluent_L2 <- as.numeric(dem_all$age_fluent_L2)
dem_all %>%
  filter(l1 == "en") %>% 
  # group_by(l1) %>%
  summarise(mean_AoA = mean(AoA_L2),
            sd_AoA = sd(AoA_L2),
            mean_fluent = mean(age_fluent_L2),
            sd_fluent = sd(age_fluent_L2),
            mean_abroad = mean(mo_ES_country),
            sd_abroad = sd(mo_ES_country),
            mean_prof = mean(DELE),
            sd_prof = sd(DELE),
            mean_use = mean(percent_l2_week),
            sd_use = sd(percent_l2_week))
#    mean_AoA sd_AoA mean_fluent sd_fluent mean_abroad sd_abroad mean_prof sd_prof mean_use sd_use
#        16.3   5.55        22.1      4.35        38.1      33.5      38.5    8.19     33.3   17.4



dem_all %>%
  filter(l1 == "en") %>% 
  group_by(l2_school) %>%
  summarise(educ_es = n())
#   l2_school  educ_es
# 1 no              26
# 2 U.              39













# ---------------------------------------------------------------------



arrange(dem_all, desc(WM_set)) %>% filter(l1 != 'ma') %>% View()

bartlett.test(WM_set ~ l1, data = dem_all %>%   filter(l1 != 'ma' & participant != 'ies09' & 
                                                         participant != 'ies16' & participant != 'aes03' &
                                                         participant != 'aes04' & participant != 'aes16' &
                                                         participant != 'aes23' & participant != 'ies05' &
                                                         participant != 'ies21' & participant != 'ies23' &
                                                         participant != 'ies33' & participant != 'aes01' &
                                                         participant != 'aes05' & participant != 'aes07' &
                                                         participant != 'aes11' & participant != 'aes14' &
                                                         participant != 'aes28' &
                                                         participant != 'ies30' & participant != 'mon20' &
                                                         participant != 'mon04' & participant != 'mon24' &
                                                         participant != 'mon21' & participant != 'mon27' & 
                                                         participant != 'mon23'))
# K-squared = 1.7481, df = 1, p-value = 0.1861

t.test(WM_set ~ l1, data = dem_all %>% filter(l1 != 'ma' & participant != 'ies09' & 
                                                participant != 'ies16' & participant != 'aes03' &
                                                participant != 'aes04' & participant != 'aes16' &
                                                participant != 'aes23' & participant != 'ies05' &
                                                participant != 'ies21' & participant != 'ies23' &
                                                participant != 'ies33' & participant != 'aes01' &
                                                participant != 'aes05' & participant != 'aes07' &
                                                participant != 'aes11' & participant != 'aes14' &
                                                participant != 'aes28' &
                                                participant != 'ies30' & participant != 'mon20' &
                                                participant != 'mon04' & participant != 'mon24' &
                                                participant != 'mon21' & participant != 'mon27' & 
                                                participant != 'mon23'), var.equal = TRUE) 
# t = -1.9516, df = 123, p-value = 0.05326


dem_all %>%
  filter(l1 != 'ma' & participant != 'ies09' & 
           participant != 'ies16' & participant != 'aes03' &
           participant != 'aes04' & participant != 'aes16' &
           participant != 'aes23' & participant != 'ies05' &
           participant != 'ies21' & participant != 'ies23' &
           participant != 'ies33' & participant != 'aes01' &
           participant != 'aes05' & participant != 'aes07' &
           participant != 'aes11' & participant != 'aes14' &
           participant != 'aes28' &
           participant != 'ies30' & participant != 'mon20' &
           participant != 'mon04' & participant != 'mon24' &
           participant != 'mon21' & participant != 'mon27' & 
           participant != 'mon23') %>%
  #filter(., group %in% c("hs", "l2")) %>%
  group_by(., l1) %>%
  summarise(#mean_perc_week_Spa = mean(percent_l2_week),
    # sd_perc_week_Spa = sd(percent_l2_week),
    # mean_DELE = mean(DELE),
    # sd_DELE = sd(DELE),
    ospan_mean = round(mean(WM_set),2),   # ospan instead of WM_set if using z scores
    ospan_sd = round(sd(WM_set),2),       # ospan instead of WM_set if using z scores
    n = length(unique(participant))) %>%
  knitr::kable()
# |l1 | ospan_mean| ospan_sd|  n|
# |:--|----------:|--------:|--:|
# |en |       8.02|     1.74| 48|
# |es |       7.08|     2.21| 24|



# -----------------------------

arrange(dem_all, desc(age)) %>% filter(l1 != 'ma') %>% View()

bartlett.test(age ~ l1, data = dem_all %>% filter(l1 != 'ma' & participant != 'mon26' & 
                                                    participant != 'mon28' & participant != 'mon27' &
                                                    participant != 'mon23' & participant != 'mon30' &
                                                    participant != 'mon24' & participant != 'mon19' &
                                                    participant != 'mon16'))
# Bartlett's K-squared = 1.7362, df = 1, p-value = 0.1876



t.test(age ~ l1, data = dem_all %>% filter(l1 != 'ma' & participant != 'mon26' & 
                                             participant != 'mon28' & participant != 'mon27' &
                                             participant != 'mon23' & participant != 'mon30' &
                                             participant != 'mon24' & participant != 'mon19' &
                                             participant != 'mon16'), var.equal = TRUE)
# t = 1.9201, df = 85, p-value = 0.0582



dem_all %>%
  filter(l1 != 'ma' & participant != 'mon26' & 
           participant != 'mon28' & participant != 'mon27' &
           participant != 'mon23' & participant != 'mon30' &
           participant != 'mon24' & participant != 'mon19' &
           participant != 'mon16') %>%
  group_by(., l1) %>%
  summarise(age_mean = mean(age),
            age_sd = sd(age),
            n = length(unique(participant))) %>%
  knitr::kable()
# |l1 | age_mean|   age_sd|  n|
# |:--|--------:|--------:|--:|
# |en | 26.75385| 4.541581| 65|
# |es | 24.45455| 5.704862| 22|







# TOSTER::TOSTtwo(m1 = 26.8, sd1 = 4.54, n1 = 65, # EN
#                 m2 = 24.7, sd2 = 4.13, n2 = 64, # MA
#                 low_eqbound_d = -0.3,
#                 high_eqbound_d = 0.3,
#                 alpha = 0.05)

