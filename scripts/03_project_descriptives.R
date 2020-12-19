# Project descriptives --------------------------------------------------------
#
# Description: get basic project descriptives for sanity checks
# Last update: 6/03/2019
#
# -----------------------------------------------------------------------------


# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------




## Linguistic questionnaire info

# Load .csv file from google drive

allinfo <- read.csv(here::here("data", "pupurri_analysis.csv"))

#glimpse(allinfo)

#unique(allinfo$participant)

allinfo$participant <- tolower(allinfo$participant)
allinfo$DELE <- as.numeric(as.character(allinfo$DELE))
allinfo$age_fluent_L2 <- as.numeric(as.character(allinfo$age_fluent_L2))
allinfo$percent_l1_week <- as.numeric(as.character(allinfo$percent_l1_week))
allinfo$percent_l2_week <- as.numeric(as.character(allinfo$pc_l2_week))
#######################################################################################################
# for groups not collapsed, see below

agg <- separate(data = allinfo,
                       col = group,
                       into = c("prof", "l1"),
                       sep = 1,
                       remove = FALSE) 

agg$l1 <- str_replace(agg$l1, "es", "en")
agg$l1 <- str_replace(agg$l1, "ms", "ma")
agg$l1 <- str_replace(agg$l1, "on", "es")

# No. participants
agg %>%
  group_by(., l1) %>%
  tally()
#   l1        n
# <chr> <int>
# 1 en       65
# 2 ma       64
# 3 es       30

# Females per group
agg %>%
  group_by(., l1, gender) %>%
  tally()
# l1    gender     n
# <chr> <fct>  <int>
# 1 en  female    48  - EN natives
# 2 en  male      17
# 3 ma  female    52  - MA natives
# 4 ma  male      12
# 5 es  female    20  - ES natives
# 6 es  male      10

# Get mean AGE as a function of group + SD
agg %>%
  group_by(l1) %>%
  summarize(min_age = min(age), 
            max_age = max(age),
            mean_age = mean(age),
            sd_age = sd(age))
# l1    min_age max_age mean_age sd_age
# 1 en       20      38     26.8   4.54
# 2 ma       18      41     24.7   4.13
# 3 es       18      45     26.2   8.82


# Get mean TIME ABROAD in months as a function of group + SD
agg %>%
  filter(., group != "mon") %>%
  group_by(., l1) %>%
  summarise(., max_abroad = max(mo_ES_country),
            min_abroad = min(mo_ES_country),
            mean_abroad = mean(mo_ES_country),
            sd_abroad = sd(mo_ES_country))
# l1    max_abroad min_abroad mean_abroad sd_abroad
# 1 en         168        2.5        38.1      33.5
# 2 ma         228        1          40.8      45.5


# Get mean DELE as a function of group + SD
agg %>%
  filter(., group %in% c("ies", "aes", "ims", "ams")) %>%
  group_by(., l1) %>%
  summarise(mean_DELE = mean(DELE),
            sd_DELE = sd(DELE))
# l1    mean_DELE sd_DELE
# 1 en       38.5    8.19
# 2 ma       39.2    7.56


# When participants started learning the L2
agg %>%
  filter(., group %in% c("ies", "aes", "ims", "ams")) %>%
  group_by(., l1) %>%
  summarise(., mean_AoA = mean(AoA_L2),
            sd_AoA = sd(AoA_L2),
            mean_fluentL2 = mean(AoA_L2),
            sd_fluentL2 = sd(AoA_L2))
# l1    mean_AoA sd_AoA mean_fluentL2 sd_fluentL2 
# 1 en      16.3   5.55          16.3        5.55
# 2 ma      18.9   3.59          18.9        3.59


# L2 use per week in %
agg %>%
  filter(., group %in% c("aes", "ies", "ams", "ims")) %>%
  group_by(., l1) %>%
  summarise(mean_l2use_week = mean(percent_l2_week),
            sd_l2use_week = sd(percent_l2_week))
# l1    mean_l2use_week sd_l2use_week
# 1 en             33.3          17.4
# 2 ma             41.6          21.7


# country of origin
agg <- agg[!is.na(agg$country),]

english <- filter(agg, l1 == 'en') # descriptives show all, not only EN speakers
table(english$country)
# count
# au bh ca ch es ir nz tw uk us 
#  3  1  2 63  0  1  2  1 25 31

prop.table(table(english$country))
# au         bh         ca         ch         es 
# 0.04615385 0.01538462 0.03076923 0.00000000 0.00000000 
# ir         nz         tw         uk         us 
# 0.01538462 0.03076923 0.00000000 0.38461538 0.47692308


#############################################################################################################

# PROFICIENCIES NOT COLLAPSED


#####
# Females per group

allinfo %>%
  group_by(., group) %>%
  tally()
#   group     n
# 1 aes      32
# 2 ams      32
# 3 ies      33
# 4 ims      32
# 5 mon      30

allinfo %>%
  filter(., gender == "female") %>%
  group_by(., group) %>%
  tally()

# group n_females
# 1 aes        24
# 2 ams        26
# 3 ies        24
# 4 ims        26
# 5 mon        20


#################################
# Get mean AGE as a function of group + SD

allinfo %>%
  group_by(., group) %>%
  summarise(min_age = min(age),
            max_age = max(age),
            mean_age = mean(age),
            sd_age = sd(age))
# group min_age max_age mean_age sd_age
# 1 aes      20      38     27.5   4.83
# 2 ams      18      41     24.8   4.37
# 3 ies      21      35     26.1   4.21
# 4 ims      19      37     24.5   3.95
# 5 mon      18      45     26.2   8.82



#################################
# Get mean TIME ABROAD in months as a function of group + SD

allinfo %>%
  filter(., group != "mon") %>%
  group_by(., group) %>%
  summarise(max_abroad = max(mo_ES_country),
            min_abroad = min(mo_ES_country),
            mean_abroad = mean(mo_ES_country),
            sd_abroad = sd(mo_ES_country))
# group max_abroad min_abroad mean_abroad sd_abroad
# 1 aes        168        7          51.3      38.8
# 2 ams        228        1          45.1      55.7
# 3 ies         72        2.5        25.3      20.9
# 4 ims        129        6          36.5      32.6




#################################
# Get mean DELE as a function of group + SD

# allinfo$DELE <- as.numeric(as.character(allinfo$DELE))

allinfo %>%
  filter(., group %in% c("ies", "aes", "ims", "ams")) %>%
  group_by(., group) %>%
  summarise(mean_DELE = mean(DELE),
            sd_DELE = sd(DELE))
# group mean_DELE sd_DELE
# 1 aes      45.4    4.26
# 2 ams      45.5    3.97
# 3 ies      31.7    4.58
# 4 ims      32.8    4.23





################################
# Get mean AoA as a function of group + SD

# allinfo$age_fluent_L2 <- as.numeric(as.character(allinfo$age_fluent_L2))

allinfo %>%
  filter(., group %in% c("ies", "aes", "ims", "ams")) %>%
  group_by(., group) %>%
  summarise(., mean_AoA = mean(AoA_L2),
            sd_AoA = sd(AoA_L2),
            mean_fluentL2 = mean(age_fluent_L2),
            sd_fluentL2 = sd(age_fluent_L2))
# group mean_AoA sd_AoA mean_fluentL2 sd_fluentL2
# 1 aes     15.1   4.35          20.8        3.73
# 2 ams     17.9   2.83          20.3        3.12
# 3 ies     17.5   6.35          23.2        4.63
# 4 ims     19.9   4.01          22.1        4.46



#################################
# Get mean L2 USE as a function of group + SD

# allinfo$pc_en_ma <- as.numeric(as.character(allinfo$pc_en_ma))
# allinfo$pc_es <- as.numeric(as.character(allinfo$pc_es))

allinfo %>%
  filter(., group %in% c("aes", "ies", "ams", "ims")) %>%
  group_by(., group) %>%
  summarise(mean_l2use_week = mean(percent_l2_week),
            sd_l2use_week = sd(percent_l2_week))
# group mean_l2use_week sd_l2use_week
# 1 aes            38.6          16.2
# 2 ams            46.6          21.8
# 3 ies            28.2          17.3
# 4 ims            36.7          20.7

