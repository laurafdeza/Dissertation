# DELE analysis ---------------------------------------------------------------



# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------


dem_all <- read_csv("./data/pupurri_analysis.csv")

### DELE

dem_all$DELE <- as.numeric(dem_all$DELE)
dem_all %>% 
  filter(., group %in% c("ies", "aes", "ims", "ams")) %>%
  group_by(., group) %>% 
  summarize(., DELE_mean = mean(DELE), 
            DELE_sd = sd(DELE), 
            n = n())

# results including all participants collected

#   group   DELE_mean  DELE_sd n
# 1 aes        45.4    4.26    32
# 2 ams        45.5    3.97    32
# 3 ies        31.7    4.58    33
# 4 ims        32.8    4.23    32


bartlett.test(DELE ~ group, data = dem_all) # all good
# Bartlett's K-squared = 0.84963, df = 4, p-value = 0.9317

# DELE tost aes vs ams (OK, not significant and between dotted lines so they are not different)
TOSTER::TOSTtwo(m1 = 45.4, sd1 = 4.26, n1 = 32, # aes
        m2 = 45.5, sd2 = 3.97, n2 = 32, # ams
        low_eqbound_d = -0.3, 
        high_eqbound_d = 0.3, 
        alpha = 0.05)


# DELE tost ies vs ims (idem)
TOSTER::TOSTtwo(m1 = 31.7, sd1 = 4.58, n1 = 33, # ies
                m2 = 32.8, sd2 = 4.23, n2 = 32, # ims
                low_eqbound_d = -0.3, 
                high_eqbound_d = 0.3, 
                alpha = 0.05)






###########################################################################

# If some participants need to be removed, filter out and repeat TOST

dem_rm <- dem_all %>%
  filter(., group %in% c("lal", "int")) %>%
  filter(DELE > 39 & ID != "LAL29" & ID != "LAL15")


dem_rm %>%
  group_by(., group) %>% 
  summarize(., DELE_mean = mean(DELE), 
            DELE_sd = sd(DELE), 
            n = n())

# group DELE_mean DELE_sd     n
# 1 int        49.0    4.20    24
# 2 lal        46.7    3.94    26


# DELE toast mon vs int (all good now)
TOSTtwo(m1 = 49.0, sd1 = 4.20, n1 = 24, # int
        m2 = 46.7, sd2 = 3.94, n2 = 26, # lal
        low_eqbound_d = -0.3, 
        high_eqbound_d = 0.3, 
        alpha = 0.05)


######################################################################################

# Graph comparing int and lal

dem_all %>%
  filter(., group != "mon") %>%
  ggplot(., aes(x = group, y = DELE)) +
  geom_boxplot()



# NOT RUN BECAUSE I DO NOT HAVE ONLY TWO GROUPS, BUT FOUR
## took from previous script, not sure why we're doing this
dem_all %>%
  filter(., group != "mon") %>%
  group_by(., ID, group) %>%
  summarize(., dele = unique(DELE)) %>%
  mutate(., groupSum = ifelse(group == "int", yes = -0.5, no = 0.5)) %>%
  lm(dele ~ groupSum, data = .) %>%
  summary




# -----------------------------------------------------------------------------
# FROM TEMPLATE SCRIPT: DUR_STRESS PROJECT, BY JOSEPH CASILLAS
# 
# dele_df2 %>%
#   group_by(., participant, Group) %>%
#   summarize(., dele = unique(dele)) %>%
#   mutate(., groupSum = ifelse(Group == "L", yes = -0.5, no = 0.5)) %>%
#   lm(dele ~ groupSum, data = .) %>%
#   summary
# 

# 
# # We administered a revised version of the Diploma de Español como
# # Lengua Extranjera (DELE) test (Eng. Diploma of Spanish as a Foreign
# # Language) to all learner participants. The learners had an
# # average DELE score of 31.18 points; however, the distribution was
# # bimodal. Thus we divided the learners into two groups: late beginning
# # learners (LB) and late advanced (LA). Specifically, the LBs’ scores
# # were 15.08 ± 0.76 standard errors (se) lower than the grand mean
# # (x̅ = 16.10, SD = 3.96). The LA group scores were significantly
# # higher by 30.17 points ± 1.51 se (t = 19.98, p < 0.001, x̅ = 46.27,
# # SD = 4.09).
