#
#
# Growth curve analysis ------------------------------------------------------
#
# - Question 1: Do visuospatial prediction abilities (continuous) influence 
#   Spanish speakers' abilities to predict verbal tense
#   based on the presence or absence (categorical) of lexical stress?
# - Question 2: Do verbal and visuospatial processing speed influence linguistic prediction?
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

# Load data
source(here::here("scripts", "00_load_libs.R"))
# source(here::here("scripts", "02_load_data.R"))
stress50 <- read_csv(here("data", "clean", "stress_50ms_final_onsetc3updated.csv"))


# Get path to saved models
gca_mods_path  <- here("mods", "vision", "gca", "cont_speed_verb")

# Load models as lists
load(paste0(gca_mods_path, "/mon_mods_psequal.Rdata"))
load(paste0(gca_mods_path, "/model_preds_psequal.Rdata"))
# load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))

# Store objects in global env
#list2env(mon_mods, globalenv())
# list2env(nested_model_comparisons, globalenv())
# list2env(model_preds, globalenv())

# -----------------------------------------------------------------------------







# Data prep -------------------------------------------------------------------

# - subset using time course
#    - We need to reduce the time course to a relevant time window that
#      that includes enough of the trajectory from before and after the
#      target syllable onset
#    - Importantly, we need to make sure that the adjusted time course
#      is centered at 200ms after the offset of the first syllable
#    - This is because the orthogonal polynomials center the time course,
#      thus the parameter estimates on the intercept and the linear slope
#      are calculated for the midpoint (0).
#    - This has an added bonus of assessing group differences at the mid
#      point (200ms after target syllable offset), which will corroborate
#      the results from the GLMMs.
#    - We can select the appropriate time course subset by selecting the
#      target syllable offset, bin 4 (200ms / 50 = 4), and keeping an
#      equal number of bins on each side:
#                     8 7 6 5 4 3 2 1 X 1 2 3 4 5 6 7 8
#                                     ^
#                     center of time course (bin 4)
#
#
# Number of bins:     1  2  3  4 5 6 7 8 9 10 11 12 13 14 15 16 17
# Actual bin number: -4 -3 -2 -1 0 1 2 3 4  5  6  7  8  9 10 11 12

wm <- read_csv("./data/clean/wm_processing_speed.csv")
vision <- read_csv("./data/clean/vision_scores_nooutliers.csv") # pred car
corsi <- read_csv("./data/clean/corsi_z_scores.csv")
wm_score <- read_csv("./data/clean/ospan_set_z_scores.csv")

vision50 <- left_join(x = stress50, y = wm, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = vision, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = corsi, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = wm_score, by = "participant", all.x=TRUE)

vision50 <- na.omit(vision50)




vision50 <- vision50 %>%
  filter(., l1 == 'es' & time_zero >= 0 & time_zero <= 8) %>%
  mutate(., #l1 = fct_relevel(l1, "es", "en", "ma"),
            stress_sum = if_else(cond == "1", -1, 1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

mon_vision <- filter(vision50, l1 == 'es') %>% select(-DELE, -percent_l2_week, 
                                                    -DELE_z, -use_z, -prof, -group)


# -----------------------------------------------------------------------------







# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + stress_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = mon_vision, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + stress_sum + ot1 | participant) +
             ot2 + (1 + stress_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + stress_sum + ot1 + ot2 | participant) +
             ot3 + (1 + stress_sum + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar   AIC   BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    9 14011 14064 -6996.6    13993                        
  # mod_ot2   14 14017 14099 -6994.6    13989  3.9025  5    0.56354  
  # mod_ot3   20 14014 14131 -6987.0    13974 15.3054  6    0.01801 *
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df   AIC   BIC  logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot3   20 14014 14131 -6987.0    13974                          
  # mod_ot4   21 13876 13999 -6916.9    13834 140.062  1  < 2.2e-16 ***
  # mod_ot5   23 13845 13980 -6899.5    13799  34.828  2  2.736e-08 ***
  # mod_ot6   26 13836 13989 -6892.1    13784  14.795  3   0.002001 ** 
  # mod_ot7   30 13828 14004 -6884.1    13768  16.119  4   0.002863 ** 
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_mon_base <- mod_ot7
  # lmer(eLog ~ 1 + ot1 + ot2 + ot3 +         
  #        (1 + stress_sum + ot1 + ot2 + ot3 | participant) +
  #        (1 + ot1 + ot2 + ot3 | target),
  #      control = lmerControl(optimizer = 'bobyqa'), #, optCtrl = list(maxfun = 3e5)
  #      data = stress_gc_subset, REML = F)    # , na.action = na.exclude 

# add stress effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_stress_0 <- update(gca_mon_base,    . ~ . + stress_sum)
gca_mon_stress_1 <- update(gca_mon_stress_0, . ~ . + ot1:stress_sum)
gca_mon_stress_2 <- update(gca_mon_stress_1, . ~ . + ot2:stress_sum)
gca_mon_stress_3 <- update(gca_mon_stress_2, . ~ . + ot3:stress_sum)

mon_stress_anova <-
  anova(gca_mon_base, gca_mon_stress_0, gca_mon_stress_1,
        gca_mon_stress_2, gca_mon_stress_3)
                 # npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base       30 13828 14004 -6884.1    13768                     
# gca_mon_stress_0   31 13830 14011 -6883.7    13768 0.6446  1     0.4220
# gca_mon_stress_1   32 13831 14019 -6883.6    13767 0.3065  1     0.5799
# gca_mon_stress_2   33 13832 14026 -6883.1    13766 1.0385  1     0.3082
# gca_mon_stress_3   34 13832 14031 -6881.8    13764 2.4517  1     0.1174


# add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_0 <- update(gca_mon_base,    . ~ . + corsi_rt)
gca_mon_corsirt_1 <- update(gca_mon_corsirt_0, . ~ . + ot1:corsi_rt)
gca_mon_corsirt_2 <- update(gca_mon_corsirt_1, . ~ . + ot2:corsi_rt)
gca_mon_corsirt_3 <- update(gca_mon_corsirt_2, . ~ . + ot3:corsi_rt)

mon_corsirt_anova <-
  anova(gca_mon_base, gca_mon_corsirt_0, gca_mon_corsirt_1,
        gca_mon_corsirt_2, gca_mon_corsirt_3)
#                   npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)
# gca_mon_base        30 13828 14004 -6884.1    13768                       
# gca_mon_corsirt_0   31 13827 14009 -6882.7    13765 2.7769  1    0.09563 .
# gca_mon_corsirt_1   32 13829 14017 -6882.6    13765 0.0633  1    0.80140  
# gca_mon_corsirt_2   33 13830 14023 -6881.8    13764 1.6762  1    0.19543  
# gca_mon_corsirt_3   34 13831 14030 -6881.6    13763 0.3912  1    0.53168

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_int_0 <- update(gca_mon_base, . ~ . + stress_sum:corsi_rt)
gca_mon_corsirt_int_1 <- update(gca_mon_corsirt_int_0,   . ~ . + ot1:stress_sum:corsi_rt)
gca_mon_corsirt_int_2 <- update(gca_mon_corsirt_int_1,   . ~ . + ot2:stress_sum:corsi_rt)
gca_mon_corsirt_int_3 <- update(gca_mon_corsirt_int_2,   . ~ . + ot3:stress_sum:corsi_rt)

mon_corsirt_int_anova <-
  anova(gca_mon_base, gca_mon_corsirt_int_0, gca_mon_corsirt_int_1,
        gca_mon_corsirt_int_2, gca_mon_corsirt_int_3)
#                       npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_base            30 13828 14004 -6884.1    13768                       
# gca_mon_corsirt_int_0   31 13830 14012 -6884.1    13768 0.0109  1    0.91683  
# gca_mon_corsirt_int_1   32 13830 14017 -6882.9    13766 2.2893  1    0.13027  
# gca_mon_corsirt_int_2   33 13828 14022 -6881.1    13762 3.6148  1    0.05727 .
# gca_mon_corsirt_int_3   34 13830 14029 -6881.1    13762 0.1195  1    0.72962



# BRANCH #1
# add verbal proc time effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_0 <- update(gca_mon_base,    . ~ . + ospan_rt)
gca_mon_ospan_1 <- update(gca_mon_ospan_0, . ~ . + ot1:ospan_rt)
gca_mon_ospan_2 <- update(gca_mon_ospan_1, . ~ . + ot2:ospan_rt)
gca_mon_ospan_3 <- update(gca_mon_ospan_2, . ~ . + ot3:ospan_rt)

mon_ospan_anova <-
  anova(gca_mon_base, gca_mon_ospan_0, gca_mon_ospan_1,
        gca_mon_ospan_2, gca_mon_ospan_3)
#                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      30 13828 14004 -6884.1    13768                     
# gca_mon_ospan_0   31 13830 14011 -6883.8    13768 0.5903  1     0.4423
# gca_mon_ospan_1   32 13832 14019 -6883.8    13768 0.0104  1     0.9186
# gca_mon_ospan_2   33 13833 14027 -6883.6    13767 0.3195  1     0.5719
# gca_mon_ospan_3   34 13835 14034 -6883.6    13767 0.1147  1     0.7348

gca_mon_ospan_i_0 <- update(gca_mon_base, . ~ . + stress_sum:ospan_rt)
gca_mon_ospan_i_1 <- update(gca_mon_ospan_i_0,   . ~ . + ot1:stress_sum:ospan_rt)
gca_mon_ospan_i_2 <- update(gca_mon_ospan_i_1,   . ~ . + ot2:stress_sum:ospan_rt)
gca_mon_ospan_i_3 <- update(gca_mon_ospan_i_2,   . ~ . + ot3:stress_sum:ospan_rt)

mon_ospan_i_anova <-
  anova(gca_mon_base, gca_mon_ospan_i_0, gca_mon_ospan_i_1,
        gca_mon_ospan_i_2, gca_mon_ospan_i_3)
#                   npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_base        30 13828 14004 -6884.1    13768                       
# gca_mon_ospan_i_0   31 13830 14011 -6883.7    13768 0.6684  1     0.4136  
# gca_mon_ospan_i_1   32 13830 14017 -6882.9    13766 1.6687  1     0.1964  
# gca_mon_ospan_i_2   33 13827 14020 -6880.5    13761 4.8970  1     0.0269 *
# gca_mon_ospan_i_3   34 13828 14027 -6879.9    13760 1.0249  1     0.3114


gca_mon_ospan_in_0 <- update(gca_mon_ospan_i_2,    . ~ . + ospan_rt:corsi_rt)
gca_mon_ospan_in_1 <- update(gca_mon_ospan_in_0,   . ~ . + ot1:ospan_rt:corsi_rt)
gca_mon_ospan_in_2 <- update(gca_mon_ospan_in_1,   . ~ . + ot2:ospan_rt:corsi_rt)
gca_mon_ospan_in_3 <- update(gca_mon_ospan_in_2,   . ~ . + ot3:ospan_rt:corsi_rt)

mon_ospan_in_anova <-
  anova(gca_mon_ospan_i_2, gca_mon_ospan_in_0, gca_mon_ospan_in_1,
        gca_mon_ospan_in_2, gca_mon_ospan_in_3)
#                    npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_ospan_i_2    33 13827 14020 -6880.5    13761                       
# gca_mon_ospan_in_0   34 13826 14025 -6878.8    13758 3.2249  1    0.07252 .
# gca_mon_ospan_in_1   35 13828 14033 -6878.8    13758 0.0136  1    0.90708  
# gca_mon_ospan_in_2   36 13830 14040 -6878.8    13758 0.1430  1    0.70531  
# gca_mon_ospan_in_3   37 13831 14048 -6878.7    13757 0.2027  1    0.65252


# add 3-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_int_0 <- update(gca_mon_ospan_i_2, . ~ . + stress_sum:ospan_rt:corsi_rt)
gca_mon_ospan_int_1 <- update(gca_mon_ospan_int_0,   . ~ . + ot1:stress_sum:ospan_rt:corsi_rt)
gca_mon_ospan_int_2 <- update(gca_mon_ospan_int_1,   . ~ . + ot2:stress_sum:ospan_rt:corsi_rt)
gca_mon_ospan_int_3 <- update(gca_mon_ospan_int_2,   . ~ . + ot3:stress_sum:ospan_rt:corsi_rt)


mon_ospan_int_anova <-
  anova(gca_mon_ospan_i_2, gca_mon_ospan_int_0, gca_mon_ospan_int_1,
        gca_mon_ospan_int_2, gca_mon_ospan_int_3)
#                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mon_ospan_i_2     33 13827 14020 -6880.5    13761                     
# gca_mon_ospan_int_0   34 13827 14026 -6879.3    13759 2.2674  1     0.1321
# gca_mon_ospan_int_1   35 13828 14034 -6879.2    13758 0.1555  1     0.6933
# gca_mon_ospan_int_2   36 13830 14041 -6878.9    13758 0.6037  1     0.4372
# gca_mon_ospan_int_3   37 13832 14049 -6878.9    13758 0.0781  1     0.7799


# BRANCH #2
# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_0 <- update(gca_mon_base,    . ~ . + car_dev)
gca_mon_car_1 <- update(gca_mon_car_0, . ~ . + ot1:car_dev)
gca_mon_car_2 <- update(gca_mon_car_1, . ~ . + ot2:car_dev)
gca_mon_car_3 <- update(gca_mon_car_2, . ~ . + ot3:car_dev)

mon_car_anova <-
  anova(gca_mon_base, gca_mon_car_0, gca_mon_car_1,
        gca_mon_car_2, gca_mon_car_3)
#               npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base    30 13828 14004 -6884.1    13768                     
# gca_mon_car_0   31 13830 14012 -6884.0    13768 0.1009  1     0.7507
# gca_mon_car_1   32 13831 14019 -6883.6    13767 0.9106  1     0.3400
# gca_mon_car_2   33 13832 14025 -6882.9    13766 1.3455  1     0.2461
# gca_mon_car_3   34 13834 14033 -6882.8    13766 0.1058  1     0.7450

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_i_0 <- update(gca_mon_base,    . ~ . + car_dev:stress_sum)
gca_mon_car_i_1 <- update(gca_mon_car_i_0, . ~ . + ot1:car_dev:stress_sum)
gca_mon_car_i_2 <- update(gca_mon_car_i_1, . ~ . + ot2:car_dev:stress_sum)
gca_mon_car_i_3 <- update(gca_mon_car_i_2, . ~ . + ot3:car_dev:stress_sum)

mon_car_i_anova <-
  anova(gca_mon_base, gca_mon_car_i_0, gca_mon_car_i_1,
        gca_mon_car_i_2, gca_mon_car_i_3)
#                 npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_base      30 13828 14004 -6884.1    13768                       
# gca_mon_car_i_0   31 13830 14012 -6884.1    13768 0.0340  1    0.85365  
# gca_mon_car_i_1   32 13829 14016 -6882.5    13765 3.1286  1    0.07693 .
# gca_mon_car_i_2   33 13831 14024 -6882.3    13765 0.4144  1    0.51975  
# gca_mon_car_i_3   34 13831 14030 -6881.6    13763 1.3208  1    0.25044


gca_mon_car_in_0 <- update(gca_mon_base,     . ~ . + car_dev:corsi_rt)
gca_mon_car_in_1 <- update(gca_mon_car_in_0, . ~ . + ot1:car_dev:corsi_rt)
gca_mon_car_in_2 <- update(gca_mon_car_in_1, . ~ . + ot2:car_dev:corsi_rt)
gca_mon_car_in_3 <- update(gca_mon_car_in_2, . ~ . + ot3:car_dev:corsi_rt)

mon_car_in_anova <-
  anova(gca_mon_base, gca_mon_car_in_0, gca_mon_car_in_1,
        gca_mon_car_in_2, gca_mon_car_in_3)
#                  npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base       30 13828 14004 -6884.1    13768                     
# gca_mon_car_in_0   31 13830 14012 -6884.1    13768 0.0310  1     0.8603
# gca_mon_car_in_1   32 13831 14019 -6883.7    13767 0.6624  1     0.4157
# gca_mon_car_in_2   33 13832 14025 -6882.7    13766 1.9848  1     0.1589
# gca_mon_car_in_3   34 13833 14032 -6882.4    13765 0.6640  1     0.4151


gca_mon_car_int_0 <- update(gca_mon_base, . ~ . + stress_sum:car_dev:corsi_rt)
gca_mon_car_int_1 <- update(gca_mon_car_int_0,   . ~ . + ot1:stress_sum:car_dev:corsi_rt)
gca_mon_car_int_2 <- update(gca_mon_car_int_1,   . ~ . + ot2:stress_sum:car_dev:corsi_rt)
gca_mon_car_int_3 <- update(gca_mon_car_int_2,   . ~ . + ot3:stress_sum:car_dev:corsi_rt)

mon_car_int_anova <-
  anova(gca_mon_base, gca_mon_car_int_0, gca_mon_car_int_1,
        gca_mon_car_int_2, gca_mon_car_int_3)
#                   npar   AIC   BIC  logLik deviance   Chisq Df Pr(>Chisq)  
# gca_mon_base        30 13828 14004 -6884.1    13768                     
# gca_mon_car_int_0   31 13830 14012 -6884.1    13768 0.0044  1     0.9469
# gca_mon_car_int_1   32 13830 14017 -6882.9    13766 2.3082  1     0.1287
# gca_mon_car_int_2   33 13832 14025 -6882.9    13766 0.1172  1     0.7320
# gca_mon_car_int_3   34 13834 14033 -6882.8    13766 0.0145  1     0.9041


car::vif(gca_mon_ospan_i_2)
# ot1                     ot2                     ot3 
# 1.491246                1.015589                1.505883 
# stress_sum:ospan_rt     ot1:stress_sum:ospan_rt ot2:stress_sum:ospan_rt 
# 1.023016                1.032439                1.015460 

}

# -----------------------------------------------------------------------------


# save models  
mod_type <- "gca_mon"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", 
              "_corsirt_0", "_corsirt_1", "_corsirt_2", 
              "_corsirt_int_0", "_corsirt_int_1", "_corsirt_int_2", 
              "_ospan_0", "_ospan_1", "_ospan_2", 
              "_ospan_i_0", "_ospan_i_1", "_ospan_i_2", 
              "_ospan_in_0", "_ospan_in_1", "_ospan_in_2", 
              "_ospan_int_0", "_ospan_int_1", "_ospan_int_2", 
              "_car_0", "_car_1", "_car_2", 
              "_car_i_0", "_car_i_1", "_car_i_2",
              "_car_in_0", "_car_in_1", "_car_in_2",
              "_car_int_0", "_car_int_1", "_car_int_2"
              )

# Store ind models in list
mon_mods_psequal <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods_psequal,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "mon_mods_psequal.Rdata"))

  


  
# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_newdata <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(corsi_rt = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan_rt = c(-1, 0, 1)))

# Get model predictions and SE
fits_mon_ospan <- predictSE(gca_mon_ospan_i_2, mon_newdata) %>%        
  as_tibble %>%
  bind_cols(mon_newdata) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_mon_ospan <- filter(fits_mon_ospan, time_zero == 4) %>%
  select(stress = stress_sum, ospan_rt, corsi_rt,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 


# Save models predictions
model_preds_psequal <- mget(c('fits_mon_ospan', 
                      "preds_mon_ospan"
                      ))

save(model_preds_psequal,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "model_preds_psequal.Rdata"))



# -----------------------------------------------------------------------------






# Save models -----------------------------------------------------------------

if(F) {
  
  # Save anova model comparisons
  # nested_model_comparisons <-
  #   mget(c("mon_stress_anova", "mon_car_anova", "mon_corsi_anova",
  #          "mon_car_int_anova", "mon_corsi_int_anova"
  #   ))
  # 
  # save(nested_model_comparisons,
  #      file = here("mods", "vision", "gca", "continuous",
  #                  "nested_model_comparisons.Rdata"))
  
  
  
  
  
  
}

# -----------------------------------------------------------------------------

