#
# Original script by Joseph
# Modified by Cristina to adapt to CrisLau Project
# Last update: 06/12/2019
# Modified by Laura to adapt to Pupurri project
# Last update: 05/08/2021
#
# Growth curve analysis ------------------------------------------------------
#
# - Question 1: Do visuospatial prediction abilities (continuous) influence 
#   English and Mandarin Chinese speakers' abilities to predict verbal tense
#   based on the presence or absence (categorical) of lexical stress 
#   at different levels of proficiency (continuous, from low intermediate 
#   to very advanced)?
#     - Hypothesis 1: No association in any population
# - Question 2: Do verbal and visuospatial WM influence linguistic prediction?
#     - Hypothesis 2: Not really.
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

# Load data
source(here::here("scripts", "02_load_data.R"))

# Get path to saved models
gca_mods_path  <- here("mods", "vision", "gca", "continuous_speed", "verb")

# Load models as lists
load(paste0(gca_mods_path, "/mon_mods.Rdata"))
load(paste0(gca_mods_path, "/en_mods.Rdata"))
load(paste0(gca_mods_path, "/ma_mods.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))

# Store objects in global env
list2env(mon_mods, globalenv())
list2env(en_mods, globalenv())
list2env(ma_mods, globalenv())
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
vision <- read_csv("./data/clean/vision_scores.csv")

vision50 <- left_join(x = stress50, y = wm, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = vision, by = "participant", all.x=TRUE)

vision50 <- na.omit(vision50)




vision50 <- vision50 %>%
  filter(., time_zero >= -4 & time_zero <= 8) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma"),
            stress_sum = if_else(cond == "1", -1, 1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

mon_vision <- filter(vision50, l1 == 'es') %>% select(-DELE, -percent_l2_week, 
                                                    -DELE_z, -use_z, -prof, -group)

en_vision <- vision50 %>%
  filter(., l1 == 'en') %>% 
  mutate(., 
         prof_std = (DELE - mean(DELE))/sd(DELE),
         use_std = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week)) %>%
  select(-DELE, -percent_l2_week, 
         -DELE_z, -use_z, -prof, -group)

ma_vision <- vision50 %>%
  filter(., l1 == 'ma') %>% 
  mutate(., 
         prof_std = (DELE - mean(DELE))/sd(DELE),
         use_std = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week)) %>%
  select(-DELE, -percent_l2_week, 
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
  #         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    9 31888 31948 -15935    31870                         
  # mod_ot2   14 31874 31967 -15923    31846 23.983  5  0.0002187 ***
  # mod_ot3   20 31867 32000 -15914    31827 19.151  6  0.0039154 ** 
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df   AIC   BIC logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot3   20 31867 32000 -15914    31827                          
  # mod_ot4   21 31743 31883 -15850    31701 126.420  1  < 2.2e-16 ***
  # mod_ot5   23 31690 31843 -15822    31644  56.489  2  5.415e-13 ***
  # mod_ot6   26 31676 31849 -15812    31624  20.036  3  0.0001669 ***
  # mod_ot7   30 31672 31872 -15806    31612  12.188  4  0.0160043 * 
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_mon_base <- mod_ot7
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +          
  #        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 + ot3 | target),
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 3e5)),
  #      data = stress_gc_subset, REML = F)    # , na.action = na.exclude 

# add stress effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_stress_0 <- update(gca_mon_base,    . ~ . + stress_sum)
gca_mon_stress_1 <- update(gca_mon_stress_0, . ~ . + ot1:stress_sum)
gca_mon_stress_2 <- update(gca_mon_stress_1, . ~ . + ot2:stress_sum)
gca_mon_stress_3 <- update(gca_mon_stress_2, . ~ . + ot3:stress_sum)

mon_stress_anova <-
  anova(gca_mon_base, gca_mon_stress_0, gca_mon_stress_1,
        gca_mon_stress_2, gca_mon_stress_3)
                 # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base       30 31672 31872 -15806    31612                     
# gca_mon_stress_0   31 31673 31880 -15806    31611 0.6879  1     0.4069
# gca_mon_stress_1   32 31673 31886 -15804    31609 2.3344  1     0.1265
# gca_mon_stress_2   33 31673 31893 -15804    31608 1.3054  1     0.2532
# gca_mon_stress_3   34 31675 31902 -15804    31607 0.3628  1     0.5470

# BRANCH #0
# add verbal proc time effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_0 <- update(gca_mon_base,    . ~ . + ospan_rt)
gca_mon_ospan_1 <- update(gca_mon_ospan_0, . ~ . + ot1:ospan_rt)
gca_mon_ospan_2 <- update(gca_mon_ospan_1, . ~ . + ot2:ospan_rt)
gca_mon_ospan_3 <- update(gca_mon_ospan_2, . ~ . + ot3:ospan_rt)

mon_ospan_anova <-
  anova(gca_mon_base, gca_mon_ospan_0, gca_mon_ospan_1,
        gca_mon_ospan_2, gca_mon_ospan_3)
#                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      30 31672 31872 -15806    31612                     
# gca_mon_ospan_0   31 31674 31880 -15806    31612 0.1533  1     0.6954
# gca_mon_ospan_1   32 31676 31889 -15806    31612 0.0311  1     0.8599
# gca_mon_ospan_2   33 31677 31897 -15806    31612 0.1430  1     0.7053
# gca_mon_ospan_3   34 31679 31905 -15805    31611 0.8029  1     0.3702


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_int_0 <- update(gca_mon_base, . ~ . + stress_sum:ospan_rt)
gca_mon_ospan_int_1 <- update(gca_mon_ospan_int_0,   . ~ . + ot1:stress_sum:ospan_rt)
gca_mon_ospan_int_2 <- update(gca_mon_ospan_int_1,   . ~ . + ot2:stress_sum:ospan_rt)
gca_mon_ospan_int_3 <- update(gca_mon_ospan_int_2,   . ~ . + ot3:stress_sum:ospan_rt)

mon_ospan_int_anova <-
  anova(gca_mon_base, gca_mon_ospan_int_0, gca_mon_ospan_int_1,
        gca_mon_ospan_int_2, gca_mon_ospan_int_3)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mon_base          30 1672 31872 -15806    31612                          
# gca_mon_ospan_int_0   31 31674 31880 -15806    31612  0.0227  1    0.88013    
# gca_mon_ospan_int_1   32 31671 31884 -15803    31607  5.0185  1    0.02508 *  
# gca_mon_ospan_int_2   33 31673 31893 -15803    31607  0.0707  1    0.79026    
# gca_mon_ospan_int_3   34 31654 31881 -15793    31586 20.6611  1  5.482e-06 ***


# BRANCH #1
# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_0 <- update(gca_mon_base,    . ~ . + car_dev)
gca_mon_car_1 <- update(gca_mon_car_0, . ~ . + ot1:car_dev)
gca_mon_car_2 <- update(gca_mon_car_1, . ~ . + ot2:car_dev)
gca_mon_car_3 <- update(gca_mon_car_2, . ~ . + ot3:car_dev)

mon_car_anova <-
  anova(gca_mon_base, gca_mon_car_0, gca_mon_car_1,
        gca_mon_car_2, gca_mon_car_3)
#               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base    30 31672 31872 -15806    31612                     
# gca_mon_car_0   31 31674 31880 -15806    31612 0.0038  1     0.9509
# gca_mon_car_1   32 31673 31887 -15805    31609 2.6477  1     0.1037
# gca_mon_car_2   33 31674 31894 -15804    31608 1.1043  1     0.2933
# gca_mon_car_3   34 31676 31903 -15804    31608 0.1512  1     0.6973    


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_int_0 <- update(gca_mon_base, . ~ . + stress_sum:car_dev)
gca_mon_car_int_1 <- update(gca_mon_car_int_0,   . ~ . + ot1:stress_sum:car_dev)
gca_mon_car_int_2 <- update(gca_mon_car_int_1,   . ~ . + ot2:stress_sum:car_dev)
gca_mon_car_int_3 <- update(gca_mon_car_int_2,   . ~ . + ot3:stress_sum:car_dev)

mon_car_int_anova <-
  anova(gca_mon_base, gca_mon_car_int_0, gca_mon_car_int_1,
        gca_mon_car_int_2, gca_mon_car_int_3)
#                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_base        30 31672 31872 -15806    31612                       
# gca_mon_car_int_0   31 31674 31880 -15806    31612 0.2137  1    0.64386  
# gca_mon_car_int_1   32 31671 31884 -15803    31607 4.7574  1    0.02917 *
# gca_mon_car_int_2   33 31673 31893 -15803    31607 0.1059  1    0.74491  
# gca_mon_car_int_3   34 31672 31899 -15802    31604 2.7405  1    0.09783 .       

# # BRANCH #2
# # add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
# gca_mon_corsi_0 <- update(gca_mon_base,    . ~ . + corsi)
# gca_mon_corsi_1 <- update(gca_mon_corsi_0, . ~ . + ot1:corsi)
# gca_mon_corsi_2 <- update(gca_mon_corsi_1, . ~ . + ot2:corsi)
# gca_mon_corsi_3 <- update(gca_mon_corsi_2, . ~ . + ot3:corsi)
# 
# mon_corsi_anova <-
#   anova(gca_mon_base, gca_mon_corsi_0, gca_mon_corsi_1,
#         gca_mon_corsi_2, gca_mon_corsi_3)
# #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# # gca_mon_base      
# 
# 
# # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
# gca_mon_corsi_int_0 <- update(gca_mon_base, . ~ . + stress_sum:corsi)
# gca_mon_corsi_int_1 <- update(gca_mon_corsi_int_0,   . ~ . + ot1:stress_sum:corsi)
# gca_mon_corsi_int_2 <- update(gca_mon_corsi_int_1,   . ~ . + ot2:stress_sum:corsi)
# gca_mon_corsi_int_3 <- update(gca_mon_corsi_int_2,   . ~ . + ot3:stress_sum:corsi)
# 
# mon_corsi_int_anova <-
#   anova(gca_mon_base, gca_mon_corsi_int_0, gca_mon_corsi_int_1,
#         gca_mon_corsi_int_2, gca_mon_corsi_int_3)
# #                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
# # gca_mon_base    


# BRANCH #3
# add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_0 <- update(gca_mon_base,    . ~ . + corsi_rt)
gca_mon_corsirt_1 <- update(gca_mon_corsirt_0, . ~ . + ot1:corsi_rt)
gca_mon_corsirt_2 <- update(gca_mon_corsirt_1, . ~ . + ot2:corsi_rt)
gca_mon_corsirt_3 <- update(gca_mon_corsirt_2, . ~ . + ot3:corsi_rt)

mon_corsirt_anova <-
  anova(gca_mon_base, gca_mon_corsirt_0, gca_mon_corsirt_1,
        gca_mon_corsirt_2, gca_mon_corsirt_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base        30 31672 31872 -15806    31612                     
# gca_mon_corsirt_0   31 31673 31880 -15806    31611 0.5603  1     0.4541
# gca_mon_corsirt_1   32 31673 31886 -15804    31609 2.2736  1     0.1316
# gca_mon_corsirt_2   33 31675 31895 -15804    31609 0.0479  1     0.8268
# gca_mon_corsirt_3   34 31675 31902 -15804    31608 1.4482  1     0.2288

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_int_0 <- update(gca_mon_base, . ~ . + stress_sum:corsi_rt)
gca_mon_corsirt_int_1 <- update(gca_mon_corsirt_int_0,   . ~ . + ot1:stress_sum:corsi_rt)
gca_mon_corsirt_int_2 <- update(gca_mon_corsirt_int_1,   . ~ . + ot2:stress_sum:corsi_rt)
gca_mon_corsirt_int_3 <- update(gca_mon_corsirt_int_2,   . ~ . + ot3:stress_sum:corsi_rt)

mon_corsirt_int_anova <-
  anova(gca_mon_base, gca_mon_corsirt_int_0, gca_mon_corsirt_int_1,
        gca_mon_corsirt_int_2, gca_mon_corsirt_int_3)
#                       npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_base            30 31672 31872 -15806    31612                       
# gca_mon_corsirt_int_0   31 31674 31880 -15806    31612 0.0479  1    0.82670  
# gca_mon_corsirt_int_1   32 31675 31889 -15806    31612 0.2606  1    0.60968  
# gca_mon_corsirt_int_2   33 31673 31893 -15803    31606 4.9447  1    0.02617 *
# gca_mon_corsirt_int_3   34 31669 31895 -15800    31601 5.9259  1    0.01492 *


summary(gca_mon_ospan_int_3)
# Estimate Std. Error t value
# (Intercept)   

}

# -----------------------------------------------------------------------------









# Random effects structure ----------------------------------------------------




# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + stress_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = en_vision, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + stress_sum + ot1 | participant) +
             ot2 + (1 + stress_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + stress_sum + ot1 + ot2 | participant) +
             ot3 + (1 + stress_sum + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    9 70261 70328 -35122    70243                          
  # mod_ot2   14 70062 70166 -35017    70034 209.472  5  < 2.2e-16 ***
  # mod_ot3   20 70044 70192 -35002    70004  30.231  6  3.553e-05 ***
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df    AIC    BIC  logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot3   20 70044 70192 -35002    70004                           
  # mod_ot4   21 69878 70035 -34918    69836 167.1765  1  < 2.2e-16 ***
  #   mod_ot5   23 69804 69975 -34879    69758  78.1364  2  < 2.2e-16 ***
  #   mod_ot6   26 69796 69989 -34872    69744  14.4768  3   0.002323 ** 
  #   mod_ot7   30 69797 70020 -34869    69737   6.7332  4   0.150677    
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_en_base <- mod_ot6
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
  #        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 | target),
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 3e5)),
  #      data = en_vision, REML = F)    # , na.action = na.exclude
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_prof_0 <- update(gca_en_base,    . ~ . + prof_std)
  gca_en_prof_1 <- update(gca_en_prof_0, . ~ . + ot1:prof_std)
  gca_en_prof_2 <- update(gca_en_prof_1, . ~ . + ot2:prof_std)
  gca_en_prof_3 <- update(gca_en_prof_2, . ~ . + ot3:prof_std)
  
  en_prof_anova <-
    anova(gca_en_base, gca_en_prof_0, gca_en_prof_1,
          gca_en_prof_2, gca_en_prof_3)
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_en_base     26 69796 69989 -34872    69744                          
  # gca_en_prof_0   27 69798 69999 -34872    69744  0.2512  1  0.6162287    
  # gca_en_prof_1   28 69788 69997 -34866    69732 11.3823  1  0.0007415 ***
  # gca_en_prof_2   29 69790 70006 -34866    69732  0.0001  1  0.9914034    
  # gca_en_prof_3   30 69792 70015 -34866    69732  0.1683  1  0.6816610 
  
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_stress_0 <- update(gca_en_prof_1,    . ~ . + stress_sum)
  gca_en_stress_1 <- update(gca_en_stress_0, . ~ . + ot1:stress_sum)
  gca_en_stress_2 <- update(gca_en_stress_1, . ~ . + ot2:stress_sum)
  gca_en_stress_3 <- update(gca_en_stress_2, . ~ . + ot3:stress_sum)
  
  en_stress_anova <-
    anova(gca_en_prof_1, gca_en_stress_0, gca_en_stress_1,
          gca_en_stress_2, gca_en_stress_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_prof_1     28 69788 69997 -34866    69732                       
  # gca_en_stress_0   29 69789 70005 -34865    69731 1.3458  1    0.24602  
  # gca_en_stress_1   30 69791 70014 -34865    69731 0.1123  1    0.73758  
  # gca_en_stress_2   31 69788 70019 -34863    69726 4.2255  1    0.03982 *
  # gca_en_stress_3   32 69790 70029 -34863    69726 0.0388  1    0.84388 

  
  
  # BRANCH #1
  # add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_car_0 <- update(gca_en_stress_2,    . ~ . + car_dev)
  gca_en_car_1 <- update(gca_en_car_0, . ~ . + ot1:car_dev)
  gca_en_car_2 <- update(gca_en_car_1, . ~ . + ot2:car_dev)
  gca_en_car_3 <- update(gca_en_car_2, . ~ . + ot3:car_dev)
  
  en_car_anova <-
    anova(gca_en_stress_2, gca_en_car_0, gca_en_car_1,
          gca_en_car_2, gca_en_car_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stress_2   31 69788 70019 -34863    69726                     
  # gca_en_car_0      32 69790 70028 -34863    69726 0.7424  1     0.3889
  # gca_en_car_1      33 69792 70037 -34863    69726 0.0558  1     0.8133
  # gca_en_car_2      34 69793 70046 -34862    69725 0.8145  1     0.3668
  # gca_en_car_3      35 69792 70053 -34861    69722 2.6240  1     0.1053   
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_car_int_0 <- update(gca_en_stress_2, . ~ . + prof_std:stress_sum:car_dev)
  gca_en_car_int_1 <- update(gca_en_car_int_0,   . ~ . + ot1:prof_std:stress_sum:car_dev)
  gca_en_car_int_2 <- update(gca_en_car_int_1,   . ~ . + ot2:prof_std:stress_sum:car_dev)
  gca_en_car_int_3 <- update(gca_en_car_int_2,   . ~ . + ot3:prof_std:stress_sum:car_dev)
  
  en_car_int_anova <-
    anova(gca_en_stress_2, gca_en_car_int_0, gca_en_car_int_1,
          gca_en_car_int_2, gca_en_car_int_3)
  #                 npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_en_stress_2    31 69788 70019 -34863    69726                        
  # gca_en_car_int_0   32 69789 70027 -34863    69725 1.3406  1   0.246931   
  # gca_en_car_int_1   33 69784 70030 -34859    69718 6.9896  1   0.008199 **
  # gca_en_car_int_2   34 69782 70035 -34857    69714 4.2484  1   0.039287 * 
  # gca_en_car_int_3   35 69780 70041 -34855    69710 3.6673  1   0.055491 .       
  
  
  # # BRANCH #2
  # # add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
  # gca_en_corsi_0 <- update(gca_en_stress_2,    . ~ . + corsi)
  # gca_en_corsi_1 <- update(gca_en_corsi_0, . ~ . + ot1:corsi)
  # gca_en_corsi_2 <- update(gca_en_corsi_1, . ~ . + ot2:corsi)
  # gca_en_corsi_3 <- update(gca_en_corsi_2, . ~ . + ot3:corsi)
  # 
  # en_corsi_anova <-
  #   anova(gca_en_stress_2, gca_en_corsi_0, gca_en_corsi_1,
  #         gca_en_corsi_2, gca_en_corsi_3)
  # #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # # gca_en_stress_2   
  # 
  # # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  # gca_en_corsi_int_0 <- update(gca_en_stress_2, . ~ . + prof_std:stress_sum:corsi)
  # gca_en_corsi_int_1 <- update(gca_en_corsi_int_0,   . ~ . + ot1:prof_std:stress_sum:corsi)
  # gca_en_corsi_int_2 <- update(gca_en_corsi_int_1,   . ~ . + ot2:prof_std:stress_sum:corsi)
  # gca_en_corsi_int_3 <- update(gca_en_corsi_int_2,   . ~ . + ot3:prof_std:stress_sum:corsi)
  # 
  # en_corsi_int_anova <-
  #   anova(gca_en_stress_2, gca_en_corsi_int_0, gca_en_corsi_int_1,
  #         gca_en_corsi_int_2, gca_en_corsi_int_3)
  # #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq) 
  # # gca_en_stress_2      
  
  
  # BRANCH #3
  # add verbal processing speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_ospan_0 <- update(gca_en_stress_2,    . ~ . + ospan_rt)
  gca_en_ospan_1 <- update(gca_en_ospan_0, . ~ . + ot1:ospan_rt)
  gca_en_ospan_2 <- update(gca_en_ospan_1, . ~ . + ot2:ospan_rt)
  gca_en_ospan_3 <- update(gca_en_ospan_2, . ~ . + ot3:ospan_rt)
  
  en_ospan_anova <-
    anova(gca_en_stress_2, gca_en_ospan_0, gca_en_ospan_1,
          gca_en_ospan_2, gca_en_ospan_3)
  #                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stress_2    31 69788 70019 -34863    69726                       
  # gca_en_ospan_0    32 69790 70028 -34863    69726 0.2728  1    0.60149  
  # gca_en_ospan_1    33 69788 70034 -34861    69722 4.0304  1    0.04469 *
  # gca_en_ospan_2    34 69790 70043 -34861    69722 0.0221  1    0.88194  
  # gca_en_ospan_3    35 69791 70051 -34860    69721 1.2523  1    0.26311
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_ospan_int_0 <- update(gca_en_ospan_1, . ~ . + stress_sum:prof_std:ospan_rt)
  gca_en_ospan_int_1 <- update(gca_en_ospan_int_0,   . ~ . + ot1:stress_sum:prof_std:ospan_rt)
  gca_en_ospan_int_2 <- update(gca_en_ospan_int_1,   . ~ . + ot2:stress_sum:prof_std:ospan_rt)
  gca_en_ospan_int_3 <- update(gca_en_ospan_int_2,   . ~ . + ot3:stress_sum:prof_std:ospan_rt)
  
  en_ospan_int_anova <-
    anova(gca_en_ospan_1, gca_en_ospan_int_0, gca_en_ospan_int_1,
          gca_en_ospan_int_2, gca_en_ospan_int_3)
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_ospan_1       33 69788 70034 -34861    69722                     
  # gca_en_ospan_int_0   34 69790 70043 -34861    69722 0.4984  1     0.4802
  # gca_en_ospan_int_1   35 69790 70050 -34860    69720 1.7644  1     0.1841
  # gca_en_ospan_int_2   36 69792 70060 -34860    69720 0.3208  1     0.5712
  # gca_en_ospan_int_3   37 69794 70069 -34860    69720 0.0126  1     0.9107
  
  
  # BRANCH #4
  # add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_corsirt_0 <- update(gca_en_stress_2,    . ~ . + corsi_rt)
  gca_en_corsirt_1 <- update(gca_en_corsirt_0, . ~ . + ot1:corsi_rt)
  gca_en_corsirt_2 <- update(gca_en_corsirt_1, . ~ . + ot2:corsi_rt)
  gca_en_corsirt_3 <- update(gca_en_corsirt_2, . ~ . + ot3:corsi_rt)
  
  en_corsirt_anova <-
    anova(gca_en_stress_2, gca_en_corsirt_0, gca_en_corsirt_1,
          gca_en_corsirt_2, gca_en_corsirt_3)
  #                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stress_2    31 69788 70019 -34863    69726                     
  # gca_en_corsirt_0   32 69790 70028 -34863    69726 0.4374  1     0.5084
  # gca_en_corsirt_1   33 69792 70038 -34863    69726 0.1440  1     0.7044
  # gca_en_corsirt_2   34 69792 70045 -34862    69724 1.6632  1     0.1972
  # gca_en_corsirt_3   35 69794 70055 -34862    69724 0.0954  1     0.7575
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_corsirt_int_0 <- update(gca_en_stress_2, . ~ . + stress_sum:prof_std:corsi_rt)
  gca_en_corsirt_int_1 <- update(gca_en_corsirt_int_0,   . ~ . + ot1:stress_sum:prof_std:corsi_rt)
  gca_en_corsirt_int_2 <- update(gca_en_corsirt_int_1,   . ~ . + ot2:stress_sum:prof_std:corsi_rt)
  gca_en_corsirt_int_3 <- update(gca_en_corsirt_int_2,   . ~ . + ot3:stress_sum:prof_std:corsi_rt)
  
  en_corsirt_int_anova <-
    anova(gca_en_stress_2, gca_en_corsirt_int_0, gca_en_corsirt_int_1,
          gca_en_corsirt_int_2, gca_en_corsirt_int_3)
  #                      npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_stress_2        31 69788 70019 -34863    69726                     
  # gca_en_corsirt_int_0   32 69790 70029 -34863    69726 0.0462  1     0.8298
  # gca_en_corsirt_int_1   33 69792 70037 -34863    69726 0.7399  1     0.3897
  # gca_en_corsirt_int_2   34 69792 70045 -34862    69724 1.7269  1     0.1888
  # gca_en_corsirt_int_3   35 69793 70054 -34862    69723 0.6204  1     0.4309
  
  
  
  summary(gca_en_)
  #               Estimate Std. Error t value
  # (Intercept)    
  


  
  
  
  
  
  
  
  
  # add use effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_use_0 <- update(gca_en_base,    . ~ . + use_std)
  gca_en_use_1 <- update(gca_en_use_0, . ~ . + ot1:use_std)
  gca_en_use_2 <- update(gca_en_use_1, . ~ . + ot2:use_std)
  gca_en_use_3 <- update(gca_en_use_2, . ~ . + ot3:use_std)
  
  en_use_anova <-
    anova(gca_en_base, gca_en_use_0, gca_en_use_1,
          gca_en_use_2, gca_en_use_3)
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_en_base     26 69796 69989 -34872    69744                          
  # gca_en_use_0   27 69797 69998 -34871    69743 0.7643  1    0.38199  
  # gca_en_use_1   28 69796 70004 -34870    69740 3.1829  1    0.07441 .
  # gca_en_use_2   29 69796 70012 -34869    69738 2.0732  1    0.14991  
  # gca_en_use_3   30 69795 70018 -34867    69735 2.9522  1    0.08576 .
  
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_stressu_0 <- update(gca_en_base,    . ~ . + stress_sum)
  gca_en_stressu_1 <- update(gca_en_stressu_0, . ~ . + ot1:stress_sum)
  gca_en_stressu_2 <- update(gca_en_stressu_1, . ~ . + ot2:stress_sum)
  gca_en_stressu_3 <- update(gca_en_stressu_2, . ~ . + ot3:stress_sum)
  
  en_stressu_anova <-
    anova(gca_en_base, gca_en_stressu_0, gca_en_stressu_1,
          gca_en_stressu_2, gca_en_stressu_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_base        26 69796 69989 -34872    69744                       
  # gca_en_stressu_0   27 69796 69997 -34871    69742 1.3545  1    0.24450  
  # gca_en_stressu_1   28 69798 70007 -34871    69742 0.1110  1    0.73904  
  # gca_en_stressu_2   29 69796 70012 -34869    69738 4.1811  1    0.04088 *
  # gca_en_stressu_3   30 69798 70021 -34869    69738 0.0510  1    0.82127  
  
  
  
  # BRANCH #1
  # add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_caru_0 <- update(gca_en_stressu_2,    . ~ . + car_dev)
  gca_en_caru_1 <- update(gca_en_caru_0, . ~ . + ot1:car_dev)
  gca_en_caru_2 <- update(gca_en_caru_1, . ~ . + ot2:car_dev)
  gca_en_caru_3 <- update(gca_en_caru_2, . ~ . + ot3:car_dev)
  
  en_caru_anova <-
    anova(gca_en_stressu_2, gca_en_caru_0, gca_en_caru_1,
          gca_en_caru_2, gca_en_caru_3)
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stressu_2   29 69796 70012 -34869    69738                     
  # gca_en_caru_0      30 69797 70021 -34869    69737 0.7352  1     0.3912
  # gca_en_caru_1      31 69799 70030 -34869    69737 0.1326  1     0.7158
  # gca_en_caru_2      32 69800 70039 -34868    69736 0.7852  1     0.3755
  # gca_en_caru_3      33 69800 70046 -34867    69734 2.5737  1     0.1087
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_caru_int_0 <- update(gca_en_stressu_2, . ~ . + use_std:stress_sum:car_dev)
  gca_en_caru_int_1 <- update(gca_en_caru_int_0,   . ~ . + ot1:use_std:stress_sum:car_dev)
  gca_en_caru_int_2 <- update(gca_en_caru_int_1,   . ~ . + ot2:use_std:stress_sum:car_dev)
  gca_en_caru_int_3 <- update(gca_en_caru_int_2,   . ~ . + ot3:use_std:stress_sum:car_dev)
  
  en_caru_int_anova <-
    anova(gca_en_stressu_2, gca_en_caru_int_0, gca_en_caru_int_1,
          gca_en_caru_int_2, gca_en_caru_int_3)
  #                 npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_en_stressu_2    29 69796 70012 -34869    69738                          
  # gca_en_caru_int_0   30 69798 70021 -34869    69738  0.0264  1   0.870844    
  # gca_en_caru_int_1   31 69792 70022 -34865    69730  8.4332  1   0.003684 ** 
  # gca_en_caru_int_2   32 69776 70014 -34856    69712 18.1044  1  2.091e-05 ***
  # gca_en_caru_int_3   33 69774 70020 -34854    69708  3.1558  1   0.075657 .  
  
  
  # # BRANCH #2
  # # add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
  # gca_en_corsi_0 <- update(gca_en_stress_2,    . ~ . + corsi)
  # gca_en_corsi_1 <- update(gca_en_corsi_0, . ~ . + ot1:corsi)
  # gca_en_corsi_2 <- update(gca_en_corsi_1, . ~ . + ot2:corsi)
  # gca_en_corsi_3 <- update(gca_en_corsi_2, . ~ . + ot3:corsi)
  # 
  # en_corsi_anova <-
  #   anova(gca_en_stress_2, gca_en_corsi_0, gca_en_corsi_1,
  #         gca_en_corsi_2, gca_en_corsi_3)
  # #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # # gca_en_stress_2   
  # 
  # # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  # gca_en_corsi_int_0 <- update(gca_en_stress_2, . ~ . + prof_std:stress_sum:corsi)
  # gca_en_corsi_int_1 <- update(gca_en_corsi_int_0,   . ~ . + ot1:prof_std:stress_sum:corsi)
  # gca_en_corsi_int_2 <- update(gca_en_corsi_int_1,   . ~ . + ot2:prof_std:stress_sum:corsi)
  # gca_en_corsi_int_3 <- update(gca_en_corsi_int_2,   . ~ . + ot3:prof_std:stress_sum:corsi)
  # 
  # en_corsi_int_anova <-
  #   anova(gca_en_stress_2, gca_en_corsi_int_0, gca_en_corsi_int_1,
  #         gca_en_corsi_int_2, gca_en_corsi_int_3)
  # #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq) 
  # # gca_en_stress_2      
  
  
  # BRANCH #3
  # add verbal processing speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_ospanu_0 <- update(gca_en_stressu_2,    . ~ . + ospan_rt)
  gca_en_ospanu_1 <- update(gca_en_ospanu_0, . ~ . + ot1:ospan_rt)
  gca_en_ospanu_2 <- update(gca_en_ospanu_1, . ~ . + ot2:ospan_rt)
  gca_en_ospanu_3 <- update(gca_en_ospanu_2, . ~ . + ot3:ospan_rt)
  
  en_ospanu_anova <-
    anova(gca_en_stressu_2, gca_en_ospanu_0, gca_en_ospanu_1,
          gca_en_ospanu_2, gca_en_ospanu_3)
  #                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stressu_2   29 69796 70012 -34869    69738                     
  # gca_en_ospanu_0    30 69798 70021 -34869    69738 0.1574  1     0.6915
  # gca_en_ospanu_1    31 69797 70028 -34868    69735 2.6049  1     0.1065
  # gca_en_ospanu_2    32 69799 70038 -34868    69735 0.0396  1     0.8423
  # gca_en_ospanu_3    33 69800 70046 -34867    69734 1.3721  1     0.2415

  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_ospanu_int_0 <- update(gca_en_stressu_2, . ~ . + stress_sum:use_std:ospan_rt)
  gca_en_ospanu_int_1 <- update(gca_en_ospanu_int_0,   . ~ . + ot1:stress_sum:use_std:ospan_rt)
  gca_en_ospanu_int_2 <- update(gca_en_ospanu_int_1,   . ~ . + ot2:stress_sum:use_std:ospan_rt)
  gca_en_ospanu_int_3 <- update(gca_en_ospanu_int_2,   . ~ . + ot3:stress_sum:use_std:ospan_rt)
  
  en_ospanu_int_anova <-
    anova(gca_en_stressu_2, gca_en_ospanu_int_0, gca_en_ospanu_int_1,
          gca_en_ospanu_int_2, gca_en_ospanu_int_3)
  #                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_stressu_2      29 69796 70012 -34869    69738                     
  # gca_en_ospanu_int_0   30 69798 70021 -34869    69738 0.4214  1     0.5162
  # gca_en_ospanu_int_1   31 69800 70030 -34869    69738 0.0547  1     0.8151
  # gca_en_ospanu_int_2   32 69801 70039 -34869    69737 0.6120  1     0.4340
  # gca_en_ospanu_int_3   33 69803 70048 -34868    69737 0.4372  1     0.5085
  
  
  # BRANCH #4
  # add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_corsirtu_0 <- update(gca_en_stressu_2,    . ~ . + corsi_rt)
  gca_en_corsirtu_1 <- update(gca_en_corsirtu_0, . ~ . + ot1:corsi_rt)
  gca_en_corsirtu_2 <- update(gca_en_corsirtu_1, . ~ . + ot2:corsi_rt)
  gca_en_corsirtu_3 <- update(gca_en_corsirtu_2, . ~ . + ot3:corsi_rt)
  
  en_corsirtu_anova <-
    anova(gca_en_stressu_2, gca_en_corsirtu_0, gca_en_corsirtu_1,
          gca_en_corsirtu_2, gca_en_corsirtu_3)
  #                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stressu_2    29 69796 70012 -34869    69738                     
  # gca_en_corsirtu_0   30 69798 70021 -34869    69738 0.4740  1     0.4912
  # gca_en_corsirtu_1   31 69799 70030 -34869    69737 0.3872  1     0.5338
  # gca_en_corsirtu_2   32 69800 70038 -34868    69736 1.7370  1     0.1875
  # gca_en_corsirtu_3   33 69801 70047 -34868    69735 0.1142  1     0.7354
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_corsirtu_int_0 <- update(gca_en_stressu_2, . ~ . + stress_sum:use_std:corsi_rt)
  gca_en_corsirtu_int_1 <- update(gca_en_corsirtu_int_0,   . ~ . + ot1:stress_sum:use_std:corsi_rt)
  gca_en_corsirtu_int_2 <- update(gca_en_corsirtu_int_1,   . ~ . + ot2:stress_sum:use_std:corsi_rt)
  gca_en_corsirtu_int_3 <- update(gca_en_corsirtu_int_2,   . ~ . + ot3:stress_sum:use_std:corsi_rt)
  
  en_corsirtu_int_anova <-
    anova(gca_en_stressu_2, gca_en_corsirtu_int_0, gca_en_corsirtu_int_1,
          gca_en_corsirtu_int_2, gca_en_corsirtu_int_3)
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stressu_2        29 69796 70012 -34869    69738                     
  # gca_en_corsirtu_int_0   30 69798 70021 -34869    69738 0.0226  1     0.8804
  # gca_en_corsirtu_int_1   31 69798 70029 -34868    69736 1.9633  1     0.1612
  # gca_en_corsirtu_int_2   32 69800 70038 -34868    69736 0.1033  1     0.7479
  # gca_en_corsirtu_int_3   33 69802 70048 -34868    69736 0.0018  1     0.9658
  
}

# -----------------------------------------------------------------------------





# Random effects structure ----------------------------------------------------


# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + stress_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = ma_vision, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + stress_sum + ot1 | participant) +
             ot2 + (1 + stress_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + stress_sum + ot1 + ot2 | participant) +
             ot3 + (1 + stress_sum + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    9 72784 72851 -36383    72766                        
  # mod_ot2   14 72672 72777 -36322    72644 121.8  5  < 2.2e-16 ***
  # mod_ot3   20 72663 72812 -36311    72623  21.5  6   0.001491 ** 
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df   AIC   BIC logLik deviance    Chisq Chi Pr(>Chisq)
  # mod_ot3   20 72663 72812 -36311    72623                          
  # mod_ot4   21 72620 72777 -36289    72578 44.3458  1  2.752e-11 ***
  # mod_ot5   23 72604 72776 -36279    72558 19.9677  2  4.614e-05 ***
  # mod_ot6   26 72607 72801 -36277    72555  3.6806  3     0.2981    
  # mod_ot7   30 72609 72833 -36275    72549  5.7561  4     0.2181  
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_ma_base <- mod_ot5
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +          
  #        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 | target),
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 3e5)),
  #      data = ma_vision, REML = F)    # , na.action = na.exclude 
  
  # add proficimacy effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_prof_0 <- update(gca_ma_base,    . ~ . + prof_std)
  gca_ma_prof_1 <- update(gca_ma_prof_0, . ~ . + ot1:prof_std)
  gca_ma_prof_2 <- update(gca_ma_prof_1, . ~ . + ot2:prof_std)
  gca_ma_prof_3 <- update(gca_ma_prof_2, . ~ . + ot3:prof_std)
  
  ma_prof_anova <-
    anova(gca_ma_base, gca_ma_prof_0, gca_ma_prof_1,
          gca_ma_prof_2, gca_ma_prof_3)
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_ma_base     23 72604 72776 -36279    72558                       
  # gca_ma_prof_0   24 72602 72782 -36277    72554 4.0220  1    0.04491 *
  # gca_ma_prof_1   25 72602 72789 -36276    72552 2.8224  1    0.09296 .
  # gca_ma_prof_2   26 72603 72797 -36275    72551 0.9656  1    0.32577  
  # gca_ma_prof_3   27 72604 72806 -36275    72550 0.8689  1    0.35126 
  
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_stress_0 <- update(gca_ma_prof_0,    . ~ . + stress_sum)
  gca_ma_stress_1 <- update(gca_ma_stress_0, . ~ . + ot1:stress_sum)
  gca_ma_stress_2 <- update(gca_ma_stress_1, . ~ . + ot2:stress_sum)
  gca_ma_stress_3 <- update(gca_ma_stress_2, . ~ . + ot3:stress_sum)
  
  ma_stress_anova <-
    anova(gca_ma_prof_0, gca_ma_stress_0, gca_ma_stress_1,
          gca_ma_stress_2, gca_ma_stress_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_prof_0     24 72602 72782 -36277    72554                       
  # gca_ma_stress_0   25 72601 72788 -36276    72551 3.2276  1     0.0724 .
  # gca_ma_stress_1   26 72603 72797 -36275    72551 0.4659  1     0.4949  
  # gca_ma_stress_2   27 72603 72804 -36274    72549 2.2476  1     0.1338  
  # gca_ma_stress_3   28 72604 72814 -36274    72548 0.0473  1     0.8279 
  
  
  # BRANCH #1
  # add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_car_0 <- update(gca_ma_prof_0,    . ~ . + car_dev)
  gca_ma_car_1 <- update(gca_ma_car_0, . ~ . + ot1:car_dev)
  gca_ma_car_2 <- update(gca_ma_car_1, . ~ . + ot2:car_dev)
  gca_ma_car_3 <- update(gca_ma_car_2, . ~ . + ot3:car_dev)
  
  ma_car_anova <-
    anova(gca_ma_prof_0, gca_ma_car_0, gca_ma_car_1,
          gca_ma_car_2, gca_ma_car_3)
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_prof_0   24 72602 72782 -36277    72554                     
  # gca_ma_car_0    25 72604 72791 -36277    72554 0.1567  1     0.6922
  # gca_ma_car_1    26 72606 72801 -36277    72554 0.0846  1     0.7712
  # gca_ma_car_2    27 72607 72809 -36277    72553 0.9138  1     0.3391
  # gca_ma_car_3    28 72608 72818 -36276    72552 1.0501  1     0.3055
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_car_int_0 <- update(gca_ma_prof_0, . ~ . + prof_std:stress_sum:car_dev)
  gca_ma_car_int_1 <- update(gca_ma_car_int_0,   . ~ . + ot1:prof_std:stress_sum:car_dev)
  gca_ma_car_int_2 <- update(gca_ma_car_int_1,   . ~ . + ot2:prof_std:stress_sum:car_dev)
  gca_ma_car_int_3 <- update(gca_ma_car_int_2,   . ~ . + ot3:prof_std:stress_sum:car_dev)
  
  ma_car_int_anova <-
    anova(gca_ma_prof_0, gca_ma_car_int_0, gca_ma_car_int_1,
          gca_ma_car_int_2, gca_ma_car_int_3)
  #                 npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_ma_prof_0      24 72602 72782 -36277    72554                     
  # gca_ma_car_int_0   25 72604 72791 -36277    72554 0.0744  1     0.7850
  # gca_ma_car_int_1   26 72606 72801 -36277    72554 0.0970  1     0.7554
  # gca_ma_car_int_2   27 72606 72808 -36276    72552 2.1930  1     0.1386
  # gca_ma_car_int_3   28 72608 72817 -36276    72552 0.2741  1     0.6006
  
  
  # # BRANCH #2
  # # add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
  # gca_ma_corsi_0 <- update(gca_ma_stress_0,    . ~ . + corsi)
  # gca_ma_corsi_1 <- update(gca_ma_corsi_0, . ~ . + ot1:corsi)
  # gca_ma_corsi_2 <- update(gca_ma_corsi_1, . ~ . + ot2:corsi)
  # gca_ma_corsi_3 <- update(gca_ma_corsi_2, . ~ . + ot3:corsi)
  # 
  # ma_corsi_anova <-
  #   anova(gca_ma_stress_0, gca_ma_corsi_0, gca_ma_corsi_1,
  #         gca_ma_corsi_2, gca_ma_corsi_3)
  # #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # # gca_ma_
  #
  #
  # # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  # gca_ma_corsi_int_0 <- update(gca_ma_stress_0, . ~ . + prof_std:stress_sum:corsi)
  # gca_ma_corsi_int_1 <- update(gca_ma_corsi_int_0,   . ~ . + ot1:prof_std:stress_sum:corsi)
  # gca_ma_corsi_int_2 <- update(gca_ma_corsi_int_1,   . ~ . + ot2:prof_std:stress_sum:corsi)
  # gca_ma_corsi_int_3 <- update(gca_ma_corsi_int_2,   . ~ . + ot3:prof_std:stress_sum:corsi)
  # 
  # ma_corsi_int_anova <-
  #   anova(gca_ma_stress_0, gca_ma_corsi_int_0, gca_ma_corsi_int_1,
  #         gca_ma_corsi_int_2, gca_ma_corsi_int_3)
  # #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq) 
  # # gca_ma_
  
  
  # BRANCH #3
  # add verbal processing speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_ospan_0 <- update(gca_ma_prof_0,    . ~ . + ospan_rt)
  gca_ma_ospan_1 <- update(gca_ma_ospan_0, . ~ . + ot1:ospan_rt)
  gca_ma_ospan_2 <- update(gca_ma_ospan_1, . ~ . + ot2:ospan_rt)
  gca_ma_ospan_3 <- update(gca_ma_ospan_2, . ~ . + ot3:ospan_rt)
  
  ma_ospan_anova <-
    anova(gca_ma_prof_0, gca_ma_ospan_0, gca_ma_ospan_1,
          gca_ma_ospan_2, gca_ma_ospan_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_prof_0    24 72602 72782 -36277    72554                     
  # gca_ma_ospan_0   25 72604 72791 -36277    72554 0.7581  1     0.3839
  # gca_ma_ospan_1   26 72604 72799 -36276    72552 1.2849  1     0.2570
  # gca_ma_ospan_2   27 72605 72807 -36275    72551 1.5260  1     0.2167
  # gca_ma_ospan_3   28 72606 72816 -36275    72550 0.6988  1     0.4032
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_ospan_int_0 <- update(gca_ma_prof_0, . ~ . + stress_sum:prof_std:ospan_rt)
  gca_ma_ospan_int_1 <- update(gca_ma_ospan_int_0,   . ~ . + ot1:stress_sum:prof_std:ospan_rt)
  gca_ma_ospan_int_2 <- update(gca_ma_ospan_int_1,   . ~ . + ot2:stress_sum:prof_std:ospan_rt)
  gca_ma_ospan_int_3 <- update(gca_ma_ospan_int_2,   . ~ . + ot3:stress_sum:prof_std:ospan_rt)
  
  ma_ospan_int_anova <-
    anova(gca_ma_prof_0, gca_ma_ospan_int_0, gca_ma_ospan_int_1,
          gca_ma_ospan_int_2, gca_ma_ospan_int_3)
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_prof_0        24 72602 72782 -36277    72554                     
  # gca_ma_ospan_int_0   25 72604 72791 -36277    72554 0.0002  1     0.9885
  # gca_ma_ospan_int_1   26 72606 72800 -36277    72554 0.5431  1     0.4611
  # gca_ma_ospan_int_2   27 72608 72810 -36277    72554 0.1479  1     0.7005
  # gca_ma_ospan_int_3   28 72610 72819 -36277    72554 0.0934  1     0.7599
  
  
  # BRANCH #4
  # add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_corsirt_0 <- update(gca_ma_prof_0,    . ~ . + corsi_rt)
  gca_ma_corsirt_1 <- update(gca_ma_corsirt_0, . ~ . + ot1:corsi_rt)
  gca_ma_corsirt_2 <- update(gca_ma_corsirt_1, . ~ . + ot2:corsi_rt)
  gca_ma_corsirt_3 <- update(gca_ma_corsirt_2, . ~ . + ot3:corsi_rt)
  
  ma_corsirt_anova <-
    anova(gca_ma_prof_0, gca_ma_corsirt_0, gca_ma_corsirt_1,
          gca_ma_corsirt_2, gca_ma_corsirt_3)
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_prof_0      24 72602 72782 -36277    72554                        
  # gca_ma_corsirt_0   25 72604 72791 -36277    72554 0.0021  1    0.96338   
  # gca_ma_corsirt_1   26 72599 72793 -36273    72547 7.5622  1    0.00596 **
  # gca_ma_corsirt_2   27 72601 72803 -36273    72547 0.0352  1    0.85127   
  # gca_ma_corsirt_3   28 72603 72812 -36273    72547 0.2289  1    0.63235  
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_corsirt_int_0 <- update(gca_ma_corsirt_1, . ~ . + stress_sum:prof_std:corsi_rt)
  gca_ma_corsirt_int_1 <- update(gca_ma_corsirt_int_0,   . ~ . + ot1:stress_sum:prof_std:corsi_rt)
  gca_ma_corsirt_int_2 <- update(gca_ma_corsirt_int_1,   . ~ . + ot2:stress_sum:prof_std:corsi_rt)
  gca_ma_corsirt_int_3 <- update(gca_ma_corsirt_int_2,   . ~ . + ot3:stress_sum:prof_std:corsi_rt)
  
  ma_corsirt_int_anova <-
    anova(gca_ma_corsirt_1, gca_ma_corsirt_int_0, gca_ma_corsirt_int_1,
          gca_ma_corsirt_int_2, gca_ma_corsirt_int_3)
  #                      npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_corsirt_1       26 72599 72793 -36273    72547                     
  # gca_ma_corsirt_int_0   27 72599 72801 -36273    72545 1.4382  1     0.2304
  # gca_ma_corsirt_int_1   28 72601 72811 -36273    72545 0.0942  1     0.7589
  # gca_ma_corsirt_int_2   29 72602 72819 -36272    72544 1.0066  1     0.3157
  # gca_ma_corsirt_int_3   30 72604 72828 -36272    72544 0.5308  1     0.4663

  
  summary(gca_ma_car_int_2)
  # Estimate Std. Error t value
  # (Intercept)   
  
  
  
  
  
  # add use effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_use_0 <- update(gca_ma_base,    . ~ . + use_std)
  gca_ma_use_1 <- update(gca_ma_use_0, . ~ . + ot1:use_std)
  gca_ma_use_2 <- update(gca_ma_use_1, . ~ . + ot2:use_std)
  gca_ma_use_3 <- update(gca_ma_use_2, . ~ . + ot3:use_std)
  
  ma_use_anova <-
    anova(gca_ma_base, gca_ma_use_0, gca_ma_use_1,
          gca_ma_use_2, gca_ma_use_3)
  #              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_ma_base    23 72604 72776 -36279    72558                       
  # gca_ma_use_0   24 72606 72786 -36279    72558 0.0077  1     0.9302  
  # gca_ma_use_1   25 72602 72789 -36276    72552 6.1287  1     0.0133 *
  # gca_ma_use_2   26 72604 72798 -36276    72552 0.2423  1     0.6226  
  # gca_ma_use_3   27 72604 72806 -36275    72550 2.1751  1     0.1403  
  
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_stressu_0 <- update(gca_ma_use_1,    . ~ . + stress_sum)
  gca_ma_stressu_1 <- update(gca_ma_stressu_0, . ~ . + ot1:stress_sum)
  gca_ma_stressu_2 <- update(gca_ma_stressu_1, . ~ . + ot2:stress_sum)
  gca_ma_stressu_3 <- update(gca_ma_stressu_2, . ~ . + ot3:stress_sum)
  
  ma_stressu_anova <-
    anova(gca_ma_use_1, gca_ma_stressu_0, gca_ma_stressu_1,
          gca_ma_stressu_2, gca_ma_stressu_3)
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_use_1       25 72602 72789 -36276    72552                       
  # gca_ma_stressu_0   26 72601 72796 -36275    72549 3.1849  1    0.07432 .
  # gca_ma_stressu_1   27 72603 72805 -36274    72549 0.4441  1    0.50513  
  # gca_ma_stressu_2   28 72602 72812 -36273    72546 2.2988  1    0.12947  
  # gca_ma_stressu_3   29 72604 72821 -36273    72546 0.0681  1    0.79414  
  
  
  
  # BRANCH #1
  # add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_caru_0 <- update(gca_ma_use_1,    . ~ . + car_dev)
  gca_ma_caru_1 <- update(gca_ma_caru_0, . ~ . + ot1:car_dev)
  gca_ma_caru_2 <- update(gca_ma_caru_1, . ~ . + ot2:car_dev)
  gca_ma_caru_3 <- update(gca_ma_caru_2, . ~ . + ot3:car_dev)
  
  ma_caru_anova <-
    anova(gca_ma_use_1, gca_ma_caru_0, gca_ma_caru_1,
          gca_ma_caru_2, gca_ma_caru_3)
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_use_1    25 72602 72789 -36276    72552                     
  # gca_ma_caru_0   26 72604 72799 -36276    72552 0.1508  1     0.6977
  # gca_ma_caru_1   27 72606 72808 -36276    72552 0.0119  1     0.9133
  # gca_ma_caru_2   28 72607 72817 -36276    72551 0.9967  1     0.3181
  # gca_ma_caru_3   29 72608 72825 -36275    72550 1.0596  1     0.3033
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_caru_int_0 <- update(gca_ma_use_1, . ~ . + use_std:stress_sum:car_dev)
  gca_ma_caru_int_1 <- update(gca_ma_caru_int_0,   . ~ . + ot1:use_std:stress_sum:car_dev)
  gca_ma_caru_int_2 <- update(gca_ma_caru_int_1,   . ~ . + ot2:use_std:stress_sum:car_dev)
  gca_ma_caru_int_3 <- update(gca_ma_caru_int_2,   . ~ . + ot3:use_std:stress_sum:car_dev)
  
  ma_caru_int_anova <-
    anova(gca_ma_use_1, gca_ma_caru_int_0, gca_ma_caru_int_1,
          gca_ma_caru_int_2, gca_ma_caru_int_3)
  #                   npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_ma_use_1        25 72602 72789 -36276    72552                       
  # gca_ma_caru_int_0   26 72604 72799 -36276    72552 0.0915  1    0.76223  
  # gca_ma_caru_int_1   27 72606 72808 -36276    72552 0.0017  1    0.96716  
  # gca_ma_caru_int_2   28 72608 72817 -36276    72552 0.6198  1    0.43113  
  # gca_ma_caru_int_3   29 72607 72823 -36274    72549 2.9869  1    0.08394 .
  
  
  # # BRANCH #2
  # # add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
  # gca_ma_corsi_0 <- update(gca_ma_stress_2,    . ~ . + corsi)
  # gca_ma_corsi_1 <- update(gca_ma_corsi_0, . ~ . + ot1:corsi)
  # gca_ma_corsi_2 <- update(gca_ma_corsi_1, . ~ . + ot2:corsi)
  # gca_ma_corsi_3 <- update(gca_ma_corsi_2, . ~ . + ot3:corsi)
  # 
  # ma_corsi_anova <-
  #   anova(gca_ma_stress_2, gca_ma_corsi_0, gca_ma_corsi_1,
  #         gca_ma_corsi_2, gca_ma_corsi_3)
  # #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # # gca_ma_stress_2   
  # 
  # # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  # gca_ma_corsi_int_0 <- update(gca_ma_stress_2, . ~ . + prof_std:stress_sum:corsi)
  # gca_ma_corsi_int_1 <- update(gca_ma_corsi_int_0,   . ~ . + ot1:prof_std:stress_sum:corsi)
  # gca_ma_corsi_int_2 <- update(gca_ma_corsi_int_1,   . ~ . + ot2:prof_std:stress_sum:corsi)
  # gca_ma_corsi_int_3 <- update(gca_ma_corsi_int_2,   . ~ . + ot3:prof_std:stress_sum:corsi)
  # 
  # ma_corsi_int_anova <-
  #   anova(gca_ma_stress_2, gca_ma_corsi_int_0, gca_ma_corsi_int_1,
  #         gca_ma_corsi_int_2, gca_ma_corsi_int_3)
  # #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq) 
  # # gca_ma_stress_2      
  
  
  # BRANCH #3
  # add verbal processing speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_ospanu_0 <- update(gca_ma_use_1,    . ~ . + ospan_rt)
  gca_ma_ospanu_1 <- update(gca_ma_ospanu_0, . ~ . + ot1:ospan_rt)
  gca_ma_ospanu_2 <- update(gca_ma_ospanu_1, . ~ . + ot2:ospan_rt)
  gca_ma_ospanu_3 <- update(gca_ma_ospanu_2, . ~ . + ot3:ospan_rt)
  
  ma_ospanu_anova <-
    anova(gca_ma_use_1, gca_ma_ospanu_0, gca_ma_ospanu_1,
          gca_ma_ospanu_2, gca_ma_ospanu_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_use_1      25 72602 72789 -36276    72552                     
  # gca_ma_ospanu_0   26 72603 72797 -36276    72551 1.2432  1     0.2649
  # gca_ma_ospanu_1   27 72604 72806 -36275    72550 0.6967  1     0.4039
  # gca_ma_ospanu_2   28 72605 72814 -36275    72549 1.2663  1     0.2605
  # gca_ma_ospanu_3   29 72606 72823 -36274    72548 0.7912  1     0.3737
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_ospanu_int_0 <- update(gca_ma_use_1, . ~ . + stress_sum:use_std:ospan_rt)
  gca_ma_ospanu_int_1 <- update(gca_ma_ospanu_int_0,   . ~ . + ot1:stress_sum:use_std:ospan_rt)
  gca_ma_ospanu_int_2 <- update(gca_ma_ospanu_int_1,   . ~ . + ot2:stress_sum:use_std:ospan_rt)
  gca_ma_ospanu_int_3 <- update(gca_ma_ospanu_int_2,   . ~ . + ot3:stress_sum:use_std:ospan_rt)
  
  ma_ospanu_int_anova <-
    anova(gca_ma_use_1, gca_ma_ospanu_int_0, gca_ma_ospanu_int_1,
          gca_ma_ospanu_int_2, gca_ma_ospanu_int_3)
  #                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_ma_use_1          25 72602 72789 -36276    72552                     
  # gca_ma_ospanu_int_0   26 72604 72799 -36276    72552 0.0568  1     0.8117
  # gca_ma_ospanu_int_1   27 72605 72807 -36275    72551 1.4026  1     0.2363
  # gca_ma_ospanu_int_2   28 72607 72816 -36275    72551 0.3333  1     0.5637
  # gca_ma_ospanu_int_3   29 72608 72825 -36275    72550 0.1433  1     0.7050
  
  
  # BRANCH #4
  # add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_corsirtu_0 <- update(gca_ma_use_1,    . ~ . + corsi_rt)
  gca_ma_corsirtu_1 <- update(gca_ma_corsirtu_0, . ~ . + ot1:corsi_rt)
  gca_ma_corsirtu_2 <- update(gca_ma_corsirtu_1, . ~ . + ot2:corsi_rt)
  gca_ma_corsirtu_3 <- update(gca_ma_corsirtu_2, . ~ . + ot3:corsi_rt)
  
  ma_corsirtu_anova <-
    anova(gca_ma_use_1, gca_ma_corsirtu_0, gca_ma_corsirtu_1,
          gca_ma_corsirtu_2, gca_ma_corsirtu_3)
  #                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_use_1        25 72602 72789 -36276    72552                       
  # gca_ma_corsirtu_0   26 72604 72799 -36276    72552 0.0013  1    0.97144  
  # gca_ma_corsirtu_1   27 72602 72803 -36274    72548 4.8129  1    0.02825 *
  # gca_ma_corsirtu_2   28 72603 72813 -36274    72547 0.0423  1    0.83706  
  # gca_ma_corsirtu_3   29 72605 72822 -36274    72547 0.2411  1    0.62343  
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_corsirtu_int_0 <- update(gca_ma_corsirtu_1, . ~ . + stress_sum:use_std:corsi_rt)
  gca_ma_corsirtu_int_1 <- update(gca_ma_corsirtu_int_0,   . ~ . + ot1:stress_sum:use_std:corsi_rt)
  gca_ma_corsirtu_int_2 <- update(gca_ma_corsirtu_int_1,   . ~ . + ot2:stress_sum:use_std:corsi_rt)
  gca_ma_corsirtu_int_3 <- update(gca_ma_corsirtu_int_2,   . ~ . + ot3:stress_sum:use_std:corsi_rt)
  
  ma_corsirtu_int_anova <-
    anova(gca_ma_corsirtu_1, gca_ma_corsirtu_int_0, gca_ma_corsirtu_int_1,
          gca_ma_corsirtu_int_2, gca_ma_corsirtu_int_3)
  #                       npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_corsirtu_1       27 72602 72803 -36274    72548                       
  # gca_ma_corsirtu_int_0   28 72601 72810 -36272    72545 2.7629  1    0.09647 .
  # gca_ma_corsirtu_int_1   29 72603 72819 -36272    72545 0.2357  1    0.62736  
  # gca_ma_corsirtu_int_2   30 72602 72826 -36271    72542 2.7791  1    0.09550 .
  # gca_ma_corsirtu_int_3   31 72601 72833 -36270    72539 2.5170  1    0.11263  
  
  
} 
  
mod_type <- "gca_mon"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_ospan_0", "_ospan_1", "_ospan_2", "_ospan_3",
              "_ospan_int_0", "_ospan_int_1", "_ospan_int_2", "_ospan_int_3",
              "_corsirt_0", "_corsirt_1", "_corsirt_2", "_corsirt_3",
              "_corsirt_int_0", "_corsirt_int_1", "_corsirt_int_2", "_corsirt_int_3",
              "_car_0", "_car_1", "_car_2", "_car_3",
              # "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3"
              # "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3"
              )

# Store ind models in list
mon_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods,
     file = here("mods", "vision", "gca", "continuous_speed", "verb",
                 "mon_mods.Rdata"))

  

mod_type <- "gca_en"
mod_spec <- c("_base", 
              "_prof_0", "_prof_1", "_prof_2", "_prof_3",
              "_use_0", "_use_1", "_use_2", "_use_3",
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_ospan_0", "_ospan_1", "_ospan_2", "_ospan_3",
              "_ospan_int_0", "_ospan_int_1", "_ospan_int_2", "_ospan_int_3",
              "_corsirt_0", "_corsirt_1", "_corsirt_2", "_corsirt_3",
              "_corsirt_int_0", "_corsirt_int_1", "_corsirt_int_2", "_corsirt_int_3",
              "_car_0", "_car_1", "_car_2", "_car_3",
              # "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              # "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3",
              "_stressu_0", "_stressu_1", "_stressu_2", "_stressu_3",
              "_ospanu_0", "_ospanu_1", "_ospanu_2", "_ospanu_3",
              "_ospanu_int_0", "_ospanu_int_1", "_ospanu_int_2", "_ospanu_int_3",
              "_corsirtu_0", "_corsirtu_1", "_corsirtu_2", "_corsirtu_3",
              "_corsirtu_int_0", "_corsirtu_int_1", "_corsirtu_int_2", "_corsirtu_int_3",
              "_caru_0", "_caru_1", "_caru_2", "_caru_3",
              "_caru_int_0", "_caru_int_1", "_caru_int_2", "_caru_int_3"
              )

# Store ind models in list
en_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(en_mods,
     file = here("mods", "vision", "gca", "continuous_speed", "verb",
                 "en_mods.Rdata"))


mod_type <- "gca_ma"
mod_spec <- c("_base", 
              "_prof_0", "_prof_1", "_prof_2", "_prof_3",
              "_use_0", "_use_1", "_use_2", "_use_3",
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_ospan_0", "_ospan_1", "_ospan_2", "_ospan_3",
              "_ospan_int_0", "_ospan_int_1", "_ospan_int_2", "_ospan_int_3",
              "_corsirt_0", "_corsirt_1", "_corsirt_2", "_corsirt_3",
              "_corsirt_int_0", "_corsirt_int_1", "_corsirt_int_2", "_corsirt_int_3",
              "_car_0", "_car_1", "_car_2", "_car_3",
              # "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              # "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3",
              "_stressu_0", "_stressu_1", "_stressu_2", "_stressu_3",
              "_ospanu_0", "_ospanu_1", "_ospanu_2", "_ospanu_3",
              "_ospanu_int_0", "_ospanu_int_1", "_ospanu_int_2", "_ospanu_int_3",
              "_corsirtu_0", "_corsirtu_1", "_corsirtu_2", "_corsirtu_3",
              "_corsirtu_int_0", "_corsirtu_int_1", "_corsirtu_int_2", "_corsirtu_int_3",
              "_caru_0", "_caru_1", "_caru_2", "_caru_3",
              "_caru_int_0", "_caru_int_1", "_caru_int_2", "_caru_int_3"
)

# Store ind models in list
ma_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(ma_mods,
     file = here("mods", "vision", "gca", "continuous_speed", "verb",
                 "ma_mods.Rdata"))
  

  
# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_ospan <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(ospan_rt = c(-1, 0, 1)))
  
mon_corsi <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(corsi_rt = c(-1, 0, 1)))

mon_car <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(car_dev = c(-1, 0, 1)))

en_ospan <- en_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan_rt = c(-1, 0, 1)))

en_corsi <- en_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(corsi_rt = c(-1, 0, 1)))

en_car <- en_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(car_dev = c(-1, 0, 1)))

ma_car <- ma_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(car_dev = c(-1, 0, 1))) 

ma_corsi <- ma_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(corsi_rt = c(-1, 0, 1)))

ma_ospan <- ma_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan_rt = c(-1, 0, 1)))


# Get model predictions and SE
fits_mon_corsirt <- predictSE(gca_mon_corsirt_int_3, mon_corsi) %>%        
  as_tibble %>%
  bind_cols(mon_corsi) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_mon_ospan <- predictSE(gca_mon_ospan_int_3, mon_ospan) %>%        
  as_tibble %>%
  bind_cols(mon_ospan) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_mon_car <- predictSE(gca_mon_car_int_1, mon_car) %>%        
  as_tibble %>%
  bind_cols(mon_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_en_ospan <- predictSE(gca_en_ospan_1, en_ospan) %>%        
  as_tibble %>%
  bind_cols(en_ospan) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_en_car <- predictSE(gca_en_car_int_2, en_car) %>%        
  as_tibble %>%
  bind_cols(en_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_en_corsi <- predictSE(gca_en_stress_2, en_corsi) %>%        
  as_tibble %>%
  bind_cols(en_corsi) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ma_corsi <- predictSE(gca_ma_corsirt_1, ma_corsi) %>%        
  as_tibble %>%
  bind_cols(ma_corsi) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ma_car <- predictSE(gca_ma_prof_0, ma_car) %>%        
  as_tibble %>%
  bind_cols(ma_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ma_ospan <- predictSE(gca_ma_prof_0, ma_ospan) %>%        
  as_tibble %>%
  bind_cols(ma_ospan) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_mon_corsi <- filter(fits_mon_corsirt, time_zero == 4) %>%
  select(stress = stress_sum, corsi_rt,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_mon_ospan <- filter(fits_mon_ospan, time_zero == 4) %>%
  select(stress = stress_sum, ospan_rt, 
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_mon_car <- filter(fits_mon_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, 
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_en_ospan <- filter(fits_en_ospan, time_zero == 4) %>%
  select(stress = stress_sum, ospan_rt, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_en_car <- filter(fits_en_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_en_corsi <- filter(fits_en_corsi, time_zero == 4) %>%
  select(stress = stress_sum, corsi_rt, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_ma_ospan <- filter(fits_ma_ospan, time_zero == 4) %>%
  select(stress = stress_sum, ospan_rt, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_ma_car <- filter(fits_ma_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_ma_corsi <- filter(fits_ma_corsi, time_zero == 4) %>%
  select(stress = stress_sum, corsi_rt, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 




# Save models predictions
model_preds <- mget(c("fits_mon_corsirt", 'fits_mon_ospan', 'fits_mon_car',
                      "fits_en_ospan", 'fits_en_car', 'fits_en_corsi',
                      "fits_ma_corsi", 'fits_ma_car', 'fits_ma_ospan',
                      "preds_mon_corsi", "preds_mon_ospan", 'preds_mon_car',
                      "preds_en_ospan", 'preds_en_car', 'preds_en_corsi',
                      "preds_ma_ospan", 'preds_ma_car', 'preds_ma_corsi'
                      ))

save(model_preds,
     file = here("mods", "vision", "gca", "continuous_speed", 'verb',
                 "model_preds.Rdata"))



# -----------------------------------------------------------------------------






# Save models -----------------------------------------------------------------

if(F) {
  
  # Save anova model comparisons
  nested_model_comparisons <-
    mget(c("mon_stress_anova", "mon_car_anova", "mon_corsi_anova",
           "mon_car_int_anova", "mon_corsi_int_anova",
           "en_stress_anova", "en_car_anova", "en_corsi_anova",
           "en_prof_anova", "en_car_int_anova", "en_corsi_int_anova",
           "ma_stress_anova", "ma_car_anova", "ma_corsi_anova",
           "ma_prof_anova", "ma_car_int_anova", "ma_corsi_int_anova"
    ))
  
  save(nested_model_comparisons,
       file = here("mods", "vision", "gca", "continuous",
                   "nested_model_comparisons.Rdata"))
  
  
  
  
  
  
}

# -----------------------------------------------------------------------------

