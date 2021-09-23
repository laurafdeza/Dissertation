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
gca_mods_path  <- here("mods", "vision", "gca", "continuous")

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

vision50 <- left_join(x = stress50, y = wm, by = "participant", all.x=TRUE)

vision50 <- na.omit(vision50)




vision50 <- vision50 %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma"),
            stress_sum = if_else(cond == "1", 1, -1)) %>%           # 1 = present, 2 = preterit        
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
  # mod_ot1    9 42756 42819 -21369    42738                         
  # mod_ot2   14 42657 42755 -21315    42629 108.82  5  < 2.2e-16 ***
  # mod_ot3   20 42549 42688 -21254    42509 120.67  6  < 2.2e-16 ***
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df   AIC   BIC logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot3   20 42549 42688 -21254    42509                          
  # mod_ot4   21 42408 42554 -21183    42366 142.536   1  < 2.2e-16 ***
  # mod_ot5   23 42196 42356 -21075    42150 216.441   2  < 2.2e-16 ***
  # mod_ot6   26 42120 42301 -21034    42068  81.740   3  < 2.2e-16 ***
  # mod_ot7   30 42114 42323 -21027    42054  13.564   4   0.008826 ** 
  
  
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
# gca_mon_base       30 42114 42323 -21027    42054                     
# gca_mon_stress_0   31 42116 42332 -21027    42054 0.5089  1     0.4756
# gca_mon_stress_1   32 42118 42340 -21027    42054 0.2444  1     0.6210
# gca_mon_stress_2   33 42120 42349 -21027    42054 0.0006  1     0.9800
# gca_mon_stress_3   34 42121 42358 -21026    42053 0.8150  1     0.3667

# BRANCH #0
# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_0 <- update(gca_mon_base,    . ~ . + ospan_rt)
gca_mon_ospan_1 <- update(gca_mon_ospan_0, . ~ . + ot1:ospan_rt)
gca_mon_ospan_2 <- update(gca_mon_ospan_1, . ~ . + ot2:ospan_rt)
gca_mon_ospan_3 <- update(gca_mon_ospan_2, . ~ . + ot3:ospan_rt)

mon_ospan_anova <-
  anova(gca_mon_base, gca_mon_ospan_0, gca_mon_ospan_1,
        gca_mon_ospan_2, gca_mon_ospan_3)
#                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      30 42114 42323 -21027    42054                     
# gca_mon_ospan_0   31 42116 42332 -21027    42054 0.1219  1     0.7270
# gca_mon_ospan_1   32 42118 42341 -21027    42054 0.1811  1     0.6704
# gca_mon_ospan_2   33 42120 42350 -21027    42054 0.2511  1     0.6163
# gca_mon_ospan_3   34 42121 42358 -21027    42053 0.4146  1     0.5197


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_int_0 <- update(gca_mon_base, . ~ . + stress_sum:ospan_rt)
gca_mon_ospan_int_1 <- update(gca_mon_ospan_int_0,   . ~ . + ot1:stress_sum:ospan_rt)
gca_mon_ospan_int_2 <- update(gca_mon_ospan_int_1,   . ~ . + ot2:stress_sum:ospan_rt)
gca_mon_ospan_int_3 <- update(gca_mon_ospan_int_2,   . ~ . + ot3:stress_sum:ospan_rt)

mon_ospan_int_anova <-
  anova(gca_mon_base, gca_mon_ospan_int_0, gca_mon_ospan_int_1,
        gca_mon_ospan_int_2, gca_mon_ospan_int_3)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mon_base          30 42114 42323 -21027    42054                        
# gca_mon_ospan_int_0   31 42116 42332 -21027    42054 0.2567  1   0.612410   
# gca_mon_ospan_int_1   32 42118 42341 -21027    42054 0.1576  1   0.691391   
# gca_mon_ospan_int_2   33 42115 42345 -21025    42049 4.5861  1   0.032231 * 
# gca_mon_ospan_int_3   34 42110 42347 -21021    42042 7.3791  1   0.006598 **


# BRANCH #1
# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_0 <- update(gca_mon_base,    . ~ . + car_dev)
gca_mon_car_1 <- update(gca_mon_car_0, . ~ . + ot1:car_dev)
gca_mon_car_2 <- update(gca_mon_car_1, . ~ . + ot2:car_dev)
gca_mon_car_3 <- update(gca_mon_car_2, . ~ . + ot3:car_dev)

mon_car_anova <-
  anova(gca_mon_base, gca_mon_car_0, gca_mon_car_1,
        gca_mon_car_2, gca_mon_car_3)
# npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base    


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_int_0 <- update(gca_mon_base, . ~ . + stress_sum:car_dev)
gca_mon_car_int_1 <- update(gca_mon_car_int_0,   . ~ . + ot1:stress_sum:car_dev)
gca_mon_car_int_2 <- update(gca_mon_car_int_1,   . ~ . + ot2:stress_sum:car_dev)
gca_mon_car_int_3 <- update(gca_mon_car_int_2,   . ~ . + ot3:stress_sum:car_dev)

mon_car_int_anova <-
  anova(gca_mon_base, gca_mon_car_int_0, gca_mon_car_int_1,
        gca_mon_car_int_2, gca_mon_car_int_3)
#                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_base        

# BRANCH #2
# add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_0 <- update(gca_mon_base,    . ~ . + corsi)
gca_mon_corsi_1 <- update(gca_mon_corsi_0, . ~ . + ot1:corsi)
gca_mon_corsi_2 <- update(gca_mon_corsi_1, . ~ . + ot2:corsi)
gca_mon_corsi_3 <- update(gca_mon_corsi_2, . ~ . + ot3:corsi)

mon_corsi_anova <-
  anova(gca_mon_base, gca_mon_corsi_0, gca_mon_corsi_1,
        gca_mon_corsi_2, gca_mon_corsi_3)
#                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_int_0 <- update(gca_mon_base, . ~ . + stress_sum:corsi)
gca_mon_corsi_int_1 <- update(gca_mon_corsi_int_0,   . ~ . + ot1:stress_sum:corsi)
gca_mon_corsi_int_2 <- update(gca_mon_corsi_int_1,   . ~ . + ot2:stress_sum:corsi)
gca_mon_corsi_int_3 <- update(gca_mon_corsi_int_2,   . ~ . + ot3:stress_sum:corsi)

mon_corsi_int_anova <-
  anova(gca_mon_base, gca_mon_corsi_int_0, gca_mon_corsi_int_1,
        gca_mon_corsi_int_2, gca_mon_corsi_int_3)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_base    


# BRANCH #4
# add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_0 <- update(gca_mon_base,    . ~ . + corsi_rt)
gca_mon_corsirt_1 <- update(gca_mon_corsirt_0, . ~ . + ot1:corsi_rt)
gca_mon_corsirt_2 <- update(gca_mon_corsirt_1, . ~ . + ot2:corsi_rt)
gca_mon_corsirt_3 <- update(gca_mon_corsirt_2, . ~ . + ot3:corsi_rt)

mon_corsirt_anova <-
  anova(gca_mon_base, gca_mon_corsirt_0, gca_mon_corsirt_1,
        gca_mon_corsirt_2, gca_mon_corsirt_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base        30 42114 42323 -21027    42054                       
# gca_mon_corsirt_0   31 42116 42332 -21027    42054 0.6660  1    0.41444  
# gca_mon_corsirt_1   32 42118 42341 -21027    42054 0.0268  1    0.86997  
# gca_mon_corsirt_2   33 42118 42348 -21026    42052 1.9012  1    0.16794  
# gca_mon_corsirt_3   34 42117 42353 -21024    42049 3.0886  1    0.07884 .


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_int_0 <- update(gca_mon_base, . ~ . + stress_sum:corsi_rt)
gca_mon_corsirt_int_1 <- update(gca_mon_corsirt_int_0,   . ~ . + ot1:stress_sum:corsi_rt)
gca_mon_corsirt_int_2 <- update(gca_mon_corsirt_int_1,   . ~ . + ot2:stress_sum:corsi_rt)
gca_mon_corsirt_int_3 <- update(gca_mon_corsirt_int_2,   . ~ . + ot3:stress_sum:corsi_rt)

mon_corsirt_int_anova <-
  anova(gca_mon_base, gca_mon_corsirt_int_0, gca_mon_corsirt_int_1,
        gca_mon_corsirt_int_2, gca_mon_corsirt_int_3)
#                       npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_base            30 42114 42323 -21027    42054                       
# gca_mon_corsirt_int_0   31 42116 42332 -21027    42054 0.1993  1    0.65530   
# gca_mon_corsirt_int_1   32 42118 42341 -21027    42054 0.4965  1    0.48105   
# gca_mon_corsirt_int_2   33 42117 42347 -21026    42051 2.4773  1    0.11550   
# gca_mon_corsirt_int_3   34 42112 42349 -21022    42044 6.9313  1    0.00847 **




summary(gca_mon_ospan_int_3)
#                           Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                1.68466    0.14549   50.93069  11.579 7.04e-16 ***
# ot1                        4.76964    0.47457   52.38551  10.050 8.04e-14 ***
# ot2                       -1.11755    0.38708   46.10240  -2.887 0.005901 ** 
# ot3                       -0.99571    0.27775   37.30254  -3.585 0.000962 ***
# stress_sum:ospan_rt       -0.02829    0.08324   29.61447  -0.340 0.736367    
# ot1:stress_sum:ospan_rt    0.04133    0.18262 6943.35910   0.226 0.820977    
# ot2:stress_sum:ospan_rt   -0.47847    0.18328 5670.67226  -2.611 0.009060 ** 
# ot3:stress_sum:ospan_rt   -0.48448    0.17807 3940.58229  -2.721 0.006543 ** 

summary(gca_mon_base)
# Estimate Std. Error t value
# (Intercept)   

summary(gca_mon_corsi_int_1)
# Estimate Std. Error t value
# (Intercept) 

summary(gca_mon_corsirt_int_3)
# Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                1.68153    0.14597   50.85627  11.520 8.69e-16 ***
# ot1                        4.76149    0.47665   52.20495   9.989 1.03e-13 ***
# ot2                       -1.11794    0.38864   46.05545  -2.877 0.006072 ** 
# ot3                       -0.99929    0.27823   36.98696  -3.592 0.000951 ***
# stress_sum:corsi_rt        0.02938    0.06356   29.01663   0.462 0.647304    
# ot1:stress_sum:corsi_rt   -0.08083    0.13436 6494.47515  -0.602 0.547463    
# ot2:stress_sum:corsi_rt    0.14256    0.13441 4800.11980   1.061 0.288906    
# ot3:stress_sum:corsi_rt   -0.33599    0.12703 2808.56623  -2.645 0.008215 ** 

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
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    9 90417 90487 -45200    90399                          
  # mod_ot2   14 90391 90499 -45182    90363  36.161  5   8.82e-07 ***
  # mod_ot3   20 90142 90296 -45051    90102 261.153  6  < 2.2e-16 ***
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df    AIC    BIC  logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot3   20 90142 90296 -45051    90102                          
  # mod_ot4   21 89774 89936 -44866    89732 369.684  1  < 2.2e-16 ***
  #   mod_ot5   23 89554 89731 -44754    89508 224.645  2  < 2.2e-16 ***
  #   mod_ot6   26 89429 89629 -44688    89377 131.175  3  < 2.2e-16 ***
  #   mod_ot7   30 89382 89613 -44661    89322  54.298  4  4.558e-11 ***
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_en_base <- mod_ot7
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
  #        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 + ot3 | target),
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
  # gca_en_base     30 89382 89613 -44661    89322                       
  # gca_en_prof_0   31 89384 89623 -44661    89322 0.3159  1    0.57409  
  # gca_en_prof_1   32 89381 89628 -44658    89317 4.9547  1    0.02602 *
  # gca_en_prof_2   33 89378 89632 -44656    89312 4.8743  1    0.02726 *
  # gca_en_prof_3   34 89379 89641 -44656    89311 0.7857  1    0.37540  
  
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_stress_0 <- update(gca_en_prof_2,    . ~ . + stress_sum)
  gca_en_stress_1 <- update(gca_en_stress_0, . ~ . + ot1:stress_sum)
  gca_en_stress_2 <- update(gca_en_stress_1, . ~ . + ot2:stress_sum)
  gca_en_stress_3 <- update(gca_en_stress_2, . ~ . + ot3:stress_sum)
  
  en_stress_anova <-
    anova(gca_en_prof_2, gca_en_stress_0, gca_en_stress_1,
          gca_en_stress_2, gca_en_stress_3)
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_prof_2     33 89378 89632 -44656    89312                       
  # gca_en_stress_0   34 89380 89642 -44656    89312 0.4687  1    0.49359  
  # gca_en_stress_1   35 89381 89651 -44656    89311 0.4470  1    0.50376  
  # gca_en_stress_2   36 89377 89654 -44652    89305 6.4153  1    0.01131 *
  # gca_en_stress_3   37 89378 89663 -44652    89304 1.1393  1    0.28580  

  
  
  # BRANCH #1
  # add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_car_0 <- update(gca_en_stress_2,    . ~ . + car_dev)
  gca_en_car_1 <- update(gca_en_car_0, . ~ . + ot1:car_dev)
  gca_en_car_2 <- update(gca_en_car_1, . ~ . + ot2:car_dev)
  gca_en_car_3 <- update(gca_en_car_2, . ~ . + ot3:car_dev)
  
  en_car_anova <-
    anova(gca_en_stress_2, gca_en_car_0, gca_en_car_1,
          gca_en_car_2, gca_en_car_3)
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stress_2   
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_car_int_0 <- update(gca_en_car_3, . ~ . + prof_std:stress_sum:car_dev)
  gca_en_car_int_1 <- update(gca_en_car_int_0,   . ~ . + ot1:prof_std:stress_sum:car_dev)
  gca_en_car_int_2 <- update(gca_en_car_int_1,   . ~ . + ot2:prof_std:stress_sum:car_dev)
  gca_en_car_int_3 <- update(gca_en_car_int_2,   . ~ . + ot3:prof_std:stress_sum:car_dev)
  
  en_car_int_anova <-
    anova(gca_en_car_3, gca_en_car_int_0, gca_en_car_int_1,
          gca_en_car_int_2, gca_en_car_int_3)
  #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_en_car_3       
  
  
  # BRANCH #2
  # add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_corsi_0 <- update(gca_en_stress_2,    . ~ . + corsi)
  gca_en_corsi_1 <- update(gca_en_corsi_0, . ~ . + ot1:corsi)
  gca_en_corsi_2 <- update(gca_en_corsi_1, . ~ . + ot2:corsi)
  gca_en_corsi_3 <- update(gca_en_corsi_2, . ~ . + ot3:corsi)
  
  en_corsi_anova <-
    anova(gca_en_stress_2, gca_en_corsi_0, gca_en_corsi_1,
          gca_en_corsi_2, gca_en_corsi_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stress_2   
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_corsi_int_0 <- update(gca_en_stress_2, . ~ . + prof_std:stress_sum:corsi)
  gca_en_corsi_int_1 <- update(gca_en_corsi_int_0,   . ~ . + ot1:prof_std:stress_sum:corsi)
  gca_en_corsi_int_2 <- update(gca_en_corsi_int_1,   . ~ . + ot2:prof_std:stress_sum:corsi)
  gca_en_corsi_int_3 <- update(gca_en_corsi_int_2,   . ~ . + ot3:prof_std:stress_sum:corsi)
  
  en_corsi_int_anova <-
    anova(gca_en_stress_2, gca_en_corsi_int_0, gca_en_corsi_int_1,
          gca_en_corsi_int_2, gca_en_corsi_int_3)
  #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq) 
  # gca_en_stress_2      
  
  
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
  # gca_en_stress_2    36 88527 88804 -44227    88455                       
  # gca_en_ospan_0    37 88527 88813 -44227    88453 1.2903  1    0.25599  
  # gca_en_ospan_1    38 88529 88822 -44227    88453 0.2684  1    0.60438  
  # gca_en_ospan_2    39 88526 88826 -44224    88448 5.3803  1    0.02037 *
  # gca_en_ospan_3    40 88527 88836 -44224    88447 0.3429  1    0.55815  
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_ospan_int_0 <- update(gca_en_ospan_2, . ~ . + stress_sum:prof_std:ospan_rt)
  gca_en_ospan_int_1 <- update(gca_en_ospan_int_0,   . ~ . + ot1:stress_sum:prof_std:ospan_rt)
  gca_en_ospan_int_2 <- update(gca_en_ospan_int_1,   . ~ . + ot2:stress_sum:prof_std:ospan_rt)
  gca_en_ospan_int_3 <- update(gca_en_ospan_int_2,   . ~ . + ot3:stress_sum:prof_std:ospan_rt)
  
  en_ospan_int_anova <-
    anova(gca_en_ospan_2, gca_en_ospan_int_0, gca_en_ospan_int_1,
          gca_en_ospan_int_2, gca_en_ospan_int_3)
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_ospan_2       39 88526 88826 -44224    88448                     
  # gca_en_ospan_int_0   40 88527 88835 -44224    88447 0.7669  1     0.3812
  # gca_en_ospan_int_1   41 88528 88844 -44223    88446 1.1359  1     0.2865
  # gca_en_ospan_int_2   42 88530 88854 -44223    88446 0.0961  1     0.7565
  # gca_en_ospan_int_3   43 88531 88863 -44223    88445 0.5003  1     0.4794     
  
  
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
  # gca_en_stress_2    36 88527 88804 -44227    88455                       
  # gca_en_corsirt_0   37 88528 88813 -44227    88454 1.0506  1     0.3054
  # gca_en_corsirt_1   38 88529 88822 -44227    88453 0.6179  1     0.4318
  # gca_en_corsirt_2   39 88531 88832 -44227    88453 0.0701  1     0.7912
  # gca_en_corsirt_3   40 88533 88841 -44226    88453 0.4617  1     0.4968
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_corsirt_int_0 <- update(gca_en_stress_2, . ~ . + stress_sum:prof_std:corsi_rt)
  gca_en_corsirt_int_1 <- update(gca_en_corsirt_int_0,   . ~ . + ot1:stress_sum:prof_std:corsi_rt)
  gca_en_corsirt_int_2 <- update(gca_en_corsirt_int_1,   . ~ . + ot2:stress_sum:prof_std:corsi_rt)
  gca_en_corsirt_int_3 <- update(gca_en_corsirt_int_2,   . ~ . + ot3:stress_sum:prof_std:corsi_rt)
  
  en_corsirt_int_anova <-
    anova(gca_en_stress_2, gca_en_corsirt_int_0, gca_en_corsirt_int_1,
          gca_en_corsirt_int_2, gca_en_corsirt_int_3)
  #                      npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_stress_2        36 88527 88804 -44227    88455                     
  # gca_en_corsirt_int_0   37 88528 88814 -44227    88454 0.4330  1     0.5105
  # gca_en_corsirt_int_1   38 88530 88823 -44227    88454 0.4293  1     0.5123
  # gca_en_corsirt_int_2   39 88531 88832 -44227    88453 0.6036  1     0.4372
  # gca_en_corsirt_int_3   40 88533 88841 -44227    88453 0.1983  1     0.6561
  
  
  summary(gca_en_ospan_2)
  #                Estimate Std. Error       df t value Pr(>|t|)    
  # (Intercept)     1.70610    0.13276 65.79324  12.851  < 2e-16 ***
  # ot1             7.06836    0.42494 71.78093  16.634  < 2e-16 ***
  # ot2            -0.91930    0.27698 52.90552  -3.319  0.00164 ** 
  # ot3            -1.80479    0.23138 42.82692  -7.800 9.37e-10 ***
  # prof_std        0.10874    0.08469 60.56105   1.284  0.20406    
  # stress_sum     -0.14351    0.09354 87.16405  -1.534  0.12858    
  # ospan_rt       -0.33650    0.19924 60.71656  -1.689  0.09637 .  
  # ot1:prof_std    0.57363    0.26764 63.83037   2.143  0.03591 *  
  # ot2:prof_std   -0.32230    0.18815 62.06239  -1.713  0.09170 .  
  # ot1:stress_sum  0.26195    0.23696 36.57400   1.105  0.27617    
  # ot2:stress_sum -0.02263    0.20905 40.96211  -0.108  0.91431    
  # ot1:ospan_rt   -0.86363    0.62955 63.98386  -1.372  0.17491    
  # ot2:ospan_rt    1.05801    0.44278 62.28462   2.389  0.01991 * 
  
  summary(gca_en_car_int_3)
  #               Estimate Std. Error t value
  # (Intercept)    
  
  summary(gca_en_stress_2)
  #               Estimate Std. Error t value
  # (Intercept)    1.73673    0.13283 65.76812  13.075  < 2e-16 ***
  # ot1             7.14714    0.42198 70.74282  16.937  < 2e-16 ***
  #   ot2            -1.01582    0.28094 54.40111  -3.616 0.000656 ***
  #   ot3            -1.80467    0.23139 42.86740  -7.799 9.33e-10 ***
  #   prof_std        0.10158    0.08661 60.57450   1.173 0.245461    
  # stress_sum     -0.14142    0.09360 86.98893  -1.511 0.134450    
  # ot1:prof_std    0.55517    0.27158 63.26195   2.044 0.045101 *  
  #   ot2:prof_std   -0.29960    0.19661 61.95229  -1.524 0.132635    
  # ot1:stress_sum  0.26403    0.23693 36.51579   1.114 0.272414    
  # ot2:stress_sum -0.02879    0.20902 40.82360  -0.138 0.891141         
  

  
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
  # mod_ot1    9 95518 95588 -47750    95500                          
  # mod_ot2   14 95477 95585 -47724    95449  51.274  5  7.599e-10 ***
  # mod_ot3   20 95314 95469 -47637    95274 175.060  6  < 2.2e-16 ***
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df   AIC   BIC logLik deviance    Chisq Chi Pr(>Chisq)
  # mod_ot3   20 95314 95469 -47637    95274                           
  # mod_ot4   21 95162 95325 -47560    95120 153.6754   1  < 2.2e-16 ***
  # mod_ot5   23 94970 95149 -47462    94924 195.9337   2  < 2.2e-16 ***
  # mod_ot6   26 94935 95137 -47442    94883  40.7541   3  7.374e-09 ***
  # mod_ot7   30 94935 95167 -47437    94875   8.5841   4    0.07238 . 
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_ma_base <- mod_ot6
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +          
  #        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 | target),
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
  # gca_ma_base     26 94935 95137 -47442    94883                       
  # gca_ma_prof_0   27 94932 95141 -47439    94878 5.6652  1     0.0173 *
  # gca_ma_prof_1   28 94931 95148 -47438    94875 2.5425  1     0.1108  
  # gca_ma_prof_2   29 94933 95158 -47438    94875 0.1014  1     0.7502  
  # gca_ma_prof_3   30 94934 95167 -47437    94874 0.5995  1     0.4388 
  
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_stress_0 <- update(gca_ma_prof_0,    . ~ . + stress_sum)
  gca_ma_stress_1 <- update(gca_ma_stress_0, . ~ . + ot1:stress_sum)
  gca_ma_stress_2 <- update(gca_ma_stress_1, . ~ . + ot2:stress_sum)
  gca_ma_stress_3 <- update(gca_ma_stress_2, . ~ . + ot3:stress_sum)
  
  ma_stress_anova <-
    anova(gca_ma_prof_0, gca_ma_stress_0, gca_ma_stress_1,
          gca_ma_stress_2, gca_ma_stress_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_prof_0     27 94932 95141 -47439    94878                       
  # gca_ma_stress_0   28 94929 95146 -47436    94873 4.8571  1    0.02753 *
  # gca_ma_stress_1   29 94931 95156 -47436    94873 0.0081  1    0.92818  
  # gca_ma_stress_2   30 94933 95166 -47436    94873 0.0214  1    0.88368  
  # gca_ma_stress_3   31 94934 95175 -47436    94872 0.4953  1    0.48159  
  
  
  # BRANCH #1
  # add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_car_0 <- update(gca_ma_stress_0,    . ~ . + car_dev)
  gca_ma_car_1 <- update(gca_ma_car_0, . ~ . + ot1:car_dev)
  gca_ma_car_2 <- update(gca_ma_car_1, . ~ . + ot2:car_dev)
  gca_ma_car_3 <- update(gca_ma_car_2, . ~ . + ot3:car_dev)
  
  ma_car_anova <-
    anova(gca_ma_stress_0, gca_ma_car_0, gca_ma_car_1,
          gca_ma_car_2, gca_ma_car_3)
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_stress_0   28 90483 90699 -45214    90427                     
  # gca_ma_car_0      29 90485 90709 -45214    90427 0.0335  1     0.8549
  # gca_ma_car_1      30 90487 90718 -45214    90427 0.0344  1     0.8529
  # gca_ma_car_2      31 90489 90728 -45213    90427 0.2442  1     0.6212
  # gca_ma_car_3      32 90490 90737 -45213    90426 0.7783  1     0.3776
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_car_int_0 <- update(gca_ma_stress_0, . ~ . + prof_std:stress_sum:car_dev)
  gca_ma_car_int_1 <- update(gca_ma_car_int_0,   . ~ . + ot1:prof_std:stress_sum:car_dev)
  gca_ma_car_int_2 <- update(gca_ma_car_int_1,   . ~ . + ot2:prof_std:stress_sum:car_dev)
  gca_ma_car_int_3 <- update(gca_ma_car_int_2,   . ~ . + ot3:prof_std:stress_sum:car_dev)
  
  ma_car_int_anova <-
    anova(gca_ma_stress_0, gca_ma_car_int_0, gca_ma_car_int_1,
          gca_ma_car_int_2, gca_ma_car_int_3)
  #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_ma_stress_0    28 90483 90699 -45214    90427                       
  # gca_ma_car_int_0   29 90485 90709 -45214    90427 0.1553  1    0.69349  
  # gca_ma_car_int_1   30 90486 90717 -45213    90426 1.4589  1    0.22711  
  # gca_ma_car_int_2   31 90484 90722 -45211    90422 4.0347  1    0.04457 *
  # gca_ma_car_int_3   32 90485 90732 -45211    90421 0.1190  1    0.73014  
  
  
  # BRANCH #2
  # add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_corsi_0 <- update(gca_ma_stress_0,    . ~ . + corsi)
  gca_ma_corsi_1 <- update(gca_ma_corsi_0, . ~ . + ot1:corsi)
  gca_ma_corsi_2 <- update(gca_ma_corsi_1, . ~ . + ot2:corsi)
  gca_ma_corsi_3 <- update(gca_ma_corsi_2, . ~ . + ot3:corsi)
  
  ma_corsi_anova <-
    anova(gca_ma_stress_0, gca_ma_corsi_0, gca_ma_corsi_1,
          gca_ma_corsi_2, gca_ma_corsi_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_stress_0   28 90483 90699 -45214    90427                     
  # gca_ma_corsi_0    29 90485 90709 -45214    90427 0.0046  1     0.9458
  # gca_ma_corsi_1    30 90485 90717 -45213    90425 1.8563  1     0.1731
  # gca_ma_corsi_2    31 90487 90726 -45212    90425 0.6863  1     0.4074
  # gca_ma_corsi_3    32 90488 90734 -45212    90424 0.8695  1     0.3511
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_corsi_int_0 <- update(gca_ma_stress_0, . ~ . + prof_std:stress_sum:corsi)
  gca_ma_corsi_int_1 <- update(gca_ma_corsi_int_0,   . ~ . + ot1:prof_std:stress_sum:corsi)
  gca_ma_corsi_int_2 <- update(gca_ma_corsi_int_1,   . ~ . + ot2:prof_std:stress_sum:corsi)
  gca_ma_corsi_int_3 <- update(gca_ma_corsi_int_2,   . ~ . + ot3:prof_std:stress_sum:corsi)
  
  ma_corsi_int_anova <-
    anova(gca_ma_stress_0, gca_ma_corsi_int_0, gca_ma_corsi_int_1,
          gca_ma_corsi_int_2, gca_ma_corsi_int_3)
  #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq) 
  # gca_ma_stress_0      28 90483 90699 -45214    90427                     
  # gca_ma_corsi_int_0   29 90485 90708 -45213    90427 0.5119  1     0.4743
  # gca_ma_corsi_int_1   30 90486 90717 -45213    90426 0.7868  1     0.3751
  # gca_ma_corsi_int_2   31 90487 90726 -45212    90425 1.2748  1     0.2589
  # gca_ma_corsi_int_3   32 90488 90735 -45212    90424 0.7072  1     0.4004
  
  
  # BRANCH #3
  # add verbal processing speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_ospan_0 <- update(gca_ma_stress_0,    . ~ . + ospan_rt)
  gca_ma_ospan_1 <- update(gca_ma_ospan_0, . ~ . + ot1:ospan_rt)
  gca_ma_ospan_2 <- update(gca_ma_ospan_1, . ~ . + ot2:ospan_rt)
  gca_ma_ospan_3 <- update(gca_ma_ospan_2, . ~ . + ot3:ospan_rt)
  
  ma_ospan_anova <-
    anova(gca_ma_stress_0, gca_ma_ospan_0, gca_ma_ospan_1,
          gca_ma_ospan_2, gca_ma_ospan_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stress_0   28 94929 95146 -47436    94873                     
  # gca_ma_ospan_0    29 94930 95155 -47436    94872 0.3554  1     0.5511
  # gca_ma_ospan_1    30 94931 95164 -47436    94871 1.1804  1     0.2773
  # gca_ma_ospan_2    31 94931 95171 -47434    94869 2.4717  1     0.1159
  # gca_ma_ospan_3    32 94932 95180 -47434    94868 0.6306  1     0.4271
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_ospan_int_0 <- update(gca_ma_stress_0, . ~ . + stress_sum:prof_std:ospan_rt)
  gca_ma_ospan_int_1 <- update(gca_ma_ospan_int_0,   . ~ . + ot1:stress_sum:prof_std:ospan_rt)
  gca_ma_ospan_int_2 <- update(gca_ma_ospan_int_1,   . ~ . + ot2:stress_sum:prof_std:ospan_rt)
  gca_ma_ospan_int_3 <- update(gca_ma_ospan_int_2,   . ~ . + ot3:stress_sum:prof_std:ospan_rt)
  
  ma_ospan_int_anova <-
    anova(gca_ma_stress_0, gca_ma_ospan_int_0, gca_ma_ospan_int_1,
          gca_ma_ospan_int_2, gca_ma_ospan_int_3)
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_stress_0      28 94929 95146 -47436    94873                     
  # gca_ma_ospan_int_0   29 94931 95156 -47436    94873 0.0340  1     0.8537
  # gca_ma_ospan_int_1   30 94931 95163 -47435    94871 2.2573  1     0.1330
  # gca_ma_ospan_int_2   31 94932 95173 -47435    94870 0.1329  1     0.7155
  # gca_ma_ospan_int_3   32 94934 95182 -47435    94870 0.4802  1     0.4883         
  
  
  # BRANCH #4
  # add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_corsirt_0 <- update(gca_ma_stress_0,    . ~ . + corsi_rt)
  gca_ma_corsirt_1 <- update(gca_ma_corsirt_0, . ~ . + ot1:corsi_rt)
  gca_ma_corsirt_2 <- update(gca_ma_corsirt_1, . ~ . + ot2:corsi_rt)
  gca_ma_corsirt_3 <- update(gca_ma_corsirt_2, . ~ . + ot3:corsi_rt)
  
  ma_corsirt_anova <-
    anova(gca_ma_stress_0, gca_ma_corsirt_0, gca_ma_corsirt_1,
          gca_ma_corsirt_2, gca_ma_corsirt_3)
  #                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stress_0    28 94929 95146 -47436    94873                     
  # gca_ma_corsirt_0   29 94931 95156 -47436    94873 0.0092  1     0.9237
  # gca_ma_corsirt_1   30 94930 95163 -47435    94870 2.6106  1     0.1061
  # gca_ma_corsirt_2   31 94930 95170 -47434    94868 2.2469  1     0.1339
  # gca_ma_corsirt_3   32 94930 95179 -47433    94866 1.5970  1     0.2063
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_corsirt_int_0 <- update(gca_ma_stress_0, . ~ . + stress_sum:prof_std:corsi_rt)
  gca_ma_corsirt_int_1 <- update(gca_ma_corsirt_int_0,   . ~ . + ot1:stress_sum:prof_std:corsi_rt)
  gca_ma_corsirt_int_2 <- update(gca_ma_corsirt_int_1,   . ~ . + ot2:stress_sum:prof_std:corsi_rt)
  gca_ma_corsirt_int_3 <- update(gca_ma_corsirt_int_2,   . ~ . + ot3:stress_sum:prof_std:corsi_rt)
  
  ma_corsirt_int_anova <-
    anova(gca_ma_stress_0, gca_ma_corsirt_int_0, gca_ma_corsirt_int_1,
          gca_ma_corsirt_int_2, gca_ma_corsirt_int_3)
  #                      npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_stress_0        28 94929 95146 -47436    94873                     
  # gca_ma_corsirt_int_0   29 94930 95155 -47436    94872 1.0005  1     0.3172
  # gca_ma_corsirt_int_1   30 94931 95164 -47436    94871 0.7891  1     0.3744
  # gca_ma_corsirt_int_2   31 94931 95171 -47434    94869 2.2060  1     0.1375
  # gca_ma_corsirt_int_3   32 94933 95181 -47434    94869 0.0514  1     0.8207
  
  
  
  
  
  summary(gca_ma_car_int_2)
  # Estimate Std. Error t value
  # (Intercept)   
  
  summary(gca_ma_stress_0)
  # Estimate Std. Error t value
  # (Intercept)  0.97319    0.09286 76.07697  10.481  < 2e-16 ***
  # ot1          4.25036    0.32049 67.91118  13.262  < 2e-16 ***
  # ot2          0.24872    0.21030 59.43678   1.183   0.2416    
  # ot3         -0.97030    0.14476 61.09732  -6.703 7.45e-09 ***
  # prof_std     0.14010    0.05549 56.98473   2.525   0.0144 *  
  # stress_sum  -0.12674    0.05603 69.56325  -2.262   0.0268 *  
  
  
} 
  
mod_type <- "gca_mon"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_ospan_0", "_ospan_1", "_ospan_2", "_ospan_3",
              "_ospan_int_0", "_ospan_int_1", "_ospan_int_2", "_ospan_int_3",
              "_corsirt_0", "_corsirt_1", "_corsirt_2", "_corsirt_3",
              "_corsirt_int_0", "_corsirt_int_1", "_corsirt_int_2", "_corsirt_int_3"
              # "_car_0", "_car_1", "_car_2", "_car_3",
              # "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              # "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              # "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3"
              )

# Store ind models in list
mon_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods,
     file = here("mods", "vision", "gca", "continuous_speed",
                 "mon_mods.Rdata"))

  

mod_type <- "gca_en"
mod_spec <- c("_base", 
              "_prof_0", "_prof_1", "_prof_2", #"_prof_3",
              "_stress_0", "_stress_1", "_stress_2", #"_stress_3",
              "_ospan_0", "_ospan_1", "_ospan_2", "_ospan_3",
              "_ospan_int_0", "_ospan_int_1", "_ospan_int_2", "_ospan_int_3",
              "_corsirt_0", "_corsirt_1", "_corsirt_2", "_corsirt_3",
              "_corsirt_int_0", "_corsirt_int_1", "_corsirt_int_2", "_corsirt_int_3"
              # "_car_0", "_car_1", "_car_2", "_car_3",
              # "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              # "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              # "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3"
              )

# Store ind models in list
en_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(en_mods,
     file = here("mods", "vision", "gca", "continuous_speed",
                 "en_mods.Rdata"))


mod_type <- "gca_ma"
mod_spec <- c("_base", 
              "_prof_0", "_prof_1", "_prof_2", "_prof_3",
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_ospan_0", "_ospan_1", "_ospan_2", "_ospan_3",
              "_ospan_int_0", "_ospan_int_1", "_ospan_int_2", "_ospan_int_3",
              "_corsirt_0", "_corsirt_1", "_corsirt_2", "_corsirt_3",
              "_corsirt_int_0", "_corsirt_int_1", "_corsirt_int_2", "_corsirt_int_3"
              # "_car_0", "_car_1", "_car_2", "_car_3",
              # "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              # "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              # "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3"
              )

# Store ind models in list
ma_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(ma_mods,
     file = here("mods", "vision", "gca", "continuous_speed",
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

en_ospan <- en_vision %>%
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

fits_en_ospan <- predictSE(gca_en_ospan_2, en_ospan) %>%        
  as_tibble %>%
  bind_cols(en_ospan) %>%
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

preds_en_ospan <- filter(fits_en_ospan, time_zero == 4) %>%
  select(stress = stress_sum, ospan_rt, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

# Save models predictions
model_preds <- mget(c("fits_mon_corsirt", 'fits_mon_ospan', 
                      "fits_en_ospan", #'fits_en_car',
                      #"fits_ma_corsi", 'fits_ma_car',
                      "preds_mon_corsi", "preds_mon_ospan",
                      # "target_offset_preds_en_corsi", "target_offset_preds_en_car",
                      # "target_offset_preds_ma_corsi", 
                      "preds_en_ospan"
                      ))

save(model_preds,
     file = here("mods", "vision", "gca", "continuous_speed",
                 "model_preds.Rdata"))











en_corsi <- en_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(corsi = c(-1, 0, 1)))

ma_car <- ma_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(car_dev = c(-1, 0, 1))) 

ma_corsi <- ma_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(corsi = c(-1, 0, 1)))


fits_en_car <- predictSE(gca_en_car_int_3, en_car) %>%        
  as_tibble %>%
  bind_cols(en_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ma_corsi <- predictSE(gca_ma_stress_0, ma_corsi) %>%        
  as_tibble %>%
  bind_cols(ma_corsi) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ma_car <- predictSE(gca_ma_car_int_2, ma_car) %>%        
  as_tibble %>%
  bind_cols(ma_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)



target_offset_preds_en_car <- filter(fits_en_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 


target_offset_preds_ma_corsi <- filter(fits_ma_corsi, time_zero == 4) %>%
  select(stress = stress_sum, corsi, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

target_offset_preds_ma_car <- filter(fits_ma_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

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

