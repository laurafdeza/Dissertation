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
source(here::here("scripts", "00_load_libs.R"))
# source(here::here("scripts", "02_load_data.R"))
stress50 <- read_csv(here("data", "clean", "stress_50ms_final_onsetc3updated.csv"))


# Get path to saved models
gca_mods_path  <- here("mods", "vision", "gca", "cont_speed_verb_equal")

# Load models as lists
load(paste0(gca_mods_path, "/mon_mods.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))

# Store objects in global env
list2env(mon_mods, globalenv())
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
vision <- read_csv("./data/clean/vision_scores_nooutliers-400.csv") # pred car
corsi <- read_csv("./data/clean/corsi_z_scores.csv")
wm_score <- read_csv("./data/clean/ospan_set_z_scores.csv")

vision50 <- left_join(x = stress50, y = wm, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = vision, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = corsi, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = wm_score, by = "participant", all.x=TRUE)

vision50 <- na.omit(vision50)




vision50 <- vision50 %>%
  filter(., l1 == 'es' & time_zero >= -4 & time_zero <= 4) %>%
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
  # mod_ot1    9 14440 14493 -7211.0    14422                        
  # mod_ot2   14 14439 14521 -7205.6    14411 10.7710  5    0.05612 .
  # mod_ot3   20 14451 14568 -7205.4    14411  0.2625  6    0.99966 
  
  mod_ot4 <- update(mod_ot1, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + ot2 + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + ot3 + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot1, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df   AIC   BIC  logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot1    9 14440 14493 -7211.0    14422                          
  # mod_ot4   10 14344 14403 -7162.0    14324 97.9045  1  < 2.2e-16 ***
  # mod_ot5   12 14328 14399 -7152.2    14304 19.6173  2  5.498e-05 ***
  # mod_ot6   16 14318 14412 -7143.0    14286 18.4595  4   0.001003 ** 
  # mod_ot7   21 14326 14449 -7141.8    14284  2.2576  5   0.812473 
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_mon_base <- mod_ot6
  # lmer(eLog ~ 1 + ot1 + ot2 +          
  #        (1 + stress_sum + ot1 | participant) +
  #        (1 + ot1 + ot2 | target),
  #      control = lmerControl(optimizer = 'bobyqa'), #, optCtrl = list(maxfun = 3e5)
  #      data = stress_gc_subset, REML = F)    # , na.action = na.exclude 

# add stress effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_stress_0 <- update(gca_mon_base,    . ~ . + stress_sum)
gca_mon_stress_1 <- update(gca_mon_stress_0, . ~ . + ot1:stress_sum)
gca_mon_stress_2 <- update(gca_mon_stress_1, . ~ . + ot2:stress_sum)

mon_stress_anova <-
  anova(gca_mon_base, gca_mon_stress_0, gca_mon_stress_1,
        gca_mon_stress_2)
                 # npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base       16 14318 14412 -7143.0    14286                     
# gca_mon_stress_0   17 14320 14419 -7142.8    14286 0.2712  1     0.6025
# gca_mon_stress_1   18 14320 14425 -7141.8    14284 2.1587  1     0.1418
# gca_mon_stress_2   19 14322 14433 -7141.8    14284 0.0054  1     0.9417

# BRANCH #0
# add verbal proc time effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_0 <- update(gca_mon_base,    . ~ . + ospan_rt)
gca_mon_ospan_1 <- update(gca_mon_ospan_0, . ~ . + ot1:ospan_rt)
gca_mon_ospan_2 <- update(gca_mon_ospan_1, . ~ . + ot2:ospan_rt)

mon_ospan_anova <-
  anova(gca_mon_base, gca_mon_ospan_0, gca_mon_ospan_1,
        gca_mon_ospan_2)
#                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      16 14318 14412 -7143.0    14286                     
# gca_mon_ospan_0   17 14319 14419 -7142.6    14285 0.7161  1     0.3974
# gca_mon_ospan_1   18 14321 14427 -7142.6    14285 0.0299  1     0.8626
# gca_mon_ospan_2   19 14322 14433 -7142.0    14284 1.1725  1     0.2789


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_int_0 <- update(gca_mon_base, . ~ . + stress_sum:ospan_rt)
gca_mon_ospan_int_1 <- update(gca_mon_ospan_int_0,   . ~ . + ot1:stress_sum:ospan_rt)
gca_mon_ospan_int_2 <- update(gca_mon_ospan_int_1,   . ~ . + ot2:stress_sum:ospan_rt)

mon_ospan_int_anova <-
  anova(gca_mon_base, gca_mon_ospan_int_0, gca_mon_ospan_int_1,
        gca_mon_ospan_int_2)
#                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mon_base          16 14318 14412 -7143.0    14286                       
# gca_mon_ospan_int_0   17 14320 14419 -7142.9    14286 0.2058  1    0.65004  
# gca_mon_ospan_int_1   18 14321 14426 -7142.4    14285 1.0012  1    0.31701  
# gca_mon_ospan_int_2   19 14317 14429 -7139.7    14279 5.4086  1    0.02004 *


# BRANCH #1
# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_0 <- update(gca_mon_base,    . ~ . + car_dev)
gca_mon_car_1 <- update(gca_mon_car_0, . ~ . + ot1:car_dev)
gca_mon_car_2 <- update(gca_mon_car_1, . ~ . + ot2:car_dev)

mon_car_anova <-
  anova(gca_mon_base, gca_mon_car_0, gca_mon_car_1,
        gca_mon_car_2)
#               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base    16 14318 14412 -7143.0    14286                     
# gca_mon_car_0   17 14319 14418 -7142.4    14285 1.1272  1     0.2884
# gca_mon_car_1   18 14319 14425 -7141.6    14283 1.5332  1     0.2156
# gca_mon_car_2   19 14320 14432 -7141.2    14282 0.8847  1     0.3469


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_int_0 <- update(gca_mon_base, . ~ . + stress_sum:car_dev)
gca_mon_car_int_1 <- update(gca_mon_car_int_0,   . ~ . + ot1:stress_sum:car_dev)
gca_mon_car_int_2 <- update(gca_mon_car_int_1,   . ~ . + ot2:stress_sum:car_dev)

mon_car_int_anova <-
  anova(gca_mon_base, gca_mon_car_int_0, gca_mon_car_int_1,
        gca_mon_car_int_2)
#                   npar   AIC   BIC  logLik deviance   Chisq Df Pr(>Chisq)  
# gca_mon_base        16 14318 14412 -7143.0    14286                     
# gca_mon_car_int_0   17 14320 14420 -7142.9    14286 0.1033  1     0.7479
# gca_mon_car_int_1   18 14321 14426 -7142.5    14285 0.8051  1     0.3696
# gca_mon_car_int_2   19 14323 14434 -7142.4    14285 0.1387  1     0.7096



# BRANCH #2
# add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_0 <- update(gca_mon_base,    . ~ . + corsi)
gca_mon_corsi_1 <- update(gca_mon_corsi_0, . ~ . + ot1:corsi)
gca_mon_corsi_2 <- update(gca_mon_corsi_1, . ~ . + ot2:corsi)

mon_corsi_anova <-
  anova(gca_mon_base, gca_mon_corsi_0, gca_mon_corsi_1,
        gca_mon_corsi_2)
#                 npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      16 14318 14412 -7143.0    14286                        
# gca_mon_corsi_0   17 14319 14419 -7142.5    14285 0.9592  1   0.327377   
# gca_mon_corsi_1   18 14319 14425 -7141.7    14283 1.5530  1   0.212697   
# gca_mon_corsi_2   19 14314 14425 -7138.0    14276 7.4109  1   0.006483 **


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_int_0 <- update(gca_mon_corsi_2, . ~ . + stress_sum:corsi)
gca_mon_corsi_int_1 <- update(gca_mon_corsi_int_0,   . ~ . + ot1:stress_sum:corsi)
gca_mon_corsi_int_2 <- update(gca_mon_corsi_int_1,   . ~ . + ot2:stress_sum:corsi)
 
mon_corsi_int_anova <-
  anova(gca_mon_corsi_2, gca_mon_corsi_int_0, gca_mon_corsi_int_1,
        gca_mon_corsi_int_2)
#                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_corsi_2       19 14314 14425 -7138.0    14276                       
# gca_mon_corsi_int_0   20 14316 14433 -7138.0    14276 0.0172  1    0.89576  
# gca_mon_corsi_int_1   21 14317 14440 -7137.3    14275 1.3864  1    0.23901  
# gca_mon_corsi_int_2   22 14316 14445 -7135.8    14272 2.9185  1    0.08757 .



# BRANCH #2.5
# add car pred effect to intercept, linear slope, quadratic, and cubic time terms to corsi
gca_mon_carcor_0 <- update(gca_mon_corsi_2,    . ~ . + car_dev)
gca_mon_carcor_1 <- update(gca_mon_carcor_0, . ~ . + ot1:car_dev)
gca_mon_carcor_2 <- update(gca_mon_carcor_1, . ~ . + ot2:car_dev)


mon_carcor_anova <-
  anova(gca_mon_corsi_2, gca_mon_carcor_0, gca_mon_carcor_1,
        gca_mon_carcor_2)
#                  npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_corsi_2    19 14314 14425 -7138.0    14276                     
# gca_mon_carcor_0   20 14314 14431 -7137.1    14274 1.7792  1     0.1822
# gca_mon_carcor_1   21 14314 14437 -7136.2    14272 1.9159  1     0.1663
# gca_mon_carcor_2   22 14315 14444 -7135.4    14271 1.5076  1     0.2195


# add 3-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_carcor_int_0 <- update(gca_mon_corsi_2, . ~ . + stress_sum:corsi:car_dev)
gca_mon_carcor_int_1 <- update(gca_mon_carcor_int_0,   . ~ . + ot1:stress_sum:corsi:car_dev)
gca_mon_carcor_int_2 <- update(gca_mon_carcor_int_1,   . ~ . + ot2:stress_sum:corsi:car_dev)


mon_carcor_int_anova <-
  anova(gca_mon_corsi_2, gca_mon_carcor_int_0, gca_mon_carcor_int_1,
        gca_mon_carcor_int_2)
#                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_corsi_2       19 14314 14425 -7138.0    14276                       
# gca_mon_carcor_int_0   20 14316 14433 -7137.9    14276 0.2456  1     0.6202
# gca_mon_carcor_int_1   21 14318 14441 -7137.9    14276 0.0059  1     0.9387
# gca_mon_carcor_int_2   22 14319 14448 -7137.3    14275 1.0920  1     0.2960




# BRANCH #3
# add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_0 <- update(gca_mon_base,    . ~ . + corsi_rt)
gca_mon_corsirt_1 <- update(gca_mon_corsirt_0, . ~ . + ot1:corsi_rt)
gca_mon_corsirt_2 <- update(gca_mon_corsirt_1, . ~ . + ot2:corsi_rt)

mon_corsirt_anova <-
  anova(gca_mon_base, gca_mon_corsirt_0, gca_mon_corsirt_1,
        gca_mon_corsirt_2)
#                   npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)
# gca_mon_base        16 14318 14412 -7143.0    14286                       
# gca_mon_corsirt_0   17 14317 14416 -7141.3    14283 3.3530  1    0.06708 .
# gca_mon_corsirt_1   18 14317 14422 -7140.3    14281 2.0050  1    0.15678  
# gca_mon_corsirt_2   19 14318 14430 -7140.1    14280 0.3365  1    0.56184

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_int_0 <- update(gca_mon_base, . ~ . + stress_sum:corsi_rt)
gca_mon_corsirt_int_1 <- update(gca_mon_corsirt_int_0,   . ~ . + ot1:stress_sum:corsi_rt)
gca_mon_corsirt_int_2 <- update(gca_mon_corsirt_int_1,   . ~ . + ot2:stress_sum:corsi_rt)

mon_corsirt_int_anova <-
  anova(gca_mon_base, gca_mon_corsirt_int_0, gca_mon_corsirt_int_1,
        gca_mon_corsirt_int_2)
#                       npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_base            16 14318 14412 -7143.0    14286                        
# gca_mon_corsirt_int_0   17 14320 14420 -7142.9    14286 0.0456  1   0.830955   
# gca_mon_corsirt_int_1   18 14322 14427 -7142.9    14286 0.0533  1   0.817493   
# gca_mon_corsirt_int_2   19 14316 14428 -7139.1    14278 7.6237  1   0.005761 **



# BRANCH #3.5
gca_mon_carcrt_0 <- update(gca_mon_corsirt_int_2,    . ~ . + car_dev)
gca_mon_carcrt_1 <- update(gca_mon_carcrt_0, . ~ . + ot1:car_dev)
gca_mon_carcrt_2 <- update(gca_mon_carcrt_1, . ~ . + ot2:car_dev)


mon_carcrt_anova <-
  anova(gca_mon_corsirt_int_2, gca_mon_carcrt_0, gca_mon_carcrt_1,
        gca_mon_carcrt_2)
#                       npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_corsirt_int_2   19 14316 14428 -7139.1    14278                     
# gca_mon_carcrt_0        20 14317 14434 -7138.5    14277 1.1863  1     0.2761
# gca_mon_carcrt_1        21 14318 14440 -7137.7    14276 1.5494  1     0.2132
# gca_mon_carcrt_2        22 14319 14448 -7137.4    14275 0.7457  1     0.3878


# add 3-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_carcrt_int_0 <- update(gca_mon_corsirt_int_2, . ~ . + stress_sum:corsi_rt:car_dev)
gca_mon_carcrt_int_1 <- update(gca_mon_carcrt_int_0,   . ~ . + ot1:stress_sum:corsi_rt:car_dev)
gca_mon_carcrt_int_2 <- update(gca_mon_carcrt_int_1,   . ~ . + ot2:stress_sum:corsi_rt:car_dev)


mon_carcrt_int_anova <-
  anova(gca_mon_corsirt_int_2, gca_mon_carcrt_int_0, gca_mon_carcrt_int_1,
        gca_mon_carcrt_int_2)
#                       npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_corsirt_int_2   19 14316 14428 -7139.1    14278                     
# gca_mon_carcrt_int_0    20 14317 14434 -7138.7    14277 0.8858  1     0.3466
# gca_mon_carcrt_int_1    21 14318 14442 -7138.2    14276 0.8554  1     0.3550
# gca_mon_carcrt_int_2    22 14320 14449 -7138.0    14276 0.4486  1     0.5030





# BRANCH #4
# add verbal wm score effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospanscore_0 <- update(gca_mon_base,    . ~ . + ospan)
gca_mon_ospanscore_1 <- update(gca_mon_ospanscore_0, . ~ . + ot1:ospan)
gca_mon_ospanscore_2 <- update(gca_mon_ospanscore_1, . ~ . + ot2:ospan)

mon_ospanscore_anova <-
  anova(gca_mon_base, gca_mon_ospanscore_0, gca_mon_ospanscore_1,
        gca_mon_ospanscore_2)
#                      npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base           16 14318 14412 -7143.0    14286                     
# gca_mon_ospanscore_0   17 14320 14419 -7142.9    14286 0.2359  1     0.6272
# gca_mon_ospanscore_1   18 14321 14426 -7142.4    14285 0.8467  1     0.3575
# gca_mon_ospanscore_2   19 14322 14433 -7141.9    14284 1.0680  1     0.3014

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospanscore_int_0 <- update(gca_mon_base, . ~ . + stress_sum:ospan)
gca_mon_ospanscore_int_1 <- update(gca_mon_ospanscore_int_0,   . ~ . + ot1:stress_sum:ospan)
gca_mon_ospanscore_int_2 <- update(gca_mon_ospanscore_int_1,   . ~ . + ot2:stress_sum:ospan)

mon_ospanscore_int_anova <-
  anova(gca_mon_base, gca_mon_ospanscore_int_0, gca_mon_ospanscore_int_1,
        gca_mon_ospanscore_int_2)
#                          npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_base               16 14318 14412 -7143.0    14286                       
# gca_mon_ospanscore_int_0   17 14319 14419 -7142.5    14285 0.9533  1    0.32889  
# gca_mon_ospanscore_int_1   18 14317 14422 -7140.4    14281 4.2806  1    0.03855 *
# gca_mon_ospanscore_int_2   19 14319 14430 -7140.4    14281 0.0002  1    0.99021  

}

# -----------------------------------------------------------------------------






  
mod_type <- "gca_mon"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", 
              "_ospan_0", "_ospan_1", "_ospan_2", 
              "_ospan_int_0", "_ospan_int_1", "_ospan_int_2", 
              "_corsirt_0", "_corsirt_1", "_corsirt_2", 
              "_corsirt_int_0", "_corsirt_int_1", "_corsirt_int_2", 
              "_car_0", "_car_1", "_car_2", 
              "_corsi_0", "_corsi_1", "_corsi_2", 
              "_car_int_0", "_car_int_1", "_car_int_2", 
              "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", 
              "_ospanscore_0", "_ospanscore_1", "_ospanscore_2", 
              "_ospanscore_int_0", "_ospanscore_int_1", "_ospanscore_int_2" ,
              '_carcrt_0', '_carcrt_1', '_carcrt_2',
              '_carcrt_int_0', '_carcrt_int_1', '_carcrt_int_2',
              '_carcor_0', '_carcor_1', '_carcor_2',
              '_carcor_int_0', '_carcor_int_1', '_carcor_int_2'
              )

# Store ind models in list
mon_mods_allequal400 <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods_allequal400,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "mon_mods_allequal400.Rdata"))

  


  
# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_ospanrt <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(ospan_rt = c(-1, 0, 1)))
  
mon_corsirt <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(corsi_rt = c(-1, 0, 1)))

mon_ospansc <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(ospan = c(-1, 0, 1)))

mon_corsisc <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(corsi = c(-1, 0, 1)))



# Get model predictions and SE
fits_mon_ospanrt <- predictSE(gca_mon_ospan_int_2, mon_ospanrt) %>%        
  as_tibble %>%
  bind_cols(mon_ospanrt) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


fits_mon_corsirt <- predictSE(gca_mon_corsirt_int_2, mon_corsirt) %>%        
  as_tibble %>%
  bind_cols(mon_corsirt) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


fits_mon_ospansc <- predictSE(gca_mon_ospanscore_int_2, mon_ospansc) %>%        
  as_tibble %>%
  bind_cols(mon_ospansc) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


fits_mon_corsisc <- predictSE(gca_mon_corsi_2, mon_corsisc) %>%        
  as_tibble %>%
  bind_cols(mon_corsisc) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at target syllable offset
preds_mon_ospanrt <- filter(fits_mon_ospanrt, time_zero == 4) %>%
  select(stress = stress_sum, ospan_rt, 
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_mon_corsirt <- filter(fits_mon_corsirt, time_zero == 4) %>%
  select(stress = stress_sum, corsi_rt,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_mon_ospansc <- filter(fits_mon_ospansc, time_zero == 4) %>%
  select(stress = stress_sum, ospan, 
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

preds_mon_corsisc <- filter(fits_mon_corsisc, time_zero == 4) %>%
  select(stress = stress_sum, corsi, 
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 






# Save models predictions
model_predsalleq400 <- mget(c("fits_mon_corsirt", 'fits_mon_ospanrt', 
                      'fits_mon_corsisc', 'fits_mon_ospansc',
                      "preds_mon_corsirt", "preds_mon_ospanrt", 
                      'preds_mon_corsisc', 'preds_mon_ospansc'
                      ))

save(model_predsalleq400,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "model_predsalleq400.Rdata"))



# -----------------------------------------------------------------------------






# Save models -----------------------------------------------------------------

if(F) {
  
  # Save anova model comparisons
  nested_model_comparisons <-
    mget(c("mon_stress_anova", "mon_car_anova", "mon_corsi_anova",
           "mon_car_int_anova", "mon_corsi_int_anova"
    ))
  
  save(nested_model_comparisons,
       file = here("mods", "vision", "gca", "continuous",
                   "nested_model_comparisons.Rdata"))
  
  
  
  
  
  
}

# -----------------------------------------------------------------------------

