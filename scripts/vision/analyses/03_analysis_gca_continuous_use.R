#
# Original script by Joseph
# Modified by Cristina to adapt to CrisLau Project
# Last update: 06/12/2019
# Modified by Laura to adapt to Pupurri project
#
# Growth curve analyisis ------------------------------------------------------
#
# - Question 1: Does visuospatial prediction abilities (continuous) influence 
#   SS, IE, AE, IM, and AM (categorical) speakers' abilities to predict verbal tense
#   based on the presence or absence (categorical) of lexical stress?
#     - No association in any population
# - Question 2: Do verbal and visuospatial WM influence linguistic prediction?
#     - Not really.
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

# Load data
source(here::here("scripts", "02_load_data.R"))

# Get path to saved models
gca_mods_path  <- here("mods", "vision", "gca", "continuous")

# Load models as lists
#load(paste0(gca_mods_path, "/ind_mods.Rdata"))
# load(paste0(gca_mods_path, "/group_mods.Rdata"))
#load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
# load(paste0(gca_mods_path, "/model_preds.Rdata"))

# Store objects in global env
# list2env(ind_mods, globalenv())
# list2env(group_mods, globalenv())
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

vision <- read_csv("./data/clean/vision_scores.csv")
corsi <- read_csv("./data/clean/corsi_z_scores.csv")

visuospatial_df <- left_join(x = vision, y = corsi, by = "participant", all.x=TRUE)

vision50 <- left_join(x = stress50, y = visuospatial_df, by = "participant", all.x=TRUE)

vision50 <- na.omit(vision50)




vision50 <- vision50 %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "mon", "aes", "ams", "ies", "ims"),
            stress_sum = if_else(cond == "1", 1, -1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

mon_vision <- filter(vision50, l1 == 'es') %>% select(-DELE, -percent_l2_week, 
                                                    -DELE_z, -use_z, -prof, -group)

l2_vision <- vision50 %>%
  filter(., l1 != 'es') %>% 
  filter(., participant != 'ies04' & participant != 'ies17' & participant != 'ies28' & participant != 'aes32') %>%
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1),
         prof_std = (DELE - mean(DELE))/sd(DELE),
         use_std = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week)) %>%
  select(-DELE, -percent_l2_week, 
         -DELE_z, -use_z, -prof, -group)




# -----------------------------------------------------------------------------






# Random effects structure ----------------------------------------------------

en_vision <- l2_vision %>%
  filter(., l1 == 'en')


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
  #         npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
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
  #           Df   AIC   BIC logLik deviance  Chisq Chi Pr(>Chisq)
  # mod_ot3   20 90142 90296 -45051    90102                          
  # mod_ot4   21 89774 89936 -44866    89732 369.684  1  < 2.2e-16 ***
  # mod_ot5   23 89554 89731 -44754    89508 224.645  2  < 2.2e-16 ***
  # mod_ot6   26 89429 89629 -44688    89377 131.175  3  < 2.2e-16 ***
  # mod_ot7   30 89382 89613 -44661    89322  54.298  4  4.558e-11 ***
  
  
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
  #      data = stress_gc_subset, REML = F)    # , na.action = na.exclude 
  
  # add L2 use effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_use_0 <- update(gca_en_base,    . ~ . + use_std)
  gca_en_use_1 <- update(gca_en_use_0, . ~ . + ot1:use_std)
  gca_en_use_2 <- update(gca_en_use_1, . ~ . + ot2:use_std)
  gca_en_use_3 <- update(gca_en_use_2, . ~ . + ot3:use_std)
  
  en_use_anova <-
    anova(gca_en_base, gca_en_use_0, gca_en_use_1,
          gca_en_use_2, gca_en_use_3)
  #              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_en_base    30 89382 89613 -44661    89322                       
  # gca_en_use_0   31 89384 89623 -44661    89322 0.5852  1    0.44428  
  # gca_en_use_1   32 89382 89629 -44659    89318 3.6959  1    0.05455 .
  # gca_en_use_2   33 89384 89638 -44659    89318 0.1008  1    0.75088  
  # gca_en_use_3   34 89382 89644 -44657    89314 4.0585  1    0.04395 *
  
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_stress_0 <- update(gca_en_use_3,    . ~ . + stress_sum)
  gca_en_stress_1 <- update(gca_en_stress_0, . ~ . + ot1:stress_sum)
  gca_en_stress_2 <- update(gca_en_stress_1, . ~ . + ot2:stress_sum)
  gca_en_stress_3 <- update(gca_en_stress_2, . ~ . + ot3:stress_sum)
  
  en_stress_anova <-
    anova(gca_en_use_3, gca_en_stress_0, gca_en_stress_1,
          gca_en_stress_2, gca_en_stress_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_use_3      34 89382 89644 -44657    89314                       
  # gca_en_stress_0   35 89383 89653 -44657    89313 0.5247  1    0.46885  
  # gca_en_stress_1   36 89385 89662 -44656    89313 0.4700  1    0.49299  
  # gca_en_stress_2   37 89381 89666 -44653    89307 6.2704  1    0.01228 *
  # gca_en_stress_3   38 89381 89674 -44653    89305 1.1365  1    0.28639 
  
  
  # BRANCH #1
  # add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_car_0 <- update(gca_en_stress_2,  . ~ . + car_dev)
  gca_en_car_1 <- update(gca_en_car_0, . ~ . + ot1:car_dev)
  gca_en_car_2 <- update(gca_en_car_1, . ~ . + ot2:car_dev)
  gca_en_car_3 <- update(gca_en_car_2, . ~ . + ot3:car_dev)
  
  en_car_anova <-
    anova(gca_en_stress_2, gca_en_car_0, gca_en_car_1,
          gca_en_car_2, gca_en_car_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stress_2   37 89381 89666 -44653    89307                       
  # gca_en_car_0      38 89382 89675 -44653    89306 0.2715  1    0.60232  
  # gca_en_car_1      39 89384 89685 -44653    89306 0.3396  1    0.56004  
  # gca_en_car_2      40 89385 89693 -44652    89305 0.9121  1    0.33956  
  # gca_en_car_3      41 89384 89700 -44651    89302 2.7163  1    0.09933 .
  
  # add 3-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_car_int_0 <- update(gca_en_stress_2, . ~ . + use_std:stress_sum:car_dev)
  gca_en_car_int_1 <- update(gca_en_car_int_0,   . ~ . + ot1:use_std:stress_sum:car_dev)
  gca_en_car_int_2 <- update(gca_en_car_int_1,   . ~ . + ot2:use_std:stress_sum:car_dev)
  gca_en_car_int_3 <- update(gca_en_car_int_2,   . ~ . + ot3:use_std:stress_sum:car_dev)
  
  en_car_int_anova <-
    anova(gca_en_stress_2, gca_en_car_int_0, gca_en_car_int_1,
          gca_en_car_int_2, gca_en_car_int_3)
  #                  npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)  
  # gca_en_stress_2    37 89381 89666 -44653    89307                          
  # gca_en_car_int_0   38 89382 89675 -44653    89306  0.8791  1     0.3485    
  # gca_en_car_int_1   39 89365 89666 -44644    89287 18.1991  1  1.990e-05 ***
  # gca_en_car_int_2   40 89366 89674 -44643    89286  1.6095  1     0.2046    
  # gca_en_car_int_3   41 89346 89662 -44632    89264 21.8736  1  2.912e-06 ***
  
  
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
  # gca_en_stress_2   37 89381 89666 -44653    89307                     
  # gca_en_corsi_0    38 89382 89675 -44653    89306 0.0840  1     0.7720
  # gca_en_corsi_1    39 89384 89685 -44653    89306 0.3094  1     0.5781
  # gca_en_corsi_2    40 89385 89694 -44653    89305 0.8651  1     0.3523
  # gca_en_corsi_3    41 89386 89702 -44652    89304 1.3957  1     0.2374
  
  # add 3-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_corsi_int_0 <- update(gca_en_stress_2, . ~ . + use_std:stress_sum:corsi)
  gca_en_corsi_int_1 <- update(gca_en_corsi_int_0,   . ~ . + ot1:use_std:stress_sum:corsi)
  gca_en_corsi_int_2 <- update(gca_en_corsi_int_1,   . ~ . + ot2:use_std:stress_sum:corsi)
  gca_en_corsi_int_3 <- update(gca_en_corsi_int_2,   . ~ . + ot3:use_std:stress_sum:corsi)
  
  en_corsi_int_anova <-
    anova(gca_en_stress_2, gca_en_corsi_int_0, gca_en_corsi_int_1,
          gca_en_corsi_int_2, gca_en_corsi_int_3)
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
  # gca_en_stress_2      37 89381 89666 -44653    89307                       
  # gca_en_corsi_int_0   38 89382 89675 -44653    89306 0.8855  1    0.34670  
  # gca_en_corsi_int_1   39 89379 89680 -44650    89301 4.7471  1    0.02935 *
  # gca_en_corsi_int_2   40 89381 89689 -44650    89301 0.0033  1    0.95443  
  # gca_en_corsi_int_3   41 89380 89696 -44649    89298 3.0244  1    0.08202 .
  
  summary(gca_en_car_int_3)
  #                                 Estimate Std. Error t value
  # (Intercept)                     1.277432   0.101375  12.601
  # ot1                             5.644590   0.339415  16.630
  # ot2                             0.286581   0.281357   1.019
  # ot3                            -1.445775   0.231671  -6.241
  # use_std                         0.074815   0.077351   0.967
  # stress_sum                     -0.005595   0.075618  -0.074
  # ot1:use_std                     0.700269   0.268076   2.612
  # ot2:use_std                    -0.062302   0.182290  -0.342
  # ot3:use_std                    -0.387157   0.173037  -2.237
  # ot1:stress_sum                  0.381584   0.235141   1.623
  # ot2:stress_sum                  0.610242   0.180742   3.376
  # use_std:stress_sum:car_dev      0.343295   0.232777   1.475
  # ot1:use_std:stress_sum:car_dev  2.261142   0.573818   3.941
  # ot2:use_std:stress_sum:car_dev  0.290952   0.552673   0.526
  # ot3:use_std:stress_sum:car_dev -2.618416   0.558030  -4.692
  summary(gca_en_corsi_int_1)
  #                              Estimate Std. Error t value
  # (Intercept)                   1.27100    0.10119  12.560
  # ot1                           5.63977    0.33711  16.730
  # ot2                           0.27401    0.27723   0.988
  # ot3                          -1.43159    0.22794  -6.280
  # use_std                       0.07063    0.07768   0.909
  # stress_sum                   -0.03780    0.07473  -0.506
  # ot1:use_std                   0.69903    0.26672   2.621
  # ot2:use_std                  -0.07259    0.18114  -0.401
  # ot3:use_std                  -0.36573    0.17385  -2.104
  # ot1:stress_sum                0.21635    0.23336   0.927
  # ot2:stress_sum                0.47686    0.17953   2.656
  # use_std:stress_sum:corsi      0.04849    0.04291   1.130
  # ot1:use_std:stress_sum:corsi  0.23200    0.10627   2.183
  
}

# -----------------------------------------------------------------------------





# Random effects structure ----------------------------------------------------

ma_vision <- l2_vision %>%
  filter(., l1 == 'ma')


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
  # mod_ot1    9 91035 91104 -45508    91017                          
  # mod_ot2   14 91001 91109 -45486    90973  44.313  5  2.001e-08 ***
  # mod_ot3   20 90846 91000 -45403    90806 166.375  6  < 2.2e-16 ***
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df   AIC   BIC logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot3   20 90846 91000 -45403    90806                           
  # mod_ot4   21 90712 90874 -45335    90670 136.4335  1  < 2.2e-16 ***
  # mod_ot5   23 90532 90709 -45243    90486 184.2964  2  < 2.2e-16 ***
  # mod_ot6   26 90490 90690 -45219    90438  47.9140  3  2.221e-10 ***
  # mod_ot7   30 90489 90720 -45215    90429   8.6034  4    0.07181 .  
  
  
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
  
  # add L2 use effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_use_0 <- update(gca_ma_base,    . ~ . + use_std)
  gca_ma_use_1 <- update(gca_ma_use_0, . ~ . + ot1:use_std)
  gca_ma_use_2 <- update(gca_ma_use_1, . ~ . + ot2:use_std)
  gca_ma_use_3 <- update(gca_ma_use_2, . ~ . + ot3:use_std)
  
  ma_use_anova <-
    anova(gca_ma_base, gca_ma_use_0, gca_ma_use_1,
          gca_ma_use_2, gca_ma_use_3)
  #              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_ma_base    26 90490 90690 -45219    90438                       
  # gca_ma_use_0   27 90492 90700 -45219    90438 0.0295  1    0.86373  
  # gca_ma_use_1   28 90490 90706 -45217    90434 3.5751  1    0.05865 .
  # gca_ma_use_2   29 90488 90712 -45215    90430 3.8287  1    0.05038 .
  # gca_ma_use_3   30 90490 90721 -45215    90430 0.2737  1    0.60087   
  
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_stress_0 <- update(gca_ma_base,    . ~ . + stress_sum)
  gca_ma_stress_1 <- update(gca_ma_stress_0, . ~ . + ot1:stress_sum)
  gca_ma_stress_2 <- update(gca_ma_stress_1, . ~ . + ot2:stress_sum)
  gca_ma_stress_3 <- update(gca_ma_stress_2, . ~ . + ot3:stress_sum)
  
  ma_stress_anova <-
    anova(gca_ma_base, gca_ma_stress_0, gca_ma_stress_1,
          gca_ma_stress_2, gca_ma_stress_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_base       26 90490 90690 -45219    90438                       
  # gca_ma_stress_0   27 90486 90695 -45216    90432 5.1776  1    0.02288 *
  # gca_ma_stress_1   28 90488 90704 -45216    90432 0.0028  1    0.95743  
  # gca_ma_stress_2   29 90490 90714 -45216    90432 0.0224  1    0.88092  
  # gca_ma_stress_3   30 90491 90723 -45216    90431 1.0253  1    0.31126  
  
  
  # BRANCH #1
  # add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_car_0 <- update(gca_ma_stress_0,  . ~ . + car_dev)
  gca_ma_car_1 <- update(gca_ma_car_0, . ~ . + ot1:car_dev)
  gca_ma_car_2 <- update(gca_ma_car_1, . ~ . + ot2:car_dev)
  gca_ma_car_3 <- update(gca_ma_car_2, . ~ . + ot3:car_dev)
  
  ma_car_anova <-
    anova(gca_ma_stress_0, gca_ma_car_0, gca_ma_car_1,
          gca_ma_car_2, gca_ma_car_3)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_stress_0   27 90486 90695 -45216    90432                     
  # gca_ma_car_0      28 90488 90704 -45216    90432 0.0261  1     0.8718
  # gca_ma_car_1      29 90490 90714 -45216    90432 0.0371  1     0.8472
  # gca_ma_car_2      30 90492 90723 -45216    90432 0.2324  1     0.6298
  # gca_ma_car_3      31 90493 90732 -45216    90431 0.7797  1     0.3772
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_car_int_0 <- update(gca_ma_stress_0, . ~ . + use_std:stress_sum:car_dev)
  gca_ma_car_int_1 <- update(gca_ma_car_int_0,   . ~ . + ot1:use_std:stress_sum:car_dev)
  gca_ma_car_int_2 <- update(gca_ma_car_int_1,   . ~ . + ot2:use_std:stress_sum:car_dev)
  gca_ma_car_int_3 <- update(gca_ma_car_int_2,   . ~ . + ot3:use_std:stress_sum:car_dev)
  
  ma_car_int_anova <-
    anova(gca_ma_stress_0, gca_ma_car_int_0, gca_ma_car_int_1,
          gca_ma_car_int_2, gca_ma_car_int_3)
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_ma_stress_0    27 90486 90695 -45216    90432                     
  # gca_ma_car_int_0   28 90488 90704 -45216    90432 0.0436  1     0.8346
  # gca_ma_car_int_1   29 90490 90714 -45216    90432 0.4497  1     0.5025
  # gca_ma_car_int_2   30 90492 90723 -45216    90432 0.0353  1     0.8510
  # gca_ma_car_int_3   31 90493 90732 -45216    90431 0.5909  1     0.4421
  
  
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
  # gca_ma_stress_0   27 90486 90695 -45216    90432                     
  # gca_ma_corsi_0    28 90488 90704 -45216    90432 0.1239  1     0.7248
  # gca_ma_corsi_1    29 90489 90712 -45215    90431 1.8316  1     0.1759
  # gca_ma_corsi_2    30 90490 90721 -45215    90430 0.6741  1     0.4116
  # gca_ma_corsi_3    31 90491 90730 -45214    90429 0.8846  1     0.3469
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_corsi_int_0 <- update(gca_ma_stress_0, . ~ . + use_std:stress_sum:corsi)
  gca_ma_corsi_int_1 <- update(gca_ma_corsi_int_0,   . ~ . + ot1:use_std:stress_sum:corsi)
  gca_ma_corsi_int_2 <- update(gca_ma_corsi_int_1,   . ~ . + ot2:use_std:stress_sum:corsi)
  gca_ma_corsi_int_3 <- update(gca_ma_corsi_int_2,   . ~ . + ot3:use_std:stress_sum:corsi)
  
  ma_corsi_int_anova <-
    anova(gca_ma_stress_0, gca_ma_corsi_int_0, gca_ma_corsi_int_1,
          gca_ma_corsi_int_2, gca_ma_corsi_int_3)
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
  # gca_ma_stress_0      27 90486 90695 -45216    90432                     
  # gca_ma_corsi_int_0   28 90488 90704 -45216    90432 0.2331  1     0.6292
  # gca_ma_corsi_int_1   29 90490 90713 -45216    90432 0.6235  1     0.4297
  # gca_ma_corsi_int_2   30 90491 90723 -45216    90431 0.2841  1     0.5940
  # gca_ma_corsi_int_3   31 90493 90732 -45216    90431 0.0401  1     0.8413
  
  summary(gca_ma_stress_0)
  # Estimate Std. Error t value
  # (Intercept)  0.98380    0.09873   9.964
  # ot1          4.26708    0.32491  13.133
  # ot2          0.30587    0.21547   1.420
  # ot3         -0.96171    0.15002  -6.410
  # stress_sum  -0.12998    0.05563  -2.336
  
} 


  

mod_type <- "gca_en"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_car_0", "_car_1", "_car_2", "_car_3",
              "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              "_use_0", "_use_1", "_use_2", "_use_3",
              "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3")

# Store ind models in list
en_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(en_mods,
     file = here("mods", "vision", "gca", "continuous", "use",
                 "en_mods.Rdata"))


mod_type <- "gca_ma"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_car_0", "_car_1", "_car_2", "_car_3",
              "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              "_use_0", "_use_1", "_use_2", "_use_3",
              "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3")

# Store ind models in list
ma_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(ma_mods,
     file = here("mods", "vision", "gca", "continuous", "use",
                 "ma_mods.Rdata"))
  
  
  
  
# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions

en_car <- en_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, car_dev) %>%
  distinct %>%
  expand_grid(., tibble(use_std = c(-1, 0, 1)))

en_corsi <- en_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, corsi) %>%
  distinct %>%
  expand_grid(., tibble(use_std = c(-1, 0, 1)))

ma_car <- ma_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, car_dev) %>%
  distinct %>%
  expand_grid(., tibble(use_std = c(-1, 0, 1)))

ma_corsi <- ma_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, corsi) %>%
  distinct %>%
  expand_grid(., tibble(use_std = c(-1, 0, 1)))

# Get model predictions and SE

fits_en_corsi <- predictSE(gca_en_corsi_int_1, en_corsi) %>%        
  as_tibble %>%
  bind_cols(en_corsi) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

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

fits_ma_car <- predictSE(gca_ma_stress_0, ma_car) %>%        
  as_tibble %>%
  bind_cols(ma_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset

target_offset_preds_en_corsi <- filter(fits_en_corsi, time_zero == 4) %>%
  select(stress = stress_sum, corsi, use_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

target_offset_preds_en_car <- filter(fits_en_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, use_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 


target_offset_preds_ma_corsi <- filter(fits_ma_corsi, time_zero == 4) %>%
  select(stress = stress_sum, corsi, use_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

target_offset_preds_ma_car <- filter(fits_ma_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, use_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

# -----------------------------------------------------------------------------






# Save models -----------------------------------------------------------------

if(F) {
  
  # Save anova model comparisons
  nested_model_comparisons <-
    mget(c(
           "en_stress_anova", "en_car_anova", "en_corsi_anova",
           "en_use_anova", "en_car_int_anova", "en_corsi_int_anova",
           "ma_stress_anova", "ma_car_anova", "ma_corsi_anova",
           "ma_use_anova", "ma_car_int_anova", "ma_corsi_int_anova"
    ))
  
  save(nested_model_comparisons,
       file = here("mods", "vision", "gca", "continuous", "use",
                   "nested_model_comparisons.Rdata"))
  
  
  
  
  # Save models predictions
  model_preds <- mget(c(
                        "fits_en_corsi", 'fits_en_car',
                        "fits_ma_corsi", 'fits_ma_car',
                        "target_offset_preds_en_corsi", "target_offset_preds_en_car",
                        "target_offset_preds_ma_corsi", "target_offset_preds_ma_car"))
  
  save(model_preds,
       file = here("mods", "vision", "gca", "continuous", "use",
                   "model_preds.Rdata"))
  
}

# -----------------------------------------------------------------------------

