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

vision <- read_csv("./data/clean/vision_scores.csv")
corsi <- read_csv("./data/clean/corsi_z_scores.csv")

visuospatial_df <- left_join(x = vision, y = corsi, by = "participant", all.x=TRUE)

vision50 <- left_join(x = stress50, y = visuospatial_df, by = "participant", all.x=TRUE)

vision50 <- na.omit(vision50)




vision50 <- vision50 %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma"),
            stress_sum = if_else(cond == "1", 1, -1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

mon_vision <- filter(vision50, l1 == 'es') %>% select(-DELE, -percent_l2_week, 
                                                    -DELE_z, -use_z, -prof, -group)

l2_vision <- vision50 %>%
  filter(., l1 != 'es') %>% 
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1),
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
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    9 29641 29700 -14811    29623                         
  # mod_ot2   14 29557 29650 -14765    29529 93.395  5  < 2.2e-16 ***
  # mod_ot3   20 29496 29628 -14728    29456 73.374  6  8.292e-14 ***
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df    AIC    BIC  logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot1    9 240042 240120 -120012   240024
  # mod_ot2   14 239367 239489 -119670   239339 684.486  5  < 2.2e-16 ***
  # mod_ot3   20 239311 239484 -119635   239271  68.534  6  8.167e-13 ***
  # mod_ot4   21 239043 239225 -119501   239001 269.485  1  < 2.2e-16 ***
  
  
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
# gca_mon_base       30 29135 29333 -14537    29075                     
# gca_mon_stress_0   31 29136 29341 -14537    29074 0.2798  1     0.5968
# gca_mon_stress_1   32 29137 29348 -14536    29073 1.4592  1     0.2271
# gca_mon_stress_2   33 29139 29357 -14536    29073 0.0328  1     0.8564
# gca_mon_stress_3   34 29140 29365 -14536    29072 0.5505  1     0.4581


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
# gca_mon_base    30 29135 29333 -14537    29075                     
# gca_mon_car_0   31 29137 29341 -14537    29075 0.0203  1     0.8867
# gca_mon_car_1   32 29138 29349 -14537    29074 0.4666  1     0.4946
# gca_mon_car_2   33 29140 29358 -14537    29074 0.4657  1     0.4950
# gca_mon_car_3   34 29141 29365 -14536    29073 0.9468  1     0.3305

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_int_0 <- update(gca_mon_base, . ~ . + stress_sum:car_dev)
gca_mon_car_int_1 <- update(gca_mon_car_int_0,   . ~ . + ot1:stress_sum:car_dev)
gca_mon_car_int_2 <- update(gca_mon_car_int_1,   . ~ . + ot2:stress_sum:car_dev)
gca_mon_car_int_3 <- update(gca_mon_car_int_2,   . ~ . + ot3:stress_sum:car_dev)

mon_car_int_anova <-
  anova(gca_mon_base, gca_mon_car_int_0, gca_mon_car_int_1,
        gca_mon_car_int_2, gca_mon_car_int_3)
#                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_full_base        30 29135 29333 -14537    29075                       
# gca_mon_car_int_0   31 29134 29339 -14536    29072 2.2782  1     0.1312  
# gca_mon_car_int_1   32 29133 29344 -14534    29069 3.6023  1     0.0577 .
# gca_mon_car_int_2   33 29133 29351 -14534    29067 1.7743  1     0.1828  
# gca_mon_car_int_3   34 29133 29358 -14533    29065 1.8393  1     0.1750  


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
# gca_mon_base      30 29135 29333 -14537    29075                     
# gca_mon_corsi_0   31 29136 29340 -14537    29074 0.9756  1     0.3233
# gca_mon_corsi_1   32 29138 29349 -14537    29074 0.0591  1     0.8080
# gca_mon_corsi_2   33 29139 29357 -14536    29073 0.7939  1     0.3729
# gca_mon_corsi_3   34 29141 29365 -14536    29073 0.0787  1     0.7790

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_int_0 <- update(gca_mon_base, . ~ . + stress_sum:corsi)
gca_mon_corsi_int_1 <- update(gca_mon_corsi_int_0,   . ~ . + ot1:stress_sum:corsi)
gca_mon_corsi_int_2 <- update(gca_mon_corsi_int_1,   . ~ . + ot2:stress_sum:corsi)
gca_mon_corsi_int_3 <- update(gca_mon_corsi_int_2,   . ~ . + ot3:stress_sum:corsi)

mon_corsi_int_anova <-
  anova(gca_mon_base, gca_mon_corsi_int_0, gca_mon_corsi_int_1,
        gca_mon_corsi_int_2, gca_mon_corsi_int_3)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_base          30 29135 29333 -14537    29075                       
# gca_mon_corsi_int_0   31 29137 29341 -14537    29075 0.0020  1    0.96453  
# gca_mon_corsi_int_1   32 29134 29346 -14535    29070 4.3148  1    0.03778 *
# gca_mon_corsi_int_2   33 29136 29354 -14535    29070 0.2038  1    0.65170  
# gca_mon_corsi_int_3   34 29135 29359 -14534    29067 3.1805  1    0.07452 .

summary(gca_mon_base)
# Estimate Std. Error t value
# (Intercept)   1.6442     0.1812   9.072
# ot1           4.2727     0.5517   7.745
# ot2          -1.3490     0.5011  -2.692
# ot3          -0.7122     0.3426  -2.079
summary(gca_mon_corsi_int_1)
# Estimate Std. Error t value
# (Intercept)           1.64545    0.18068   9.107
# ot1                   4.28638    0.55098   7.780
# ot2                  -1.34921    0.50133  -2.691
# ot3                  -0.70364    0.34364  -2.048
# stress_sum:corsi      0.01145    0.07841   0.146
# ot1:stress_sum:corsi  0.31797    0.15227   2.088

}

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
  #      data = stress_gc_subset, REML = F)    # , na.action = na.exclude 
  
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
  # gca_en_stress_2   36 89377 89654 -44652    89305                       
  # gca_en_car_0      37 89379 89664 -44652    89305 0.1490  1    0.69948  
  # gca_en_car_1      38 89381 89673 -44652    89305 0.0519  1    0.81973  
  # gca_en_car_2      39 89382 89682 -44652    89304 0.9988  1    0.31760  
  # gca_en_car_3      40 89379 89688 -44650    89299 4.0381  1    0.04448 *
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_car_int_0 <- update(gca_en_car_3, . ~ . + prof_std:stress_sum:car_dev)
  gca_en_car_int_1 <- update(gca_en_car_int_0,   . ~ . + ot1:prof_std:stress_sum:car_dev)
  gca_en_car_int_2 <- update(gca_en_car_int_1,   . ~ . + ot2:prof_std:stress_sum:car_dev)
  gca_en_car_int_3 <- update(gca_en_car_int_2,   . ~ . + ot3:prof_std:stress_sum:car_dev)
  
  en_car_int_anova <-
    anova(gca_en_car_3, gca_en_car_int_0, gca_en_car_int_1,
          gca_en_car_int_2, gca_en_car_int_3)
  #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_en_car_3       40 89379 89688 -44650    89299                        
  # gca_en_car_int_0   41 89378 89694 -44648    89296 3.7205  1   0.053748 . 
  # gca_en_car_int_1   42 89373 89696 -44644    89289 7.2424  1   0.007120 **
  # gca_en_car_int_2   43 89374 89706 -44644    89288 0.1062  1   0.744464   
  # gca_en_car_int_3   44 89370 89709 -44641    89282 6.7112  1   0.009581 **
  
  
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
  # gca_en_stress_2   36 89377 89654 -44652    89305                     
  # gca_en_corsi_0    37 89379 89664 -44652    89305 0.1060  1     0.7448
  # gca_en_corsi_1    38 89380 89673 -44652    89304 0.3695  1     0.5433
  # gca_en_corsi_2    39 89382 89682 -44652    89304 0.6930  1     0.4051
  # gca_en_corsi_3    40 89382 89691 -44651    89302 1.3230  1     0.2501
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_corsi_int_0 <- update(gca_en_stress_2, . ~ . + prof_std:stress_sum:corsi)
  gca_en_corsi_int_1 <- update(gca_en_corsi_int_0,   . ~ . + ot1:prof_std:stress_sum:corsi)
  gca_en_corsi_int_2 <- update(gca_en_corsi_int_1,   . ~ . + ot2:prof_std:stress_sum:corsi)
  gca_en_corsi_int_3 <- update(gca_en_corsi_int_2,   . ~ . + ot3:prof_std:stress_sum:corsi)
  
  en_corsi_int_anova <-
    anova(gca_en_stress_2, gca_en_corsi_int_0, gca_en_corsi_int_1,
          gca_en_corsi_int_2, gca_en_corsi_int_3)
  #                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq) 
  # gca_en_stress_2      36 89377 89654 -44652    89305                     
  # gca_en_corsi_int_0   37 89378 89664 -44652    89304 0.3781  1     0.5386
  # gca_en_corsi_int_1   38 89380 89673 -44652    89304 0.2526  1     0.6152
  # gca_en_corsi_int_2   39 89379 89680 -44651    89301 2.6225  1     0.1054
  # gca_en_corsi_int_3   40 89379 89688 -44650    89299 1.9862  1     0.1587
  
  summary(gca_en_car_int_3)
  # Estimate Std. Error t value
  # (Intercept)   1.272282   0.099738  12.756
  # ot1                              5.535813   0.338411  16.358
  # ot2                              0.268889   0.270682   0.993
  # ot3                             -1.420396   0.226300  -6.277
  # prof_std                         0.088586   0.058920   1.504
  # stress_sum                      -0.034351   0.074103  -0.464
  # car_dev                         -0.304621   0.279246  -1.091
  # ot1:prof_std                     0.444011   0.211643   2.098
  # ot2:prof_std                    -0.328686   0.143594  -2.289
  # ot1:stress_sum                   0.211281   0.233730   0.904
  # ot2:stress_sum                   0.504341   0.177028   2.849
  # ot1:car_dev                     -0.373806   0.998171  -0.374
  # ot2:car_dev                      0.629074   0.637154   0.987
  # ot3:car_dev                      1.404928   0.634363   2.215
  # prof_std:stress_sum:car_dev      0.346966   0.145858   2.379
  # ot1:prof_std:stress_sum:car_dev  0.905886   0.358009   2.530
  # ot2:prof_std:stress_sum:car_dev -0.007305   0.349868  -0.021
  # ot3:prof_std:stress_sum:car_dev -0.935750   0.360172  -2.598
  summary(gca_en_stress_2)
  # Estimate Std. Error t value
  # (Intercept)           1.25761    0.09927  12.669
  # ot1             5.52541    0.33505  16.491
  # ot2             0.28947    0.26966   1.073
  # ot3            -1.36833    0.22847  -5.989
  # prof_std        0.09103    0.05914   1.539
  # stress_sum     -0.03643    0.07484  -0.487
  # ot1:prof_std    0.44794    0.21190   2.114
  # ot2:prof_std   -0.33260    0.14404  -2.309
  # ot1:stress_sum  0.20675    0.23327   0.886
  # ot2:stress_sum  0.48724    0.17797   2.738
  
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
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
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
  #           Df    AIC    BIC  logLik deviance   Chisq Chi Pr(>Chisq)
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
  
  # add proficimacy effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_prof_0 <- update(gca_ma_base,    . ~ . + prof_std)
  gca_ma_prof_1 <- update(gca_ma_prof_0, . ~ . + ot1:prof_std)
  gca_ma_prof_2 <- update(gca_ma_prof_1, . ~ . + ot2:prof_std)
  gca_ma_prof_3 <- update(gca_ma_prof_2, . ~ . + ot3:prof_std)
  
  ma_prof_anova <-
    anova(gca_ma_base, gca_ma_prof_0, gca_ma_prof_1,
          gca_ma_prof_2, gca_ma_prof_3)
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_ma_base     26 90490 90690 -45219    90438                       
  # gca_ma_prof_0   27 90486 90695 -45216    90432 5.2713  1    0.02168 *
  # gca_ma_prof_1   28 90486 90701 -45215    90430 2.7687  1    0.09612 .
  # gca_ma_prof_2   29 90487 90711 -45215    90429 0.2154  1    0.64258  
  # gca_ma_prof_3   30 90489 90720 -45214    90429 0.4164  1    0.51874  
  
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_stress_0 <- update(gca_ma_prof_0,    . ~ . + stress_sum)
  gca_ma_stress_1 <- update(gca_ma_stress_0, . ~ . + ot1:stress_sum)
  gca_ma_stress_2 <- update(gca_ma_stress_1, . ~ . + ot2:stress_sum)
  gca_ma_stress_3 <- update(gca_ma_stress_2, . ~ . + ot3:stress_sum)
  
  ma_stress_anova <-
    anova(gca_ma_prof_0, gca_ma_stress_0, gca_ma_stress_1,
          gca_ma_stress_2, gca_ma_stress_3)
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_prof_0     27 90486 90695 -45216    90432                       
  # gca_ma_stress_0   28 90483 90699 -45214    90427 5.2262  1    0.02225 *
  # gca_ma_stress_1   29 90485 90709 -45214    90427 0.0066  1    0.93518  
  # gca_ma_stress_2   30 90487 90718 -45214    90427 0.0236  1    0.87779  
  # gca_ma_stress_3   31 90488 90727 -45213    90426 0.9928  1    0.31905  
  
  
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
  
  summary(gca_ma_car_int_2)
  # Estimate Std. Error t value
  # (Intercept)   0.98297    0.09497  10.350
  # ot1                              4.27198    0.32432  13.172
  # ot2                              0.30751    0.21389   1.438
  # ot3                             -0.95600    0.14978  -6.382
  # prof_std                         0.14404    0.05932   2.428
  # stress_sum                      -0.12947    0.05555  -2.331
  # prof_std:stress_sum:car_dev     -0.09301    0.15706  -0.592
  # ot1:prof_std:stress_sum:car_dev -0.51841    0.38511  -1.346
  summary(gca_ma_stress_0)
  # Estimate Std. Error t value
  # (Intercept)  0.98253    0.09501  10.341
  # ot1          4.26435    0.32529  13.109
  # ot2          0.30544    0.21540   1.418
  # ot3         -0.96040    0.15039  -6.386
  # prof_std     0.14454    0.05930   2.437
  # stress_sum  -0.13061    0.05567  -2.346
  
} 
  
mod_type <- "gca_mon"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_car_0", "_car_1", "_car_2", "_car_3",
              "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3")

# Store ind models in list
mon_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods,
     file = here("mods", "vision", "gca", "continuous",
                 "mon_mods.Rdata"))

  

mod_type <- "gca_en"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_car_0", "_car_1", "_car_2", "_car_3",
              "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              "_prof_0", "_prof_1", "_prof_2", "_prof_3",
              "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3")

# Store ind models in list
en_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(en_mods,
     file = here("mods", "vision", "gca", "continuous",
                 "en_mods.Rdata"))


mod_type <- "gca_ma"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_car_0", "_car_1", "_car_2", "_car_3",
              "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              "_prof_0", "_prof_1", "_prof_2", "_prof_3",
              "_car_int_0", "_car_int_1", "_car_int_2", "_car_int_3",
              "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3")

# Store ind models in list
ma_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(ma_mods,
     file = here("mods", "vision", "gca", "continuous",
                 "ma_mods.Rdata"))
  
  
  
  
# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_car <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  # mutate(l1_sum = as.character(l1_sum)) %>% 
  expand_grid(., tibble(car_dev = c(-1, 0, 1)))
  
mon_corsi <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  # mutate(l1_sum = as.character(l1_sum)) %>% 
  expand_grid(., tibble(corsi = c(-1, 0, 1)))

en_car <- en_vision %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(car_dev = c(-1, 0, 1)))

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

# Get model predictions and SE
fits_mon_corsi <- predictSE(gca_mon_corsi_int_1, mon_corsi) %>%        
  as_tibble %>%
  bind_cols(mon_corsi) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_mon_car <- predictSE(gca_mon_base, mon_car) %>%        
  as_tibble %>%
  bind_cols(mon_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


fits_en_corsi <- predictSE(gca_en_stress_2, en_corsi) %>%        
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

fits_ma_car <- predictSE(gca_ma_car_int_2, ma_car) %>%        
  as_tibble %>%
  bind_cols(ma_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
target_offset_preds_mon_corsi <- filter(fits_mon_corsi, time_zero == 4) %>%
  select(stress = stress_sum, corsi,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

target_offset_preds_mon_car <- filter(fits_mon_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

target_offset_preds_en_corsi <- filter(fits_en_corsi, time_zero == 4) %>%
  select(stress = stress_sum, corsi, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

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
  
  
  
  
  # Save models predictions
  model_preds <- mget(c("fits_mon_corsi", 'fits_mon_car', 
                        "fits_en_corsi", 'fits_en_car',
                        "fits_ma_corsi", 'fits_ma_car',
                        "target_offset_preds_mon_corsi", "target_offset_preds_mon_car",
                        "target_offset_preds_en_corsi", "target_offset_preds_en_car",
                        "target_offset_preds_ma_corsi", "target_offset_preds_ma_car"))
  
  save(model_preds,
       file = here("mods", "vision", "gca", "continuous",
                   "model_preds.Rdata"))
  
}

# -----------------------------------------------------------------------------

