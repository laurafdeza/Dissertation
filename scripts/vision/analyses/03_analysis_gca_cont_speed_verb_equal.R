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
vision <- read_csv("./data/clean/vision_scores.csv") # pred car
corsi <- read_csv("./data/clean/corsi_z_scores.csv")

vision50 <- left_join(x = stress50, y = wm, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = vision, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = corsi, by = "participant", all.x=TRUE)

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
  #         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    9 22390 22447 -11186    22372                          
  # mod_ot2   14 22378 22466 -11175    22350 22.6051  5  0.0004016 ***
  # mod_ot3   20 22387 22513 -11173    22347  2.9911  6  0.8099630    
  
  mod_ot4 <- update(mod_ot2, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + ot3 + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot2, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df   AIC   BIC logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot2   14 22378 22466 -11175    22350                          
  # mod_ot4   15 22288 22383 -11129    22258 91.0894  1  < 2.2e-16 ***
  # mod_ot5   17 22271 22378 -11118    22237 21.4141  2  2.239e-05 ***
  # mod_ot6   20 22272 22398 -11116    22232  4.6564  3     0.1988    
  # mod_ot7   25 22279 22437 -11115    22229  3.1346  5     0.6792  
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_mon_base <- mod_ot5
  # lmer(eLog ~ 1 + (ot1 + ot2) +          
  #        (1 + stress_sum + (ot1 + ot2) | participant) +
  #        (1 + ot1 | target),
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 3e5)),
  #      data = stress_gc_subset, REML = F)    # , na.action = na.exclude 

# add stress effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_stress_0 <- update(gca_mon_base,    . ~ . + stress_sum)
gca_mon_stress_1 <- update(gca_mon_stress_0, . ~ . + ot1:stress_sum)
gca_mon_stress_2 <- update(gca_mon_stress_1, . ~ . + ot2:stress_sum)

mon_stress_anova <-
  anova(gca_mon_base, gca_mon_stress_0, gca_mon_stress_1,
        gca_mon_stress_2)
                 # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base       17 22271 22378 -11118    22237                     
# gca_mon_stress_0   18 22273 22386 -11118    22237 0.0152  1     0.9017
# gca_mon_stress_1   19 22274 22394 -11118    22236 1.2699  1     0.2598
# gca_mon_stress_2   20 22276 22402 -11118    22236 0.2534  1     0.6147

# BRANCH #0
# add verbal proc time effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_0 <- update(gca_mon_base,    . ~ . + ospan_rt)
gca_mon_ospan_1 <- update(gca_mon_ospan_0, . ~ . + ot1:ospan_rt)
gca_mon_ospan_2 <- update(gca_mon_ospan_1, . ~ . + ot2:ospan_rt)

mon_ospan_anova <-
  anova(gca_mon_base, gca_mon_ospan_0, gca_mon_ospan_1,
        gca_mon_ospan_2)
#                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      17 22271 22378 -11118    22237                     
# gca_mon_ospan_0   18 22273 22386 -11118    22237 0.1179  1     0.7314
# gca_mon_ospan_1   19 22275 22395 -11118    22237 0.0061  1     0.9376
# gca_mon_ospan_2   20 22277 22403 -11118    22237 0.3223  1     0.5702


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_int_0 <- update(gca_mon_base, . ~ . + stress_sum:ospan_rt)
gca_mon_ospan_int_1 <- update(gca_mon_ospan_int_0,   . ~ . + ot1:stress_sum:ospan_rt)
gca_mon_ospan_int_2 <- update(gca_mon_ospan_int_1,   . ~ . + ot2:stress_sum:ospan_rt)

mon_ospan_int_anova <-
  anova(gca_mon_base, gca_mon_ospan_int_0, gca_mon_ospan_int_1,
        gca_mon_ospan_int_2)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mon_base          17 22271 22378 -11118    22237                        
# gca_mon_ospan_int_0   18 22273 22386 -11118    22237 0.1752  1   0.675490   
# gca_mon_ospan_int_1   19 22272 22392 -11117    22234 2.7569  1   0.096833 . 
# gca_mon_ospan_int_2   20 22265 22391 -11113    22225 9.0521  1   0.002624 **


# BRANCH #1
# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_0 <- update(gca_mon_base,    . ~ . + car_dev)
gca_mon_car_1 <- update(gca_mon_car_0, . ~ . + ot1:car_dev)
gca_mon_car_2 <- update(gca_mon_car_1, . ~ . + ot2:car_dev)

mon_car_anova <-
  anova(gca_mon_base, gca_mon_car_0, gca_mon_car_1,
        gca_mon_car_2)
#               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base    17 22271 22378 -11118    22237                          
# gca_mon_car_0   18 22273 22386 -11118    22237  0.0326  1  0.8567766    
# gca_mon_car_1   19 22286 22405 -11124    22248  0.0000  1  1.0000000    
# gca_mon_car_2   20 22276 22402 -11118    22236 11.5099  1  0.0006922 ***


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_int_0 <- update(gca_mon_car_2, . ~ . + stress_sum:car_dev)
gca_mon_car_int_1 <- update(gca_mon_car_int_0,   . ~ . + ot1:stress_sum:car_dev)
gca_mon_car_int_2 <- update(gca_mon_car_int_1,   . ~ . + ot2:stress_sum:car_dev)

mon_car_int_anova <-
  anova(gca_mon_car_2, gca_mon_car_int_0, gca_mon_car_int_1,
        gca_mon_car_int_2)
#                   npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)  
# gca_mon_car_2       20 22276 22402 -11118    22236                          
# gca_mon_car_int_0   21 22278 22410 -11118    22236  0.1911  1  0.6619903    
# gca_mon_car_int_1   22 22290 22428 -11123    22246  0.0000  1  1.0000000    
# gca_mon_car_int_2   23 22278 22422 -11116    22232 13.9467  1  0.0001881 ***

# BRANCH #2
# add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_0 <- update(gca_mon_base,    . ~ . + corsi)
gca_mon_corsi_1 <- update(gca_mon_corsi_0, . ~ . + ot1:corsi)
gca_mon_corsi_2 <- update(gca_mon_corsi_1, . ~ . + ot2:corsi)

mon_corsi_anova <-
  anova(gca_mon_base, gca_mon_corsi_0, gca_mon_corsi_1,
        gca_mon_corsi_2)
#                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      


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

mon_corsirt_anova <-
  anova(gca_mon_base, gca_mon_corsirt_0, gca_mon_corsirt_1,
        gca_mon_corsirt_2)
#                   npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)
# gca_mon_base        17 22271 22378 -11118    22237                          
# gca_mon_corsirt_0   18 22284 22397 -11124    22248  0.0000  1  1.0000000    
# gca_mon_corsirt_1   19 22273 22393 -11118    22235 12.2818  1  0.0004574 ***
# gca_mon_corsirt_2   20 22275 22401 -11117    22235  0.6368  1  0.4248541    

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_int_0 <- update(gca_mon_corsirt_1, . ~ . + stress_sum:corsi_rt)
gca_mon_corsirt_int_1 <- update(gca_mon_corsirt_int_0,   . ~ . + ot1:stress_sum:corsi_rt)
gca_mon_corsirt_int_2 <- update(gca_mon_corsirt_int_1,   . ~ . + ot2:stress_sum:corsi_rt)

mon_corsirt_int_anova <-
  anova(gca_mon_corsirt_1, gca_mon_corsirt_int_0, gca_mon_corsirt_int_1,
        gca_mon_corsirt_int_2)
#                       npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_corsirt_1       19 22273 22393 -11118    22235                       
# gca_mon_corsirt_int_0   20 22275 22401 -11118    22235 0.0313  1    0.85964  
# gca_mon_corsirt_int_1   21 22276 22408 -11117    22234 1.9184  1    0.16603  
# gca_mon_corsirt_int_2   22 22272 22411 -11114    22228 5.4551  1    0.01951 *



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
              # "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3",
              "_car_int_0", "_car_int_1", "_car_int_2" 
              # "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", "_corsi_int_3"
              )

# Store ind models in list
mon_mods_equal <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods_equal,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "mon_mods_equal.Rdata"))

  


  
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






# Save models predictions
model_preds <- mget(c("fits_mon_corsirt", 'fits_mon_ospan', 'fits_mon_car',
                      "preds_mon_corsi", "preds_mon_ospan", 'preds_mon_car'
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
           "mon_car_int_anova", "mon_corsi_int_anova"
    ))
  
  save(nested_model_comparisons,
       file = here("mods", "vision", "gca", "continuous",
                   "nested_model_comparisons.Rdata"))
  
  
  
  
  
  
}

# -----------------------------------------------------------------------------

