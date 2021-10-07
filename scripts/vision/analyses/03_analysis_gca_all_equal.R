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
  #         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    9 15237 15290 -7609.3    15219                        
  # mod_ot2   14 15232 15315 -7602.1    15204 14.3462  5    0.01355 *
  # mod_ot3   20 15239 15358 -7599.7    15199  4.9433  6    0.55110    
  
  mod_ot4 <- update(mod_ot2, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + ot3 + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot2, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           Df   AIC   BIC logLik deviance   Chisq Chi Pr(>Chisq)
  # mod_ot2   14 15232 15315 -7602.1    15204                           
  # mod_ot4   15 15116 15205 -7543.3    15086 117.7395  1    < 2e-16 ***
  # mod_ot5   17 15102 15202 -7534.0    15068  18.5875  2    9.2e-05 ***
  # mod_ot6   20 15103 15222 -7531.6    15063   4.7474  3     0.1913    
  # mod_ot7   25 15110 15258 -7530.0    15060   3.1294  5     0.6800    
  
  
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
                 # npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base       17 15102 15202 -7534.0    15068                     
# gca_mon_stress_0   18 15104 15210 -7533.9    15068 0.2303  1     0.6313
# gca_mon_stress_1   19 15103 15216 -7532.6    15065 2.5649  1     0.1093
# gca_mon_stress_2   20 15105 15223 -7532.5    15065 0.1090  1     0.7413

# BRANCH #0
# add verbal proc time effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_0 <- update(gca_mon_base,    . ~ . + ospan_rt)
gca_mon_ospan_1 <- update(gca_mon_ospan_0, . ~ . + ot1:ospan_rt)
gca_mon_ospan_2 <- update(gca_mon_ospan_1, . ~ . + ot2:ospan_rt)

mon_ospan_anova <-
  anova(gca_mon_base, gca_mon_ospan_0, gca_mon_ospan_1,
        gca_mon_ospan_2)
#                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      17 15102 15202 -7534.0    15068                     
# gca_mon_ospan_0   18 15103 15210 -7533.6    15067 0.6800  1     0.4096
# gca_mon_ospan_1   19 15105 15218 -7533.6    15067 0.0583  1     0.8092
# gca_mon_ospan_2   20 15106 15225 -7533.2    15066 0.7430  1     0.3887


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospan_int_0 <- update(gca_mon_base, . ~ . + stress_sum:ospan_rt)
gca_mon_ospan_int_1 <- update(gca_mon_ospan_int_0,   . ~ . + ot1:stress_sum:ospan_rt)
gca_mon_ospan_int_2 <- update(gca_mon_ospan_int_1,   . ~ . + ot2:stress_sum:ospan_rt)

mon_ospan_int_anova <-
  anova(gca_mon_base, gca_mon_ospan_int_0, gca_mon_ospan_int_1,
        gca_mon_ospan_int_2)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mon_base          17 15102 15202 -7534.0    15068                        
# gca_mon_ospan_int_0   18 15103 15210 -7533.7    15067 0.4923  1   0.482924   
# gca_mon_ospan_int_1   19 15104 15217 -7533.2    15066 0.9496  1   0.329817   
# gca_mon_ospan_int_2   20 15100 15218 -7529.8    15060 6.9563  1   0.008352 **


# BRANCH #1
# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_0 <- update(gca_mon_base,    . ~ . + car_dev)
gca_mon_car_1 <- update(gca_mon_car_0, . ~ . + ot1:car_dev)
gca_mon_car_2 <- update(gca_mon_car_1, . ~ . + ot2:car_dev)

mon_car_anova <-
  anova(gca_mon_base, gca_mon_car_0, gca_mon_car_1,
        gca_mon_car_2)
#               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base    17 15102 15202 -7534.0    15068                     
# gca_mon_car_0   18 15104 15210 -7533.9    15068 0.1521  1     0.6966
# gca_mon_car_1   19 15106 15218 -7533.9    15068 0.0633  1     0.8014
# gca_mon_car_2   20 15108 15226 -7533.8    15068 0.0394  1     0.8428


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_int_0 <- update(gca_mon_base, . ~ . + stress_sum:car_dev)
gca_mon_car_int_1 <- update(gca_mon_car_int_0,   . ~ . + ot1:stress_sum:car_dev)
gca_mon_car_int_2 <- update(gca_mon_car_int_1,   . ~ . + ot2:stress_sum:car_dev)

mon_car_int_anova <-
  anova(gca_mon_base, gca_mon_car_int_0, gca_mon_car_int_1,
        gca_mon_car_int_2)
#                   npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)  
# gca_mon_base        17 15102 15202 -7534.0    15068                     
# gca_mon_car_int_0   18 15104 15210 -7534.0    15068 0.0113  1     0.9155
# gca_mon_car_int_1   19 15105 15218 -7533.7    15067 0.5652  1     0.4522
# gca_mon_car_int_2   20 15107 15225 -7533.6    15067 0.2103  1     0.6466

# BRANCH #2
# add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_0 <- update(gca_mon_base,    . ~ . + corsi)
gca_mon_corsi_1 <- update(gca_mon_corsi_0, . ~ . + ot1:corsi)
gca_mon_corsi_2 <- update(gca_mon_corsi_1, . ~ . + ot2:corsi)

mon_corsi_anova <-
  anova(gca_mon_base, gca_mon_corsi_0, gca_mon_corsi_1,
        gca_mon_corsi_2)
#                 npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base      17 15102 15202 -7534.0    15068                       
# gca_mon_corsi_0   18 15102 15209 -7533.2    15066 1.4445  1    0.22941  
# gca_mon_corsi_1   19 15104 15217 -7533.1    15066 0.2387  1    0.62514  
# gca_mon_corsi_2   20 15100 15218 -7530.0    15060 6.2299  1    0.01256 *


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_int_0 <- update(gca_mon_corsi_2, . ~ . + stress_sum:corsi)
gca_mon_corsi_int_1 <- update(gca_mon_corsi_int_0,   . ~ . + ot1:stress_sum:corsi)
gca_mon_corsi_int_2 <- update(gca_mon_corsi_int_1,   . ~ . + ot2:stress_sum:corsi)
 
mon_corsi_int_anova <-
  anova(gca_mon_corsi_2, gca_mon_corsi_int_0, gca_mon_corsi_int_1,
        gca_mon_corsi_int_2)
#                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_corsi_2       20 15100 15218 -7530.0    15060                     
# gca_mon_corsi_int_0   21 15102 15226 -7530.0    15060 0.0005  1     0.9826
# gca_mon_corsi_int_1   22 15103 15233 -7529.6    15059 0.8358  1     0.3606
# gca_mon_corsi_int_2   23 15103 15239 -7528.3    15057 2.5995  1     0.1069   


# BRANCH #3
# add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_0 <- update(gca_mon_base,    . ~ . + corsi_rt)
gca_mon_corsirt_1 <- update(gca_mon_corsirt_0, . ~ . + ot1:corsi_rt)
gca_mon_corsirt_2 <- update(gca_mon_corsirt_1, . ~ . + ot2:corsi_rt)

mon_corsirt_anova <-
  anova(gca_mon_base, gca_mon_corsirt_0, gca_mon_corsirt_1,
        gca_mon_corsirt_2)
#                   npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)
# gca_mon_base        17 15102 15202 -7534.0    15068                       
# gca_mon_corsirt_0   18 15100 15207 -7532.1    15064 3.8097  1    0.05096 .
# gca_mon_corsirt_1   19 15100 15212 -7530.8    15062 2.4805  1    0.11526  
# gca_mon_corsirt_2   20 15102 15220 -7530.8    15062 0.1238  1    0.72495  

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsirt_int_0 <- update(gca_mon_base, . ~ . + stress_sum:corsi_rt)
gca_mon_corsirt_int_1 <- update(gca_mon_corsirt_int_0,   . ~ . + ot1:stress_sum:corsi_rt)
gca_mon_corsirt_int_2 <- update(gca_mon_corsirt_int_1,   . ~ . + ot2:stress_sum:corsi_rt)

mon_corsirt_int_anova <-
  anova(gca_mon_base, gca_mon_corsirt_int_0, gca_mon_corsirt_int_1,
        gca_mon_corsirt_int_2)
#                       npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_base            17 15102 15202 -7534.0    15068                        
# gca_mon_corsirt_int_0   18 15104 15210 -7533.9    15068 0.1278  1   0.720707   
# gca_mon_corsirt_int_1   19 15106 15218 -7533.9    15068 0.0128  1   0.909785   
# gca_mon_corsirt_int_2   20 15099 15217 -7529.6    15059 8.6361  1   0.003296 **





# BRANCH #4
# add verbal wm score effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospanscore_0 <- update(gca_mon_base,    . ~ . + ospan)
gca_mon_ospanscore_1 <- update(gca_mon_ospanscore_0, . ~ . + ot1:ospan)
gca_mon_ospanscore_2 <- update(gca_mon_ospanscore_1, . ~ . + ot2:ospan)

mon_ospanscore_anova <-
  anova(gca_mon_base, gca_mon_ospanscore_0, gca_mon_ospanscore_1,
        gca_mon_ospanscore_2)
#                      npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base           17 15102 15202 -7534.0    15068                     
# gca_mon_ospanscore_0   18 15104 15210 -7533.9    15068 0.0949  1     0.7581
# gca_mon_ospanscore_1   19 15106 15218 -7533.9    15068 0.0274  1     0.8686
# gca_mon_ospanscore_2   20 15108 15226 -7533.8    15068 0.1239  1     0.7248

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_ospanscore_int_0 <- update(gca_mon_base, . ~ . + stress_sum:ospan)
gca_mon_ospanscore_int_1 <- update(gca_mon_ospanscore_int_0,   . ~ . + ot1:stress_sum:ospan)
gca_mon_ospanscore_int_2 <- update(gca_mon_ospanscore_int_1,   . ~ . + ot2:stress_sum:ospan)

mon_ospanscore_int_anova <-
  anova(gca_mon_base, gca_mon_ospanscore_int_0, gca_mon_ospanscore_int_1,
        gca_mon_ospanscore_int_2)
#                          npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq) 
# gca_mon_base               17 15102 15202 -7534.0    15068                        
# gca_mon_ospanscore_int_0   18 15104 15210 -7533.9    15068 0.1144  1    0.73524  
# gca_mon_ospanscore_int_1   19 15103 15215 -7532.4    15065 2.9939  1    0.08358 .
# gca_mon_ospanscore_int_2   20 15105 15223 -7532.4    15065 0.0980  1    0.75428  

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
              "_ospanscore_int_0", "_ospanscore_int_1", "_ospanscore_int_2" 
              )

# Store ind models in list
mon_mods_allset <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods_allset,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "mon_mods_allset.Rdata"))

  


  
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

