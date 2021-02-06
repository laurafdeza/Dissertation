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
gca_mods_path  <- here("mods", "vision", "gca")

# Load models as lists
#load(paste0(gca_mods_path, "/ind_mods.Rdata"))
load(paste0(gca_mods_path, "/group_mods.Rdata"))
#load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))

# Store objects in global env
# list2env(ind_mods, globalenv())
list2env(group_mods, globalenv())
# list2env(nested_model_comparisons, globalenv())
list2env(model_preds, globalenv())

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




stress_gc_subset <- vision50 %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "mon", "aes", "ams", "ies", "ims"),
            stress_sum = if_else(cond == "1", 1, -1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------







# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + stress_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = stress_gc_subset, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + stress_sum + ot1 | participant) +
             ot2 + (1 + stress_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + stress_sum + ot1 + ot2 | participant) +
             ot3 + (1 + stress_sum + ot1 + ot2 + ot3 | participant))
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3, mod_ot4)
  
## All participants
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
gca_full_mod_base <-
  lmer(eLog ~ 1 + stress_sum * car_dev * corsi * (ot1 + ot2 + ot3) +          
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 3e5)),
       data = stress_gc_subset, REML = F)    # , na.action = na.exclude 

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_group_0 <- update(gca_full_mod_base,    . ~ . + group)
gca_full_mod_group_1 <- update(gca_full_mod_group_0, . ~ . + ot1:group)
gca_full_mod_group_2 <- update(gca_full_mod_group_1, . ~ . + ot2:group)
gca_full_mod_group_3 <- update(gca_full_mod_group_2, . ~ . + ot3:group)

full_group_anova <-
  anova(gca_full_mod_base, gca_full_mod_group_0, gca_full_mod_group_1,
        gca_full_mod_group_2, gca_full_mod_group_3)
#                        Df    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
# gca_full_mod_base      49 219731 220152 -109817   219633                          
# gca_full_mod_group_0   53 219713 220168 -109804   219607 26.0053  4  3.157e-05 ***
# gca_full_mod_group_1   57 219712 220202 -109799   219598  8.9010  4    0.06362 .  
# gca_full_mod_group_2   61 219692 220215 -109785   219570 28.5872  4  9.481e-06 ***
# gca_full_mod_group_3   65 219694 220252 -109782   219564  5.6758  4    0.22470    

# add 3-way int to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_int_0 <- update(gca_full_mod_group_2, . ~ . + stress_sum:car_dev:group:corsi)
gca_full_mod_int_1 <- update(gca_full_mod_int_0,   . ~ . + ot1:stress_sum:car_dev:group:corsi)
gca_full_mod_int_2 <- update(gca_full_mod_int_1,   . ~ . + ot2:stress_sum:car_dev:group:corsi)
gca_full_mod_int_3 <- update(gca_full_mod_int_2,   . ~ . + ot3:stress_sum:car_dev:group:corsi)

full_int_anova <-
  anova(gca_full_mod_group_2, gca_full_mod_int_0, gca_full_mod_int_1,
        gca_full_mod_int_2, gca_full_mod_int_3)
#                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_full_mod_group_2   61 219692 220215 -109785   219570                         
# gca_full_mod_int_0     65 219698 220256 -109784   219568  1.3627  4   0.850658   
# gca_full_mod_int_1     69 219691 220283 -109776   219553 15.4050  4   0.003931 **
# gca_full_mod_int_2     73 219694 220321 -109774   219548  4.4122  4   0.353085   
# gca_full_mod_int_3     77 219697 220358 -109771   219543  5.8860  4   0.207827     



# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "aes"))
gca_full_mod_int_1_aes <- update(gca_full_mod_int_1)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ams"))
gca_full_mod_int_1_ams <- update(gca_full_mod_int_1)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ies"))
gca_full_mod_int_1_ies <- update(gca_full_mod_int_1)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ims"))
gca_full_mod_int_1_ims <- update(gca_full_mod_int_1)

#stress_gc_subset %<>% mutate(., group = fct_relevel(group, "mon"))

}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, stress_sum,
                car_dev, corsi) %>%
  distinct

# Get model predictions and SE
fits_all <- predictSE(gca_full_mod, new_dat_all) %>%        
  as_tibble %>%
  bind_cols(new_dat_all) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, SS = "mon", AE = "aes", IE = "ies", AM = "ams", IM = "ims"))

# Filter preds at target offset
target_offset_preds <- filter(fits_all, time_zero == 4) %>%
  select(group, stress = stress_sum, corsi, motion = car_dev,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

# -----------------------------------------------------------------------------






# Save models -----------------------------------------------------------------

if(F) {

mod_type <- "gca_full_mod"
mod_spec <- c("_base", "_group_0", "_group_1", "_group_2", "_group_3",
              '_int_0', '_int_1', '_int_2', "_int_3",
              '_int_1_ams', '_int_1_ims', '_int_1_aes', '_int_1_ies')
  
full_mods <- mget(c(paste0(mod_type, mod_spec)))
  
save(full_mods,
       file = here("mods", "vision", "gca",
                   "full_mods.Rdata"))
  
  
  
  
# Save anova model comparisons
nested_model_comparisons <-
  mget(c(#"mon_cond_anova", "mon_car_anova", "mon_corsi_anova",
         
         # "aes_cond_anova", "aes_car_anova", "aes_corsi_anova",
         # 
         # "ies_cond_anova", "ies_car_anova", "ies_corsi_anova",
         # 
         # "ams_cond_anova", "ams_car_anova", "ams_corsi_anova",
         # 
         # "ims_cond_anova", "ims_car_anova", "ims_corsi_anova",
         
         "full_group_anova", "full_int_anova"))

save(nested_model_comparisons,
     file = here("mods", "vision", "gca",
                 "nested_model_comparisons.Rdata"))









# Save models predictions
model_preds <- mget(c("fits_all", "target_offset_preds"))

save(model_preds,
     file = here("mods", "vision", "gca",
                 "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------




# -----------------------------------------------------------------------------

# Full model adding 1 by 1
# Model fit not continuous either

gca_full_mod_base <-
  lmer(eLog ~ 1 + group * (ot1 + ot2 + ot3) +          
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 3e5)),
       data = stress_gc_subset, REML = F) 

gca_full_mod_stress_0 <- update(gca_full_mod_base,    . ~ . + stress_sum)
gca_full_mod_stress_1 <- update(gca_full_mod_stress_0, . ~ . + ot1:stress_sum)
gca_full_mod_stress_2 <- update(gca_full_mod_stress_1, . ~ . + ot2:stress_sum)
gca_full_mod_stress_3 <- update(gca_full_mod_stress_2, . ~ . + ot3:stress_sum)

full_stress_anova <-
  anova(gca_full_mod_base, gca_full_mod_stress_0, gca_full_mod_stress_1,
        gca_full_mod_stress_2, gca_full_mod_stress_3)
#                       npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
# gca_full_mod_base       37 219713 220030 -109819   219639                          
# gca_full_mod_stress_0   38 219715 220041 -109819   219639  0.0179  1    0.89361    
# gca_full_mod_stress_1   39 219717 220052 -109819   219639  0.1171  1    0.73218    
# gca_full_mod_stress_2   40 219698 220042 -109809   219618 20.5763  1   5.73e-06 ***
# gca_full_mod_stress_3   41 219697 220049 -109807   219615  3.2966  1    0.06942 .  

gca_full_mod_corsi_0 <- update(gca_full_mod_stress_2,    . ~ . + corsi)
gca_full_mod_corsi_1 <- update(gca_full_mod_corsi_0, . ~ . + ot1:corsi)
gca_full_mod_corsi_2 <- update(gca_full_mod_corsi_1, . ~ . + ot2:corsi)
gca_full_mod_corsi_3 <- update(gca_full_mod_corsi_2, . ~ . + ot3:corsi)

full_corsi_anova <-
  anova(gca_full_mod_stress_2, gca_full_mod_corsi_0, gca_full_mod_corsi_1,
        gca_full_mod_corsi_2, gca_full_mod_corsi_3)
#                       npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_full_mod_stress_2   40 219698 220042 -109809   219618                       
# gca_full_mod_corsi_0    41 219699 220051 -109809   219617 0.7701  1    0.38020  
# gca_full_mod_corsi_1    42 219701 220062 -109809   219617 0.0321  1    0.85772  
# gca_full_mod_corsi_2    43 219700 220069 -109807   219614 3.0281  1    0.08184 .
# gca_full_mod_corsi_3    44 219697 220074 -109804   219609 5.5764  1    0.01820 *

gca_full_mod_car_0 <- update(gca_full_mod_corsi_3,  . ~ . + car_dev)
gca_full_mod_car_1 <- update(gca_full_mod_car_0, . ~ . + ot1:car_dev)
gca_full_mod_car_2 <- update(gca_full_mod_car_1, . ~ . + ot2:car_dev)
gca_full_mod_car_3 <- update(gca_full_mod_car_2, . ~ . + ot3:car_dev)

full_car_anova <-
  anova(gca_full_mod_corsi_3, gca_full_mod_car_0, gca_full_mod_car_1,
        gca_full_mod_car_2, gca_full_mod_car_3)
#                      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_full_mod_corsi_3   44 219697 220074 -109804   219609                     
# gca_full_mod_car_0     45 219699 220085 -109804   219609 0.0670  1     0.7957
# gca_full_mod_car_1     46 219700 220095 -109804   219608 0.8865  1     0.3464
# gca_full_mod_car_2     47 219701 220105 -109804   219607 0.7503  1     0.3864
# gca_full_mod_car_3     48 219703 220115 -109804   219607 0.0008  1     0.9779

new_dat_all1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, stress_sum, corsi) %>%
  distinct

# Get model predictions and SE
fits_all1 <- predictSE(gca_full_mod_corsi_3, new_dat_all1) %>%        
  as_tibble %>%
  bind_cols(new_dat_all1) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, SS = "mon", AE = "aes", IE = "ies", AM = "ams", IM = "ims"))







# Individual models -----------------------------------------------------------



############# Individual models not run


#
# only mon
#

gca_mod_mon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +          
         (1 + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 3e5)),    # 2e4
       REML = F, na.action = na.exclude,
       data = filter(stress_gc_subset, group == 'mon')) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + stress_sum) # singular
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0, . ~ . + ot1:stress_sum) # singular
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1, . ~ . + ot2:stress_sum) # singular
gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2, . ~ . + ot3:stress_sum) # singular

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2, gca_mod_mon_cond_3) # 
#                      Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_base     25 29464 29629 -14707    29414                     
# gca_mod_mon_cond_0   26 29465 29636 -14706    29413 1.0536  1     0.3047
# gca_mod_mon_cond_1   27 29466 29644 -14706    29412 0.5960  1     0.4401
# gca_mod_mon_cond_2   28 29468 29652 -14706    29412 0.5675  1     0.4513
# gca_mod_mon_cond_3   29 29470 29661 -14706    29412 0.1315  1     0.7169


# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_car_0 <- update(gca_mod_mon_cond_3,. ~ . + car_dev) # singular
gca_mod_mon_car_1 <- update(gca_mod_mon_car_0, . ~ . + ot1:car_dev) # singular
gca_mod_mon_car_2 <- update(gca_mod_mon_car_1, . ~ . + ot2:car_dev) # singular
gca_mod_mon_car_3 <- update(gca_mod_mon_car_2, . ~ . + ot3:car_dev) # singular

mon_car_anova <-
  anova(gca_mod_mon_cond_3, gca_mod_mon_car_0, gca_mod_mon_car_1,
        gca_mod_mon_car_2, gca_mod_mon_car_3)
#                      Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_cond_3   29 29470 29661 -14706    29412                     
# gca_mod_mon_car_0    30 29471 29669 -14706    29411 0.1713  1     0.6790
# gca_mod_mon_car_1    31 29473 29678 -14705    29411 0.5067  1     0.4766
# gca_mod_mon_car_2    32 29475 29686 -14705    29411 0.0042  1     0.9485
# gca_mod_mon_car_3    33 29477 29694 -14705    29411 0.2041  1     0.6514


# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_wm_0 <- update(gca_mod_mon_car_3,   . ~ . + WM_set) # singular
gca_mod_mon_wm_1 <- update(gca_mod_mon_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_mon_wm_2 <- update(gca_mod_mon_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_mon_wm_3 <- update(gca_mod_mon_wm_2, . ~ . + ot3:WM_set) # singular

mon_wm_anova <-
  anova(gca_mod_mon_car_3, gca_mod_mon_wm_0, gca_mod_mon_wm_1,
        gca_mod_mon_wm_2, gca_mod_mon_wm_3)
#                     Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_car_3   33 29477 29694 -14705    29411                     
# gca_mod_mon_wm_0    34 29479 29703 -14705    29411 0.1063  1     0.7444
# gca_mod_mon_wm_1    35 29481 29712 -14705    29411 0.0101  1     0.9201
# gca_mod_mon_wm_2    36 29482 29719 -14705    29410 0.8890  1     0.3457
# gca_mod_mon_wm_3    37 29484 29728 -14705    29410 0.0069  1     0.9338


# add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_corsi_0 <- update(gca_mod_mon_wm_3,   . ~ . + corsi_pt) # singular
gca_mod_mon_corsi_1 <- update(gca_mod_mon_corsi_0, . ~ . + ot1:corsi_pt) # singular
gca_mod_mon_corsi_2 <- update(gca_mod_mon_corsi_1, . ~ . + ot2:corsi_pt) # singular
gca_mod_mon_corsi_3 <- update(gca_mod_mon_corsi_2, . ~ . + ot3:corsi_pt) # singular

mon_corsi_anova <-
  anova(gca_mod_mon_wm_3, gca_mod_mon_corsi_0, gca_mod_mon_corsi_1,
        gca_mod_mon_corsi_2, gca_mod_mon_corsi_3)
#                       Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_wm_3      37 29484 29728 -14705    29410                       
# gca_mod_mon_corsi_0   38 29480 29730 -14702    29404 5.9864  1    0.01442 *
# gca_mod_mon_corsi_1   39 29482 29739 -14702    29404 0.1318  1    0.71661  
# gca_mod_mon_corsi_2   40 29483 29747 -14701    29403 0.7601  1    0.38329  
# gca_mod_mon_corsi_3   41 29484 29755 -14701    29402 0.5118  1    0.47437  




#######

# only aes

#######

gca_mod_aes_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
       REML = F,
       data = filter(stress_gc_subset, group == "aes")) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_cond_0 <- update(gca_mod_aes_base,   . ~ . + stress_sum) # singular
gca_mod_aes_cond_1 <- update(gca_mod_aes_cond_0, . ~ . + ot1:stress_sum) # singular
gca_mod_aes_cond_2 <- update(gca_mod_aes_cond_1, . ~ . + ot2:stress_sum) # singular
gca_mod_aes_cond_3 <- update(gca_mod_aes_cond_2, . ~ . + ot3:stress_sum) # singular

aes_cond_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_cond_0, gca_mod_aes_cond_1,
        gca_mod_aes_cond_2, gca_mod_aes_cond_3) # 
#                      Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_aes_base     25 47600 47777 -23775    47550                       
# gca_mod_aes_cond_0   26 47601 47785 -23775    47549 1.1557  1    0.28235  
# gca_mod_aes_cond_1   27 47603 47794 -23775    47549 0.0476  1    0.82734  
# gca_mod_aes_cond_2   28 47599 47797 -23772    47543 6.0072  1    0.01425 *
# gca_mod_aes_cond_3   29 47601 47806 -23772    47543 0.1151  1    0.73442  

# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_car_0 <- update(gca_mod_aes_cond_3,   . ~ . + car_dev) # singular
gca_mod_aes_car_1 <- update(gca_mod_aes_car_0, . ~ . + ot1:car_dev) # singular
gca_mod_aes_car_2 <- update(gca_mod_aes_car_1, . ~ . + ot2:car_dev) # singular
gca_mod_aes_car_3 <- update(gca_mod_aes_car_2, . ~ . + ot3:car_dev) # singular

aes_car_anova <-
  anova(gca_mod_aes_cond_3, gca_mod_aes_car_0, gca_mod_aes_car_1,
        gca_mod_aes_car_2, gca_mod_aes_car_3)
#                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_aes_cond_3   29 47601 47806 -23772    47543                     
# gca_mod_aes_car_0    30 47603 47815 -23772    47543 0.0027  1     0.9584
# gca_mod_aes_car_1    31 47605 47824 -23771    47543 0.3124  1     0.5762
# gca_mod_aes_car_2    32 47606 47832 -23771    47542 0.7113  1     0.3990
# gca_mod_aes_car_3    33 47607 47840 -23770    47541 1.1637  1     0.2807


# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_wm_0 <- update(gca_mod_aes_car_3,   . ~ . + WM_set) # singular
gca_mod_aes_wm_1 <- update(gca_mod_aes_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_aes_wm_2 <- update(gca_mod_aes_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_aes_wm_3 <- update(gca_mod_aes_wm_2, . ~ . + ot3:WM_set) # singular

aes_wm_anova <-
  anova(gca_mod_aes_car_3, gca_mod_aes_wm_0, gca_mod_aes_wm_1,
        gca_mod_aes_wm_2, gca_mod_aes_wm_3)
#                     Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_aes_car_3   33 47607 47840 -23770    47541                     
# gca_mod_aes_wm_0    34 47609 47849 -23770    47541 0.0512  1     0.8209
# gca_mod_aes_wm_1    35 47611 47858 -23770    47541 0.0419  1     0.8378
# gca_mod_aes_wm_2    36 47611 47866 -23770    47539 1.3110  1     0.2522
# gca_mod_aes_wm_3    37 47613 47875 -23770    47539 0.0007  1     0.9785


# add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_corsi_0 <- update(gca_mod_aes_wm_3,   . ~ . + corsi_pt) # singular
gca_mod_aes_corsi_1 <- update(gca_mod_aes_corsi_0, . ~ . + ot1:corsi_pt) # singular
gca_mod_aes_corsi_2 <- update(gca_mod_aes_corsi_1, . ~ . + ot2:corsi_pt) # singular
gca_mod_aes_corsi_3 <- update(gca_mod_aes_corsi_2, . ~ . + ot3:corsi_pt) # singular

aes_corsi_anova <-
  anova(gca_mod_aes_wm_3, gca_mod_aes_corsi_0, gca_mod_aes_corsi_1,
        gca_mod_aes_corsi_2, gca_mod_aes_corsi_3)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_aes_wm_3      37 47613 47875 -23770    47539                     
# gca_mod_aes_corsi_0   38 47615 47884 -23770    47539 0.1972  1     0.6570
# gca_mod_aes_corsi_1   39 47617 47893 -23770    47539 0.2230  1     0.6368
# gca_mod_aes_corsi_2   40 47619 47902 -23770    47539 0.0865  1     0.7687
# gca_mod_aes_corsi_3   41 47621 47911 -23770    47539 0.0471  1     0.8282



#
# only ies
#

gca_mod_ies_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ies")) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_cond_0 <- update(gca_mod_ies_base,   . ~ . + stress_sum) # singular
gca_mod_ies_cond_1 <- update(gca_mod_ies_cond_0,   . ~ . + ot1:stress_sum) # singular
gca_mod_ies_cond_2 <- update(gca_mod_ies_cond_1,   . ~ . + ot2:stress_sum) # singular
gca_mod_ies_cond_3 <- update(gca_mod_ies_cond_2,   . ~ . + ot3:stress_sum) # singular

ies_cond_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_cond_0, gca_mod_ies_cond_1,    
        gca_mod_ies_cond_2, gca_mod_ies_cond_3)
#                    Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ies_base   25 49516 49693 -24733    49466                       
# gca_mod_ies_cond_0   26 49517 49702 -24733    49465 0.2707  1    0.60287  
# gca_mod_ies_cond_1   27 49519 49710 -24732    49465 0.3949  1    0.52976  
# gca_mod_ies_cond_2   28 49521 49719 -24732    49465 0.0717  1    0.78882  
# gca_mod_ies_cond_3   29 49520 49725 -24731    49462 3.3299  1    0.06803 .

# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_car_0 <- update(gca_mod_ies_cond_3,   . ~ . + car_dev) # singular
gca_mod_ies_car_1 <- update(gca_mod_ies_car_0, . ~ . + ot1:car_dev) # singular
gca_mod_ies_car_2 <- update(gca_mod_ies_car_1, . ~ . + ot2:car_dev) # singular
gca_mod_ies_car_3 <- update(gca_mod_ies_car_2, . ~ . + ot3:car_dev) # singular

ies_car_anova <-
  anova(gca_mod_ies_cond_3, gca_mod_ies_car_0, gca_mod_ies_car_1,
        gca_mod_ies_car_2, gca_mod_ies_car_3)
#                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ies_cond_3   29 49520 49725 -24731    49462                       
# gca_mod_ies_car_0    30 49521 49734 -24731    49461 0.2958  1     0.5865  
# gca_mod_ies_car_1    31 49523 49743 -24730    49461 0.2261  1     0.6344  
# gca_mod_ies_car_2    32 49521 49748 -24728    49457 4.4796  1     0.0343 *
# gca_mod_ies_car_3    33 49522 49756 -24728    49456 0.2254  1     0.6349  



# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_wm_0 <- update(gca_mod_ies_car_3,   . ~ . + WM_set) # singular
gca_mod_ies_wm_1 <- update(gca_mod_ies_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_ies_wm_2 <- update(gca_mod_ies_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_ies_wm_3 <- update(gca_mod_ies_wm_2, . ~ . + ot3:WM_set) # singular

ies_wm_anova <-
  anova(gca_mod_ies_car_3, gca_mod_ies_wm_0, gca_mod_ies_wm_1,
        gca_mod_ies_wm_2, gca_mod_ies_wm_3)
#                     Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_car_3   33 49522 49756 -24728    49456                     
# gca_mod_ies_wm_0    34 49523 49765 -24728    49455 0.9087  1     0.3405
# gca_mod_ies_wm_1    35 49525 49774 -24728    49455 0.0037  1     0.9514
# gca_mod_ies_wm_2    36 49526 49781 -24727    49454 1.9726  1     0.1602
# gca_mod_ies_wm_3    37 49527 49790 -24727    49453 0.1186  1     0.7306

# add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_corsi_0 <- update(gca_mod_ies_wm_3,   . ~ . + corsi_pt) # singular
gca_mod_ies_corsi_1 <- update(gca_mod_ies_corsi_0, . ~ . + ot1:corsi_pt) # singular
gca_mod_ies_corsi_2 <- update(gca_mod_ies_corsi_1, . ~ . + ot2:corsi_pt) # singular
gca_mod_ies_corsi_3 <- update(gca_mod_ies_corsi_2, . ~ . + ot3:corsi_pt) # singular

ies_corsi_anova <-
  anova(gca_mod_ies_wm_3, gca_mod_ies_corsi_0, gca_mod_ies_corsi_1,
        gca_mod_ies_corsi_2, gca_mod_ies_corsi_3)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ies_wm_3      37 49527 49790 -24727    49453                       
# gca_mod_ies_corsi_0   38 49529 49798 -24726    49453 0.3799  1    0.53766  
# gca_mod_ies_corsi_1   39 49529 49806 -24726    49451 1.7607  1    0.18454  
# gca_mod_ies_corsi_2   40 49531 49814 -24725    49451 0.7055  1    0.40095  
# gca_mod_ies_corsi_3   41 49527 49818 -24722    49445 5.6728  1    0.01723 *



#
# only ams
#

gca_mod_ams_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ams")) # singular

# add cond effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_cond_0 <- update(gca_mod_ams_base,   . ~ . + stress_sum) # singular
gca_mod_ams_cond_1 <- update(gca_mod_ams_cond_0,   . ~ . + ot1:stress_sum) # singular
gca_mod_ams_cond_2 <- update(gca_mod_ams_cond_1,   . ~ . + ot2:stress_sum) # singular
gca_mod_ams_cond_3 <- update(gca_mod_ams_cond_2,   . ~ . + ot3:stress_sum) # singular

ams_cond_anova <-
  anova(gca_mod_ams_base, gca_mod_ams_cond_0, gca_mod_ams_cond_1,
        gca_mod_ams_cond_2, gca_mod_ams_cond_3)
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_base     25 45266 45441 -22608    45216                       
# gca_mod_ams_cond_0   26 45268 45450 -22608    45216 0.0006  1    0.98044  
# gca_mod_ams_cond_1   27 45270 45459 -22608    45216 0.0073  1    0.93197  
# gca_mod_ams_cond_2   28 45268 45464 -22606    45212 4.2239  1    0.03986 *
# gca_mod_ams_cond_3   29 45270 45473 -22606    45212 0.0003  1    0.98634


# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_car_0 <- update(gca_mod_ams_cond_3,   . ~ . + car_dev) # singular
gca_mod_ams_car_1 <- update(gca_mod_ams_car_0, . ~ . + ot1:car_dev) # singular
gca_mod_ams_car_2 <- update(gca_mod_ams_car_1, . ~ . + ot2:car_dev) # singular
gca_mod_ams_car_3 <- update(gca_mod_ams_car_2, . ~ . + ot3:car_dev) # singular

ams_car_anova <-
  anova(gca_mod_ams_cond_3, gca_mod_ams_car_0, gca_mod_ams_car_1,
        gca_mod_ams_car_2, gca_mod_ams_car_3)
#                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ams_cond_3   29 45270 45473 -22606    45212                     
# gca_mod_ams_car_0    30 45271 45481 -22606    45211 0.4839  1     0.4866
# gca_mod_ams_car_1    31 45272 45489 -22605    45210 1.1675  1     0.2799
# gca_mod_ams_car_2    32 45273 45497 -22604    45209 1.5987  1     0.2061
# gca_mod_ams_car_3    33 45274 45506 -22604    45208 0.1201  1     0.7289


# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_wm_0 <- update(gca_mod_ams_car_3,   . ~ . + WM_set) # singular
gca_mod_ams_wm_1 <- update(gca_mod_ams_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_ams_wm_2 <- update(gca_mod_ams_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_ams_wm_3 <- update(gca_mod_ams_wm_2, . ~ . + ot3:WM_set) # singular

ams_wm_anova <-
  anova(gca_mod_ams_car_3, gca_mod_ams_wm_0, gca_mod_ams_wm_1,
        gca_mod_ams_wm_2, gca_mod_ams_wm_3)
#                     Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_car_3   33 45274 45506 -22604    45208                       
# gca_mod_ams_wm_0    34 45276 45514 -22604    45208 0.1149  1    0.73460  
# gca_mod_ams_wm_1    35 45274 45519 -22602    45204 4.5085  1    0.03373 *
# gca_mod_ams_wm_2    36 45274 45526 -22601    45202 1.4752  1    0.22453  
# gca_mod_ams_wm_3    37 45276 45535 -22601    45202 0.1148  1    0.73471  


# add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_corsi_0 <- update(gca_mod_ams_wm_3,   . ~ . + corsi_pt) # singular
gca_mod_ams_corsi_1 <- update(gca_mod_ams_corsi_0, . ~ . + ot1:corsi_pt) # singular
gca_mod_ams_corsi_2 <- update(gca_mod_ams_corsi_1, . ~ . + ot2:corsi_pt) # singular
gca_mod_ams_corsi_3 <- update(gca_mod_ams_corsi_2, . ~ . + ot3:corsi_pt) # singular

ams_corsi_anova <-
  anova(gca_mod_ams_wm_3, gca_mod_ams_corsi_0, gca_mod_ams_corsi_1,
        gca_mod_ams_corsi_2, gca_mod_ams_corsi_3)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ams_wm_3      37 45276 45535 -22601    45202                     
# gca_mod_ams_corsi_0   38 45277 45543 -22600    45201 1.2731  1     0.2592
# gca_mod_ams_corsi_1   39 45279 45552 -22600    45201 0.0417  1     0.8382
# gca_mod_ams_corsi_2   40 45280 45560 -22600    45200 0.6377  1     0.4246
# gca_mod_ams_corsi_3   41 45281 45568 -22600    45199 0.9876  1     0.3203



#
# only ims
#

gca_mod_ims_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ims")) # singular

# add cond effect to imsercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_cond_0 <- update(gca_mod_ims_base,   . ~ . + stress_sum) # singular
gca_mod_ims_cond_1 <- update(gca_mod_ims_cond_0,   . ~ . + ot1:stress_sum) # singular
gca_mod_ims_cond_2 <- update(gca_mod_ims_cond_1,   . ~ . + ot2:stress_sum) # singular
gca_mod_ims_cond_3 <- update(gca_mod_ims_cond_2,   . ~ . + ot3:stress_sum) # singular

ims_cond_anova <-
  anova(gca_mod_ims_base, gca_mod_ims_cond_0, gca_mod_ims_cond_1, # none singular, and none significant
        gca_mod_ims_cond_2, gca_mod_ims_cond_3)
#                      Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ims_base     25 47176 47352 -23563    47126                     
# gca_mod_ims_cond_0   26 47178 47361 -23563    47126 0.3962  1     0.5291
# gca_mod_ims_cond_1   27 47179 47369 -23562    47125 1.1522  1     0.2831
# gca_mod_ims_cond_2   28 47179 47376 -23562    47123 1.8059  1     0.1790
# gca_mod_ims_cond_3   29 47179 47383 -23561    47121 1.6463  1     0.1995


# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_car_0 <- update(gca_mod_ims_cond_3,   . ~ . + car_dev) # singular
gca_mod_ims_car_1 <- update(gca_mod_ims_car_0, . ~ . + ot1:car_dev) # singular
gca_mod_ims_car_2 <- update(gca_mod_ims_car_1, . ~ . + ot2:car_dev) # singular
gca_mod_ims_car_3 <- update(gca_mod_ims_car_2, . ~ . + ot3:car_dev) # singular

ims_car_anova <-
  anova(gca_mod_ims_cond_3, gca_mod_ims_car_0, gca_mod_ims_car_1,
        gca_mod_ims_car_2, gca_mod_ims_car_3)
#                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ims_cond_3   29 47179 47383 -23561    47121                       
# gca_mod_ims_car_0    30 47177 47389 -23559    47117 3.9619  1    0.04654 *
# gca_mod_ims_car_1    31 47173 47391 -23556    47111 6.3780  1    0.01155 *
# gca_mod_ims_car_2    32 47175 47400 -23556    47111 0.0044  1    0.94728  
# gca_mod_ims_car_3    33 47176 47408 -23555    47110 1.1026  1    0.29369  


# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_wm_0 <- update(gca_mod_ims_car_3,   . ~ . + WM_set) # singular
gca_mod_ims_wm_1 <- update(gca_mod_ims_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_ims_wm_2 <- update(gca_mod_ims_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_ims_wm_3 <- update(gca_mod_ims_wm_2, . ~ . + ot3:WM_set) # singular

ims_wm_anova <-
  anova(gca_mod_ims_car_3, gca_mod_ims_wm_0, gca_mod_ims_wm_1,
        gca_mod_ims_wm_2, gca_mod_ims_wm_3)
#                     Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ims_car_3   33 47176 47408 -23555    47110                        
# gca_mod_ims_wm_0    34 47178 47417 -23555    47110 0.4237  1   0.515094   
# gca_mod_ims_wm_1    35 47179 47426 -23555    47109 0.0888  1   0.765761   
# gca_mod_ims_wm_2    36 47172 47425 -23550    47100 9.6190  1   0.001926 **
# gca_mod_ims_wm_3    37 47174 47434 -23550    47100 0.2083  1   0.648107 

# add visuospatial WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_corsi_0 <- update(gca_mod_ims_wm_3,   . ~ . + corsi_pt) # singular
gca_mod_ims_corsi_1 <- update(gca_mod_ims_corsi_0, . ~ . + ot1:corsi_pt) # singular
gca_mod_ims_corsi_2 <- update(gca_mod_ims_corsi_1, . ~ . + ot2:corsi_pt) # singular
gca_mod_ims_corsi_3 <- update(gca_mod_ims_corsi_2, . ~ . + ot3:corsi_pt) # singular

ims_corsi_anova <-
  anova(gca_mod_ims_wm_3, gca_mod_ims_corsi_0, gca_mod_ims_corsi_1,
        gca_mod_ims_corsi_2, gca_mod_ims_corsi_3)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ims_wm_3      37 47174 47434 -23550    47100                     
# gca_mod_ims_corsi_0   38 47175 47442 -23550    47099 0.6494  1     0.4203
# gca_mod_ims_corsi_1   39 47175 47449 -23548    47097 1.9123  1     0.1667
# gca_mod_ims_corsi_2   40 47177 47458 -23548    47097 0.2083  1     0.6481
# gca_mod_ims_corsi_3   41 47178 47466 -23548    47096 0.9455  1     0.3309


mod_start <- "gca_mod_"
mod_spec <- c("_base", "_cond_0", "_cond_1", "_cond_2", "_cond_3", 
              "_wm_0", "_wm_1", "_wm_2", "_wm_3",
              "_car_0", "_car_1", "_car_2", "_car_3",
              "_corsi_0", "_corsi_1", "_corsi_2", "_corsi_3")

# Store ind models in list
ind_mods <- mget(c(paste0(mod_start, "mon", mod_spec),
                   paste0(mod_start, "aes", mod_spec),
                   paste0(mod_start, "ies", mod_spec),
                   paste0(mod_start, "ams", mod_spec),
                   paste0(mod_start, "ims", mod_spec)
))

save(ind_mods,
     file = here("mods", "vision", "gca",
                 "individual_mods.Rdata"))


# -----------------------------------------------------------------------------
