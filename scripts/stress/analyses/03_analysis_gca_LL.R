#
# Original script by Joseph
# Modified by Cristina to adapt to CrisLau Project
# Then by Laura 
# Last update: 09/02/2020
#
# Growth curve analyisis ------------------------------------------------------
#
# - Question 1: Are the groups different from each other in when they begin
#   to fixate on the target?
#     - test 5 groups at each level of 'condition'
#     - hypothesis: SS has steeper slope for both conditions
# - Question 2: W/in groups, is the there a difference between
#   oxytone/paroxytone items?
#     - test oxytone vs. paroxytone for each group
#     - hypothesis: steeper slope/earlier break in oxytone condition
# - Question 3: Does verbal WM mediate fixation on the target?
#     - compare WM and fixations across times across groups
#     - hypothesis: higher WM helps in fixating on target earlier
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# Load data
source(here::here("scripts", "02_load_data.R"))
stress50 <- read_csv(here("data", "clean", "stress_50ms_final_onsetc3updated.csv"))

# Get path to saved models
gca_mods_path  <- here("mods", "stress", "gca", "LL_changes")

# Load models as lists
load(paste0(gca_mods_path, "/gca_mon_mods.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_mods.Rdata"))
load(paste0(gca_mods_path, "/nested_model_l2.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))

# # Store objects in global env
# list2env(gca_mon_mods, globalenv())
# list2env(gca_l2_mods, globalenv())
# list2env(gca_l2_mods_sum, globalenv())
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



# stress_50 <- na.omit(stress50)


stress_gc_subset <- stress50 %>%
  # select(., -WM_set) %>%
  filter(., time_zero >= -4 & time_zero <= 4) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma")#,
         #   condition_sum = if_else(cond == "1", -1, 1)
         ) %>%       # 1 = present (now -1), 2 = past (now 1)
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")



# -----------------------------------------------------------------------------







#################### MONOLINGUAL SPEAKERS ########################################

# Build up random effects to test time terms
if(F){
  
  mon_data <- filter(stress_gc_subset, l1 == 'es') %>% select(-DELE, -percent_l2_week)
  
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),   # , optCtrl=list(maxfun=2e5)
         data = mon_data, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + ot1 | participant) +
             ot2 + (1 + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + ot1 + ot2 | participant) +
             ot3 + (1 + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    6 24003 24041 -11995    23991                          
  # mod_ot2   10 23981 24045 -11980    23961 29.7212  4  5.578e-06 ***
  # mod_ot3   15 23988 24084 -11979    23958  2.6159  5      0.759

  
  mod_ot0 <- update(mod_ot2, . ~ . + (1 | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + ot1 | target) +
                       + (1 + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + ot1 + ot2 | target) +
                       + ot3 + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot2, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot2    10 23981 24045 -11980    23961                          
  # mod_ot0    11 23889 23959 -11934    23867 93.7183  1  < 2.2e-16 ***
  # mod_ot1a   13 23872 23955 -11923    23846 20.7574  2  3.109e-05 ***
  # mod_ot2a   16 23874 23976 -11921    23842  4.4507  3     0.2167    
  # mod_ot3a   21 23882 24015 -11920    23840  2.4904  5     0.7779   
  
}



# Individual model MON -----------------------------------------------------------

gca_mod_mon_base <- mod_ot1a
  # lmer(eLog ~ 1 + (ot1 + ot2) +
  #        (1 + (ot1 + ot2) | participant) +
  #        (1 + ot1 | target),
  #      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
  #      REML = F,
  #      data = filter(mon_data)) 

## add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
#gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + condition_sum) 
#gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0,   . ~ . + ot1:condition_sum) 
#gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1,   . ~ . + ot2:condition_sum) 
#gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2,   . ~ . + ot3:condition_sum) 

#mon_cond_anova <-
#  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
#        gca_mod_mon_cond_2, gca_mod_mon_cond_3)
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base     


gca_mod_mon_final <- gca_mod_mon_base


mod_type <- "gca_mod_mon"
mod_spec <- c('_base')
              

# Store ind models in list
gca_mon_mods <- mget(c(paste0(mod_type, mod_spec)))

save(gca_mon_mods,
     file = here("mods", "stress", "gca", "LL_changes",
                 "gca_mon_mods.Rdata"))

 

#################### L2 SPEAKERS ########################################

l2_data <- stress_gc_subset%>%
  filter(., l1 != 'es') %>% 
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1))

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){

  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),   # , optCtrl=list(maxfun=2e5)
         data = l2_data, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + ot1 | participant) +
             ot2 + (1 + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + ot1 + ot2 | participant) +
             ot3 + (1 + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  
  #         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    6 100945 100991 -50466   100933                         
  # mod_ot2   10 100873 100951 -50426   100853 79.681  4     <2e-16 ***
  # mod_ot3   15 100878 100995 -50424   100848  4.383  5     0.4957    
  
  
  
  mod_ot0 <- update(mod_ot2, . ~ . + (1 + DELE_z + use_z | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 + DELE_z + use_z | target) + 
                       (1 + DELE_z + use_z + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + DELE_z + use_z + ot1 | target) +
                       + (1 + DELE_z + use_z + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + DELE_z + use_z + ot1 + ot2 | target) +
                       + ot3 + (1 + DELE_z + use_z + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot2, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot2    10 100873 100951 -50426   100853                           
  # mod_ot0    16 100732 100857 -50350   100700 152.5634  6  < 2.2e-16 ***
  # mod_ot1a   20 100706 100862 -50333   100666  34.2228  4  6.708e-07 ***
  # mod_ot2a   25 100704 100899 -50327   100654  11.7772  5    0.03797 *  
  # mod_ot3a   32 100714 100963 -50325   100650   4.6444  7    0.70327    
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_l2_mod_base <- mod_ot2a
   #lmer(eLog ~ 1 + (ot1 + ot2) +         
   #       (1 + DELE_z + use_z + (ot1 + ot2) | participant) +
   #       (1 + ot1 + ot2 | target), 
   #     control = lmerControl(optimizer = 'bobyqa',
   #                           optCtrl = list(maxfun = 2e4)),
   #     data = l2_data, REML = F)

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_l2_mod_l1_0 <- update(gca_l2_mod_base, . ~ . + l1_sum) 
gca_l2_mod_l1_1 <- update(gca_l2_mod_l1_0, . ~ . + ot1:l1_sum) 
gca_l2_mod_l1_2 <- update(gca_l2_mod_l1_1, . ~ . + ot2:l1_sum) 

l2_l1_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_l1_0, gca_l2_mod_l1_1,
        gca_l2_mod_l1_2)
#                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_l2_mod_base   25 100704 100899 -50327   100654                     
# gca_l2_mod_l1_0   26 100707 100910 -50328   100655 0.0000  1     1.0000
# gca_l2_mod_l1_1   27 100716 100927 -50331   100662 0.0000  1     1.0000
# gca_l2_mod_l1_2   28 100717 100936 -50331   100661 0.5951  1     0.4404


# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_dele_0 <- update(gca_l2_mod_base,   . ~ . + DELE_z) 
gca_l2_mod_dele_1 <- update(gca_l2_mod_dele_0, . ~ . + ot1:DELE_z) 
gca_l2_mod_dele_2 <- update(gca_l2_mod_dele_1, . ~ . + ot2:DELE_z)

l2_dele_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_dele_0, gca_l2_mod_dele_1,
        gca_l2_mod_dele_2)
#                     npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_l2_mod_base     25 100704 100899 -50327   100654                          
# gca_l2_mod_dele_0   26 100712 100915 -50330   100660  0.0000  1  1.0000000    
# gca_l2_mod_dele_1   27 100699 100910 -50323   100645 15.0735  1  0.0001034 ***
# gca_l2_mod_dele_2   28 100701 100919 -50322   100645  0.7087  1  0.3998799 



# add L2 use effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_use_0 <- update(gca_l2_mod_dele_1,  . ~ . + use_z) 
gca_l2_mod_use_1 <- update(gca_l2_mod_use_0, . ~ . + ot1:use_z) 
gca_l2_mod_use_2 <- update(gca_l2_mod_use_1, . ~ . + ot2:use_z)

l2_use_anova <-
  anova(gca_l2_mod_dele_1, gca_l2_mod_use_0, gca_l2_mod_use_1,
        gca_l2_mod_use_2)
#                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_l2_mod_dele_1   27 100699 100910 -50323   100645                       
# gca_l2_mod_use_0    28 100700 100918 -50322   100644 1.3387  1    0.24726  
# gca_l2_mod_use_1    29 100700 100926 -50321   100642 2.2794  1    0.13111  
# gca_l2_mod_use_2    30 100698 100932 -50319   100638 3.9104  1    0.04799 *



# add interaction

gca_l2_mod_int_0 <- update(gca_l2_mod_use_2,     . ~ . + l1_sum:DELE_z:use_z) 
gca_l2_mod_int_1 <- update(gca_l2_mod_int_0, . ~ . + ot1:l1_sum:DELE_z:use_z) 
gca_l2_mod_int_2 <- update(gca_l2_mod_int_1, . ~ . + ot2:l1_sum:DELE_z:use_z)

l2_int_anova <-
  anova(gca_l2_mod_use_2, gca_l2_mod_int_0, gca_l2_mod_int_1,
        gca_l2_mod_int_2)
#                      npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)   
# gca_l2_mod_use_2   30 100698 100932 -50319   100638                     
# gca_l2_mod_int_0   31 100699 100941 -50319   100637 0.2564  1     0.6126
# gca_l2_mod_int_1   32 100700 100949 -50318   100636 1.5717  1     0.2100
# gca_l2_mod_int_2   33 100711 100968 -50322   100645 0.0000  1     1.0000

gca_l2_mod_final <- gca_l2_mod_use_2

summary(gca_l2_mod_final)
#             Estimate Std. Error t value
# (Intercept)  0.06372    0.06967   0.915
# ot1          1.35457    0.14908   9.086
# ot2          0.67896    0.07619   8.912
# DELE_z       0.10416    0.05693   1.830
# use_z       -0.06765    0.06533  -1.036
# ot1:DELE_z   0.23245    0.11182   2.079
# ot1:use_z    0.21189    0.12212   1.735
# ot2:use_z    0.16292    0.09139   1.783


mod_type <- "gca_l2_mod"
mod_spec <- c('_base', 
              "_l1_0", "_l1_1", "_l1_2", 
              "_dele_0", "_dele_1", "_dele_2", 
              "_use_0", "_use_1", "_use_2", 
              "_int_0", "_int_1", "_int_2",
              '_final') 

# Store ind models in list
gca_l2_mods <- mget(c(paste0(mod_type, mod_spec)))

save(gca_l2_mods,
     file = here("mods", "stress", "gca", "LL_changes",
                 "gca_l2_mods.Rdata"))


# Save anova model comparisons
nested_model_comparisons <-
  mget(c('l2_l1_anova', 'l2_dele_anova', 
         'l2_use_anova', 'l2_int_anova'
  ))

save(nested_model_comparisons,
     file = here("mods", "stress", "gca", "LL_changes",
                 "nested_model_comparisons_l2.Rdata"))


}

# -----------------------------------------------------------------------------







# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_mon <- mon_data %>%
  dplyr::select(time_zero, ot1:ot2) %>% 
  distinct
  
# Get model predictions and SE
fits_all_mon <- predictSE(gca_mon_mods$gca_mod_mon_base, new_dat_mon) %>%  
  as_tibble %>%
  bind_cols(new_dat_mon) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at first syllable offset
target_offset_preds_mon <- filter(fits_all_mon, time_zero == 4) %>%
  select(elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))


# -----------------------------------------------------------------------------


new_l2_dat <- l2_data %>%
  dplyr::select(l1_sum, time_zero, ot1:ot2) %>%
  distinct %>%
  # mutate(l1_sum = as.character(l1_sum)) %>% 
  expand_grid(., tibble(DELE_z = c(-1, 0, 1))) %>%
  expand_grid(., tibble(use_z = c(-1, 0, 1)))

fits_all_l2 <- predictSE(gca_l2_mod_final, new_l2_dat) %>%
  as_tibble %>%
  bind_cols(new_l2_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at first syllable offset
target_offset_preds_l2 <- filter(fits_all_l2, time_zero == 4) %>% #
  select(l1 = l1_sum, DELE = DELE_z, `L2 use` = use_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))

 


model_preds <- mget(c("fits_all_mon", "fits_all_l2", 
  "target_offset_preds_mon", "target_offset_preds_l2"))

save(model_preds,
     file = here("mods", "stress", "gca", "LL_changes",
                 "model_preds.Rdata"))
