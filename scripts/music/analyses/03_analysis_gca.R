#
# Original script by Joseph
# Modified by Cristina to adapt to CrisLau Project
# Last updat: 06/12/2019
#
# Growth curve analyisis ------------------------------------------------------
#
# - Question 1: Are the groups different from each other in when they begin
#   to fixate on the target?
#     - test 3 groups at each level of 'condition'
#     - hypothesis: SS has steeper slope for both conditions
# - Question 2: W/in groups, is the there a difference between
#   oxytone/paroxytone items?
#     - test oxytone vs. paroxytone for each group
#     - hypothesis: steeper slope/earlier break in oxytone condition
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

# Load data
source(here::here("scripts", "02_load_data.R"))

# Get path to saved models
gca_mods_path  <- here("reports", "mods", "gca")

# Load models as lists
#load(paste0(gca_mods_path, "/ind_mods.Rdata"))
load(paste0(gca_mods_path, "/full_mods.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
#load(paste0(gca_mods_path, "/model_preds.Rdata"))

# Store objects in global env
list2env(ind_mods, globalenv())
list2env(full_mods, globalenv())
list2env(nested_model_comparisons, globalenv())
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

stress_gc_subset <- stress50 %>%
  filter(.,participant != "mon04" & participant != "mon05" &
           participant != "mon06" &
           participant != "mon11" & participant != "mon016" & 
           participant != "mon18" & participant != "mon20" & 
           participant != "mon21" & participant != "mon23" &
           participant != "mon24" & participant != "mon27" &
           participant != "ies05" & participant != "ies09" & 
           participant != 'ies16' & participant != 'ies21' &
           participant != "ies23" & participant != "ies30" &
           participant != "ies33" & participant != "aes01" &
           participant != 'aes03' & participant != 'aes04' &
           participant != "aes05" & participant != "aes07" &
           participant != "aes11" & participant != "aes14" &
           participant != 'aes16' & participant != 'aes23',
           time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "mon", "aes", "ams", "ies", "ims"),
            condition_sum = if_else(cond == "1", 1, -1)) %>%       # 1 = present, 2 = past
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------







# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){

  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + condition_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = stress_gc_subset, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + condition_sum + ot1 | participant) +
             ot2 + (1 + condition_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + condition_sum + ot1 + ot2 | participant) +
             ot3 + (1 + condition_sum + ot1 + ot2 + ot3 | participant))
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3, mod_ot4)
  
# We retain the most complex model: mod_ot4
  ## All participants
  #         Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
  
  
  ## WM homogeneous
  
  
  
  
  
}

# -----------------------------------------------------------------------------





# Individual models -----------------------------------------------------------

#
# only mon
#

gca_mod_mon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
       REML = F,
       data = filter(stress_gc_subset, group == "mon")) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + condition_sum) # singular
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2,   . ~ . + ot3:condition_sum) # singular

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2, gca_mod_mon_cond_3) # nada es sig

# WM homog 
# Error in anova.merMod(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,  : 
#                         models were not all fitted to the same size of dataset

#
# only aes
#

gca_mod_aes_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "aes")) # singular

# add coda effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_cond_0 <- update(gca_mod_aes_base,   . ~ . + condition_sum) # singular
gca_mod_aes_cond_1 <- update(gca_mod_aes_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_aes_cond_2 <- update(gca_mod_aes_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_aes_cond_3 <- update(gca_mod_aes_cond_2,   . ~ . + ot3:condition_sum) # singular

aes_cond_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_cond_0, gca_mod_aes_cond_1,
        gca_mod_aes_cond_2, gca_mod_aes_cond_3) 

# WM homog: none significant
#
# only ies
#

gca_mod_ies_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ies")) # singular

# add cond effect to iesercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_cond_0 <- update(gca_mod_ies_base,   . ~ . + condition_sum) # singular
gca_mod_ies_cond_1 <- update(gca_mod_ies_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ies_cond_2 <- update(gca_mod_ies_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ies_cond_3 <- update(gca_mod_ies_cond_2,   . ~ . + ot3:condition_sum) # singular

ies_cond_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_cond_0, gca_mod_ies_cond_1,    # WM homog: no sig
        gca_mod_ies_cond_2, gca_mod_ies_cond_3)


#
# only ams
#

gca_mod_ams_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ams")) # singular

# add cond effect to amsercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_cond_0 <- update(gca_mod_ams_base,   . ~ . + condition_sum) # singular
gca_mod_ams_cond_1 <- update(gca_mod_ams_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ams_cond_2 <- update(gca_mod_ams_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ams_cond_3 <- update(gca_mod_ams_cond_2,   . ~ . + ot3:condition_sum) # singular

ams_cond_anova <-
  anova(gca_mod_ams_base, gca_mod_ams_cond_0, gca_mod_ams_cond_1,
        gca_mod_ams_cond_2, gca_mod_ams_cond_3) # WM homog: none significant


#
# only ims
#

gca_mod_ims_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ims")) # singular

# add cond effect to imsercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_cond_0 <- update(gca_mod_ims_base,   . ~ . + condition_sum) # singular
gca_mod_ims_cond_1 <- update(gca_mod_ims_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ims_cond_2 <- update(gca_mod_ims_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ims_cond_3 <- update(gca_mod_ims_cond_2,   . ~ . + ot3:condition_sum) # singular

ims_cond_anova <-
  anova(gca_mod_ims_base, gca_mod_ims_cond_0, gca_mod_ims_cond_1, # none singular, and none significant
        gca_mod_ims_cond_2, gca_mod_ims_cond_3) # WM hom: none significant


# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_full_mod_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +          # El modelo original es (ot) * condition, but idk why and diff not significant
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 2e4)),
       data = stress_gc_subset, REML = F) # singular

# 
# gca_full_mod_base1 <-
#   lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + condition_sum +
#          (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
#          (1 + ot1 + ot2 + ot3 | target),
#        control = lmerControl(optimizer = 'bobyqa',
#                              optCtrl = list(maxfun = 2e4)),
#        data = stress_gc_subset, REML = F) # singular
# # 
# gca_full_mod_base2 <-
#   lmer(eLog ~ 1 + (ot1 + ot2 + ot3) * condition_sum +       # Este modelo es el original: por qué interacción?
#          (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
#          (1 + ot1 + ot2 + ot3 | target),
#        control = lmerControl(optimizer = 'bobyqa',
#                              optCtrl = list(maxfun = 2e4)),
#        data = stress_gc_subset, REML = F) # singular
# 
# anova(gca_full_mod_base, gca_full_mod_base1, gca_full_mod_base2) # WM hom: none significant

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_group_0 <- update(gca_full_mod_base,    . ~ . + group) # singular
gca_full_mod_group_1 <- update(gca_full_mod_group_0, . ~ . + ot1:group) # singular
gca_full_mod_group_2 <- update(gca_full_mod_group_1, . ~ . + ot2:group) # singular
gca_full_mod_group_3 <- update(gca_full_mod_group_2, . ~ . + ot3:group) # singular

full_group_anova <-
  anova(gca_full_mod_base, gca_full_mod_group_0, gca_full_mod_group_1,
        gca_full_mod_group_2, gca_full_mod_group_3)
## All
#                      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    

## WM homog
#                      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_full_mod_base    30 191702 191956 -95821   191642                              
# gca_full_mod_group_0 34 191684 191973 -95808   191616 25.3292      4  4.320e-05 ***
# gca_full_mod_group_1 38 191686 192008 -95805   191610  6.1765      4     0.1864    
# gca_full_mod_group_2 42 191689 192045 -95802   191605  5.0043      4     0.2869    
# gca_full_mod_group_3 46 191661 192051 -95784   191569 36.2443      4  2.577e-07 ***
################################

# add condition effect to intercept, linear slope, quadratic, and cubic time terms

## only run w/ wm homog
gca_full_mod_cond_0 <- update(gca_full_mod_group_3, . ~ . + condition_sum) # singular
gca_full_mod_cond_1 <- update(gca_full_mod_cond_0, . ~ . + ot1:condition_sum) # singular
gca_full_mod_cond_2 <- update(gca_full_mod_cond_1, . ~ . + ot2:condition_sum) # singular
gca_full_mod_cond_3 <- update(gca_full_mod_cond_2, . ~ . + ot3:condition_sum) # singular

full_cond_anova <-
  anova(gca_full_mod_group_3, gca_full_mod_cond_0, gca_full_mod_cond_1,
        gca_full_mod_cond_2, gca_full_mod_cond_3)
# WM homog
#                      Df    AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_full_mod_group_2 46 191661 192051 -95784   191569                         
# gca_full_mod_cond_0  47 191662 192061 -95784   191568 0.7111      1     0.3991
# gca_full_mod_cond_1  48 191664 192071 -95784   191568 0.1660      1     0.6837
# gca_full_mod_cond_2  49 191666 192082 -95784   191568 0.0107      1     0.9176
# gca_full_mod_cond_3  50 191666 192090 -95783   191566 1.6603      1     0.1976

################################

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_int_0 <- update(gca_full_mod_group_3, . ~ . + condition_sum:group) # singular        
gca_full_mod_int_1 <- update(gca_full_mod_int_0,   . ~ . + ot1:condition_sum:group) # singular
gca_full_mod_int_2 <- update(gca_full_mod_int_1,   . ~ . + ot2:condition_sum:group) # singular
gca_full_mod_int_3 <- update(gca_full_mod_int_2,   . ~ . + ot3:condition_sum:group) # singular

full_int_anova <-
  anova(gca_full_mod_group_3, gca_full_mod_int_0, gca_full_mod_int_1,
        gca_full_mod_int_2, gca_full_mod_int_3)
## All
#                      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    

## WM homg
#                      Df    AIC    BIC logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_full_mod_group_3 46 191661 192051 -95784   191569                             
# gca_full_mod_int_0   51 191668 192101 -95783   191566  2.4344      5    0.78635  
# gca_full_mod_int_1   56 191669 192144 -95779   191557  9.1347      5    0.10381  
# gca_full_mod_int_2   61 191674 192192 -95776   191552  4.9361      5    0.42372  
# gca_full_mod_int_3   66 191673 192233 -95771   191541 10.9550      5    0.05228 .
---

# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ams"))
gca_full_mod_gr3_relevel <- update(gca_full_mod_group_3)

}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, condition_sum) %>%
  distinct

# Get model predictions and SE
fits_all <- predictSE(gca_full_mod_group_3, new_dat_all) %>%        #change depending on significance
  as_tibble %>%
  bind_cols(new_dat_all) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, SS = "mon", AE = "aes", IE = "ies", AM = "ams", IM = "ims"))

# Filter preds at target offset
target_offset_preds <- filter(fits_all, time_zero == 4) %>%
  select(group, cond = condition_sum,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {
# Build model names programatically
mod_type <- "gca_mod_"
mod_spec <- c("_base", "_cond_0",
              "_cond_1", "_cond_2", "_cond_3", "_int0", "_int1", "_int2",
              "_int3")

# Store ind models in list
ind_mods <- mget(c(paste0(mod_type, "mon", mod_spec),
                   paste0(mod_type, "aes", mod_spec),
                   paste0(mod_type, "ies", mod_spec),
                   paste0(mod_type, "ams", mod_spec),
                   paste0(mod_type, "ims", mod_spec)
                   ))

save(ind_mods,
     file = here("reports", "mods", "gca",
                 "ind_mods.Rdata"))

# Store full (ot1, ot2, ot3, group, coda, cond) models in list
full_mods <- mget(c(
  "gca_full_mod_base", "gca_full_mod_group_0", "gca_full_mod_group_1",
  "gca_full_mod_group_2", "gca_full_mod_group_3", "gca_full_mod_int_0",
  "gca_full_mod_int_1", "gca_full_mod_int_2", "gca_full_mod_int_3",
  "gca_full_mod_int_relevel"))

save(full_mods,
     file = here("reports", "mods", "gca",
                 "full_mods.Rdata"))

# Save anova model comparisons
nested_model_comparisons <-
  mget(c("mon_cond_anova", #"mon_int_anova",
         "aes_cond_anova", #"aes_int_anova",
         "ies_cond_anova", #"ies_int_anova",
         "ams_cond_anova", #"ams_int_anova",
         "ims_cond_anova", #"ims_int_anova",
         "full_group_anova", "full_int_anova"))

save(nested_model_comparisons,
     file = here("reports", "mods", "gca",
                 "nested_model_comparisons.Rdata"))

# Save models predictions
model_preds <- mget(c("fits_all", "target_offset_preds"))

save(model_preds,
     file = here("reports", "mods", "gca",
                 "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------

