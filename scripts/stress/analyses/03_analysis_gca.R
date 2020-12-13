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

# Get path to saved models
gca_mods_path  <- here("mods", "stress", "gca")

# Load models as lists
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

wm <- read_csv("./data/clean/ospan_set_z_scores.csv")

stress_50 <- merge(x = stress50, y = wm, by = "participant", all.x=TRUE)

stress_50 <- na.omit(stress_50)



stress_gc_subset <- stress_50 %>%
  select(., -WM_set) %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
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
         control = lmerControl(optimizer = 'bobyqa'),   # , optCtrl=list(maxfun=2e5)
         data = stress_gc_subset, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + condition_sum + ot1 | participant) +
             ot2 + (1 + condition_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + condition_sum + ot1 + ot2 | participant) +
             ot3 + (1 + condition_sum + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  
# We retain the most complex model: mod_ot4
  #   npar (Df?)    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    9 240042 240120 -120012   240024                          
  # mod_ot2   14 239367 239489 -119670   239339 684.486  5  < 2.2e-16 ***
  # mod_ot3   20 239311 239484 -119635   239271  68.534  6  8.167e-13 ***
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  mod_ot0 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + ot1 | target) +
                       + (1 + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + ot1 + ot2 | target) +
                       + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot3    20 239311 239484 -119635   239271                          
  # mod_ot0    21 239043 239225 -119501   239001 269.485  1  < 2.2e-16 ***
  # mod_ot1a   23 238962 239162 -119458   238916  84.849  2  < 2.2e-16 ***
  # mod_ot2a   26 238958 239183 -119453   238906  10.359  3    0.01575 *  
  # mod_ot3a   30 238932 239192 -119436   238872  34.373  4  6.248e-07 ***
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_full_mod_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +         
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target), 
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 2e4)),
       data = stress_gc_subset, REML = F)

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_group_0 <- update(gca_full_mod_base,   . ~ . + group) 
gca_full_mod_group_1 <- update(gca_full_mod_group_0, . ~ . + ot1:group) 
gca_full_mod_group_2 <- update(gca_full_mod_group_1, . ~ . + ot2:group) 
gca_full_mod_group_3 <- update(gca_full_mod_group_2, . ~ . + ot3:group)

full_group_anova <-
  anova(gca_full_mod_base, gca_full_mod_group_0, gca_full_mod_group_1,
        gca_full_mod_group_2, gca_full_mod_group_3)
#                        Df    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
# gca_full_mod_base      30 238877 239137 -119408   238817                          
# gca_full_mod_group_0   34 238856 239151 -119394   238788 28.8172  4  8.515e-06 ***
# gca_full_mod_group_1   38 238852 239182 -119388   238776 11.5783  4    0.02078 *  
# gca_full_mod_group_2   42 238830 239194 -119373   238746 30.0158  4  4.858e-06 ***
# gca_full_mod_group_3   46 238829 239228 -119368   238737  9.2664  4    0.05478 .   


# add stress effect to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_stress_0 <- update(gca_full_mod_group_2,   . ~ . + condition_sum) 
gca_full_mod_stress_1 <- update(gca_full_mod_stress_0, . ~ . + ot1:condition_sum) 
gca_full_mod_stress_2 <- update(gca_full_mod_stress_1, . ~ . + ot2:condition_sum) # singular
gca_full_mod_stress_3 <- update(gca_full_mod_stress_2, . ~ . + ot3:condition_sum)

full_stress_anova <-
  anova(gca_full_mod_group_2, gca_full_mod_stress_0, gca_full_mod_stress_1,
        gca_full_mod_stress_2, gca_full_mod_stress_3)
#                       npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
# gca_full_mod_group_2    42 238830 239194 -119373   238746                       
# gca_full_mod_stress_0   43 238832 239205 -119373   238746 0.2368  1    0.62650  
# gca_full_mod_stress_1   44 238832 239214 -119372   238744 1.6426  1    0.19997  
# gca_full_mod_stress_2   45 238830 239220 -119370   238740 4.7536  1    0.02924 *
# gca_full_mod_stress_3   46 238830 239229 -119369   238738 1.5118  1    0.21887  

gca_full_mod_wm_0 <- update(gca_full_mod_stress_2, . ~ . + ospan) 
gca_full_mod_wm_1 <- update(gca_full_mod_wm_0, . ~ . + ot1:ospan) 
gca_full_mod_wm_2 <- update(gca_full_mod_wm_1, . ~ . + ot2:ospan) 
gca_full_mod_wm_3 <- update(gca_full_mod_wm_2, . ~ . + ot3:ospan)

full_wm_anova <-
  anova(gca_full_mod_stress_2, gca_full_mod_wm_0, gca_full_mod_wm_1,
        gca_full_mod_wm_2, gca_full_mod_wm_3)
#                       npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_full_mod_stress_2   45 238830 239220 -119370   238740                     
# gca_full_mod_wm_0       46 238830 239229 -119369   238738 1.5110  1     0.2190
# gca_full_mod_wm_1       47 238832 239240 -119369   238738 0.0182  1     0.8925
# gca_full_mod_wm_2       48 238833 239249 -119368   238737 1.3959  1     0.2374
# gca_full_mod_wm_3       49 238834 239259 -119368   238736 0.1762  1     0.6747


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_int_0 <- update(gca_full_mod_stress_2, . ~ . + group:condition_sum)  
gca_full_mod_int_1 <- update(gca_full_mod_int_0, . ~ . + ot1:group:condition_sum) # singular
gca_full_mod_int_2 <- update(gca_full_mod_int_1, . ~ . + ot2:group:condition_sum) 
gca_full_mod_int_3 <- update(gca_full_mod_int_2, . ~ . + ot3:group:condition_sum) #singular

full_int_anova <- anova(gca_full_mod_stress_2, gca_full_mod_int_0, gca_full_mod_int_1, 
                         gca_full_mod_int_2, gca_full_mod_int_3)
#                       npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)   
# gca_full_mod_stress_2   45 238830 239220 -119370   238740                         
# gca_full_mod_int_0      49 238834 239259 -119368   238736  3.7668  4   0.438483   
# gca_full_mod_int_1      53 238837 239297 -119366   238731  4.6380  4   0.326497   
# gca_full_mod_int_2      57 238828 239322 -119357   238714 16.9387  4   0.001987 **
# gca_full_mod_int_3      62 238837 239374 -119356   238713  1.4529  5   0.918437 




# when WM included for interaction, not converging, bc of that and lack of main interaction, WM removed
gca_full_mod_int_a <- update(gca_full_mod_stress_2, . ~ . + group:condition_sum:ospan)      
gca_full_mod_int_b <- update(gca_full_mod_int_a, . ~ . + ot1:group:condition_sum:ospan) # singular
gca_full_mod_int_c <- update(gca_full_mod_int_b, . ~ . + ot2:group:condition_sum:ospan, 
                             control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 4e5))) 
gca_full_mod_int_d <- update(gca_full_mod_int_c, . ~ . + ot3:group:condition_sum:ospan) #singular

full_int_anova <- anova(gca_full_mod_stress_2, gca_full_mod_int_a, gca_full_mod_int_b, 
                        gca_full_mod_int_c, gca_full_mod_int_d)





# ---

summary(gca_full_mod_int_2) # mon reference

# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ams"))
gca_full_mod_int_2_ams <- update(gca_full_mod_int_2)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "mon"))
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ims"))
gca_full_mod_int_2_ims <- update(gca_full_mod_int_2)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "mon"))
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "aes"))
gca_full_mod_int_2_aes <- update(gca_full_mod_int_2)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "mon"))
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ies"))
gca_full_mod_int_2_ies <- update(gca_full_mod_int_2) 


}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, condition_sum) %>%
  distinct

# Get model predictions and SE
fits_all <- predictSE(gca_full_mod_int_2, new_dat_all) %>%  
  as_tibble %>%
  bind_cols(new_dat_all) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, SS = "mon", AE = "aes", IE = "ies", AM = "ams", IM = "ims"),
         group = fct_relevel(group, "SS", "AE", "AM", "IE", "IM"))

# Filter preds at target offset
target_offset_preds <- filter(fits_all, time_zero == 4) %>%
  select(group, stress = condition_sum,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {
# Build model names programatically
mod_type <- "gca_full_mod"
mod_spec <- c("_base", "_group_0", "_group_1", "_group_2", "_group_3",
              "_stress_0", "_stress_1", "_stress_2", "_stress_3",
              "_int_0", "_int_1", "_int_2", "_int_3",
              "_int_2_ams", "_int_2_ims", "_int_2_aes", "_int_2_ies")
  
# Store ind models in list
gca_full_mods <- mget(c(paste0(mod_type, mod_spec)))    # int_2 final model
  
save(gca_full_mods,
     file = here("mods", "stress", "gca",
                   "gca_full_mods.Rdata"))
  

# # final model
# saveRDS(gca_full_mod_int_2, file = here("mods", "stress", "gca", "final_model.rds"))

# Save anova model comparisons
nested_model_comparisons <-
  mget(c(#"mon_cond_anova", "mon_wm_anova", "mon_int_anova",
  #        "aes_cond_anova", "aes_wm_anova", "aes_int_anova",
  #        "ies_cond_anova", "ies_wm_anova", "ies_int_anova",
  #        "ams_cond_anova", "ams_wm_anova", "ams_int_anova",
  #        "ims_cond_anova", "ims_wm_anova", "ims_int_anova",
         "full_group_anova", "full_stress_anova", "full_int_anova"
         ))

save(nested_model_comparisons,
     file = here("mods", "stress", "gca",
                 "nested_model_comparisons.Rdata"))

# Save models predictions
model_preds <- mget(c("fits_all", "target_offset_preds"))

save(model_preds,
     file = here("mods", "stress", "gca",
                 "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------





### NOT UPDATED


################################

# # add 3-way int to intercept, linear slope, quadratic, and cubic time terms
# gca_full_mod_wm_0 <- update(gca_full_mod_group_3, . ~ . + WM_set:condition_sum:group) # singular        
# gca_full_mod_wm_1 <- update(gca_full_mod_wm_0,   . ~ . + ot1:WM_set:condition_sum:group) # singular
# gca_full_mod_wm_2 <- update(gca_full_mod_wm_1,   . ~ . + ot2:WM_set:condition_sum:group) # singular
# gca_full_mod_wm_3 <- update(gca_full_mod_wm_2,   . ~ . + ot3:WM_set:condition_sum:group) # singular
# 
# full_wm_anova <-
#   anova(gca_full_mod_group_3, gca_full_mod_wm_0, gca_full_mod_wm_1,
#         gca_full_mod_wm_2, gca_full_mod_wm_3)
# #                      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# # 
# 
# 
# mod_spec <- c("_wm_0", "_wm_1", "_wm_2", "_wm_3")
# 
# # Store ind models in list
# full_mods_wm <- mget(c(paste0(mod_type, mod_spec)))
# 
# save(full_mods_wm,
#      file = here("mods", "stress", "gca",
#                  "full_mods_wm.Rdata"))







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

# save(gca_mod_mon_base,
#      file = here("mods", "stress", "gca", "gca_mod_mon_base.Rdata"))

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + condition_sum) # singular
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2,   . ~ . + ot3:condition_sum) # singular

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2, gca_mod_mon_cond_3)
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base     30 44237 44447 -22088    44177                     
# gca_mod_mon_cond_0   31 44239 44456 -22088    44177 0.2026  1     0.6526
# gca_mod_mon_cond_1   32 44240 44464 -22088    44176 0.5423  1     0.4615
# gca_mod_mon_cond_2   33 44242 44473 -22088    44176 0.0079  1     0.9290
# gca_mod_mon_cond_3   34 44244 44482 -22088    44176 0.0194  1     0.8894


# mod_type <- "gca_mod_"
# mod_spec <- c('_base', "_cond_0", "_cond_1", "_cond_2", "_cond_3")
# 
# # Store ind models in list
# mon_cond_mods <- mget(c(paste0(mod_type, "mon", mod_spec)))
# 
# save(mon_cond_mods,
#      file = here("mods", "stress", "gca",
#                  "mon_cond_mods.Rdata"))





# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_wm_0 <- update(gca_mod_mon_base, . ~ . + WM_set) # singular
gca_mod_mon_wm_1 <- update(gca_mod_mon_wm_0,   . ~ . + ot1:WM_set) # singular
gca_mod_mon_wm_2 <- update(gca_mod_mon_wm_1,   . ~ . + ot2:WM_set) # singular
gca_mod_mon_wm_3 <- update(gca_mod_mon_wm_2,   . ~ . + ot3:WM_set) # singular

mon_wm_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_wm_0, gca_mod_mon_wm_1,
        gca_mod_mon_wm_2, gca_mod_mon_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base   


# mod_type <- "gca_mod_"
# mod_spec <- c("_wm_0", "_wm_1", "_wm_2", "_wm_3")
# 
# # Store ind models in list
# mon_wm_mods <- mget(c(paste0(mod_type, "mon", mod_spec)))
# 
# save(mon_wm_mods,
#      file = here("mods", "stress", "gca",
#                  "mon_wm_mods.Rdata"))





# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_int_0 <- update(gca_mod_mon_base, . ~ . + condition_sum:WM_set) # singular
gca_mod_mon_int_1 <- update(gca_mod_mon_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_mon_int_2 <- update(gca_mod_mon_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_mon_int_3 <- update(gca_mod_mon_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

mon_int_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_int_0, gca_mod_mon_int_1,
        gca_mod_mon_int_2, gca_mod_mon_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_mon_base  30 


# mod_type <- "gca_mod_"
# mod_spec <- c("_int_0", "_int_1", "_int_2", "_int_3")
# 
# # Store ind models in list
# mon_int_mods <- mget(c(paste0(mod_type, "mon", mod_spec)))
# 
# save(mon_int_mods,
#      file = here("mods", "stress", "gca",
#                  "mon_int_mods.Rdata"))





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
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_aes_base     30 30 47567 47779 -23753    47507                       
# gca_mod_aes_cond_0   31 47568 47787 -23753    47506 0.9101  1    0.34008  
# gca_mod_aes_cond_1   32 47570 47796 -23753    47506 0.0271  1    0.86917  
# gca_mod_aes_cond_2   33 47565 47799 -23750    47499 6.2033  1    0.01275 *
# gca_mod_aes_cond_3   34 47567 47808 -23750    47499 0.1145  1    0.73506  


# mod_type <- "gca_mod_"
# mod_spec <- c('_base', "_cond_0", "_cond_1", "_cond_2", "_cond_3")
# 
# # Store ind models in list
# aes_cond_mods <- mget(c(paste0(mod_type, "aes", mod_spec)))
# 
# save(aes_cond_mods,
#      file = here("mods", "stress", "gca",
#                  "aes_cond_mods.Rdata"))


# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_wm_0 <- update(gca_mod_aes_cond_3, . ~ . + WM_set) # singular
gca_mod_aes_wm_1 <- update(gca_mod_aes_wm_0,   . ~ . + ot1:WM_set) # singular
gca_mod_aes_wm_2 <- update(gca_mod_aes_wm_1,   . ~ . + ot2:WM_set) # singular
gca_mod_aes_wm_3 <- update(gca_mod_aes_wm_2,   . ~ . + ot3:WM_set) # singular

aes_wm_anova <-
  anova(gca_mod_aes_cond_3, gca_mod_aes_wm_0, gca_mod_aes_wm_1,
        gca_mod_aes_wm_2, gca_mod_aes_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_aes_cond_3 34    


# mod_type <- "gca_mod_"
# mod_spec <- c("_wm_0", "_wm_1", "_wm_2", "_wm_3")
# 
# # Store ind models in list
# aes_wm_mods <- mget(c(paste0(mod_type, "aes", mod_spec)))
# 
# save(aes_wm_mods,
#      file = here("mods", "stress", "gca",
#                  "aes_wm_mods.Rdata"))





# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_int_0 <- update(gca_mod_aes_cond_3, . ~ . + condition_sum:WM_set) # singular
gca_mod_aes_int_1 <- update(gca_mod_aes_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_aes_int_2 <- update(gca_mod_aes_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_aes_int_3 <- update(gca_mod_aes_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

aes_int_anova <-
  anova(gca_mod_aes_cond_3, gca_mod_aes_int_0, gca_mod_aes_int_1,
        gca_mod_aes_int_2, gca_mod_aes_int_3)
#                     Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_aes_cond_3  34 


# mod_type <- "gca_mod_"
# mod_spec <- c("_int_0", "_int_1", "_int_2", "_int_3")
# 
# # Store ind models in list
# aes_int_mods <- mget(c(paste0(mod_type, "aes", mod_spec)))
# 
# save(aes_int_mods,
#      file = here("mods", "stress", "gca",
#                  "aes_int_mods.Rdata"))





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
  anova(gca_mod_ies_base, gca_mod_ies_cond_0, gca_mod_ies_cond_1,    
        gca_mod_ies_cond_2, gca_mod_ies_cond_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_base   30 49428 49641 -24684    49368                       
# gca_mod_ies_cond_0   31 49430 49650 -24684    49368 0.1038  1    0.74729  
# gca_mod_ies_cond_1   32 49432 49659 -24684    49368 0.3824  1    0.53632  
# gca_mod_ies_cond_2   33 49434 49668 -24684    49368 0.0885  1    0.76606  
# gca_mod_ies_cond_3   34 49432 49673 -24682    49364 3.4647  1    0.06269 .
# 
# save(ies_cond_mods,
#      file = here("mods", "stress", "gca",
#                  "ies_cond_mods.Rdata"))




# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_wm_0 <- update(gca_mod_ies_base, . ~ . + WM_set) # singular
gca_mod_ies_wm_1 <- update(gca_mod_ies_wm_0,   . ~ . + ot1:WM_set) # singular
gca_mod_ies_wm_2 <- update(gca_mod_ies_wm_1,   . ~ . + ot2:WM_set) # singular
gca_mod_ies_wm_3 <- update(gca_mod_ies_wm_2,   . ~ . + ot3:WM_set) # singular

ies_wm_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_wm_0, gca_mod_ies_wm_1,
        gca_mod_ies_wm_2, gca_mod_ies_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_base   30 


# mod_type <- "gca_mod_"
# mod_spec <- c("_wm_0", "_wm_1", "_wm_2", "_wm_3")
# 
# # Store ind models in list
# ies_wm_mods <- mget(c(paste0(mod_type, "ies", mod_spec)))
# 
# save(ies_wm_mods,
#      file = here("mods", "stress", "gca",
#                  "ies_wm_mods.Rdata"))




# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_int_0 <- update(gca_mod_ies_base, . ~ . + condition_sum:WM_set) # singular
gca_mod_ies_int_1 <- update(gca_mod_ies_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_ies_int_2 <- update(gca_mod_ies_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_ies_int_3 <- update(gca_mod_ies_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

ies_int_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_int_0, gca_mod_ies_int_1,
        gca_mod_ies_int_2, gca_mod_ies_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ies_base  30 

# mod_type <- "gca_mod_"
# mod_spec <- c("_int_0", "_int_1", "_int_2", "_int_3")
# 
# # Store ind models in list
# ies_int_mods <- mget(c(paste0(mod_type, "ies", mod_spec)))
# 
# save(ies_int_mods,
#      file = here("mods", "stress", "gca",
#                  "ies_int_mods.Rdata"))




#
# only ams
#

gca_mod_ams_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ams"))

# add cond effect to amsercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_cond_0 <- update(gca_mod_ams_base,   . ~ . + condition_sum)
gca_mod_ams_cond_1 <- update(gca_mod_ams_cond_0,   . ~ . + ot1:condition_sum)
gca_mod_ams_cond_2 <- update(gca_mod_ams_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ams_cond_3 <- update(gca_mod_ams_cond_2,   . ~ . + ot3:condition_sum) # singular

ams_cond_anova <-
  anova(gca_mod_ams_base, gca_mod_ams_cond_0, gca_mod_ams_cond_1,
        gca_mod_ams_cond_2, gca_mod_ams_cond_3) 
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_base     30 48201 48413 -24070    48141                       
# gca_mod_ams_cond_0   31 48203 48422 -24070    48141 0.0002  1    0.99002  
# gca_mod_ams_cond_1   32 48205 48431 -24070    48141 0.0573  1    0.81080  
# gca_mod_ams_cond_2   33 48201 48434 -24068    48135 5.7881  1    0.01613 *
# gca_mod_ams_cond_3   34 48203 48443 -24068    48135 0.0116  1    0.91414

# mod_type <- "gca_mod_"
# mod_spec <- c("_base", "_cond_0", "_cond_1", "_cond_2", "_cond_3")
# 
# # Store ind models in list
# ams_cond_mods <- mget(c(paste0(mod_type, "ams", mod_spec)))
# 
# save(ams_cond_mods,
#      file = here("mods", "stress", "gca",
#                  "ams_cond_mods.Rdata"))





# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_wm_0 <- update(gca_mod_ams_cond_3, . ~ . + WM_set) # singular
gca_mod_ams_wm_1 <- update(gca_mod_ams_wm_0,   . ~ . + ot1:WM_set) # singular
gca_mod_ams_wm_2 <- update(gca_mod_ams_wm_1,   . ~ . + ot2:WM_set) # singular
gca_mod_ams_wm_3 <- update(gca_mod_ams_wm_2,   . ~ . + ot3:WM_set) # singular

ams_wm_anova <-
  anova(gca_mod_ams_cond_3, gca_mod_ams_wm_0, gca_mod_ams_wm_1,
        gca_mod_ams_wm_2, gca_mod_ams_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_cond_3 34 


# mod_type <- "gca_mod_"
# mod_spec <- c("_wm_0", "_wm_1", "_wm_2", "_wm_3")
# 
# # Store ind models in list
# ams_wm_mods <- mget(c(paste0(mod_type, "ams", mod_spec)))
# 
# save(ams_wm_mods,
#      file = here("mods", "stress", "gca",
#                  "ams_wm_mods.Rdata"))




# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_int_0 <- update(gca_mod_ams_cond_3, . ~ . + condition_sum:WM_set) # singular
gca_mod_ams_int_1 <- update(gca_mod_ams_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_ams_int_2 <- update(gca_mod_ams_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_ams_int_3 <- update(gca_mod_ams_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

ams_int_anova <-
  anova(gca_mod_ams_cond_3, gca_mod_ams_int_0, gca_mod_ams_int_1,
        gca_mod_ams_int_2, gca_mod_ams_int_3)
#                     Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ams_cond_3  34 

# mod_type <- "gca_mod_"
# mod_spec <- c("_int_0", "_int_1", "_int_2", "_int_3")
# 
# # Store ind models in list
# ams_int_mods <- mget(c(paste0(mod_type, "ams", mod_spec)))
# 
# save(ams_int_mods,
#      file = here("mods", "stress", "gca",
#                  "ams_int_mods.Rdata"))




#
# only ims
#

gca_mod_ims_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ims")) 

# add cond effect to imsercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_cond_0 <- update(gca_mod_ims_base,   . ~ . + condition_sum) 
gca_mod_ims_cond_1 <- update(gca_mod_ims_cond_0,   . ~ . + ot1:condition_sum)
gca_mod_ims_cond_2 <- update(gca_mod_ims_cond_1,   . ~ . + ot2:condition_sum)
gca_mod_ims_cond_3 <- update(gca_mod_ims_cond_2,   . ~ . + ot3:condition_sum)

ims_cond_anova <-
  anova(gca_mod_ims_base, gca_mod_ims_cond_0, gca_mod_ims_cond_1, 
        gca_mod_ims_cond_2, gca_mod_ims_cond_3) 
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_base   30 48559 48771 -24250    48499                     
# gca_mod_ims_cond_0   31 48561 48780 -24250    48499 0.1040  1     0.7471
# gca_mod_ims_cond_1   32 48562 48788 -24249    48498 1.5108  1     0.2190
# gca_mod_ims_cond_2   33 48562 48796 -24248    48496 1.4883  1     0.2225
# gca_mod_ims_cond_3   34 48563 48803 -24247    48495 1.6916  1     0.1934


# mod_type <- "gca_mod_"
# mod_spec <- c("_base", "_cond_0", "_cond_1", "_cond_2", "_cond_3")
# 
# # Store ind models in list
# ims_cond_mods <- mget(c(paste0(mod_type, "ims", mod_spec)))
# 
# save(ims_cond_mods,
#      file = here("mods", "stress", "gca",
#                  "ims_cond_mods.Rdata"))



# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_wm_0 <- update(gca_mod_ims_base, . ~ . + WM_set) # singular
gca_mod_ims_wm_1 <- update(gca_mod_ims_wm_0,   . ~ . + ot1:WM_set) # singular
gca_mod_ims_wm_2 <- update(gca_mod_ims_wm_1,   . ~ . + ot2:WM_set) # singular
gca_mod_ims_wm_3 <- update(gca_mod_ims_wm_2,   . ~ . + ot3:WM_set) # singular

ims_wm_anova <-
  anova(gca_mod_ims_base, gca_mod_ims_wm_0, gca_mod_ims_wm_1,
        gca_mod_ims_wm_2, gca_mod_ims_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_base   30 

# mod_type <- "gca_mod_"
# mod_spec <- c("_wm_0", "_wm_1", "_wm_2", "_wm_3")

# # Store ind models in list
# ims_wm_mods <- mget(c(paste0(mod_type, "ims", mod_spec)))
# 
# save(ims_wm_mods,
#      file = here("mods", "stress", "gca",
#                  "ims_wm_mods.Rdata"))


# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_int_0 <- update(gca_mod_ims_wm_3, . ~ . + condition_sum:WM_set) # singular
gca_mod_ims_int_1 <- update(gca_mod_ims_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_ims_int_2 <- update(gca_mod_ims_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_ims_int_3 <- update(gca_mod_ims_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

ims_int_anova <-
  anova(gca_mod_ims_wm_3, gca_mod_ims_int_0, gca_mod_ims_int_1,
        gca_mod_ims_int_2, gca_mod_ims_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ims_wm_3  34


mod_type <- "gca_mod_"
mod_group <- c("mon", "aes", "ies", "ams", "ims")
mod_spec <- c("_base", "_cond_0", "_cond_1", "_cond_2", "_cond_3")

# Store ind models in list
group_mods <- mget(c(paste0(mod_type, mod_group, mod_spec)))

save(group_mods,
     file = here("mods", "stress", "gca",
                 "group_mods.Rdata"))


group_anovas <- mget(c(paste0(mod_group, "_cond_anova")))

save(group_anovas,
     file = here("mods", "stress", "gca",
                 "group_mod_anovas.Rdata"))


# -----------------------------------------------------------------------------
