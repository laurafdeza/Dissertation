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

# source(here::here("scripts", "00_load_libs.R"))

# Load data
source(here::here("scripts", "02_load_data.R"))

# Get path to saved models
gca_mods_path  <- here("mods", "stress", "gca")

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
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "mon", "aes", "ams", "ies", "ims"),
            condition_sum = if_else(cond == "Present", 1, -1)) %>%       # 1 = present, 2 = past
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
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3, mod_ot4)
  
# We retain the most complex model: mod_ot4
  #         Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
  # mod_ot1  9 235958 236036 -117970   235940                             
  # mod_ot2 14 235784 235905 -117878   235756 184.42      5  < 2.2e-16 ***
  # mod_ot3 20 235260 235434 -117610   235220 535.35      6  < 2.2e-16 ***
  # mod_ot4 21 234752 234934 -117355   234710 510.32      1  < 2.2e-16 ***
  # ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
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
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base   30 42753 42963 -21346    42693                         
# gca_mod_mon_cond_0 31 42755 42972 -21346    42693 0.4824      1     0.4873
# gca_mod_mon_cond_1 32 42756 42980 -21346    42692 0.1857      1     0.6665
# gca_mod_mon_cond_2 33 42758 42989 -21346    42692 0.0000      1     0.9998
# gca_mod_mon_cond_3 34 42760 42998 -21346    42692 0.0628      1     0.8021


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
# gca_mod_mon_base   30 42753 42963 -21346    42693                         
# gca_mod_mon_wm_0   31 42755 42972 -21346    42693 0.4829      1     0.4871
# gca_mod_mon_wm_1   32 42756 42980 -21346    42692 0.4634      1     0.4960
# gca_mod_mon_wm_2   33 42758 42989 -21346    42692 0.3297      1     0.5658
# gca_mod_mon_wm_3   34 42759 42997 -21346    42691 0.4147      1     0.5196


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
# gca_mod_mon_base  30 42753 42963 -21346    42693                         
# gca_mod_mon_int_0 31 42755 42972 -21346    42693 0.5023      1     0.4785
# gca_mod_mon_int_1 32 42756 42980 -21346    42692 0.1085      1     0.7418
# gca_mod_mon_int_2 33 42757 42988 -21346    42691 1.3188      1     0.2508
# gca_mod_mon_int_3 34 42759 42997 -21346    42691 0.1465      1     0.7019


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
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_aes_base   30 45808 46020 -22874    45748                           
# gca_mod_aes_cond_0 31 45809 46028 -22873    45747 0.9178      1    0.33805  
# gca_mod_aes_cond_1 32 45810 46037 -22873    45746 0.1099      1    0.74027  
# gca_mod_aes_cond_2 33 45812 46045 -22873    45746 0.5609      1    0.45389  
# gca_mod_aes_cond_3 34 45809 46049 -22870    45741 4.8460      1    0.02771 *


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
# gca_mod_aes_cond_3 34 45809 46049 -22870    45741                         
# gca_mod_aes_wm_0   35 45811 46058 -22870    45741 0.0132      1     0.9086
# gca_mod_aes_wm_1   36 45813 46068 -22870    45741 0.0013      1     0.9707
# gca_mod_aes_wm_2   37 45815 46076 -22870    45741 0.5454      1     0.4602
# gca_mod_aes_wm_3   38 45814 46083 -22869    45738 2.4517      1     0.1174   


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
# gca_mod_aes_cond_3  34 45809 46049 -22870    45741                         
# gca_mod_aes_int_0   35 45809 46056 -22870    45739 2.0545      1     0.1518
# gca_mod_aes_int_1   36 45809 46064 -22869    45737 1.5953      1     0.2066
# gca_mod_aes_int_2   37 45811 46073 -22869    45737 0.0204      1     0.8863
# gca_mod_aes_int_3   38 45812 46081 -22868    45736 0.9360      1     0.3333


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
# gca_mod_ies_base   30 48125 48338 -24033    48065                         
# gca_mod_ies_cond_0 31 48127 48347 -24033    48065 0.1617      1     0.6876
# gca_mod_ies_cond_1 32 48129 48356 -24033    48065 0.0155      1     0.9010
# gca_mod_ies_cond_2 33 48129 48363 -24032    48063 2.0231      1     0.1549
# gca_mod_ies_cond_3 34 48131 48372 -24031    48063 0.5790      1     0.4467


# mod_type <- "gca_mod_"
# mod_spec <- c('_base', "_cond_0", "_cond_1", "_cond_2", "_cond_3")
# 
# # Store ind models in list
# ies_cond_mods <- mget(c(paste0(mod_type, "ies", mod_spec)))
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
# gca_mod_ies_base   30 48125 48338 -24033    48065                           
# gca_mod_ies_wm_0   31 48127 48346 -24032    48065 0.7025      1    0.40196  
# gca_mod_ies_wm_1   32 48125 48352 -24031    48061 3.2136      1    0.07303 .
# gca_mod_ies_wm_2   33 48126 48360 -24030    48060 1.3777      1    0.24049  
# gca_mod_ies_wm_3   34 48128 48369 -24030    48060 0.1095      1    0.74066


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
# gca_mod_ies_base  30 48125 48338 -24033    48065                         
# gca_mod_ies_int_0 31 48127 48346 -24032    48065 0.5965      1     0.4399
# gca_mod_ies_int_1 32 48128 48354 -24032    48064 1.1650      1     0.2804
# gca_mod_ies_int_2 33 48129 48363 -24031    48063 0.7127      1     0.3985
# gca_mod_ies_int_3 34 48131 48372 -24031    48063 0.2568      1     0.6124


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
       data = filter(stress_gc_subset, group == "ams")) # singular

# add cond effect to amsercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_cond_0 <- update(gca_mod_ams_base,   . ~ . + condition_sum) # singular
gca_mod_ams_cond_1 <- update(gca_mod_ams_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ams_cond_2 <- update(gca_mod_ams_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ams_cond_3 <- update(gca_mod_ams_cond_2,   . ~ . + ot3:condition_sum) # singular

ams_cond_anova <-
  anova(gca_mod_ams_base, gca_mod_ams_cond_0, gca_mod_ams_cond_1,
        gca_mod_ams_cond_2, gca_mod_ams_cond_3) 
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_base   30 46872 47084 -23406    46812                           
# gca_mod_ams_cond_0 31 46874 47093 -23406    46812 0.2855      1    0.59313  
# gca_mod_ams_cond_1 32 46876 47102 -23406    46812 0.1743      1    0.67634  
# gca_mod_ams_cond_2 33 46877 47110 -23406    46811 0.4518      1    0.50148  
# gca_mod_ams_cond_3 34 46875 47115 -23404    46807 4.1374      1    0.04195 *


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
# gca_mod_ams_cond_3 34 46875 47115 -23404    46807                         
# gca_mod_ams_wm_0   35 46876 47123 -23403    46806 0.9389      1     0.3326
# gca_mod_ams_wm_1   36 46876 47130 -23402    46804 2.1466      1     0.1429
# gca_mod_ams_wm_2   37 46877 47139 -23402    46803 0.7272      1     0.3938
# gca_mod_ams_wm_3   38 46878 47146 -23401    46802 1.5855      1     0.2080


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
# gca_mod_ams_cond_3  34 46875 47115 -23404    46807                            
# gca_mod_ams_int_0   35 46876 47123 -23403    46806 1.2657      1   0.260582   
# gca_mod_ams_int_1   36 46875 47130 -23402    46803 2.4152      1   0.120166   
# gca_mod_ams_int_2   37 46870 47132 -23398    46796 7.1744      1   0.007395 **
# gca_mod_ams_int_3   38 46869 47138 -23397    46793 2.7948      1   0.094571 . 


mod_type <- "gca_mod_"
mod_spec <- c("_int_0", "_int_1", "_int_2", "_int_3")

# Store ind models in list
ams_int_mods <- mget(c(paste0(mod_type, "ams", mod_spec)))

save(ams_int_mods,
     file = here("mods", "stress", "gca",
                 "ams_int_mods.Rdata"))




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
  anova(gca_mod_ims_base, gca_mod_ims_cond_0, gca_mod_ims_cond_1, 
        gca_mod_ims_cond_2, gca_mod_ims_cond_3) 
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_base   30 47683 47895 -23812    47623                         
# gca_mod_ims_cond_0 31 47685 47904 -23811    47623 0.6415      1     0.4232
# gca_mod_ims_cond_1 32 47686 47912 -23811    47622 0.1620      1     0.6873
# gca_mod_ims_cond_2 33 47687 47920 -23810    47621 1.4141      1     0.2344
# gca_mod_ims_cond_3 34 47689 47929 -23810    47621 0.0662      1     0.7969


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
# gca_mod_ims_base   30 47683 47895 -23812    47623                            
# gca_mod_ims_wm_0   31 47685 47904 -23811    47623 0.6301      1   0.427336   
# gca_mod_ims_wm_1   32 47685 47911 -23811    47621 1.3248      1   0.249738   
# gca_mod_ims_wm_2   33 47686 47919 -23810    47620 1.1515      1   0.283230   
# gca_mod_ims_wm_3   34 47681 47921 -23807    47613 6.8058      1   0.009086 **


mod_type <- "gca_mod_"
mod_spec <- c("_wm_0", "_wm_1", "_wm_2", "_wm_3")

# Store ind models in list
ims_wm_mods <- mget(c(paste0(mod_type, "ims", mod_spec)))

save(ims_wm_mods,
     file = here("mods", "stress", "gca",
                 "ims_wm_mods.Rdata"))


# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_int_0 <- update(gca_mod_ims_wm_3, . ~ . + condition_sum:WM_set) # singular
gca_mod_ims_int_1 <- update(gca_mod_ims_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_ims_int_2 <- update(gca_mod_ims_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_ims_int_3 <- update(gca_mod_ims_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

ims_int_anova <-
  anova(gca_mod_ims_wm_3, gca_mod_ims_int_0, gca_mod_ims_int_1,
        gca_mod_ims_int_2, gca_mod_ims_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ims_wm_3  34 47681 47921 -23807    47613                           
# gca_mod_ims_int_0 35 47683 47930 -23806    47613 0.4897      1    0.48405  
# gca_mod_ims_int_1 36 47685 47939 -23806    47613 0.2054      1    0.65043  
# gca_mod_ims_int_2 37 47684 47945 -23805    47610 2.9785      1    0.08438 .
# gca_mod_ims_int_3 38 47686 47954 -23805    47610 0.0155      1    0.90097  



# mod_type <- "gca_mod_"
# mod_spec <- c("_int_0", "_int_1", "_int_2", "_int_3")
# 
# # Store ind models in list
# ims_int_mods <- mget(c(paste0(mod_type, "ims", mod_spec)))
# 
# save(ims_int_mods,
#      file = here("mods", "stress", "gca",
#                  "ims_int_mods.Rdata"))


# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_full_mod_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) * condition_sum +         
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 2e4)),
       data = stress_gc_subset, REML = F)

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_group_0 <- update(gca_full_mod_base,   . ~ . + group) 
gca_full_mod_group_1 <- update(gca_full_mod_group_0, . ~ . + ot1:group) 
gca_full_mod_group_2 <- update(gca_full_mod_group_1, . ~ . + ot2:group) 
gca_full_mod_group_3 <- update(gca_full_mod_group_2, . ~ . + ot3:group) # singular

full_group_anova <-
  anova(gca_full_mod_base, gca_full_mod_group_0, gca_full_mod_group_1,
        gca_full_mod_group_2, gca_full_mod_group_3)
#                      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_full_mod_base    34 232118 232412 -116025   232050                              
# gca_full_mod_group_0 38 232099 232428 -116012   232023 26.4471      4  2.571e-05 ***
# gca_full_mod_group_1 42 232100 232464 -116008   232016  7.5214      4    0.11077    
# gca_full_mod_group_2 46 232096 232495 -116002   232004 11.3494      4    0.02291 *  
# gca_full_mod_group_3 50 232063 232496 -115981   231963 41.5468      4  2.071e-08 ***


mod_type <- "gca_full_mod"
mod_spec <- c("_base", "_group_0", "_group_1", "_group_2", "_group_3")

# Store ind models in list
full_mods_group <- mget(c(paste0(mod_type, mod_spec)))

save(full_mods_group,
     file = here("mods", "stress", "gca",
                 "full_mods_group.Rdata"))


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_int_0 <- update(gca_full_mod_group_3, . ~ . + group:condition_sum) # singular
gca_full_mod_int_1 <- update(gca_full_mod_int_0, . ~ . + ot1:group:condition_sum) # singular
gca_full_mod_int_2 <- update(gca_full_mod_int_1, . ~ . + ot2:group:condition_sum) # singular
gca_full_mod_int_3 <- update(gca_full_mod_int_2, . ~ . + ot3:group:condition_sum) #singular

full_int_anova <- anova(gca_full_mod_group_3, gca_full_mod_int_0, gca_full_mod_int_1, 
                         gca_full_mod_int_2, gca_full_mod_int_3)
#                      Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_full_mod_group_3 50 232063 232496 -115981   231963                           
# gca_full_mod_int_0   54 232068 232536 -115980   231960 2.6354      4    0.62056  
# gca_full_mod_int_1   58 232067 232570 -115975   231951 9.2215      4    0.05579 .
# gca_full_mod_int_2   62 232070 232607 -115973   231946 5.2637      4    0.26129  
# gca_full_mod_int_3   66 232068 232640 -115968   231936 9.4575      4    0.05063 .


mod_spec <- c("_int_0", "_int_1", "_int_2", "_int_3")

# Store ind models in list
full_mods_int <- mget(c(paste0(mod_type, mod_spec)))

save(full_mods_int,
     file = here("mods", "stress", "gca",
                 "full_mods_int.Rdata"))

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

# ---

summary(gca_full_mod_group_3) # mon reference

# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ams"))
gca_full_mod_group_3_ams <- update(gca_full_mod_group_3) # singular

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ims"))
gca_full_mod_group_3_ims <- update(gca_full_mod_group_3) # singular

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "aes"))
gca_full_mod_group_3_aes <- update(gca_full_mod_group_3) # singular

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ies"))
gca_full_mod_group_3_ies <- update(gca_full_mod_group_3) # singular

mod_spec <- c("_group_3_ams", "_group_3_ims", "_group_3_aes", "_group_3_ies")

# Store ind models in list
full_mods_refactor <- mget(c(paste0(mod_type, mod_spec)))

save(full_mods_refactor,
     file = here("mods", "stress", "gca",
                 "full_mods_refactor.Rdata"))

full_mod_group_relevel_anova <- anova(gca_full_mod_group_3, gca_full_mod_group_3_ams, gca_full_mod_group_3_ims,
      gca_full_mod_group_3_aes, gca_full_mod_group_3_ies)
#                          Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
# gca_full_mod_group_3     50 232063 232496 -115981   231963                                # mon
# gca_full_mod_group_3_ams 50 232063 232496 -115981   231963     0      0     <2e-16 ***
# gca_full_mod_group_3_ims 50 232063 232496 -115981   231963     0      0          1    
# gca_full_mod_group_3_aes 50 232063 232496 -115981   231963     0      0     <2e-16 ***
# gca_full_mod_group_3_ies 50 232063 232496 -115981   231963     0      0          1   

# We keep gca_full_mod_int_3 (mon reference) as final full model

}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, condition_sum) %>%
  distinct

# Get model predictions and SE
fits_all <- predictSE(gca_full_mod_group_3, new_dat_all) %>%  
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
# # Build model names programatically
# mod_type <- "gca_mod_"
# mod_spec <- c("_base", "_cond_0",
#               "_cond_1", "_cond_2", "_cond_3", 
#               "_wm_0", "_wm_1", "_wm_2", "_wm_3",
#               "_int_0", "_int_1", "_int_2",
#               "_int_3")
# 
# # Store ind models in list
# ind_mods <- mget(c(paste0(mod_type, "mon", mod_spec),
#                    paste0(mod_type, "aes", mod_spec),
#                    paste0(mod_type, "ies", mod_spec),
#                    paste0(mod_type, "ams", mod_spec),
#                    paste0(mod_type, "ims", mod_spec)
#                    ))
# 
# save(ind_mods,
#      file = here("mods", "stress", "gca",
#                  "ind_mods.Rdata"))
# 
# # Store full (ot1, ot2, ot3, group, coda, cond) models in list
# full_mods <- mget(c(
#   "gca_full_mod_base", "gca_full_mod_base1", "gca_full_mod_base2",
#   "gca_full_mod_base3", "gca_full_mod_base4", "gca_full_mod_group_0",
#   "gca_full_mod_group_1", "gca_full_mod_group_2", "gca_full_mod_group_3",
#   "gca_full_mod_int_0", "gca_full_mod_int_1", "gca_full_mod_int_2",
#   "gca_full_mod_int_3", "gca_full_mod_int_0_ams", "gca_full_mod_int_0_ims",
#   "gca_full_mod_int_0_aes", "gca_full_mod_int_0_ies"))
#   
#  
# 
# save(full_mods,
#      file = here("mods", "stress", "gca",
#                  "full_mods.Rdata"))
# 
# # final model
# save(gca_full_mod_int_0, file = here("mods", "stress", "gca", "final_model.Rdata"))

# Save anova model comparisons
nested_model_comparisons <-
  mget(c(#"mon_cond_anova", "mon_wm_anova", "mon_int_anova",
  #        "aes_cond_anova", "aes_wm_anova", "aes_int_anova",
  #        "ies_cond_anova", "ies_wm_anova", "ies_int_anova",
  #        "ams_cond_anova", "ams_wm_anova", "ams_int_anova",
  #        "ims_cond_anova", "ims_wm_anova", "ims_int_anova",
         "full_group_anova", "full_int_anova", #"full_wm_anova",
         'full_mod_group_relevel_anova'))

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

