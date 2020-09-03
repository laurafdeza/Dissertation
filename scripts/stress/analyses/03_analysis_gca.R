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
  #         Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
  # mod_ot1  9 235270 235348 -117626   235252                             
  # mod_ot2 14 234469 234590 -117220   234441 810.75      5  < 2.2e-16 ***
  # mod_ot3 20 233688 233861 -116824   233648 792.81      6  < 2.2e-16 ***
  # mod_ot4 21 232813 232995 -116385   232771 877.46      1  < 2.2e-16 ***
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

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + condition_sum) # singular
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2,   . ~ . + ot3:condition_sum) # singular

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2, gca_mod_mon_cond_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base   30 41993 42203 -20967    41933                         
# gca_mod_mon_cond_0 31 41995 42212 -20966    41933 0.5922      1     0.4416
# gca_mod_mon_cond_1 32 41997 42221 -20966    41933 0.0193      1     0.8895
# gca_mod_mon_cond_2 33 41999 42230 -20966    41933 0.0047      1     0.9456
# gca_mod_mon_cond_3 34 42001 42238 -20966    41933 0.1931      1     0.6603

# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_wm_0 <- update(gca_mod_mon_cond_3, . ~ . + WM_set) # singular
gca_mod_mon_wm_1 <- update(gca_mod_mon_wm_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_mon_wm_2 <- update(gca_mod_mon_wm_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_mon_wm_3 <- update(gca_mod_mon_wm_2,   . ~ . + ot3:condition_sum) # singular

mon_wm_anova <-
  anova(gca_mod_mon_cond_3, gca_mod_mon_wm_0, gca_mod_mon_wm_1,
        gca_mod_mon_wm_2, gca_mod_mon_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_cond_3 34 42001 42238 -20966    41933                         
# gca_mod_mon_wm_0   35 42002 42247 -20966    41932 0.6328      1     0.4263
# gca_mod_mon_wm_1   35 42002 42247 -20966    41932 0.0000      0     1.0000
# gca_mod_mon_wm_2   35 42002 42247 -20966    41932 0.0000      0     1.0000
# gca_mod_mon_wm_3   35 42002 42247 -20966    41932 0.0000      0     1.0000

# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_int_0 <- update(gca_mod_mon_wm_3, . ~ . + condition_sum:WM_set) # singular
gca_mod_mon_int_1 <- update(gca_mod_mon_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_mon_int_2 <- update(gca_mod_mon_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_mon_int_3 <- update(gca_mod_mon_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

mon_int_anova <-
  anova(gca_mod_mon_wm_3, gca_mod_mon_int_0, gca_mod_mon_int_1,
        gca_mod_mon_int_2, gca_mod_mon_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_mon_wm_3  35 42002 42247 -20966    41932                           
# gca_mod_mon_int_0 36 41998 42250 -20963    41926 6.1973      1    0.01279 *
# gca_mod_mon_int_1 37 42000 42259 -20963    41926 0.0158      1    0.89982  
# gca_mod_mon_int_2 38 42000 42266 -20962    41924 2.0180      1    0.15545  
# gca_mod_mon_int_3 39 42001 42274 -20962    41923 0.2935      1    0.58797  






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
# gca_mod_aes_base   30 44963 45175 -22451    44903                           
# gca_mod_aes_cond_0 31 44964 45183 -22451    44902 0.7533      1    0.38545  
# gca_mod_aes_cond_1 32 44965 45191 -22450    44901 1.3402      1    0.24700  
# gca_mod_aes_cond_2 33 44967 45200 -22450    44901 0.0784      1    0.77948  
# gca_mod_aes_cond_3 34 44965 45205 -22448    44897 3.9651      1    0.04645 *

# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_wm_0 <- update(gca_mod_aes_cond_3, . ~ . + WM_set) # singular
gca_mod_aes_wm_1 <- update(gca_mod_aes_wm_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_aes_wm_2 <- update(gca_mod_aes_wm_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_aes_wm_3 <- update(gca_mod_aes_wm_2,   . ~ . + ot3:condition_sum) # singular

aes_wm_anova <-
  anova(gca_mod_aes_cond_3, gca_mod_aes_wm_0, gca_mod_aes_wm_1,
        gca_mod_aes_wm_2, gca_mod_aes_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_aes_cond_3 34 44965 45205 -22448    44897                         
# gca_mod_aes_wm_0   35 44966 45214 -22448    44896 0.0825      1      0.774
# gca_mod_aes_wm_1   35 44966 45214 -22448    44896 0.0000      0      1.000
# gca_mod_aes_wm_2   35 44966 45214 -22448    44896 0.0000      0      1.000
# gca_mod_aes_wm_3   35 44966 45214 -22448    44896 0.0000      0      1.000

# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_int_0 <- update(gca_mod_aes_wm_3, . ~ . + condition_sum:WM_set) # singular
gca_mod_aes_int_1 <- update(gca_mod_aes_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_aes_int_2 <- update(gca_mod_aes_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_aes_int_3 <- update(gca_mod_aes_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

aes_int_anova <-
  anova(gca_mod_aes_wm_3, gca_mod_aes_int_0, gca_mod_aes_int_1,
        gca_mod_aes_int_2, gca_mod_aes_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_aes_wm_3  35 44966 45214 -22448    44896                           
# gca_mod_aes_int_0 36 44966 45221 -22447    44894 2.0135      1    0.15590  
# gca_mod_aes_int_1 37 44966 45227 -22446    44892 2.7414      1    0.09778 .
# gca_mod_aes_int_2 38 44968 45236 -22446    44892 0.0001      1    0.99313  
# gca_mod_aes_int_3 39 44969 45245 -22446    44891 0.6389      1    0.42411 




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
# gca_mod_ies_base   30 47472 47685 -23706    47412                         
# gca_mod_ies_cond_0 31 47474 47694 -23706    47412 0.1599      1     0.6892
# gca_mod_ies_cond_1 32 47476 47703 -23706    47412 0.0297      1     0.8632
# gca_mod_ies_cond_2 33 47476 47710 -23705    47410 2.2726      1     0.1317
# gca_mod_ies_cond_3 34 47478 47719 -23705    47410 0.0117      1     0.9140


# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_wm_0 <- update(gca_mod_ies_cond_3, . ~ . + WM_set) # singular
gca_mod_ies_wm_1 <- update(gca_mod_ies_wm_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ies_wm_2 <- update(gca_mod_ies_wm_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ies_wm_3 <- update(gca_mod_ies_wm_2,   . ~ . + ot3:condition_sum) # singular

ies_wm_anova <-
  anova(gca_mod_ies_cond_3, gca_mod_ies_wm_0, gca_mod_ies_wm_1,
        gca_mod_ies_wm_2, gca_mod_ies_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_cond_3 34 47478 47719 -23705    47410                        
# gca_mod_ies_wm_0   35 47479 47727 -23705    47409 0.586      1      0.444
# gca_mod_ies_wm_1   35 47479 47727 -23705    47409 0.000      0      1.000
# gca_mod_ies_wm_2   35 47479 47727 -23705    47409 0.000      0      1.000
# gca_mod_ies_wm_3   35 47479 47727 -23705    47409 0.000      0      1.000

# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_int_0 <- update(gca_mod_ies_wm_3, . ~ . + condition_sum:WM_set) # singular
gca_mod_ies_int_1 <- update(gca_mod_ies_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_ies_int_2 <- update(gca_mod_ies_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_ies_int_3 <- update(gca_mod_ies_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

ies_int_anova <-
  anova(gca_mod_ies_wm_3, gca_mod_ies_int_0, gca_mod_ies_int_1,
        gca_mod_ies_int_2, gca_mod_ies_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ies_wm_3  35 47479 47727 -23705    47409                              
# gca_mod_ies_int_0 36 47478 47733 -23703    47406  3.2702      1    0.07055 .  
# gca_mod_ies_int_1 37 47478 47740 -23702    47404  2.4199      1    0.11981    
# gca_mod_ies_int_2 38 47463 47732 -23693    47387 17.0424      1  3.655e-05 ***
# gca_mod_ies_int_3 39 47458 47735 -23690    47380  6.2786      1    0.01222 *  




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
# gca_mod_ams_base   30 46196 46408 -23068    46136                         
# gca_mod_ams_cond_0 31 46198 46416 -23068    46136 0.3503      1     0.5539
# gca_mod_ams_cond_1 32 46198 46424 -23067    46134 1.4069      1     0.2356
# gca_mod_ams_cond_2 33 46200 46433 -23067    46134 0.0002      1     0.9900
# gca_mod_ams_cond_3 34 46200 46440 -23066    46132 1.8740      1     0.1710


# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_wm_0 <- update(gca_mod_ams_cond_3, . ~ . + WM_set) # singular
gca_mod_ams_wm_1 <- update(gca_mod_ams_wm_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ams_wm_2 <- update(gca_mod_ams_wm_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ams_wm_3 <- update(gca_mod_ams_wm_2,   . ~ . + ot3:condition_sum) # singular

ams_wm_anova <-
  anova(gca_mod_ams_cond_3, gca_mod_ams_wm_0, gca_mod_ams_wm_1,
        gca_mod_ams_wm_2, gca_mod_ams_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_cond_3 34 46200 46440 -23066    46132                        
# gca_mod_ams_wm_0   35 46202 46449 -23066    46132 0.645      1     0.4219
# gca_mod_ams_wm_1   35 46202 46449 -23066    46132 0.000      0     1.0000
# gca_mod_ams_wm_2   35 46202 46449 -23066    46132 0.000      0     1.0000
# gca_mod_ams_wm_3   35 46202 46449 -23066    46132 0.000      0     1.0000


# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_int_0 <- update(gca_mod_ams_wm_3, . ~ . + condition_sum:WM_set) # singular
gca_mod_ams_int_1 <- update(gca_mod_ams_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_ams_int_2 <- update(gca_mod_ams_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_ams_int_3 <- update(gca_mod_ams_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

ams_int_anova <-
  anova(gca_mod_ams_wm_3, gca_mod_ams_int_0, gca_mod_ams_int_1,
        gca_mod_ams_int_2, gca_mod_ams_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ams_wm_3  35 46202 46449 -23066    46132                            
# gca_mod_ams_int_0 36 46203 46457 -23065    46131 0.8302      1   0.362209   
# gca_mod_ams_int_1 37 46204 46465 -23065    46130 1.0670      1   0.301626   
# gca_mod_ams_int_2 38 46198 46467 -23061    46122 7.5275      1   0.006076 **
# gca_mod_ams_int_3 39 46200 46476 -23061    46122 0.1096      1   0.740616   







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
# gca_mod_ims_base   30 47123 47334 -23531    47063                         
# gca_mod_ims_cond_0 31 47124 47343 -23531    47062 0.8159      1     0.3664
# gca_mod_ims_cond_1 32 47126 47352 -23531    47062 0.0511      1     0.8212
# gca_mod_ims_cond_2 33 47127 47360 -23530    47061 1.1927      1     0.2748
# gca_mod_ims_cond_3 34 47129 47369 -23530    47061 0.0045      1     0.9464


# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_wm_0 <- update(gca_mod_ims_cond_3, . ~ . + WM_set) # singular
gca_mod_ims_wm_1 <- update(gca_mod_ims_wm_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ims_wm_2 <- update(gca_mod_ims_wm_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ims_wm_3 <- update(gca_mod_ims_wm_2,   . ~ . + ot3:condition_sum) # singular

ims_wm_anova <-
  anova(gca_mod_ims_cond_3, gca_mod_ims_wm_0, gca_mod_ims_wm_1,
        gca_mod_ims_wm_2, gca_mod_ims_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_cond_3 34 47129 47369 -23530    47061                         
# gca_mod_ims_wm_0   35 47129 47377 -23530    47059 1.0615      1     0.3029
# gca_mod_ims_wm_1   35 47129 47377 -23530    47059 0.0000      0     1.0000
# gca_mod_ims_wm_2   35 47129 47377 -23530    47059 0.0000      0     1.0000
# gca_mod_ims_wm_3   35 47129 47377 -23530    47059 0.0000      0     1.0000


# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_int_0 <- update(gca_mod_ims_wm_3, . ~ . + condition_sum:WM_set) # singular
gca_mod_ims_int_1 <- update(gca_mod_ims_int_0,   . ~ . + ot1:condition_sum:WM_set) # singular
gca_mod_ims_int_2 <- update(gca_mod_ims_int_1,   . ~ . + ot2:condition_sum:WM_set) # singular
gca_mod_ims_int_3 <- update(gca_mod_ims_int_2,   . ~ . + ot3:condition_sum:WM_set) # singular

ims_int_anova <-
  anova(gca_mod_ims_wm_3, gca_mod_ims_int_0, gca_mod_ims_int_1,
        gca_mod_ims_int_2, gca_mod_ims_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ims_wm_3  35 47129 47377 -23530    47059                         
# gca_mod_ims_int_0 36 47131 47385 -23530    47059 0.3421      1     0.5586
# gca_mod_ims_int_1 37 47133 47394 -23530    47059 0.0193      1     0.8897
# gca_mod_ims_int_2 38 47133 47401 -23529    47057 1.9838      1     0.1590
# gca_mod_ims_int_3 39 47134 47410 -23528    47056 0.6766      1     0.4108



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
       data = stress_gc_subset, REML = F) # singular

gca_full_mod_base1 <- update(gca_full_mod_base,   . ~ . + condition_sum:WM_set) # singular
gca_full_mod_base2 <- update(gca_full_mod_base1,  . ~ . + ot3:condition_sum) # singular
gca_full_mod_base3 <- update(gca_full_mod_base2,  . ~ . + ot2:condition_sum:WM_set) # singular
gca_full_mod_base4 <- update(gca_full_mod_base3,  . ~ . + ot3:condition_sum:WM_set) # singular

anova(gca_full_mod_base, gca_full_mod_base1, gca_full_mod_base2,
      gca_full_mod_base3, gca_full_mod_base4)
#                    Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# gca_full_mod_base  30 228578 228838 -114259   228518                            
# gca_full_mod_base1 31 228578 228846 -114258   228516 2.0142      1   0.155837   
# gca_full_mod_base2 32 228577 228854 -114257   228513 2.5509      1   0.110233   
# gca_full_mod_base3 33 228579 228865 -114257   228513 0.0832      1   0.772959   
# gca_full_mod_base4 34 228572 228866 -114252   228504 9.2979      1   0.002294 **


# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_group_0 <- update(gca_full_mod_base4,    . ~ . + group) # singular
gca_full_mod_group_1 <- update(gca_full_mod_group_0, . ~ . + ot1:group) # singular
gca_full_mod_group_2 <- update(gca_full_mod_group_1, . ~ . + ot2:group) # singular
gca_full_mod_group_3 <- update(gca_full_mod_group_2, . ~ . + ot3:group) # singular

full_group_anova <-
  anova(gca_full_mod_base4, gca_full_mod_group_0, gca_full_mod_group_1,
        gca_full_mod_group_2, gca_full_mod_group_3)
#                      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_full_mod_base4   34 228572 228866 -114252   228504                              
# gca_full_mod_group_0 38 228552 228881 -114238   228476 28.1017      4  1.189e-05 ***
# gca_full_mod_group_1 42 228552 228916 -114234   228468  7.3807      4     0.1171    
# gca_full_mod_group_2 46 228554 228953 -114231   228462  5.9537      4     0.2026    
# gca_full_mod_group_3 50 228515 228949 -114208   228415 46.8047      4  1.675e-09 ***


################################

# add 3-way int to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_int_0 <- update(gca_full_mod_group_3, . ~ . + WM_set:condition_sum:group) # singular        
gca_full_mod_int_1 <- update(gca_full_mod_int_0,   . ~ . + ot1:WM_set:condition_sum:group) # singular
gca_full_mod_int_2 <- update(gca_full_mod_int_1,   . ~ . + ot2:WM_set:condition_sum:group) # singular
gca_full_mod_int_3 <- update(gca_full_mod_int_2,   . ~ . + ot3:WM_set:condition_sum:group) # singular

full_int_anova <-
  anova(gca_full_mod_group_3, gca_full_mod_int_0, gca_full_mod_int_1,
        gca_full_mod_int_2, gca_full_mod_int_3)
#                      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_full_mod_group_3 50 228515 228949 -114208   228415                             
# gca_full_mod_int_0   54 228521 228988 -114206   228413  2.9440      4   0.567247   
# gca_full_mod_int_1   59 228514 229025 -114198   228396 16.8671      5   0.004759 **
# gca_full_mod_int_2   63 228516 229061 -114195   228390  5.9976      4   0.199325   
# gca_full_mod_int_3   67 228518 229098 -114192   228384  5.9228      4   0.204985   

# ---

summary(gca_full_mod_int_0) # mon reference   -    FINAL MODEL
# Fixed effects:
#                                 Estimate Std. Error t value
# (Intercept)                    2.5907306  0.1591699  16.277
# ot1                            4.2883790  0.4918307   8.719
# ot2                           -2.7427766  0.4120458  -6.656
# ot3                            0.2595204  0.2658637   0.976
# groupaes                      -0.4734012  0.1938739  -2.442
# groupams                      -0.6687256  0.1938830  -3.449
# groupies                      -0.7231428  0.1924898  -3.757
# groupims                      -1.1617195  0.1938813  -5.992
# condition_sum:WM_set          -0.0007732  0.0141744  -0.055
# ot3:condition_sum              0.4956811  0.2555710   1.940
# ot1:groupaes                   2.7396911  0.5450331   5.027
# ot1:groupams                   1.2414774  0.5451114   2.277
# ot1:groupies                   1.6122155  0.5412879   2.978
# ot1:groupims                   1.6363756  0.5451111   3.002
# ot2:groupaes                   0.1686496  0.4063359   0.415
# ot2:groupams                   0.8416685  0.4064497   2.071
# ot2:groupies                   0.7713272  0.4037018   1.911
# ot2:groupims                   0.9564437  0.4064500   2.353
# ot3:groupaes                  -1.8122507  0.2987619  -6.066
# ot3:groupams                  -1.6702070  0.2989170  -5.588
# ot3:groupies                  -1.8091552  0.2970251  -6.091
# ot3:groupims                  -1.5288318  0.2989155  -5.115
# ot2:condition_sum:WM_set       0.0055887  0.0194531   0.287
# ot3:condition_sum:WM_set      -0.0861903  0.0282826  -3.047
# groupaes:condition_sum:WM_set -0.0224302  0.0150556  -1.490
# groupams:condition_sum:WM_set -0.0136383  0.0157737  -0.865
# groupies:condition_sum:WM_set -0.0060240  0.0154352  -0.390
# groupims:condition_sum:WM_set -0.0089286  0.0160987  -0.555

#save(gca_full_mod_int_0, file = here("mods", "stress", "gca", "full_model_mon.Rdata"))  


# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ams"))
gca_full_mod_int_0_ams <- update(gca_full_mod_int_0)

summary(gca_full_mod_int_0_ams)
# Fixed effects:
#                                Estimate Std. Error t value
# (Intercept)                    1.922006   0.155354  12.372
# ot1                            5.529861   0.482085  11.471
# ot2                           -1.901104   0.405616  -4.687
# ot3                           -1.410683   0.260507  -5.415
# groupmon                       0.668726   0.193883   3.449
# groupaes                       0.195324   0.190749   1.024
# groupies                      -0.054417   0.189343  -0.287
# groupims                      -0.492994   0.190757  -2.584
# condition_sum:WM_set          -0.014411   0.011423  -1.262
# ot3:condition_sum              0.495683   0.255571   1.940
# ot1:groupmon                  -1.241482   0.545112  -2.277
# ot1:groupaes                   1.498213   0.536251   2.794
# ot1:groupies                   0.370738   0.532441   0.696
# ot1:groupims                   0.394896   0.536331   0.736
# ot2:groupmon                  -0.841669   0.406449  -2.071
# ot2:groupaes                  -0.673019   0.399804  -1.683
# ot2:groupies                  -0.070341   0.397123  -0.177
# ot2:groupims                   0.114775   0.399919   0.287
# ot3:groupmon                   1.670209   0.298917   5.588
# ot3:groupaes                  -0.142043   0.293997  -0.483
# ot3:groupies                  -0.138948   0.292227  -0.475
# ot3:groupims                   0.141376   0.294151   0.481
# ot2:condition_sum:WM_set       0.005589   0.019453   0.287
# ot3:condition_sum:WM_set      -0.086190   0.028283  -3.047
# groupmon:condition_sum:WM_set  0.013638   0.015774   0.865
# groupaes:condition_sum:WM_set -0.008792   0.012915  -0.681
# groupies:condition_sum:WM_set  0.007614   0.013393   0.569
# groupims:condition_sum:WM_set  0.004710   0.014193   0.332


stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ims"))
gca_full_mod_int_0_ims <- update(gca_full_mod_int_0)

summary(gca_full_mod_int_0_ims)
# Fixed effects:
#                                Estimate Std. Error t value
# (Intercept)                    1.429012   0.155351   9.199
# ot1                            5.924759   0.482084  12.290
# ot2                           -1.786328   0.405613  -4.404
# ot3                           -1.269307   0.260504  -4.872
# groupams                       0.492994   0.190758   2.584
# groupmon                       1.161720   0.193881   5.992
# groupaes                       0.688318   0.190748   3.609
# groupies                       0.438577   0.189342   2.316
# condition_sum:WM_set          -0.009702   0.012056  -0.805
# ot3:condition_sum              0.495683   0.255571   1.940
# ot1:groupams                  -0.394897   0.536330  -0.736
# ot1:groupmon                  -1.636377   0.545111  -3.002
# ot1:groupaes                   1.103316   0.536251   2.057
# ot1:groupies                  -0.024159   0.532440  -0.045
# ot2:groupams                  -0.114776   0.399920  -0.287
# ot2:groupmon                  -0.956443   0.406450  -2.353
# ot2:groupaes                  -0.787794   0.399805  -1.970
# ot2:groupies                  -0.185117   0.397123  -0.466
# ot3:groupams                  -0.141375   0.294151  -0.481
# ot3:groupmon                   1.528832   0.298915   5.115
# ot3:groupaes                  -0.283419   0.293995  -0.964
# ot3:groupies                  -0.280323   0.292225  -0.959
# ot2:condition_sum:WM_set       0.005589   0.019453   0.287
# ot3:condition_sum:WM_set      -0.086190   0.028283  -3.047
# groupams:condition_sum:WM_set -0.004710   0.014193  -0.332
# groupmon:condition_sum:WM_set  0.008929   0.016099   0.555
# groupaes:condition_sum:WM_set -0.013502   0.013354  -1.011
# groupies:condition_sum:WM_set  0.002905   0.013806   0.210

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "aes"))
gca_full_mod_int_0_aes <- update(gca_full_mod_int_0)

summary(gca_full_mod_int_0_aes)
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                    2.117330   0.155342  13.630
# ot1                            7.028073   0.481998  14.581
# ot2                           -2.574124   0.405503  -6.348
# ot3                           -1.552728   0.260329  -5.964
# groupims                      -0.688318   0.190748  -3.609
# groupams                      -0.195324   0.190749  -1.024
# groupmon                       0.473401   0.193874   2.442
# groupies                      -0.249742   0.189334  -1.319
# condition_sum:WM_set          -0.023203   0.009913  -2.341
# ot3:condition_sum              0.495683   0.255571   1.940
# ot1:groupims                  -1.103316   0.536251  -2.057
# ot1:groupams                  -1.498213   0.536250  -2.794
# ot1:groupmon                  -2.739694   0.545032  -5.027
# ot1:groupies                  -1.127475   0.532361  -2.118
# ot2:groupims                   0.787794   0.399804   1.970
# ot2:groupams                   0.673019   0.399804   1.683
# ot2:groupmon                  -0.168651   0.406335  -0.415
# ot2:groupies                   0.602678   0.397008   1.518
# ot3:groupims                   0.283419   0.293995   0.964
# ot3:groupams                   0.142043   0.293996   0.483
# ot3:groupmon                   1.812252   0.298762   6.066
# ot3:groupies                   0.003095   0.292071   0.011
# ot2:condition_sum:WM_set       0.005589   0.019453   0.287
# ot3:condition_sum:WM_set      -0.086190   0.028283  -3.047
# groupims:condition_sum:WM_set  0.013502   0.013354   1.011
# groupams:condition_sum:WM_set  0.008792   0.012915   0.681
# groupmon:condition_sum:WM_set  0.022430   0.015056   1.490
# groupies:condition_sum:WM_set  0.016406   0.012468   1.316

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ies"))
gca_full_mod_int_0_ies <- update(gca_full_mod_int_0)

summary(gca_full_mod_int_0_ies)
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                    1.867588   0.153617  12.157
# ot1                            5.900599   0.477769  12.350
# ot2                           -1.971447   0.402866  -4.894
# ot3                           -1.549633   0.258332  -5.999
# groupaes                       0.249742   0.189334   1.319
# groupims                      -0.438577   0.189342  -2.316
# groupams                       0.054417   0.189344   0.287
# groupmon                       0.723143   0.192490   3.757
# condition_sum:WM_set          -0.006797   0.010800  -0.629
# ot3:condition_sum              0.495682   0.255571   1.940
# ot1:groupaes                   1.127474   0.532362   2.118
# ot1:groupims                   0.024157   0.532440   0.045
# ot1:groupams                  -0.370738   0.532441  -0.696
# ot1:groupmon                  -1.612222   0.541288  -2.978
# ot2:groupaes                  -0.602678   0.397009  -1.518
# ot2:groupims                   0.185116   0.397123   0.466
# ot2:groupams                   0.070341   0.397123   0.177
# ot2:groupmon                  -0.771328   0.403702  -1.911
# ot3:groupaes                  -0.003095   0.292071  -0.011
# ot3:groupims                   0.280324   0.292225   0.959
# ot3:groupams                   0.138948   0.292227   0.475
# ot3:groupmon                   1.809158   0.297025   6.091
# ot2:condition_sum:WM_set       0.005589   0.019453   0.287
# ot3:condition_sum:WM_set      -0.086190   0.028283  -3.047
# groupaes:condition_sum:WM_set -0.016406   0.012468  -1.316
# groupims:condition_sum:WM_set -0.002905   0.013806  -0.210
# groupams:condition_sum:WM_set -0.007614   0.013393  -0.569
# groupmon:condition_sum:WM_set  0.006024   0.015435   0.390

}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, condition_sum, WM_set) %>%
  distinct

# Get model predictions and SE
fits_all <- predictSE(gca_full_mod_int_0, new_dat_all) %>%  
  as_tibble %>%
  bind_cols(new_dat_all) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, SS = "mon", AE = "aes", IE = "ies", AM = "ams", IM = "ims"))

# Filter preds at target offset
target_offset_preds <- filter(fits_all, time_zero == 4) %>%
  select(group, cond = condition_sum, WM_set,
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
              "_cond_1", "_cond_2", "_cond_3", 
              "_wm_0", "_wm_1", "_wm_2", "_wm_3",
              "_int_0", "_int_1", "_int_2",
              "_int_3")

# Store ind models in list
ind_mods <- mget(c(#paste0(mod_type, "mon", mod_spec),
                   #paste0(mod_type, "aes", mod_spec),
                   #paste0(mod_type, "ies", mod_spec),
                   #paste0(mod_type, "ams", mod_spec),
                   paste0(mod_type, "ims", mod_spec)
                   ))

save(ind_mods,
     file = here("mods", "stress", "gca",
                 "ind_mods.Rdata"))

# Store full (ot1, ot2, ot3, group, coda, cond) models in list
full_mods <- mget(c(
  "gca_full_mod_base", "gca_full_mod_base1", "gca_full_mod_base2",
  "gca_full_mod_base3", "gca_full_mod_base4", "gca_full_mod_group_0",
  "gca_full_mod_group_1", "gca_full_mod_group_2", "gca_full_mod_group_3",
  "gca_full_mod_int_0", "gca_full_mod_int_1", "gca_full_mod_int_2",
  "gca_full_mod_int_3", "gca_full_mod_int_0_ams", "gca_full_mod_int_0_ims",
  "gca_full_mod_int_0_aes", "gca_full_mod_int_0_ies"))
  
 

save(full_mods,
     file = here("mods", "stress", "gca",
                 "full_mods.Rdata"))

# final model
save(gca_full_mod_int_0, file = here("mods", "stress", "gca", "final_model.Rdata"))

# Save anova model comparisons
nested_model_comparisons <-
  mget(c(#"mon_cond_anova", "mon_wm_anova", "mon_int_anova",
         #"aes_cond_anova", "aes_wm_anova", "aes_int_anova",
         #"ies_cond_anova", "ies_wm_anova", "ies_int_anova",
         #"ams_cond_anova", "ams_wm_anova", "ams_int_anova",
         "ims_cond_anova", "ims_wm_anova", "ims_int_anova",
         "full_group_anova", "full_int_anova"))

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

