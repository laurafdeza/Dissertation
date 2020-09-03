#
# Original script by Joseph
# Modified by Cristina to adapt to CrisLau Project
# Last updat: 06/12/2019
#
# GLMMs -----------------------------------------------------------------------
#
# - Target fixation as a function of group, stress (condition), and coda
#   at the offset of the first syllable (time_zero == 20)
# - This model builds on the t-test analyses by looking for between group
#   differences in target fixation at the offfset of first syllable
#
# -----------------------------------------------------------------------------



# Load data and models --------------------------------------------------------

source(here::here("scripts", "02_load_data.R"))

prop_0_mod_0     <- readRDS(here("mods", "glmm", "0_prop_0_mod_0.rds"))
prop_0_mod_group <- readRDS(here("mods", "glmm", "1_prop_0_mod_group.rds"))
prop_0_mod_cond  <- readRDS(here("mods", "glmm", "2_prop_0_mod_cond.rds"))
prop_0_mod_coda  <- readRDS(here("mods", "glmm", "3_prop_0_mod_coda.rds"))
prop_0_mod_int1  <- readRDS(here("mods", "glmm", "4_prop_0_mod_int1.rds"))
prop_0_mod_int2  <- readRDS(here("mods", "glmm", "5_prop_0_mod_int2.rds"))
prop_0_mod_int3  <- readRDS(here("mods", "glmm", "6_prop_0_mod_int3.rds"))
prop_0_mod_full  <- readRDS(here("mods", "glmm", "7_prop_0_mod_full.rds"))
prop_0_mod_final <- readRDS(here("mods", "glmm", "8_prop_0_mod_final.rds"))
# -----------------------------------------------------------------------------



# Data prep -------------------------------------------------------------------

# Filter time course to offset of 1st syllable (time_zero == 20, to account for 200 ms to launch saccade)
# Create sum coded fixed factors (condition)


df_stress <- stress10 %>%
  filter(., time_zero == 20) %>%
  mutate(., condition_sum = if_else(cond == "1", 1, -1)) #   1 = present

# -----------------------------------------------------------------------------




# Random effects building -----------------------------------------------------

if(F) {

prop_0_ranefA <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 | participant),
                    data = df_stress, family = 'binomial',                 
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_ranefB <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefA, prop_0_ranefB, refit = F) # keep intercept for target bc significant
#               Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# prop_0_ranefA  2 32959 32971 -16478    32955                             
# prop_0_ranefB  3 32531 32549 -16263    32525 430.36      1  < 2.2e-16 ***

prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 + condition_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep slope for condition bc significant
#               Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# prop_0_ranefB  3 32531 32549 -16263    32525                             
# prop_0_ranefC  5 31381 31411 -15686    31371 1153.7      2  < 2.2e-16 ***

 }

# -----------------------------------------------------------------------------










# Test fixed effects ----------------------------------------------------------

if(F) {
prop_0_mod_0 <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 + condition_sum | participant) +   
                    (1 | target),
                    data = df_stress, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- update(prop_0_mod_0,     . ~ . + group)
prop_0_mod_cond  <- update(prop_0_mod_group, . ~ . + condition_sum)
prop_0_mod_wm    <- update(prop_0_mod_cond,  . ~ . + WM_set)
prop_0_mod_int1  <- update(prop_0_mod_wm,    . ~ . + group:condition_sum)
prop_0_mod_int2  <- update(prop_0_mod_int1,  . ~ . + group:WM_set)

anova(prop_0_mod_0, prop_0_mod_group, test = "Chisq")    # main effect of group bc significant
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# prop_0_mod_0      5 31381 31411 -15686    31371                             
# prop_0_mod_group  9 31338 31390 -15660    31320 51.885      4  1.458e-10 ***

anova(prop_0_mod_group, prop_0_mod_cond, test = "Chisq") # no effect of condition bc not significant
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# prop_0_mod_group  9 31338 31390 -15660    31320                           
# prop_0_mod_cond  10 31336 31395 -15658    31316 3.1279      1    0.07696 .

anova(prop_0_mod_group, prop_0_mod_int1, test = "Chisq") # no interaction
#                  Df   AIC   BIC logLik deviance Chisq Chi Df Pr(>Chisq)
# prop_0_mod_group  9 31338 31390 -15660    31320                        
# prop_0_mod_int1  15 31340 31428 -15655    31310 9.211      6     0.1621

anova(prop_0_mod_group, prop_0_mod_int2, test = "Chisq") # no interaction
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# prop_0_mod_group  9 31338 31390 -15660    31320                         
# prop_0_mod_int2  19 31343 31454 -15653    31305 14.264     10     0.1613

prop_0_mod_int3  <- update(prop_0_mod_group,  . ~ . + condition_sum:WM_set) # failed to converge if prop_0_mod_int2 updated instead of _group
prop_0_mod_int4  <- update(prop_0_mod_int3,  . ~ . + group:condition_sum:WM_set)

anova(prop_0_mod_group, prop_0_mod_int3, test = "Chisq")
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# prop_0_mod_group  9 31338 31390 -15660    31320                           
# prop_0_mod_int3  10 31335 31393 -15657    31315 4.7093      1       0.03 *

anova(prop_0_mod_int3, prop_0_mod_int4, test = "Chisq")
#                 Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# prop_0_mod_int3 10 31335 31393 -15657    31315                         
# prop_0_mod_int4 14 31335 31417 -15654    31307 7.7437      4     0.1014

summary(prop_0_mod_int3) # Adv EN reference
# Fixed effects:
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.294472   0.125432  -2.348   0.0189 *  
# groupams              0.285887   0.155665   1.837   0.0663 .  
# groupies              0.201884   0.156261   1.292   0.1964    
# groupims              0.002413   0.156981   0.015   0.9877    
# groupmon              1.076471   0.159714   6.740 1.58e-11 ***
# condition_sum:WM_set -0.016300   0.007513  -2.170   0.0300 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

df_stress$group <- factor(df_stress$group, levels = c("mon", "aes", "ies", "ams", "ims"))

prop_0_mod_int3_mon <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        
                    group + condition_sum:WM_set +
                    (1 + condition_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

summary(prop_0_mod_int3_mon) # mon reference
#                         Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)           0.781990   0.128431   6.089 1.14e-09 ***
#   groupaes             -1.076456   0.159692  -6.741 1.57e-11 ***
#   groupies             -0.874581   0.157428  -5.555 2.77e-08 ***
#   groupams             -0.790570   0.159168  -4.967 6.80e-07 ***
#   groupims             -1.074047   0.158598  -6.772 1.27e-11 ***
#   condition_sum:WM_set -0.016300   0.007513  -2.170     0.03 *  

df_stress$group <- factor(df_stress$group, levels = c("ies", "mon", "aes", "ams", "ims"))

prop_0_mod_int3_ies <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        
                               group + condition_sum:WM_set +
                               (1 + condition_sum | participant) +
                               (1 | target),
                             data = df_stress, family = 'binomial',
                             control = glmerControl(optimizer = 'bobyqa'))

summary(prop_0_mod_int3_ies) # Int EN reference
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.092581   0.123522  -0.750    0.454    
# groupmon              0.874574   0.157447   5.555 2.78e-08 ***
# groupaes             -0.201886   0.156250  -1.292    0.196    
# groupams              0.084001   0.155737   0.539    0.590    
# groupims             -0.199471   0.154232  -1.293    0.196    
# condition_sum:WM_set -0.016300   0.007513  -2.170    0.030 * 

df_stress$group <- factor(df_stress$group, levels = c("ams", "ims", "ies", "aes", "mon"))

prop_0_mod_int3_ams <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        
                               group + condition_sum:WM_set +
                               (1 + condition_sum | participant) +
                               (1 | target),
                             data = df_stress, family = 'binomial',
                             control = glmerControl(optimizer = 'bobyqa'))

summary(prop_0_mod_int3_ams) # Adv MA reference
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.008583   0.124845  -0.069   0.9452    
# groupims             -0.283471   0.156502  -1.811   0.0701 .  
# groupies             -0.084008   0.155739  -0.539   0.5896    
# groupaes             -0.285877   0.155661  -1.837   0.0663 .  
# groupmon              0.790581   0.159173   4.967 6.81e-07 ***
# condition_sum:WM_set -0.016300   0.007513  -2.170   0.0300 *  

df_stress$group <- factor(df_stress$group, levels = c("ims", "ams", "ies", "aes", "mon"))

prop_0_mod_int3_ims <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        
                               group + condition_sum:WM_set +
                               (1 + condition_sum | participant) +
                               (1 | target),
                             data = df_stress, family = 'binomial',
                             control = glmerControl(optimizer = 'bobyqa'))

summary(prop_0_mod_int3_ims) # Int MA reference
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.292055   0.124864  -2.339   0.0193 *  
# groupams              0.283473   0.156537   1.811   0.0702 .  
# groupies              0.199465   0.154293   1.293   0.1961    
# groupaes             -0.002407   0.157023  -0.015   0.9878    
# groupmon              1.074046   0.158644   6.770 1.29e-11 ***
# condition_sum:WM_set -0.016300   0.007513  -2.170   0.0300 *  


prop_0_mod_final <- prop_0_mod_int3_mon


}

# -----------------------------------------------------------------------------








# Save models -----------------------------------------------------------------

if(F) {

saveRDS(prop_0_mod_0, here("mods", "stress", 
                            "glmm", "0_prop_0_mod_0.rds"))
saveRDS(prop_0_mod_group, here("mods", "stress",
                               "glmm", "1_prop_0_mod_group.rds"))
saveRDS(prop_0_mod_cond, here("mods", "stress",
                              "glmm", "2_prop_0_mod_cond.rds"))
saveRDS(prop_0_mod_wm, here("mods", "stress",
                            "glmm", "3_prop_0_mod_wm.rds"))
saveRDS(prop_0_mod_int1, here("mods", "stress",
                              "glmm", "4_prop_0_mod_int1.rds"))
saveRDS(prop_0_mod_int2, here("mods", "stress",
                              "glmm", "5_prop_0_mod_int2.rds"))
saveRDS(prop_0_mod_int3, here("mods", "stress", 
                              "glmm", "6_prop_0_mod_int3.rds"))
saveRDS(prop_0_mod_int4, here("mods", "stress", 
                              "glmm", "7_prop_0_mod_int4.rds"))
saveRDS(prop_0_mod_int3_mon, here("mods", "stress", 
                              "glmm", "8_prop_0_mod_int3_mon.rds"))
saveRDS(prop_0_mod_int3_ies, here("mods", "stress", 
                              "glmm", "9_prop_0_mod_int3_ies.rds"))
saveRDS(prop_0_mod_int3_ams, here("mods", "stress", 
                              "glmm", "10_prop_0_mod_int3_ams.rds"))
saveRDS(prop_0_mod_int3_ims, here("mods", "stress", 
                              "glmm", "11_prop_0_mod_int3_ims.rds"))
saveRDS(prop_0_mod_final, here("mods", "stress",
                               "glmm", "12_prop_0_mod_final.rds"))
}

# -----------------------------------------------------------------------------









# Model descriptives ----------------------------------------------------------
# 
MuMIn::r.squaredGLMM(prop_0_mod_final) 
#                   R2m       R2c
# theoretical 0.1345652 0.7399663
# delta       0.1274131 0.7006375

summary(prop_0_mod_final)
#                         Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)           0.781990   0.128431   6.089 1.14e-09 ***
#   groupaes             -1.076456   0.159692  -6.741 1.57e-11 ***
#   groupies             -0.874581   0.157428  -5.555 2.77e-08 ***
#   groupams             -0.790570   0.159168  -4.967 6.80e-07 ***
#   groupims             -1.074047   0.158598  -6.772 1.27e-11 ***
#   condition_sum:WM_set -0.016300   0.007513  -2.170     0.03 *  
  
confint(prop_0_mod_final, method = "Wald")
#                            2.5 %       97.5 %
# (Intercept)           0.53026893  1.033710373
# groupaes             -1.38944642 -0.763466265
# groupies             -1.18313461 -0.566026779
# groupams             -1.10253406 -0.478605652
# groupims             -1.38489304 -0.763200152
# condition_sum:WM_set -0.03102512 -0.001574841
