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

source(here::here("scripts", "00_load_libs.R"))
source(here::here("scripts", "02_load_data.R"))

prop_0_mod_0     <- readRDS(here("mods", 'stress', "glmm", "0_prop_0_mod_0.rds"))
prop_0_mod_group <- readRDS(here("mods", 'stress', "glmm", "1_prop_0_mod_group.rds")) # aes as reference
prop_0_mod_cond  <- readRDS(here("mods", 'stress', "glmm", "2_prop_0_mod_cond.rds"))
prop_0_mod_wm    <- readRDS(here("mods", 'stress', "glmm", "3_prop_0_mod_wm.rds"))
prop_0_mod_int1  <- readRDS(here("mods", 'stress', "glmm", "4_prop_0_mod_int1.rds"))
prop_0_mod_int2  <- readRDS(here("mods", 'stress', "glmm", "5_prop_0_mod_int2.rds"))
prop_0_mod_int3  <- readRDS(here("mods", 'stress', "glmm", "6_prop_0_mod_int3.rds"))
prop_0_mod_int4  <- readRDS(here("mods", 'stress', "glmm", "7_prop_0_mod_int4.rds"))
prop_0_mod_group_mon <- readRDS(here("mods", 'stress', "glmm", "8_prop_0_mod_group_mon.rds"))
prop_0_mod_group_ies <- readRDS(here("mods", 'stress', "glmm", "9_prop_0_mod_group_ies.rds"))
prop_0_mod_group_ams <- readRDS(here("mods", 'stress', "glmm", "10_prop_0_mod_group_ams.rds"))
prop_0_mod_group_ims <- readRDS(here("mods", 'stress', "glmm", "11_prop_0_mod_group_ims.rds"))
prop_0_mod_final <- readRDS(here("mods", 'stress', "glmm", "12_prop_0_mod_final.rds"))

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
# prop_0_ranefA  2 27835 27847 -13916    27831                             
# prop_0_ranefB  3 27294 27312 -13644    27288 543.32      1  < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 + condition_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep slope for condition bc significant
#               Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# prop_0_ranefB  3 27294 27312 -13644    27288                             
# prop_0_ranefC  5 26240 26269 -13115    26230 1057.7      2  < 2.2e-16 ***

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

anova(prop_0_mod_0, prop_0_mod_group, test = "Chisq")    # main effect of group bc significant
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# prop_0_mod_0      5 26240 26269 -13115    26230                             
# prop_0_mod_group  9 26192 26245 -13087    26174 56.088      4  1.922e-11 ***

anova(prop_0_mod_group, prop_0_mod_cond, test = "Chisq") # no effect of condition bc not significant
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# prop_0_mod_group  9 26192 26245 -13087    26174                           
# prop_0_mod_cond  10 26191 26250 -13086    26171 2.9634      1    0.08517 .

anova(prop_0_mod_group, prop_0_mod_wm, test = "Chisq") # no effect of verbal WM
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# prop_0_mod_group  9 26192 26245 -13087    26174                         
# prop_0_mod_wm    11 26193 26257 -13086    26171 2.9802      2     0.2253

anova(prop_0_mod_group, prop_0_mod_int1, test = "Chisq") # no interaction group x condition
#                  Df   AIC   BIC logLik deviance Chisq Chi Df Pr(>Chisq)
# prop_0_mod_group  9 26192 26245 -13087    26174                         
# prop_0_mod_int1  15 26199 26287 -13085    26169 4.7696      6     0.5737

prop_0_mod_int2  <- update(prop_0_mod_group,  . ~ . + group:WM_set) # start here with interactions because it was not converging if the int were just added to previous model
prop_0_mod_int3  <- update(prop_0_mod_group,  . ~ . + condition_sum:WM_set) 
prop_0_mod_int4  <- update(prop_0_mod_int3,  . ~ . + group:condition_sum:WM_set)

anova(prop_0_mod_group, prop_0_mod_int2, test = "Chisq") # no interaction group x wm
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# prop_0_mod_group  9 26192 26245 -13087    26174                         
# prop_0_mod_int2  14 26198 26280 -13085    26170 3.8718      5      0.568

anova(prop_0_mod_group, prop_0_mod_int3, test = "Chisq") # no interaction condition x wm
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# prop_0_mod_group  9 26192 26245 -13087    26174                           
# prop_0_mod_int3  10 26191 26249 -13085    26171 3.5546      1    0.05938 .

anova(prop_0_mod_group, prop_0_mod_int4, test = "Chisq") # no 3-way interaction
#                 Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# prop_0_mod_int3 9 26192 26245 -13087    26174                         
# prop_0_mod_int4  14 26197 26278 -13084    26169 5.4048      5     0.3685

summary(prop_0_mod_group) # Adv EN reference
# Fixed effects:
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           1.2368     0.1713   7.221 5.16e-13 ***
# groupams             -0.3046     0.2189  -1.391  0.16414    
# groupies             -0.4047     0.2172  -1.863  0.06249 .  
# groupims             -0.6760     0.2174  -3.109  0.00188 ** 
# groupmon              1.0022     0.2271   4.413 1.02e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

df_stress$group <- factor(df_stress$group, levels = c("mon", "aes", "ies", "ams", "ims"))

prop_0_mod_group_mon <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        
                    group + 
                    (1 + condition_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

summary(prop_0_mod_group_mon) # mon reference
#                         Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)           2.2390     0.1803  12.420  < 2e-16 ***
#   groupaes             -1.0022     0.2271  -4.413 1.02e-05 ***
#   groupies             -1.4069     0.2242  -6.275 3.49e-10 ***
#   groupams             -1.3068     0.2259  -5.786 7.21e-09 ***
#   groupims             -1.6783     0.2256  -7.438 1.02e-13 ***

df_stress$group <- factor(df_stress$group, levels = c("ies", "mon", "aes", "ams", "ims"))

prop_0_mod_group_ies <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        
                               group +
                               (1 + condition_sum | participant) +
                               (1 | target),
                             data = df_stress, family = 'binomial',
                             control = glmerControl(optimizer = 'bobyqa'))

summary(prop_0_mod_group_ies) # Int EN reference
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           0.8321     0.1670   4.983 6.27e-07 ***
# groupmon              1.4069     0.2241   6.277 3.45e-10 ***
# groupaes              0.4047     0.2172   1.863   0.0624 .  
# groupams              0.1001     0.2156   0.464   0.6425    
# groupims             -0.2714     0.2153  -1.261   0.2075    

df_stress$group <- factor(df_stress$group, levels = c("ams", "ims", "ies", "aes", "mon"))

prop_0_mod_group_ams <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        
                               group + 
                               (1 + condition_sum | participant) +
                               (1 | target),
                             data = df_stress, family = 'binomial',
                             control = glmerControl(optimizer = 'bobyqa'))

summary(prop_0_mod_group_ams) # Adv MA reference
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           0.9322     0.1692   5.511 3.57e-08 ***
# groupims             -0.3714     0.2170  -1.712    0.087 .  
# groupies             -0.1001     0.2157  -0.464    0.643    
# groupaes              0.3046     0.2190   1.391    0.164    
# groupmon              1.3068     0.2259   5.784 7.28e-09 ***

df_stress$group <- factor(df_stress$group, levels = c("ims", "ams", "ies", "aes", "mon"))

prop_0_mod_group_ims <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        
                               group +
                               (1 + condition_sum | participant) +
                               (1 | target),
                             data = df_stress, family = 'binomial',
                             control = glmerControl(optimizer = 'bobyqa'))

summary(prop_0_mod_group_ims) # Int MA reference
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           0.5607     0.1687   3.324 0.000887 ***
# groupams              0.3714     0.2169   1.712 0.086870 .  
# groupies              0.2714     0.2153   1.260 0.207549    
# groupaes              0.6760     0.2174   3.110 0.001873 ** 
# groupmon              1.6783     0.2256   7.438 1.02e-13 ***


prop_0_mod_final <- prop_0_mod_group_mon


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
saveRDS(prop_0_mod_group_mon, here("mods", "stress", 
                              "glmm", "8_prop_0_mod_group_mon.rds"))
saveRDS(prop_0_mod_group_ies, here("mods", "stress", 
                              "glmm", "9_prop_0_mod_group_ies.rds"))
saveRDS(prop_0_mod_group_ams, here("mods", "stress", 
                              "glmm", "10_prop_0_mod_group_ams.rds"))
saveRDS(prop_0_mod_group_ims, here("mods", "stress", 
                              "glmm", "11_prop_0_mod_group_ims.rds"))
saveRDS(prop_0_mod_final, here("mods", "stress",
                               "glmm", "12_prop_0_mod_final.rds"))
}

# -----------------------------------------------------------------------------









# Model descriptives ----------------------------------------------------------
# 
MuMIn::r.squaredGLMM(prop_0_mod_final) 
#                   R2m       R2c
# theoretical 0.1644337 0.8341553
# delta       0.1528810 0.7755495

summary(prop_0_mod_final)
#                         Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)           2.2390     0.1803  12.420  < 2e-16 ***
#   groupaes             -1.0022     0.2271  -4.413 1.02e-05 ***
#   groupies             -1.4069     0.2242  -6.275 3.49e-10 ***
#   groupams             -1.3068     0.2259  -5.786 7.21e-09 ***
#   groupims             -1.6783     0.2256  -7.438 1.02e-13 ***
  
confint(prop_0_mod_final, method = "Wald")
#                            2.5 %       97.5 %
# (Intercept)           1.885679  2.5923575
# groupaes             -1.447381 -0.5571072
# groupies             -1.846331 -0.9674771
# groupams             -1.749521 -0.8641394
# groupims             -2.120526 -1.2360267



# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
prop_0_mod_final %>%
  group_by(., group, condition_sum, participant) %>%
  #  summarise(., meanFix = mean(target_prop)) %>%
  ggplot(., aes(x = group, y = mean(target_prop), 
                dodge = condition_sum, color = condition_sum, 
                group = interaction(group, condition_sum))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Group') + 
  #  scale_x_discrete(labels = c("mon", "aes", "ies", "ams", "ims")) +         
  ggtitle('Mean target fixations as a function of group\nand target type') +    
  #  scale_color_brewer(palette = "Set1", name = '', labels = c('Present', 'Preterit')) +       
  theme_bw(base_size = 16, base_family = 'Times') -> stress_rel_target_fix

# Graph to check the effects of WM

prop_0_mod_final %>%
#  group_by(., group, condition_sum, participant) %>%
  ggplot(., aes(x = group, y = mean(target_prop), 
                dodge = WM_set, color = WM_set,
                group = interaction(group, WM_set))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Group') + 
  #  scale_x_discrete(labels = c("mon", "aes", "ies", "ams", "ims")) +       
  ggtitle('Mean target fixations as a function of group and verbal WM') + 
  # scale_color_brewer(palette = "Set1", name = '', labels = c('Present', 'Preterit')) + 
  theme_bw(base_size = 16, base_family = 'Times') -> stress_rel_coda_fix


