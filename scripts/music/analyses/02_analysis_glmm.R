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

prop_0_mod_0     <- readRDS(here("reports", "mods", "glmm", "0_prop_0_mod_0.rds"))
prop_0_mod_group <- readRDS(here("reports", "mods", "glmm", "1_prop_0_mod_group.rds"))
prop_0_mod_cond  <- readRDS(here("reports", "mods", "glmm", "2_prop_0_mod_cond.rds"))
prop_0_mod_coda  <- readRDS(here("reports", "mods", "glmm", "3_prop_0_mod_coda.rds"))
prop_0_mod_int1  <- readRDS(here("reports", "mods", "glmm", "4_prop_0_mod_int1.rds"))
prop_0_mod_int2  <- readRDS(here("reports", "mods", "glmm", "5_prop_0_mod_int2.rds"))
prop_0_mod_int3  <- readRDS(here("reports", "mods", "glmm", "6_prop_0_mod_int3.rds"))
prop_0_mod_full  <- readRDS(here("reports", "mods", "glmm", "7_prop_0_mod_full.rds"))
prop_0_mod_final <- readRDS(here("reports", "mods", "glmm", "8_prop_0_mod_final.rds"))
# -----------------------------------------------------------------------------



# Data prep -------------------------------------------------------------------

# Get data
# Filter time course to offset of 1st syllable (time_zero == 20)
# Create sum coded fixed factors (condition)


df_stress <- stress10 %>%
  filter(., time_zero == 20) %>%
  mutate(., condition_sum = if_else(cond == "1", 1, -1)) #        1 = present


# -----------------------------------------------------------------------------

# Adding verbal WM, pitch and rhythm as covariates


# Add verbal working memory 
###### and phonotactic frequency as a covariate    ->    NOT ADDED

dem <- read_csv(here("data", "pupurri_analysis.csv"))
dem <- dem %>%
  select(., participant, WM_set)

dem$participant <- tolower(dem$participant)




pitch <- read_csv(here("data", "clean", "pitch.csv"))

pitch <- pitch %>%
  select(., -X1) %>%
  rename(., participant = subject_id) %>%
  mutate(participant = str_replace(participant, "ae", "aes")) %>%
  mutate(participant = str_replace(participant, "ie", "ies")) %>%
  mutate(participant = str_replace(participant, "am", "ams")) %>%
  mutate(participant = str_replace(participant, "im", "ims")) %>%
  mutate(participant = str_replace(participant, "mo", "mon"))





rhythm <- read_csv(here("data", "clean", "rhythm.csv"))

rhythm <- rhythm %>%
  select(., -X1) %>%
  rename(., participant = subject_id) %>%
  mutate(participant = str_replace(participant, "ae", "aes")) %>%
  mutate(participant = str_replace(participant, "ie", "ies")) %>%
  mutate(participant = str_replace(participant, "am", "ams")) %>%
  mutate(participant = str_replace(participant, "im", "ims")) %>%
  mutate(participant = str_replace(participant, "mo", "mon"))




# Add verbal wm score to eyetracking data frame
df_stress <- merge(x = df_stress, y = dem, by = "participant", all.x=TRUE)
df_stress <- merge(x = df_stress, y = pitch, by = "participant", all.x=TRUE)
df_stress <- merge(x = df_stress, y = rhythm, by = "participant", all.x=TRUE)


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
# prop_0_ranefA  2 24283 24295 -12140    24279                             
# prop_0_ranefB  3 23664 23682 -11829    23658 620.75      1  < 2.2e-16 ***

prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 + condition_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep slope for condition bc significant
#               Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
# prop_0_ranefB  3 23664 23682 -11829    23658                             
# prop_0_ranefC  5 22393 22422 -11192    22383 1274.9      2  < 2.2e-16 ***

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
prop_0_mod_int1  <- update(prop_0_mod_cond,  . ~ . + group:condition_sum)


anova(prop_0_mod_0, prop_0_mod_group, test = "Chisq")    # main effect of group bc significant
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)   
# prop_0_mod_0      5 22393 22422 -11192    22383                            
# prop_0_mod_group  9 22383 22436 -11182    22365 18.331      4   0.001063 **
anova(prop_0_mod_group, prop_0_mod_cond, test = "Chisq") # no effect of condition bc not significant
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# prop_0_mod_group  9 22383 22436 -11182    22365                         
# prop_0_mod_cond  10 22384 22442 -11182    22364 1.3647      1     0.2427
anova(prop_0_mod_group, prop_0_mod_int1, test = "Chisq") # no interaction
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# prop_0_mod_group  9 22383 22436 -11182    22365                         
# prop_0_mod_int1  14 22390 22472 -11181    22362 3.1652      5     0.6745
summary(prop_0_mod_group) 
# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   2.0975     0.2318   9.048  < 2e-16 ***
# groupams     -0.5335     0.3008  -1.774  0.07610 .  
# groupies     -0.4760     0.2993  -1.590  0.11174    
# groupims     -0.8539     0.3000  -2.846  0.00442 ** 
# groupmon      0.3604     0.3105   1.161  0.24572    



df_stress$group <- factor(df_stress$group, levels = c("mon", "aes", "ies", "ams", "ims"))

prop_0_mod_final <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        
                    group +
                    (1 + condition_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))


summary(prop_0_mod_final)
# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   2.4579     0.2401  10.238  < 2e-16 ***
# groupaes     -0.3604     0.3107  -1.160  0.24600    
# groupies     -0.8364     0.3059  -2.734  0.00625 ** 
# groupams     -0.8940     0.3080  -2.903  0.00370 ** 
# groupims     -1.2144     0.3084  -3.937 8.25e-05 ***


}

# -----------------------------------------------------------------------------


# Add covariates to model

prop_0_mod_wm <- update(prop_0_mod_final,     . ~ . + WM_set)
prop_0_mod_int2 <- update(prop_0_mod_final,     . ~ . + group:WM_set)

anova(prop_0_mod_final, prop_0_mod_wm, test = "Chisq")
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# prop_0_mod_final  9 22383 22436 -11182    22365                         
# prop_0_mod_wm    10 22385 22443 -11182    22365 0.3064      1     0.5799

anova(prop_0_mod_final, prop_0_mod_int2, test = "Chisq")
#                  Df   AIC   BIC logLik deviance Chisq Chi Df Pr(>Chisq)
# prop_0_mod_final  9 22383 22436 -11182    22365                        
# prop_0_mod_int2  14 22392 22473 -11182    22364 1.475      5     0.9159


prop_mod_DO_down <- update(prop_0_mod_final,     . ~ . + DO_down)
prop_mod_DO_up <- update(prop_mod_DO_down, . ~ . + DO_up)
prop_mod_MI_down <- update(prop_mod_DO_up,     . ~ . + MI_down)
prop_mod_MI_up <- update(prop_mod_MI_down, . ~ . + MI_up)
prop_mod_SOL_down <- update(prop_mod_MI_up,     . ~ . + SOL_down)
prop_mod_SOL_up <- update(prop_mod_SOL_down, . ~ . + SOL_up)
prop_mod_int3 <- update(prop_mod_SOL_up,     . ~ . + group:DO_down)
prop_mod_int4 <- update(prop_mod_int3,     . ~ . + group:DO_up)
prop_mod_int5 <- update(prop_mod_int4,     . ~ . + group:MI_down) # fails to converge
prop_mod_int6 <- update(prop_mod_int5,     . ~ . + group:MI_up)

prop_mod_int7 <- update(prop_mod_int6,     . ~ . + group:SOL_down)

prop_mod_int8 <- update(prop_mod_int7,     . ~ . + group:SOL_up)

anova(prop_0_mod_final, prop_mod_DO_down, test = "Chisq")
anova(prop_0_mod_final, prop_0_mod_int2, test = "Chisq")














# Save models -----------------------------------------------------------------

if(F) {

saveRDS(prop_0_mod_0, here("reports", "mods",
                            "glmm", "0_prop_0_mod_0.rds"))
saveRDS(prop_0_mod_group, here("reports", "mods",
                               "glmm", "1_prop_0_mod_group.rds"))
saveRDS(prop_0_mod_cond, here("reports", "mods",
                              "glmm", "2_prop_0_mod_cond.rds"))
# saveRDS(prop_0_mod_coda, here("models", "stress", "1_adv_int_mon",
#                               "glmm", "3_prop_0_mod_coda.rds"))
# saveRDS(prop_0_mod_int1, here("2CrisLau", "models", "stress", "1_adv_int_mon",
#                               "glmm", "4_prop_0_mod_int1.rds"))
saveRDS(prop_0_mod_int2, here("reports", "mods",
                            "glmm", "5_prop_0_mod_int1.rds"))
# saveRDS(prop_0_mod_int3, here("2CrisLau", "models", "stress", "1_adv_int_mon",
#                               "glmm", "6_prop_0_mod_int3.rds"))
# saveRDS(prop_0_mod_full, here("2CrisLau", "models", "stress", "1_adv_int_mon",
#                               "glmm", "7_prop_0_mod_full.rds"))
saveRDS(prop_0_mod_final, here("reports", "mods",
                               "glmm", "8_prop_0_mod_final.rds"))
}

# -----------------------------------------------------------------------------









# Model descriptives ----------------------------------------------------------
# 
MuMIn::r.squaredGLMM(prop_0_mod_final) #original script said full instead of final
#                   R2m       R2c
# theoretical 0.1213596 0.7382825
# delta       0.1148699 0.6988027

summary(prop_0_mod_final)
confint(prop_0_mod_final, method = "Wald")

# Fixed effects:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.295702   0.126700  -2.334   0.0196 *  
# groupams     0.286862   0.156219   1.836   0.0663 .  
# groupies     0.203980   0.157013   1.299   0.1939    
# groupims     0.004274   0.157878   0.027   0.9784    
# groupmon     1.079088   0.160794   6.711 1.93e-11 ***


#                   2.5 %      97.5 %
# (Intercept) -0.54402887 -0.04737495
# groupams    -0.01932276  0.59304580
# groupies    -0.10376047  0.51172048
# groupims    -0.30516002  0.31370847
# groupmon     0.76393733  1.39423920







# Relevel to test other group vs the others
df_stress$group <- factor(df_stress$group, levels = c("aes", "ies", "ams", "ims", "mon"))

summary(glmer(cbind(targetCount, 10 - targetCount) ~ 1 +      #change according to results above
        group + coda_sum +
        (1 + condition_sum * coda_sum | participant) +
        (1 | target),
        data = df_stress, family = 'binomial',
        control = glmerControl(optimizer = 'bobyqa')))

# Fixed effects

# -----------------------------------------------------------------------------



##############################################################################
## This has been copied from an old script, needs to be adapted to this one ##
##############################################################################


# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
stress_df_0_prop %>%
#  filter(., exp == "stress_related") %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix, 
                dodge = condition, color = condition,
                group = interaction(group, condition))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Group') + 
  scale_x_discrete(labels = c("mon", "aes", "ies", "ams", "ims")) +         # reorder?
  ggtitle('Mean target fixations as a function of group and target type') +    
  scale_color_brewer(palette = "Set1", name = '', labels = c('Present', 'Preterit')) +       #reorder
  theme_bw(base_size = 16, base_family = 'Times') -> stress_rel_target_fix





### WM ACCOUNTED

dem <- read_csv(here("data", "pupurri_analysis.csv"))
dem <- dem %>%
  select(., participant, WM_set)
dem$participant <- tolower(dem$participant)

stress10wm <- merge(x = stress10, y = dem[ , c("participant", "WM_set")], by = "participant", all.x=TRUE)

df_stress <- stress10wm %>%
  filter(., time_zero == 20) %>%
  mutate(., condition_sum = if_else(cond == "1", 1, -1))

# RANDOM EFFECTS -----------------------------------------------------------

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

prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                         (1 + condition_sum | participant) +
                         (1 | target),
                       data = df_stress, family = 'binomial', 
                       control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep slope for condition bc significant



prop_0_ranefD <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                         (1 + condition_sum + WM_set | participant) +
                         (1 | target),
                       data = df_stress, family = 'binomial', 
                       control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefC, prop_0_ranefD, refit = F) # Don't keep slope for WM

prop_0_ranefE <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                         (1 + condition_sum * WM_set | participant) +
                         (1 | target), 
                       data = df_stress, family = 'binomial',
                       control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefC, prop_0_ranefE, refit = F) # Don't keep interaction


# FIXED EFFECTS -------------------------------------------------------------

prop_0_mod_0 <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                        (1 + condition_sum | participant) +   # change random effects according to results from above
                        (1 | target),
                      data = df_stress, family = 'binomial', 
                      control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- update(prop_0_mod_0,     . ~ . + group)
prop_0_mod_cond  <- update(prop_0_mod_group, . ~ . + condition_sum)
prop_0_mod_wm  <- update(prop_0_mod_cond,  . ~ . + WM_set)
prop_0_mod_int1  <- update(prop_0_mod_wm,  . ~ . + group:condition_sum)
prop_0_mod_int2  <- update(prop_0_mod_int1,  . ~ . + group:WM_set)
prop_0_mod_int3  <- update(prop_0_mod_int2,  . ~ . + WM_set:condition_sum)
prop_0_mod_full  <- update(prop_0_mod_int3,  . ~ . + group:WM_set:condition_sum)


anova(prop_0_mod_0, prop_0_mod_group, test = "Chisq")    # main effect of group bc significant
anova(prop_0_mod_group, prop_0_mod_cond, test = "Chisq") # no effect of condition bc not significant
anova(prop_0_mod_group, prop_0_mod_wm, test = "Chisq") # no effect of wm
anova(prop_0_mod_group, prop_0_mod_int1, test = "Chisq")  # no group x condition interaction
anova(prop_0_mod_group, prop_0_mod_int2, test = "Chisq")  # no group x wm interaction
anova(prop_0_mod_group, prop_0_mod_int3, test = "Chisq")  # no condition x wm interaction
anova(prop_0_mod_group, prop_0_mod_full, test = "Chisq")  # no three way interaction

summary(prop_0_mod_group) # aes as reference
# Fixed effects:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.295703   0.126733  -2.333   0.0196 *  
# groupams     0.286862   0.156258   1.836   0.0664 .  
# groupies     0.203980   0.157043   1.299   0.1940    
# groupims     0.004274   0.157911   0.027   0.9784    
# groupmon     1.079087   0.160839   6.709 1.96e-11 ***

df_stress$group <- factor(df_stress$group, levels = c("mon", "aes", "ies", "ams", "ims"))

prop_0_mod_final <- glmer(cbind(target_count, 10 - target_count) ~ 1 +        # modify according to anova results above
                            group +
                            (1 + condition_sum | participant) +
                            (1 | target),
                          data = df_stress, family = 'binomial',
                          control = glmerControl(optimizer = 'bobyqa'))
summary(prop_0_mod_final) # mon as reference
#   Fixed effects:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)   0.7834     0.1287   6.087 1.15e-09 ***
#   groupaes     -1.0791     0.1609  -6.708 1.97e-11 ***
#   groupies     -0.8751     0.1579  -5.541 3.00e-08 ***
#   groupams     -0.7922     0.1601  -4.949 7.45e-07 ***
#   groupims     -1.0748     0.1592  -6.752 1.46e-11 ***