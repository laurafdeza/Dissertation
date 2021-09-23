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
stress50 <- read_csv(here("data", "clean", "stress_50ms_allparticipants.csv"))

# Get path to saved models
gca_ss_mods_path  <- here("mods", "stress", "gca")
gca_en_mods_path  <- here("mods", "use_prof", "gca")

# Load models as lists
load(paste0(gca_ss_mods_path, "/gca_mon_mods.Rdata"))
load(paste0(gca_en_mods_path, "/gca_en_mods.Rdata"))
# load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
load(paste0(gca_ss_mods_path, "/model_preds.Rdata")) # and then select "fits_all_mon" and  "target_offset_preds_mon"
# # load(paste0(gca_mods_path, "/model_preds_en.Rdata"))

# Store objects in global env
list2env(gca_mon_mods, globalenv())
list2env(gca_en_mods, globalenv())
# list2env(gca_l2_mods_sum, globalenv())
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



# stress_50 <- na.omit(stress50)


stress_gc_subset <- stress50 %>%
  # select(., -WM_set) %>%
  filter(., time_zero >= -4 & time_zero <= 12 & l1 != 'ma') %>%
  select(., -prof) %>%
  mutate(., #l1 = fct_relevel(l1, "es", "en", "ma"),
            condition_sum = if_else(cond == "1", -1, 1)) %>%       # 1 = present, 2 = past
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")



# -----------------------------------------------------------------------------







#################### MONOLINGUAL SPEAKERS ########################################

# Build up random effects to test time terms
if(F){
  
  mon_data <- filter(stress_gc_subset, l1 == 'es') 
  
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + condition_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),   # , optCtrl=list(maxfun=2e5)
         data = mon_data, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + condition_sum + ot1 | participant) +
             ot2 + (1 + condition_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + condition_sum + ot1 + ot2 | participant) +
             ot3 + (1 + condition_sum + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    9 44274 44337 -22128    44256                         
  # mod_ot2   14 44165 44263 -22069    44137 118.87  5  < 2.2e-16 ***
  # mod_ot3   20 44058 44198 -22009    44018 119.48  6  < 2.2e-16 ***

  
  mod_ot0 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + ot1 | target) +
                       + (1 + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + ot1 + ot2 | target) +
                       + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot3    20 44058 44198 -22009    44018                          
  # mod_ot0    21 43919 44066 -21938    43877 141.042  1    < 2e-16 ***
  # mod_ot1a   23 43705 43866 -21830    43659 217.748  2    < 2e-16 ***
  # mod_ot2a   26 43621 43803 -21785    43569  89.723  3    < 2e-16 ***
  # mod_ot3a   30 43619 43829 -21779    43559  10.426  4    0.03384 *  
  
}



# Individual model MON -----------------------------------------------------------

gca_mod_mon_base <- mod_ot3a
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
  #        (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 + ot3 | target),
  #      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
  #      REML = F,
  #      data = filter(mon_data)) 

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + condition_sum) 
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0,   . ~ . + ot1:condition_sum) 
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1,   . ~ . + ot2:condition_sum) 
gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2,   . ~ . + ot3:condition_sum) 

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2, gca_mod_mon_cond_3)
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base     30 43619 43829 -21779    43559                     
# gca_mod_mon_cond_0   31 43621 43838 -21779    43559 0.2552  1     0.6134
# gca_mod_mon_cond_1   32 43622 43846 -21779    43558 0.8792  1     0.3484
# gca_mod_mon_cond_2   33 43624 43855 -21779    43558 0.0309  1     0.8605
# gca_mod_mon_cond_3   34 43625 43863 -21779    43557 0.3546  1     0.5515


gca_mod_mon_final <- gca_mod_mon_base


# ---------------------------------------------------------------------------






#################### L2 SPEAKERS ########################################

l2_data <- stress_gc_subset%>%
  filter(., l1 == 'en') %>% 
  mutate(., use_z = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week),
         DELE_z = (DELE - mean(DELE))/sd(DELE))

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){

  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + condition_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),   # , optCtrl=list(maxfun=2e5)
         data = l2_data, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + condition_sum + ot1 | participant) +
             ot2 + (1 + condition_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + condition_sum + ot1 + ot2 | participant) +
             ot3 + (1 + condition_sum + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  
  #           npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    9 96113 96183 -48047    96095                          
  # mod_ot2   14 96083 96191 -48027    96055  40.004  5  1.491e-07 ***
  # mod_ot3   20 95830 95985 -47895    95790 265.072  6  < 2.2e-16 ***
  
  
  
  mod_ot0 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + ot1 | target) +
                       + (1 + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + ot1 + ot2 | target) +
                       + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot3    20 95830 95985 -47895    95790                          
  # mod_ot0    21 95509 95673 -47734    95467 322.054  1  < 2.2e-16 ***
  # mod_ot1a   23 95284 95463 -47619    95238 229.683  2  < 2.2e-16 ***
  # mod_ot2a   26 95137 95339 -47543    95085 152.618  3  < 2.2e-16 ***
  # mod_ot3a   30 95089 95322 -47514    95029  56.447  4  1.616e-11 ***
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_en_mod_base <- mod_ot3a
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +         
  #        (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 + ot3 | target), 
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 2e4)),
  #      data = l2_data, REML = F)


# add stress effect to intercept, linear slope, quadratic, and cubic time terms

gca_en_mod_stress_0 <- update(gca_en_mod_base,   . ~ . + condition_sum) 
gca_en_mod_stress_1 <- update(gca_en_mod_stress_0, . ~ . + ot1:condition_sum) 
gca_en_mod_stress_2 <- update(gca_en_mod_stress_1, . ~ . + ot2:condition_sum) 
gca_en_mod_stress_3 <- update(gca_en_mod_stress_2, . ~ . + ot3:condition_sum)

en_stress_anova <-
  anova(gca_en_mod_base, gca_en_mod_stress_0, gca_en_mod_stress_1,
        gca_en_mod_stress_2, gca_en_mod_stress_3)
#                     npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)    
# gca_en_mod_base       30 95089 95322 -47514    95029                       
# gca_en_mod_stress_0   31 95091 95332 -47514    95029 0.0659  1     0.7973  
# gca_en_mod_stress_1   32 95092 95341 -47514    95028 0.4974  1     0.4806  
# gca_en_mod_stress_2   33 95088 95345 -47511    95022 5.9643  1     0.0146 *
# gca_en_mod_stress_3   34 95088 95352 -47510    95020 2.1182  1     0.1456  



# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms

gca_en_mod_dele_0 <- update(gca_en_mod_stress_2,   . ~ . + DELE_z) 
gca_en_mod_dele_1 <- update(gca_en_mod_dele_0, . ~ . + ot1:DELE_z) 
gca_en_mod_dele_2 <- update(gca_en_mod_dele_1, . ~ . + ot2:DELE_z)
gca_en_mod_dele_3 <- update(gca_en_mod_dele_2, . ~ . + ot3:DELE_z)

en_dele_anova <-
  anova(gca_en_mod_stress_2, gca_en_mod_dele_0, gca_en_mod_dele_1,
        gca_en_mod_dele_2, gca_en_mod_dele_3)
#                     npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_en_mod_stress_2   33 95088 95345 -47511    95022                       
# gca_en_mod_dele_0     34 95090 95354 -47511    95022 0.1282  1    0.72029  
# gca_en_mod_dele_1     35 95087 95359 -47509    95017 4.6084  1    0.03181 *
# gca_en_mod_dele_2     36 95086 95365 -47507    95014 3.8899  1    0.04858 *
# gca_en_mod_dele_3     37 95086 95373 -47506    95012 1.6546  1    0.19834  



# add en use effect to intercept, linear slope, quadratic, and cubic time terms

gca_en_mod_use_0 <- update(gca_en_mod_dele_2,  . ~ . + use_z) 
gca_en_mod_use_1 <- update(gca_en_mod_use_0, . ~ . + ot1:use_z) 
gca_en_mod_use_2 <- update(gca_en_mod_use_1, . ~ . + ot2:use_z)
gca_en_mod_use_3 <- update(gca_en_mod_use_2, . ~ . + ot3:use_z)

en_use_anova <-
  anova(gca_en_mod_dele_2, gca_en_mod_use_0, gca_en_mod_use_1,
        gca_en_mod_use_2, gca_en_mod_use_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_en_mod_dele_2   36 95086 95365 -47507    95014                       
# gca_en_mod_use_0    37 95087 95374 -47506    95013 0.9609  1    0.32697  
# gca_en_mod_use_1    38 95088 95383 -47506    95012 0.6941  1    0.40478  
# gca_en_mod_use_2    39 95089 95392 -47506    95011 0.5591  1    0.45462  
# gca_en_mod_use_3    40 95085 95396 -47503    95005 5.9072  1    0.01508 *



# add interaction

gca_en_mod_int_0 <- update(gca_en_mod_use_3,   . ~ . + use_z:condition_sum:DELE_z) 
gca_en_mod_int_1 <- update(gca_en_mod_int_0, . ~ . + ot1:use_z:condition_sum:DELE_z) 
gca_en_mod_int_2 <- update(gca_en_mod_int_1, . ~ . + ot2:use_z:condition_sum:DELE_z)
gca_en_mod_int_3 <- update(gca_en_mod_int_2, . ~ . + ot3:use_z:condition_sum:DELE_z)

en_int_anova <-
  anova(gca_en_mod_use_3, gca_en_mod_int_0, gca_en_mod_int_1,
        gca_en_mod_int_2, gca_en_mod_int_3)
#                       npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_en_mod_use_3   40 95085 95396 -47503    95005                     
# gca_en_mod_int_0   41 95086 95405 -47502    95004 1.4667  1     0.2259
# gca_en_mod_int_1   42 95088 95414 -47502    95004 0.0732  1     0.7867
# gca_en_mod_int_2   43 95089 95424 -47502    95003 0.4179  1     0.5180
# gca_en_mod_int_3   44 95091 95433 -47502    95003 0.2777  1     0.5982


gca_en_mod_final <- gca_en_mod_use_3

summary(gca_en_mod_final)
#                   Estimate Std. Error t value
# (Intercept)       1.25875    0.09220  13.652
# ot1                5.53887    0.31597  17.530
# ot2                0.28544    0.27272   1.047
# ot3               -1.29942    0.21926  -5.926
# condition_sum     -0.02302    0.07058  -0.326
# DELE_z             0.09316    0.06174   1.509
# use_z              0.01694    0.06558   0.258
# ot1:condition_sum -0.20997    0.21668  -0.969
# ot2:condition_sum -0.46575    0.17295  -2.693
# ot1:DELE_z         0.37150    0.21944   1.693
# ot2:DELE_z        -0.32132    0.14898  -2.157
# ot1:use_z          0.35992    0.22773   1.580
# ot2:use_z          0.08984    0.14925   0.602
# ot3:use_z         -0.35469    0.14181  -2.501 

# save Correlation matrix
vc <- vcov(gca_en_mod_final)

# diagonal matrix of standard deviations associated with vcov
S <- sqrt(diag(diag(vc), nrow(vc), nrow(vc)))

# convert vc to a correlation matrix
solve(S) %*% vc %*% solve(S) 

# detect multicollinearity in a regression model: 
# if Variance inflation factor (VIF) score high, problematic (e.g., VIF = 9, usually >5)
# the VIF for a regression model variable is equal to the ratio 
# of the overall model variance to the variance of a model 
# that includes only that single independent variable.
car::vif(gca_en_mod_final)
# ot1               ot2               ot3 
# 1.069002          1.522301          1.436685 
# condition_sum            DELE_z             use_z 
# 1.053702          1.340362          1.511513 
# ot1:condition_sum ot2:condition_sum        ot1:DELE_z 
# 1.071644          1.025330          1.342828 
# ot2:DELE_z         ot1:use_z         ot2:use_z 
# 1.124697          1.452632          1.134742 
# ot3:use_z 
# 1.170925 
eigen(vc)
# An eigenvector is the direction of the line (vertical, horizontal, 45 degrees etc.).
# An eigenvalue is a number, telling you how much variance there is in the data in that direction, 
# how spread out the data is on the line. 
# The eigenvector with the highest eigenvalue is therefore the principal component.
# In a regression analysis, if there is one or more small eigenvalues of the X'X, matrix
# it implies there are near-linear relationships (dependencies) among regressors 
# and we may have multicollinearity issues in the analysis.
performance::check_collinearity(gca_en_mod_final)
# Low Correlation
# 
# Term               VIF Increased SE Tolerance
# ot1               1.07         1.03      0.94
# ot2               1.52         1.23      0.66
# ot3               1.44         1.20      0.70
# condition_sum     1.05         1.03      0.95
# DELE_z            1.34         1.16      0.75
# use_z             1.51         1.23      0.66
# ot1:condition_sum 1.07         1.04      0.93
# ot2:condition_sum 1.03         1.01      0.98
# ot1:DELE_z        1.34         1.16      0.74
# ot2:DELE_z        1.12         1.06      0.89
# ot1:use_z         1.45         1.21      0.69
# ot2:use_z         1.13         1.07      0.88
# ot3:use_z         1.17         1.08      0.85

coll_info <- performance::check_collinearity(gca_en_mod_final)
plot(coll_info)


# # autocorrelation of residuals
# lmtest::dwtest(l2_data$DELE_z ~ l2_data$use_z)
# # Durbin-Watson test [not robust for anomalies in the data]
# # 
# # data:  l2_data$DELE_z ~ l2_data$use_z
# # DW = 0.0034153, p-value < 2.2e-16
# # alternative hypothesis: true autocorrelation is greater than 0
# 
# 
# acf(residuals(gca_en_mod_final, retype="normalized"))
# pacf(residuals(gca_en_mod_final, retype="normalized"))
# #Each spike that rises above or falls below the dashed lines is considered to be statistically significant. 
# #This means the spike has a value that is significantly different from zero. 
# #If a spike is significantly different from zero, that is evidence of autocorrelation. 
# #A spike that’s close to zero is evidence against autocorrelation.
# # Also, the smaller the AIC value the better

plot(residuals(gca_en_mod_final, retype="normalized")~l2_data$DELE_z)
plot(residuals(gca_en_mod_final, retype="normalized")~l2_data$use_z)

# Breusch-Godfrey test for autocorrelation
mod <- lm(l2_data$DELE_z ~ l2_data$use_z)
res = mod$res 
n = length(res) 
mod2 = lm(res[-n] ~ res[-1]) 
summary(mod2) # If the coefficient for res[-1] is significant, you have evidence of autocorrelation in the residuals.

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -1.77544 -0.00111 -0.00003  0.00117  2.36499 
# 
# Coefficients:
#              Estimate Std. Error  t value Pr(>|t|)    
# (Intercept) 8.286e-05  4.198e-04    0.197    0.844    
# res[-1]     9.983e-01  4.415e-04 2261.342   <2e-16 ***
# ---
# Signif. codes:  
# 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05556 on 17510 degrees of freedom
# Multiple R-squared:  0.9966,	Adjusted R-squared:  0.9966 
# F-statistic: 5.114e+06 on 1 and 17510 DF,  p-value: < 2.2e-16





mod_type <- "gca_en_mod"
mod_spec <- c('_base', 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3", 
              "_dele_0", "_dele_1", "_dele_2", "_dele_3", 
              "_use_0", "_use_1", "_use_2", "_use_3", 
              "_int_0", "_int_1", "_int_2", "_int_3",
              '_final') 

# Store ind models in list
gca_en_mods <- mget(c(paste0(mod_type, mod_spec)))

save(gca_en_mods,
     file = here("mods", "use_prof", "gca",
                 "gca_en_mods.Rdata"))




# Save anova model comparisons
en_nested_model_comparisons <-
  mget(c(
         'en_stress_anova', 
         'en_dele_anova', 
         'en_use_anova', 'en_int_anova'
  ))

save(en_nested_model_comparisons,
     file = here("mods", "use_prof", "gca",
                 "nested_model_comparisons_en.Rdata"))


}

# -----------------------------------------------------------------------------








# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
en_dat <- l2_data %>%
  dplyr::select(time_zero, ot1:ot3, condition_sum) %>%
  distinct %>%
  expand_grid(., tibble(DELE_z = c(-1, 0, 1))) %>%
  expand_grid(., tibble(use_z = c(-1, 0, 1)))
  



# Get model predictions and SE
fits_all_en <- predictSE(gca_en_mod_final, en_dat) %>%
  as_tibble %>%
  bind_cols(en_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at target syllable offset
target_offset_preds_en <- filter(fits_all_en, time_zero == 4) %>% #
  select(Stress = condition_sum, DELE = DELE_z, 
         `Spanish use` = use_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))


# Save models predictions
model_preds_en <- mget(c(
  "fits_all_en", 
  "target_offset_preds_en"))

save(model_preds_en,
     file = here("mods", "use_prof", "gca",
                 "model_preds_en.Rdata"))
