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
  filter(., time_zero >= -4 & time_zero <= 8 & l1 != 'ma') %>%
  select(., -prof) %>%
  mutate(., #l1 = fct_relevel(l1, "es", "en", "ma"),
            condition_sum = if_else(cond == "1", -1, 1)) %>%       # 1 = present (paroxytone), 2 = preterite (oxytone)
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

  mod_ot0 <-
    lmer(eLog ~ 1 + 
           (1 + condition_sum | participant),
         control = lmerControl(optimizer = 'bobyqa'),   # , optCtrl=list(maxfun=2e5)
         data = l2_data, weights = 1/wts, REML = F)
  
  mod_ot1 <-
    update(mod_ot0, . ~ . -(1 + condition_sum | participant) +
              ot1 + (1 + condition_sum + ot1 | participant))
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + condition_sum + ot1 | participant) +
             ot2 + (1 + condition_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + condition_sum + ot1 + ot2 | participant) +
             ot3 + (1 + condition_sum + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot0, mod_ot1, mod_ot2, mod_ot3)
  
  #         npar   AIC   BIC logLik deviance    Chisq Df Pr(>Chisq)    
  # mod_ot0    5 76024 76062 -38007    76014                           
  # mod_ot1    9 74693 74761 -37338    74675 1338.732  4  < 2.2e-16 ***
  # mod_ot2   14 74497 74603 -37235    74469  205.794  5  < 2.2e-16 ***
  # mod_ot3   20 74481 74631 -37221    74441   28.227  6  8.514e-05 ***
  
  
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                       + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                       + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #          npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot3    20 74481 74631 -37221    74441                           
  # mod_ot4   21 74344 74502 -37151    74302 138.9667  1  < 2.2e-16 ***
  # mod_ot5   23 74273 74445 -37113    74227  75.5274  2  < 2.2e-16 ***
  # mod_ot6   26 74266 74461 -37107    74214  12.5629  3   0.005684 ** 
  # mod_ot7   30 74268 74494 -37104    74208   5.7216  4   0.220926 
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_en_mod_base <- mod_ot6
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +         
  #        (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 | target), 
  #      control = lmerControl(optimizer = 'bobyqa'), #, optCtrl = list(maxfun = 2e5)
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
# gca_en_mod_base       26 74266 74461 -37107    74214                       
# gca_en_mod_stress_0   27 74268 74471 -37107    74214 0.0770  1    0.78146  
# gca_en_mod_stress_1   28 74270 74480 -37107    74214 0.0007  1    0.97916  
# gca_en_mod_stress_2   29 74267 74484 -37104    74209 5.3626  1    0.02057 *
# gca_en_mod_stress_3   30 74269 74494 -37104    74209 0.1317  1    0.71664  



# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms

gca_en_mod_dele_0 <- update(gca_en_mod_stress_2,   . ~ . + DELE_z) 
gca_en_mod_dele_1 <- update(gca_en_mod_dele_0, . ~ . + ot1:DELE_z) 
gca_en_mod_dele_2 <- update(gca_en_mod_dele_1, . ~ . + ot2:DELE_z)
gca_en_mod_dele_3 <- update(gca_en_mod_dele_2, . ~ . + ot3:DELE_z)

en_dele_anova <-
  anova(gca_en_mod_stress_2, gca_en_mod_dele_0, gca_en_mod_dele_1,
        gca_en_mod_dele_2, gca_en_mod_dele_3)
#                     npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_en_mod_stress_2   29 74267 74484 -37104    74209                         
# gca_en_mod_dele_0     30 74269 74494 -37104    74209  0.0909  1   0.763026   
# gca_en_mod_dele_1     31 74260 74493 -37099    74198 10.1466  1   0.001446 **
# gca_en_mod_dele_2     32 74262 74503 -37099    74198  0.1349  1   0.713408   
# gca_en_mod_dele_3     33 74264 74512 -37099    74198  0.0300  1   0.862587 



# add en use effect to intercept, linear slope, quadratic, and cubic time terms

gca_en_mod_use_0 <- update(gca_en_mod_dele_1,  . ~ . + use_z) 
gca_en_mod_use_1 <- update(gca_en_mod_use_0, . ~ . + ot1:use_z) 
gca_en_mod_use_2 <- update(gca_en_mod_use_1, . ~ . + ot2:use_z)
gca_en_mod_use_3 <- update(gca_en_mod_use_2, . ~ . + ot3:use_z)

en_use_anova <-
  anova(gca_en_mod_dele_1, gca_en_mod_use_0, gca_en_mod_use_1,
        gca_en_mod_use_2, gca_en_mod_use_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_en_mod_dele_1   31 74260 74493 -37099    74198                       
# gca_en_mod_use_0    32 74262 74502 -37099    74198 0.9186  1    0.33785  
# gca_en_mod_use_1    33 74264 74511 -37099    74198 0.0692  1    0.79245  
# gca_en_mod_use_2    34 74261 74517 -37097    74193 4.0842  1    0.04329 *
# gca_en_mod_use_3    35 74260 74523 -37095    74190 3.5691  1    0.05886 .



# add interaction

gca_en_mod_int_0 <- update(gca_en_mod_use_2,   . ~ . + use_z:condition_sum:DELE_z) 
gca_en_mod_int_1 <- update(gca_en_mod_int_0, . ~ . + ot1:use_z:condition_sum:DELE_z) 
gca_en_mod_int_2 <- update(gca_en_mod_int_1, . ~ . + ot2:use_z:condition_sum:DELE_z)
gca_en_mod_int_3 <- update(gca_en_mod_int_2, . ~ . + ot3:use_z:condition_sum:DELE_z)

en_int_anova <-
  anova(gca_en_mod_use_2, gca_en_mod_int_0, gca_en_mod_int_1,
        gca_en_mod_int_2, gca_en_mod_int_3)
#                       npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_en_mod_use_2   34 74261 74517 -37097    74193                        
# gca_en_mod_int_0   35 74263 74525 -37096    74193 0.6668  1   0.414159   
# gca_en_mod_int_1   36 74264 74534 -37096    74192 0.6259  1   0.428882   
# gca_en_mod_int_2   37 74266 74543 -37096    74192 0.5780  1   0.447101   
# gca_en_mod_int_3   38 74260 74545 -37092    74184 7.7543  1   0.005358 **


gca_en_mod_final <- gca_en_mod_int_3

summary(gca_en_mod_final)
#                                  Estimate Std. Error t value
# (Intercept)                     0.7480763  0.0887181   8.432
# ot1                             4.0140242  0.2717134  14.773
# ot2                             1.2471908  0.1761377   7.081
# ot3                            -0.2323991  0.1281251  -1.814
# condition_sum                   0.0583371  0.0768372   0.759
# DELE_z                          0.1001021  0.0620580   1.613
# use_z                          -0.0244102  0.0636151  -0.384
# ot1:condition_sum              -0.0008418  0.2029645  -0.004
# ot2:condition_sum              -0.2788809  0.1242150  -2.245
# ot1:DELE_z                      0.5648338  0.1873369   3.015
# ot1:use_z                       0.1972897  0.1963398   1.005
# ot2:use_z                       0.3102034  0.1456208   2.130
# condition_sum:DELE_z:use_z     -0.0447141  0.0476480  -0.938
# ot1:condition_sum:DELE_z:use_z  0.0462346  0.0966282   0.478
# ot2:condition_sum:DELE_z:use_z -0.0828551  0.0891732  -0.929
# ot3:condition_sum:DELE_z:use_z -0.2504781  0.0893679  -2.803


# save Correlation matrix
vc <- vcov(gca_en_mods_verb$gca_en_mod_final)

# diagonal matrix of standard deviations associated with vcov
S <- sqrt(diag(diag(vc), nrow(vc), nrow(vc)))

# convert vc to a correlation matrix
solve(S) %*% vc %*% solve(S) 
#               [,6]         [,7]          [,10]

# [6,]  1.0000000000 -0.303218564    0.327288752
# [7,] -0.3032185636  1.000000000   -0.109179929
# [10,]  0.3272887522 -0.109179929    1.000000000
# [11,] -0.1006372540  0.377055665   -0.285676303
# [12,]  0.0013351496  0.221964559   -0.024054224

#             [,11]        [,12]      
# [6,] -0.100637254  0.001335150  
# [7,]  0.377055665  0.221964559 
# [10,] -0.285676303 -0.024054224  
# [11,]  1.000000000  0.324792484 
# [12,]  0.324792484  1.000000000 



# detect multicollinearity in a regression model: 
# if Variance inflation factor (VIF) score high, problematic (e.g., VIF = 9, usually >5)
# the VIF for a regression model variable is equal to the ratio 
# of the overall model variance to the variance of a model 
# that includes only that single independent variable.
car::vif(gca_en_mod_final)
# ot1                            ot2                            ot3 
# 1.029197                       1.134310                       1.115692 
# condition_sum                         DELE_z                          use_z 
# 1.076962                       1.236494                       1.304176 
# ot1:condition_sum              ot2:condition_sum                     ot1:DELE_z 
# 1.029077                       1.061443                       1.223921 
# ot1:use_z                      ot2:use_z     condition_sum:DELE_z:use_z 
# 1.371885                       1.143594                       1.056934 
# ot1:condition_sum:DELE_z:use_z ot2:condition_sum:DELE_z:use_z ot3:condition_sum:DELE_z:use_z 
# 1.039666                       1.070603                       1.014135 
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
# ot1               1.03         1.01      0.97
# ot2               1.13         1.07      0.88
# ot3               1.12         1.06      0.90
# condition_sum     1.08         1.04      0.93
# DELE_z            1.24         1.11      0.81
# use_z             1.30         1.14      0.77
# ot1:condition_sum 1.03         1.01      0.97
# ot2:condition_sum 1.06         1.03      0.94
# ot1:DELE_z        1.22         1.11      0.82
# ot1:use_z         1.37         1.17      0.73
# ot2:use_z         1.14         1.07      0.87
# condition_sum:DELE_z:use_z     1.06         1.03      0.95
# ot1:condition_sum:DELE_z:use_z 1.04         1.02      0.96
# ot2:condition_sum:DELE_z:use_z 1.07         1.03      0.93
# ot3:condition_sum:DELE_z:use_z 1.01         1.01      0.99
coll_info <- performance::check_collinearity(gca_en_mod_final)

collinearity <- plot(coll_info) 


ggsave('collinearity_verb',
       plot = collinearity, dpi = 600, device = "png",
       path = here("figs", "use_prof", "gca_verb"),
       height = 5, width = 28, units = 'in')
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
# -1.77469 -0.00145 -0.00004  0.00152  2.36372 
# 
# Coefficients:
#              Estimate Std. Error  t value Pr(>|t|)    
# (Intercept) 0.0001079  0.0005469    0.197    0.844    
# res[-1]     0.9977895  0.0005751 1734.936   <2e-16 ***
# ---
# Signif. codes:  
# 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06341 on 13439 degrees of freedom
# Multiple R-squared:  0.9956,	Adjusted R-squared:  0.9956 
# F-statistic: 3.01e+06 on 1 and 13439 DF,  p-value: < 2.2e-16





mod_type <- "gca_en_mod"
mod_spec <- c('_base', 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3", 
              "_dele_0", "_dele_1", "_dele_2", "_dele_3", 
              "_use_0", "_use_1", "_use_2", "_use_3", 
              "_int_0", "_int_1", "_int_2", "_int_3",
              '_final') 

# Store ind models in list
gca_en_mods_verb <- mget(c(paste0(mod_type, mod_spec)))

save(gca_en_mods_verb,
     file = here("mods", "use_prof", "gca",
                 "gca_en_mods_verb.Rdata"))




# Save anova model comparisons
en_nested_model_comparisons_verb <-
  mget(c(
         'en_stress_anova', 
         'en_dele_anova', 
         'en_use_anova', 'en_int_anova'
  ))

save(en_nested_model_comparisons_verb,
     file = here("mods", "use_prof", "gca",
                 "nested_model_comparisons_en_verb.Rdata"))


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
model_preds_en_verb <- mget(c(
  "fits_all_en", 
  "target_offset_preds_en"))

save(model_preds_en_verb,
     file = here("mods", "use_prof", "gca",
                 "model_preds_en_verb.Rdata"))
