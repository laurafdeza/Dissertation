#
# Original script by Joseph
# Modified by Laura for pupurri project
# Last update: 09/02/2020
#
# Growth curve analysis ------------------------------------------------------
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
gca_mods_path  <- here("mods", "use_prof", "gca")

# Load models as lists
load(paste0(gca_mods_path, "/gca_mod_mon_simple.Rdata"))
load(paste0(gca_mods_path, "/gca_en_mods_simple.Rdata"))
load(paste0(gca_ss_mods_path, "/model_preds_simple.Rdata")) 


# Store objects in global env
list2env(gca_mod_mon_simple, globalenv())
list2env(gca_en_mods_simple, globalenv())
list2env(model_preds_simple, globalenv())

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
  filter(., time_zero >= -4 & time_zero <= 4 & l1 != 'ma') %>%
  select(., -prof) %>%
  mutate(., #l1 = fct_relevel(l1, "es", "en", "ma"),
            condition_sum = if_else(cond == "1", -1, 1)) %>%       # 1 = present (paroxytone), 2 = preterite (oxytone)
  poly_add_columns(., time_zero, degree = 2, prefix = "ot")



# -----------------------------------------------------------------------------







#################### MONOLINGUAL SPEAKERS ########################################

# Build up random effects to test time terms
if(F){
  
  mon_data <- filter(stress_gc_subset, l1 == 'es') 
  
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + condition_sum | participant),
         control = lmerControl(optimizer = 'bobyqa'),   # , optCtrl=list(maxfun=2e5)
         data = mon_data, weights = 1/wts, REML = F)
  
  
  
  mod_ot2 <-
    update(mod_ot1, . ~ . + ot2) 
  
  
  anova(mod_ot1, mod_ot2)
  #         npar    AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    6 24010 24048 -11999    23998                         
  # mod_ot2    7 23993 24037 -11989    23979 19.441  1  1.037e-05 ***  
  
  
  
  mod_ot0 <- update(mod_ot2, . ~ . + (1 | target))
  
  anova(mod_ot2, mod_ot0)
  #          npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot2    7 23993 24037 -11989    23979                        
  # mod_ot0    8 23897 23948 -11940    23881 97.65  1  < 2.2e-16 ***

  
}



# Individual model MON -----------------------------------------------------------

gca_mod_mon_base <- mod_ot0
  # lmer(eLog ~ 1 + (ot1 + ot2) +
  #        (1 + condition_sum + (ot1 + ot2) | participant) +
  #        (1 + ot1 | target),
  #      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
  #      REML = F,
  #      data = filter(mon_data)) 

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + condition_sum) 
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0,   . ~ . + ot1:condition_sum) 
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1,   . ~ . + ot2:condition_sum) 

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2)
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base     8 23897 23948 -11940    23881          
# gca_mod_mon_cond_0    9 23899 23956 -11940    23881 0.0014  1   0.9704
# gca_mod_mon_cond_1   10 23899 23962 -11939    23879 2.3519  1   0.1251
# gca_mod_mon_cond_2   11 23900 23970 -11939    23878 0.1346  1   0.7137


gca_mod_mon_final <- gca_mod_mon_base


mod_type <- "gca_mod_mon"
mod_spec <- c('_base',
              '_cond_0', '_cond_1', '_cond_2',
              '_final')
# Store ind models in list
gca_mod_mon_simple <- mget(c(paste0(mod_type, mod_spec)))
save(gca_mod_mon_simple,
     file = here("mods", "use_prof", "gca",
                 "gca_mod_mon_simple.Rdata"))


summary(gca_mod_mon_final)
confint(gca_mod_mon_final)
# 2.5 %    97.5 %
  # .sig01       0.4238570 0.7736134
# .sig02       0.2068132 0.4735688
# .sig03      -0.4745975 0.7423405
# .sig04       0.1040267 0.3697595
# .sigma       2.9692408 3.0987791
# (Intercept)  0.6790108 0.6926235
# ot1          2.4168149 2.4328652
# ot2          0.5079351 0.5240516



car::vif(gca_mod_mon_final)
# ot1      ot2 
# 1.01593 1.01593 

performance::r2_nakagawa(gca_mod_mon_final, by_group = FALSE, tolerance = 1e-05)
#Conditional R2: 0.107
#Marginal R2: 0.066
#Warning message:
#  Random slopes not present as fixed effects. This artificially inflates the conditional random effect variances.
performance::r2_nakagawa(gca_mod_mon_cond_2, by_group = FALSE, tolerance = 1e-05)
# Conditional R2: 0.113
# Marginal R2: 0.066


rsq::rsq(gca_mod_mon_final)
# $model
#[1] 0.1239552

#$fixed
#[1] 0.05886097

#$random
#[1] 0.06509421

rsq::rsq(gca_mod_mon_final, adj = TRUE)
#$model
#[1] 0.1235485

#$fixed
#[1] 0.05842405

#$random
#[1] 0.06512443




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
  
  anova(mod_ot0, mod_ot1)
  
  #         npar   AIC   BIC logLik deviance    Chisq Df Pr(>Chisq)    
  # mod_ot0    5 52614 52649 -26302    52604                         
  # mod_ot1    9 52379 52443 -26180    52361 242.73  4  < 2.2e-16 ***
  
  mod_ot4 <- update(mod_ot0, . ~ . + (1 | target)) # singular if base mod_ot1
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + ot1 + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                       + ot2 + (1 + ot1 + ot2 | target))
  
  anova(mod_ot0, mod_ot4, mod_ot5, mod_ot6)
  #         npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot0    5 52614 52649 -26302    52604                          
  # mod_ot4    6 52555 52598 -26271    52543  60.971  1  5.792e-15 ***
  # mod_ot5    9 52280 52344 -26131    52262 280.576  3  < 2.2e-16 ***
  # mod_ot6   13 52233 52325 -26103    52207  55.533  4  2.512e-11 ***
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_en_mod_base <- update(mod_ot4, . ~ . + ot1 + ot2) # singular if more complex ranef


# add stress effect to intercept, linear slope, quadratic, and cubic time terms

gca_en_mod_stress_0 <- update(gca_en_mod_base,   . ~ . + condition_sum) 
gca_en_mod_stress_1 <- update(gca_en_mod_stress_0, . ~ . + ot1:condition_sum) 
gca_en_mod_stress_2 <- update(gca_en_mod_stress_1, . ~ . + ot2:condition_sum) 

en_stress_anova <-
  anova(gca_en_mod_base, gca_en_mod_stress_0, gca_en_mod_stress_1,
        gca_en_mod_stress_2)
#                     npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
# gca_en_mod_base        8 52263 52320 -26124    52247                        
# gca_en_mod_stress_0    9 52264 52328 -26123    52246 0.9764  1   0.323099   
# gca_en_mod_stress_1   10 52257 52328 -26118    52237 9.1602  1   0.002473 **
# gca_en_mod_stress_2   11 52258 52337 -26118    52236 0.5969  1   0.439781                


# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms

gca_en_mod_dele_0 <- update(gca_en_mod_stress_1,   . ~ . + DELE_z) 
gca_en_mod_dele_1 <- update(gca_en_mod_dele_0, . ~ . + ot1:DELE_z) 
gca_en_mod_dele_2 <- update(gca_en_mod_dele_1, . ~ . + ot2:DELE_z)

en_dele_anova <-
  anova(gca_en_mod_stress_1, gca_en_mod_dele_0, gca_en_mod_dele_1,
        gca_en_mod_dele_2)
#                     npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)   
# gca_en_mod_stress_1   10 52257 52328 -26118    52237                         
# gca_en_mod_dele_0     11 52259 52337 -26118    52237  0.0108  1   0.917296   
# gca_en_mod_dele_1     12 52251 52337 -26114    52227 10.0495  1   0.001524 **
# gca_en_mod_dele_2     13 52253 52346 -26113    52227  0.0138  1   0.906650 



# add en use effect to intercept, linear slope, quadratic, and cubic time terms

gca_en_mod_use_0 <- update(gca_en_mod_dele_1,  . ~ . + use_z) 
gca_en_mod_use_1 <- update(gca_en_mod_use_0, . ~ . + ot1:use_z) 
gca_en_mod_use_2 <- update(gca_en_mod_use_1, . ~ . + ot2:use_z)

en_use_anova <-
  anova(gca_en_mod_dele_1, gca_en_mod_use_0, gca_en_mod_use_1,
        gca_en_mod_use_2)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_en_mod_dele_1   12 52251 52337 -26114    52227                        
# gca_en_mod_use_0    13 52251 52344 -26112    52225 2.0395  1   0.153262   
# gca_en_mod_use_1    14 52252 52352 -26112    52224 0.4623  1   0.496538   
# gca_en_mod_use_2    15 52248 52355 -26109    52218 6.6751  1   0.009777 **                       


# add interaction

gca_en_mod_int_0 <- update(gca_en_mod_use_2,   . ~ . + use_z:condition_sum:DELE_z) 
gca_en_mod_int_1 <- update(gca_en_mod_int_0, . ~ . + ot1:use_z:condition_sum:DELE_z) 
gca_en_mod_int_2 <- update(gca_en_mod_int_1, . ~ . + ot2:use_z:condition_sum:DELE_z)

en_int_anova <-
  anova(gca_en_mod_use_2, gca_en_mod_int_0, gca_en_mod_int_1,
        gca_en_mod_int_2)
#                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_en_mod_use_2   15 52248 52355 -26109    52218                     
# gca_en_mod_int_0   16 52249 52363 -26108    52217 1.1133  1     0.2914
# gca_en_mod_int_1   17 52250 52371 -26108    52216 0.9115  1     0.3397
# gca_en_mod_int_2   18 52250 52378 -26107    52214 1.8093  1     0.1786


gca_en_mod_final <- gca_en_mod_use_2

summary(gca_en_mod_final)
#                                Estimate Std. Error t value
# (Intercept)                     0.07026    0.08645   0.813
# ot1                             1.46896    0.09409  15.612
# ot2                             0.66833    0.09423   7.092
# condition_sum                   0.07445    0.08301   0.897
# DELE_z                          0.02204    0.05844   0.377
# use_z                          -0.07715    0.05832  -1.323
# ot1:condition_sum               0.29063    0.09388   3.096
# ot1:DELE_z                      0.31130    0.09955   3.127
# ot1:use_z                      -0.07203    0.09606  -0.750
# ot2:use_z                       0.23924    0.09255   2.585

confint(gca_en_mod_final)
#                         2.5 %      97.5 %
# .sig01             0.28676190  0.47995627
# .sig02            -0.37741052  0.33776464
# .sig03             0.29894485  0.49287118
# .sig04             0.27597175  0.51421434
# .sigma             3.09545607  3.18656155
# (Intercept)        0.06971405  0.07989671
# ot1                1.46793478  1.48263542
# ot2                0.66769266  0.67920699
# condition_sum      0.07392547  0.08361477
# DELE_z             0.02166521  0.02854774
# use_z             -0.07751982 -0.07067952
# ot1:condition_sum  0.29002537  0.30114593
# ot1:DELE_z         0.31066796  0.32241908
# ot1:use_z         -0.07264262 -0.06126149
# ot2:use_z          0.23863895  0.24962341

# save Correlation matrix
vc <- vcov(gca_en_mod_final)

# diagonal matrix of standard deviations associated with vcov
S <- sqrt(diag(diag(vc), nrow(vc), nrow(vc)))

# convert vc to a correlation matrix
solve(S) %*% vc %*% solve(S) 
#              [,5]          [,6]         [,8]         [,9]         [,10]
#[5,]   1.000000000  
#[6,]  -0.306193771   1.000000000  
#[8,]  -0.068555808   0.010357519  1.000000000 
#[9,]   0.011925281  -0.075898402 -0.272056147  1.000000000 
#[10,]  0.016066384   0.026912433 -0.046288272 -0.027172118  1.0000000000   


# detect multicollinearity in a regression model: 
# if Variance inflation factor (VIF) score high, problematic (e.g., VIF = 9, usually >5)
# the VIF for a regression model variable is equal to the ratio 
# of the overall model variance to the variance of a model 
# that includes only that single independent variable.
car::vif(gca_en_mod_final)
#        ot1               ot2     condition_sum 
#   1.007652          1.006508          1.001696 
#     DELE_z             use_z ot1:condition_sum 
#   1.110559          1.112600          1.004033 
# ot1:DELE_z         ot1:use_z         ot2:use_z 
#   1.091605          1.094047          1.008253  


eigen(vc)
# An eigenvector is the direction of the line (vertical, horizontal, 45 degrees etc.).
# An eigenvalue is a number, telling you how much variance there is in the data in that direction, 
# how spread out the data is on the line. 
# The eigenvector with the highest eigenvalue is therefore the principal component.
# In a regression analysis, if there is one or more small eigenvalues of the X'X, matrix
# it implies there are near-linear relationships (dependencies) among regressors 
# and we may have multicollinearity issues in the analysis.

#$values
#[1] 0.012273967 0.009552314 0.009064854 0.008705354
#[5] 0.008117338 0.007422325 0.006881458 0.006660462
#[9] 0.004416302 0.002333985

#$vectors
#[,1]         [,2]         [,3]         [,4]
#[1,] -0.02394634  0.144338070 -0.030017224 -0.134834372
#[2,] -0.03925290 -0.576326522  0.412484137  0.583373237
#[3,] -0.07183820  0.651580529  0.383545705  0.077576249
#[4,]  0.02363884 -0.026767104 -0.096745152  0.064717365
#[5,]  0.04329073 -0.007961807  0.016046556  0.006914341
#[6,] -0.04065422 -0.016096576 -0.002634913 -0.038274257
#[7,] -0.07811259 -0.023195315  0.805048206 -0.380884208
#[8,] -0.74438067  0.104899854 -0.145171493  0.001271710
#[9,]  0.65265895  0.191418328  0.003174776  0.042289458
#[10,]  0.04971474 -0.416215244 -0.055108808 -0.694919725
#[,5]          [,6]         [,7]        [,8]
#[1,]  0.006850713  0.9611005513  0.183638999 -0.04320999
#[2,]  0.263722889  0.1487424458  0.100021265 -0.23056633
#[3,]  0.637604448 -0.0843793339  0.003062000 -0.05564865
#4,]  0.097728170 -0.1350617708  0.866152864  0.45505301
#[5,]  0.032819217 -0.0007609131 -0.008563023  0.02683591
#[6,]  0.009366044 -0.0059557034 -0.034499588  0.05763354
#[7,] -0.421167550 -0.0417260175  0.115150047  0.08714174
#[8,] -0.165988802 -0.1172541742  0.288916822 -0.53437248
#[9,] -0.198068319 -0.0941783555  0.316747905 -0.61714557
#[10,]  0.518635525 -0.0663675169  0.087846793 -0.23820992
#[,9]        [,10]
#[1,]  0.002594259 -0.012386888
#[2,]  0.030256669 -0.009210050
#[3,]  0.018653492  0.006580297
#[4,]  0.009180032  0.002446832
#[5,] -0.710416827 -0.700852326
#[6,]  0.699630522 -0.708871192
#[7,] -0.013235904  0.005782501
#[8,] -0.032081975 -0.049599441
#[9,]  0.056705219 -0.055523421
# [10,] -0.007121701  0.021064448



performance::r2_nakagawa(gca_en_mod_final, by_group = FALSE, tolerance = 1e-05)
# Conditional R2: 0.071
# Marginal R2: 0.031

performance::check_collinearity(gca_en_mod_final)
# Low Correlation
# 
# Term               VIF Increased SE Tolerance
# ot1               1.01         1.00      0.99
# ot2               1.01         1.00      0.99
# condition_sum     1.00         1.00      1.00
# DELE_z            1.11         1.05      0.90
# use_z             1.11         1.05      0.90
# ot1:condition_sum 1.00         1.00      1.00
# ot1:DELE_z        1.09         1.04      0.92
# ot1:use_z         1.09         1.05      0.91
# ot2:use_z         1.01         1.00      0.99

coll_info <- performance::check_collinearity(gca_en_mod_final)

collinearity <- plot(coll_info) 


ggsave('collinearity_simple',
       plot = collinearity, dpi = 600, device = "png",
       path = here("figs", "use_prof", "gca_verb"),
       height = 5, width = 28, units = 'in')

rsq::rsq(gca_en_mod_final)
# $model
#[1] 0.1043436
#$fixed
#[1] 0.04770466
#$random
#[1] 0.05663898
rsq::rsq(gca_en_mod_final, adj = TRUE)
#$model
#[1] 0.1034765

#$fixed
#[1] 0.04678269

#$random
#[1] 0.05669381


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
# -1.77325 -0.00209 -0.00006  0.00219  2.36125
# 
# Coefficients:
#              Estimate Std. Error  t value Pr(>|t|)    
# (Intercept) 0.0001558  0.0007899    0.197    0.844 
# res[-1]     0.9968069  0.0008306 1200.146   <2e-16 ***
# ---
# Signif. codes:  
# 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07619 on 9303 degrees of freedom
# Multiple R-squared:  0.9936,	Adjusted R-squared:  0.9936 
# F-statistic: 1.44e+06 on 1 and 9303 DF,  p-value: < 2.2e-16


mod_type <- "gca_en_mod"
mod_spec <- c('_base', 
              "_stress_0", "_stress_1", "_stress_2", 
              "_dele_0", "_dele_1", "_dele_2", 
              "_use_0", "_use_1", "_use_2",
              "_int_0", "_int_1", "_int_2", 
              '_final') 

# Store ind models in list
gca_en_mods_simple <- mget(c(paste0(mod_type, mod_spec)))

save(gca_en_mods_simple,
     file = here("mods", "use_prof", "gca",
                 "gca_en_mods_simple.Rdata"))




# Save anova model comparisons
en_nested_model_comparisons_simple <-
  mget(c(
         'en_stress_anova', 
         'en_dele_anova', 
         'en_use_anova', 'en_int_anova'
  ))

save(en_nested_model_comparisons_simple,
     file = here("mods", "use_prof", "gca",
                 "nested_model_comparisons_en_simple.Rdata"))


}

# -----------------------------------------------------------------------------








# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_dat <- mon_data %>%
  dplyr::select(time_zero, ot1:ot2, condition_sum) %>%
  distinct 

en_dat <- l2_data %>%
  dplyr::select(time_zero, ot1:ot2, condition_sum) %>%
  distinct %>%
  expand_grid(., tibble(DELE_z = c(-1, 0, 1))) %>%
  expand_grid(., tibble(use_z = c(-1, 0, 1)))


# Get model predictions and SE
fits_all_mon <- predictSE(gca_mod_mon_final, mon_dat) %>%
  as_tibble %>%
  bind_cols(mon_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_all_en <- predictSE(gca_en_mod_final, en_dat) %>%
  as_tibble %>%
  bind_cols(en_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at target syllable offset
target_offset_preds_mon <- filter(fits_all_mon, time_zero == 4) %>% 
  select(Stress = condition_sum,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))

target_offset_preds_en <- filter(fits_all_en, time_zero == 4) %>% 
  select(Stress = condition_sum, DELE = DELE_z, 
         `Spanish use` = use_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))


# Save models predictions
model_preds_simple <- mget(c(
  "fits_all_en", "fits_all_mon",
  "target_offset_preds_en", "target_offset_preds_mon"))

save(model_preds_simple,
     file = here("mods", "use_prof", "gca",
                 "model_preds_simple.Rdata"))
