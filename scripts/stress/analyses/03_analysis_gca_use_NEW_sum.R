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
load(paste0(gca_mods_path, "/gca_mon_mods.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_mods.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_mods_sum.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))
# load(paste0(gca_mods_path, "/model_preds_mon.Rdata"))

# Store objects in global env
list2env(gca_mon_mods, globalenv())
list2env(gca_l2_mods, globalenv())
list2env(gca_l2_mods_sum, globalenv())
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



# stress_50 <- na.omit(stress50)


stress_gc_subset <- stress50 %>%
  # select(., -WM_set) %>%
  filter(., time_zero >= -4 & time_zero <= 6) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma"),
            condition_sum = if_else(cond == "1", 1, -1)) %>%       # 1 = present, 2 = past
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")



# -----------------------------------------------------------------------------







#################### MONOLINGUAL SPEAKERS ########################################

# Build up random effects to test time terms
if(F){
  
  mon_data <- filter(stress_gc_subset, l1 == 'es') %>% select(-DELE, -percent_l2_week)
  
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


mod_type <- "gca_mod_mon"
mod_spec <- c('_base', 
              "_cond_0", "_cond_1", "_cond_2", "_cond_3",
              '_final')
              # , 
              # "_wm_0", "_wm_1", "_wm_2", "_wm_3",
              # "_int_0", "_int_1", "_int_2", "_int_3",
              # '_full')

# Store ind models in list
gca_mon_mods <- mget(c(paste0(mod_type, mod_spec)))

save(gca_mon_mods,
     file = here("mods", "stress", "gca",
                 "gca_mon_mods_sum.Rdata"))

# # add WM effect to intercept, linear slope, quadratic, and cubic time terms
# gca_mod_mon_wm_0 <- update(gca_mod_mon_base, . ~ . + ospan) 
# gca_mod_mon_wm_1 <- update(gca_mod_mon_wm_0,   . ~ . + ot1:ospan) 
# gca_mod_mon_wm_2 <- update(gca_mod_mon_wm_1,   . ~ . + ot2:ospan) 
# gca_mod_mon_wm_3 <- update(gca_mod_mon_wm_2,   . ~ . + ot3:ospan) 
# 
# mon_wm_anova <-
#   anova(gca_mod_mon_base, gca_mod_mon_wm_0, gca_mod_mon_wm_1,
#         gca_mod_mon_wm_2, gca_mod_mon_wm_3)
# #                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# # gca_mod_mon_base   
# 
# 
# # add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
# gca_mod_mon_int_0 <- update(gca_mod_mon_base, . ~ . + condition_sum:ospan) 
# gca_mod_mon_int_1 <- update(gca_mod_mon_int_0,   . ~ . + ot1:condition_sum:ospan) 
# gca_mod_mon_int_2 <- update(gca_mod_mon_int_1,   . ~ . + ot2:condition_sum:ospan) 
# gca_mod_mon_int_3 <- update(gca_mod_mon_int_2,   . ~ . + ot3:condition_sum:ospan) 
# 
# mon_int_anova <-
#   anova(gca_mod_mon_base, gca_mod_mon_int_0, gca_mod_mon_int_1,
#         gca_mod_mon_int_2, gca_mod_mon_int_3)
# #                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# # gca_mod_mon_base  







#################### L2 SPEAKERS ########################################

l2_data <- stress_gc_subset%>%
  filter(., l1 != 'es') %>% 
  filter(., participant != 'ies04' & participant != 'ies17' & participant != 'ies28' & participant != 'aes32') %>%
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1))

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
  
  #           npar(Df?)    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1         9 185950 186026 -92966   185932                          
  # mod_ot2   14 185872 185990 -92922   185844  87.736  5  < 2.2e-16 ***
  # mod_ot3   20 185459 185628 -92710   185419 424.968  6  < 2.2e-16 ***  
  
  
  
  mod_ot0 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + ot1 | target) +
                       + (1 + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + ot1 + ot2 | target) +
                       + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot3    20 185459 185628 -92710   185419                          
  # mod_ot0    21 184996 185173 -92477   184954 465.356  1  < 2.2e-16 ***
  # mod_ot1a   23 184615 184809 -92284   184569 385.320  2  < 2.2e-16 ***
  # mod_ot2a   26 184485 184704 -92217   184433 135.484  3  < 2.2e-16 ***
  # mod_ot3a   30 184435 184688 -92188   184375  57.797  4  8.418e-12 ***
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_l2_mod_base <- mod_ot3a
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +         
  #        (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 + ot3 | target), 
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 2e4)),
  #      data = l2_data, REML = F)

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_l2_mod_l1_0 <- update(gca_l2_mod_base, . ~ . + l1_sum) 
gca_l2_mod_l1_1 <- update(gca_l2_mod_l1_0, . ~ . + ot1:l1_sum) 
gca_l2_mod_l1_2 <- update(gca_l2_mod_l1_1, . ~ . + ot2:l1_sum) 
gca_l2_mod_l1_3 <- update(gca_l2_mod_l1_2, . ~ . + ot3:l1_sum)

l2_l1_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_l1_0, gca_l2_mod_l1_1,
        gca_l2_mod_l1_2, gca_l2_mod_l1_3)
#                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_l2_mod_base   30 184435 184688 -92188   184375                        
# gca_l2_mod_l1_0   31 184437 184698 -92188   184375 0.2224  1   0.637211   
# gca_l2_mod_l1_1   32 184431 184701 -92183   184367 8.3263  1   0.003908 **
# gca_l2_mod_l1_2   33 184433 184711 -92183   184367 0.3552  1   0.551194   
# gca_l2_mod_l1_3   34 184431 184718 -92182   184363 3.1126  1   0.077688 . 


# add stress effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_stress_0 <- update(gca_l2_mod_l1_1,   . ~ . + condition_sum) 
gca_l2_mod_stress_1 <- update(gca_l2_mod_stress_0, . ~ . + ot1:condition_sum) 
gca_l2_mod_stress_2 <- update(gca_l2_mod_stress_1, . ~ . + ot2:condition_sum) 
gca_l2_mod_stress_3 <- update(gca_l2_mod_stress_2, . ~ . + ot3:condition_sum)

l2_stress_anova <-
  anova(gca_l2_mod_l1_1, gca_l2_mod_stress_0, gca_l2_mod_stress_1,
        gca_l2_mod_stress_2, gca_l2_mod_stress_3)
#                     npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)    
# gca_l2_mod_l1_1       32 184431 184701 -92183   184367                       
# gca_l2_mod_stress_0   33 184429 184707 -92182   184363 3.7177  1    0.05384 .
# gca_l2_mod_stress_1   34 184431 184717 -92181   184363 0.3422  1    0.55855  
# gca_l2_mod_stress_2   35 184427 184722 -92179   184357 5.3341  1    0.02091 *
# gca_l2_mod_stress_3   36 184429 184732 -92178   184357 0.8045  1    0.36974  


### ROUTE 1 (proficiency alone)

# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_dele_0 <- update(gca_l2_mod_stress_2,   . ~ . + DELE_z) 
gca_l2_mod_dele_1 <- update(gca_l2_mod_dele_0, . ~ . + ot1:DELE_z) 
gca_l2_mod_dele_2 <- update(gca_l2_mod_dele_1, . ~ . + ot2:DELE_z)
gca_l2_mod_dele_3 <- update(gca_l2_mod_dele_2, . ~ . + ot3:DELE_z)

l2_dele_anova <-
  anova(gca_l2_mod_stress_2, gca_l2_mod_dele_0, gca_l2_mod_dele_1,
        gca_l2_mod_dele_2, gca_l2_mod_dele_3)
#                     npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_l2_mod_stress_2   35 184427 184722 -92179   184357                        
# gca_l2_mod_dele_0     36 184424 184728 -92176   184352 5.0895  1   0.024071 * 
# gca_l2_mod_dele_1     37 184418 184730 -92172   184344 8.2057  1   0.004176 **
# gca_l2_mod_dele_2     38 184419 184739 -92171   184343 1.3448  1   0.246196   
# gca_l2_mod_dele_3     39 184420 184749 -92171   184342 0.8998  1   0.342824   

# add interaction

gca_l2_mod_dele_int_0 <- update(gca_l2_mod_dele_1,   . ~ . + l1_sum:condition_sum:DELE_z) 
gca_l2_mod_dele_int_1 <- update(gca_l2_mod_dele_int_0, . ~ . + ot1:l1_sum:condition_sum:DELE_z) 
gca_l2_mod_dele_int_2 <- update(gca_l2_mod_dele_int_1, . ~ . + ot2:l1_sum:condition_sum:DELE_z)
gca_l2_mod_dele_int_3 <- update(gca_l2_mod_dele_int_2, . ~ . + ot3:l1_sum:condition_sum:DELE_z)

l2_dele_int_anova <-
  anova(gca_l2_mod_dele_1, gca_l2_mod_dele_int_0, gca_l2_mod_dele_int_1,
        gca_l2_mod_dele_int_2, gca_l2_mod_dele_int_3)
#                       npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_l2_mod_dele_1       37 184418 184730 -92172   184344                     
# gca_l2_mod_dele_int_0   38 184419 184740 -92172   184343 0.7050  1     0.4011
# gca_l2_mod_dele_int_1   39 184419 184748 -92171   184341 2.0175  1     0.1555
# gca_l2_mod_dele_int_2   40 184420 184757 -92170   184340 1.8167  1     0.1777
# gca_l2_mod_dele_int_3   41 184421 184767 -92170   184339 0.5044  1     0.4776

gca_l2_mod_dele_final <- gca_l2_mod_dele_1

summary(gca_l2_mod_dele_final)
#                   Estimate Std. Error t value
# (Intercept)        1.16866    0.07939  14.720
# ot1                5.00528    0.26127  19.157
# ot2                0.25565    0.20006   1.278
# ot3               -1.15544    0.15603  -7.405
# l1_sum            -0.09161    0.04484  -2.043
# condition_sum     -0.07411    0.05540  -1.338
# DELE_z             0.19586    0.05646   3.469
# ot1:l1_sum        -0.47295    0.15154  -3.121
# ot1:condition_sum  0.13821    0.17255   0.801
# ot2:condition_sum  0.30803    0.12363   2.492
# ot1:DELE_z         0.56237    0.19068   2.949



### ROUTE 2 (L2 use alone)

# add L2 use effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_use_0 <- update(gca_l2_mod_stress_2,  . ~ . + use_z) 
gca_l2_mod_use_1 <- update(gca_l2_mod_use_0, . ~ . + ot1:use_z) 
gca_l2_mod_use_2 <- update(gca_l2_mod_use_1, . ~ . + ot2:use_z)
gca_l2_mod_use_3 <- update(gca_l2_mod_use_2, . ~ . + ot3:use_z)

l2_use_anova <-
  anova(gca_l2_mod_stress_2, gca_l2_mod_use_0, gca_l2_mod_use_1,
        gca_l2_mod_use_2, gca_l2_mod_use_3)
#                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_l2_mod_stress_2   35 184427 184722 -92179   184357                       
# gca_l2_mod_use_0      36 184429 184732 -92178   184357 0.6689  1    0.41342  
# gca_l2_mod_use_1      37 184425 184737 -92175   184351 6.0747  1    0.01371 *
# gca_l2_mod_use_2      38 184424 184744 -92174   184348 2.4578  1    0.11694  
# gca_l2_mod_use_3      39 184424 184753 -92173   184346 1.9997  1    0.15733  

# add interaction

gca_l2_mod_use_int_0 <- update(gca_l2_mod_use_1,     . ~ . + l1_sum:condition_sum:use_z) 
gca_l2_mod_use_int_1 <- update(gca_l2_mod_use_int_0, . ~ . + ot1:l1_sum:condition_sum:use_z) 
gca_l2_mod_use_int_2 <- update(gca_l2_mod_use_int_1, . ~ . + ot2:l1_sum:condition_sum:use_z)
gca_l2_mod_use_int_3 <- update(gca_l2_mod_use_int_2, . ~ . + ot3:l1_sum:condition_sum:use_z)

l2_use_int_anova <-
  anova(gca_l2_mod_use_1, gca_l2_mod_use_int_0, gca_l2_mod_use_int_1,
        gca_l2_mod_use_int_2, gca_l2_mod_use_int_3)
#                      npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)   
# gca_l2_mod_use_1       37 184425 184737 -92175   184351                          
# gca_l2_mod_use_int_0   38 184426 184746 -92175   184350  0.5185  1  0.4714918    
# gca_l2_mod_use_int_1   39 184428 184756 -92175   184350  0.4172  1  0.5183491    
# gca_l2_mod_use_int_2   40 184428 184765 -92174   184348  2.1487  1  0.1426932    
# gca_l2_mod_use_int_3   41 184418 184763 -92168   184336 11.9486  1  0.0005469 ***

gca_l2_mod_use_final <- gca_l2_mod_use_int_3

summary(gca_l2_mod_use_final)
#             Estimate Std. Error t value
# (Intercept)  1.11998    0.08164  13.718
# ot1                             5.01338    0.26526  18.900
# ot2                             0.25831    0.19798   1.305
# ot3                            -1.15299    0.15484  -7.446
# l1_sum                         -0.09488    0.04762  -1.992
# condition_sum                  -0.07287    0.05539  -1.316
# use_z                           0.03420    0.06618   0.517
# ot1:l1_sum                     -0.53704    0.15554  -3.453
# ot1:condition_sum               0.12107    0.17310   0.699
# ot2:condition_sum               0.29629    0.12311   2.407
# ot1:use_z                       0.54530    0.21575   2.528
# l1_sum:condition_sum:use_z     -0.02249    0.03460  -0.650
# ot1:l1_sum:condition_sum:use_z  0.03811    0.08538   0.446
# ot2:l1_sum:condition_sum:use_z -0.15421    0.08298  -1.858
# ot3:l1_sum:condition_sum:use_z -0.29314    0.08466  -3.463


mod_type <- "gca_l2_mod"
mod_spec <- c('_base', 
              "_l1_0", "_l1_1", "_l1_2", "_l1_3", 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3", 
              "_dele_0", "_dele_1", "_dele_2", "_dele_3", 
              "_dele_int_0", "_dele_int_1", "_dele_int_2", "_dele_int_3",
              "_use_0", "_use_1", "_use_2", "_use_3", 
              "_use_int_0", "_use_int_1", "_use_int_2", "_use_int_3",
              '_dele_final', '_use_final') 

# Store ind models in list
gca_l2_mods <- mget(c(paste0(mod_type, mod_spec)))

save(gca_l2_mods,
     file = here("mods", "stress", "gca",
                 "gca_l2_mods_sum.Rdata"))


# Save anova model comparisons
nested_model_comparisons <-
  mget(c("mon_cond_anova", #"mon_wm_anova", "mon_int_anova",
         'l2_l1_anova', 'l2_stress_anova', 
         'l2_dele_anova', 'l2_dele_int_anova',
         'l2_use_anova', 'l2_use_int_anova'#,
         #'l2_wm_anova', 'l2_wm_int_anova'
  ))

save(nested_model_comparisons,
     file = here("mods", "stress", "gca",
                 "nested_model_comparisons_l2_sum.Rdata"))


}

# -----------------------------------------------------------------------------





##            MON NOT MODIFIED - L2 yes


# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_mon <- mon_data %>%
  dplyr::select(time_zero, ot1:ot3, condition_sum) %>% #, ospan
  distinct
  

# write_csv(dele_dat_l2, here::here('dele_dat_l2.csv'))
# dele_dat_l2 <- read_csv(here::here('dele_dat_l2.csv'))

# Get model predictions and SE
fits_all_mon <- predictSE(gca_mod_mon_final, new_dat_mon) %>%  
  as_tibble %>%
  bind_cols(new_dat_mon) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at target syllable offset
target_offset_preds_mon <- filter(fits_all_mon, time_zero == 4) %>%
  select(stress = condition_sum, #ospan,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))




# Save models predictions
model_preds <- mget(c(#"fits_all_mon", "fits_all_l2_dele", "fits_all_l2_use", "fits_all_l2_wm",
                      #"target_offset_preds_mon", "target_offset_preds_l2_wm",
                      "target_offset_preds_l2_dele",
                      "target_offset_preds_l2_use"))

save(model_preds,
     file = here("mods", "stress", "gca",
                 "model_preds_l2.Rdata"))

# -----------------------------------------------------------------------------


dele_dat_l2 <- l2_data %>%
  dplyr::select(l1_sum, time_zero, ot1:ot3, condition_sum) %>%
  distinct %>%
  # mutate(l1_sum = as.character(l1_sum)) %>% 
  expand_grid(., tibble(DELE_z = c(-1, 0, 1)))

use_dat_l2 <- l2_data %>%
  dplyr::select(l1_sum, time_zero, ot1:ot3, condition_sum) %>%
  distinct %>%
  # mutate(l1_sum = as.character(l1_sum)) %>% 
  expand_grid(., tibble(use_z = c(-1, 0, 1)))

# dele_dat_l2 <- dele_dat_l2[ which(dele_dat_l2$l1!='es'), ]
# write_csv(dele_dat_l2, here::here('dele_dat_l2.csv'))
# dele_dat_l2 <- read_csv(here::here('dele_dat_l2.csv'))


fits_all_l2_dele <- predictSE(gca_l2_mod_dele_final, dele_dat_l2) %>%
  as_tibble %>%
  bind_cols(dele_dat_l2) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)
# ,
#          l1 = fct_recode(l1, EN = "en", MA = "ma"),
#          l1 = fct_relevel(l1, "EN", "MA"))

fits_all_l2_use <- predictSE(gca_l2_mod_use_final, use_dat_l2) %>%
  as_tibble %>%
  bind_cols(use_dat_l2) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at target syllable offset
target_offset_preds_dele <- filter(fits_all_l2_dele, time_zero == 4) %>% #
  select(l1 = l1_sum, stress = condition_sum, DELE = DELE_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))

 
target_offset_preds_use <- filter(fits_all_l2_use, time_zero == 4) %>% #
  select(l1 = l1_sum, stress = condition_sum, use = use_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))

  

model_preds <- mget(c("fits_all_mon", 
  "fits_all_l2_dele", "fits_all_l2_use", #"fits_all_l2_wm",
  "target_offset_preds_mon",
  #"target_offset_preds_l2_wm",
  "target_offset_preds_dele",
  "target_offset_preds_use"))

save(model_preds,
     file = here("mods", "stress", "gca",
                 "model_preds_l2_sum.Rdata"))
