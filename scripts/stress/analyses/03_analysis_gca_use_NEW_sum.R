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
  filter(., time_zero >= -4 & time_zero <= 12) %>%
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
  #          npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)    
  # mod_ot3    20 44058 44198 -22009    44018                          
  # mod_ot0    21 43919 44066 -21938    43877 141.042  1    < 2e-16 ***
  # mod_ot1a   23 43705 43866 -21830    43659 217.748  2    < 2e-16 ***
  # mod_ot2a   26 43621 43803 -21785    43569  89.723  3    < 2e-16 ***
  # mod_ot3a   30 43619 43829 -21779    43559  10.426  4    0.03384 *  
  
}



# Individual model MON -----------------------------------------------------------

gca_mod_mon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
       REML = F,
       data = filter(mon_data)) 

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + condition_sum) 
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0,   . ~ . + ot1:condition_sum) 
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1,   . ~ . + ot2:condition_sum) 
gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2,   . ~ . + ot3:condition_sum) 

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2, gca_mod_mon_cond_3)
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base     30 42753 42963 -21346    42693                     
# gca_mod_mon_cond_0   31 42755 42972 -21346    42693 0.4824  1     0.4873
# gca_mod_mon_cond_1   32 42756 42980 -21346    42692 0.1857  1     0.6665
# gca_mod_mon_cond_2   33 42758 42989 -21346    42692 0.0000  1     0.9998
# gca_mod_mon_cond_3   34 42760 42998 -21346    42692 0.0628  1     0.8021


# add WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_wm_0 <- update(gca_mod_mon_base, . ~ . + ospan) 
gca_mod_mon_wm_1 <- update(gca_mod_mon_wm_0,   . ~ . + ot1:ospan) 
gca_mod_mon_wm_2 <- update(gca_mod_mon_wm_1,   . ~ . + ot2:ospan) 
gca_mod_mon_wm_3 <- update(gca_mod_mon_wm_2,   . ~ . + ot3:ospan) 

mon_wm_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_wm_0, gca_mod_mon_wm_1,
        gca_mod_mon_wm_2, gca_mod_mon_wm_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base   30 42753 42963 -21346    42693                     
# gca_mod_mon_wm_0   31 42755 42972 -21346    42693 0.4829  1     0.4871
# gca_mod_mon_wm_1   32 42756 42980 -21346    42692 0.4634  1     0.4960
# gca_mod_mon_wm_2   33 42758 42989 -21346    42692 0.3297  1     0.5658
# gca_mod_mon_wm_3   34 42759 42997 -21346    42691 0.4147  1     0.5196


# add condition x WM int to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_int_0 <- update(gca_mod_mon_base, . ~ . + condition_sum:ospan) 
gca_mod_mon_int_1 <- update(gca_mod_mon_int_0,   . ~ . + ot1:condition_sum:ospan) 
gca_mod_mon_int_2 <- update(gca_mod_mon_int_1,   . ~ . + ot2:condition_sum:ospan) 
gca_mod_mon_int_3 <- update(gca_mod_mon_int_2,   . ~ . + ot3:condition_sum:ospan) 

mon_int_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_int_0, gca_mod_mon_int_1,
        gca_mod_mon_int_2, gca_mod_mon_int_3)
#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_mon_base  30   42753 42963 -21346    42693                       
# gca_mod_mon_int_0   31 42749 42966 -21344    42687 5.6757  1     0.0172 *
# gca_mod_mon_int_1   32 42749 42973 -21343    42685 2.1568  1     0.1419  
# gca_mod_mon_int_2   33 42749 42980 -21341    42683 2.4539  1     0.1172  
# gca_mod_mon_int_3   34 42750 42989 -21341    42682 0.3066  1     0.5798 

gca_mod_mon_full <- gca_mod_mon_int_0


mod_type <- "gca_mod_mon"
mod_spec <- c('_base', 
              "_cond_0", "_cond_1", "_cond_2", "_cond_3", 
              "_wm_0", "_wm_1", "_wm_2", "_wm_3",
              "_int_0", "_int_1", "_int_2", "_int_3",
              '_full')

# Store ind models in list
gca_mon_mods <- mget(c(paste0(mod_type, mod_spec)))

save(gca_mon_mods,
     file = here("mods", "stress", "gca",
                 "gca_mon_mods.Rdata"))





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
gca_l2_mod_base <- 
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +         
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target), 
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 2e4)),
       data = l2_data, REML = F)

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_l2_mod_l1_0 <- update(gca_l2_mod_base, . ~ . + l1_sum) 
gca_l2_mod_l1_1 <- update(gca_l2_mod_l1_0, . ~ . + ot1:l1_sum) 
gca_l2_mod_l1_2 <- update(gca_l2_mod_l1_1, . ~ . + ot2:l1_sum) 
gca_l2_mod_l1_3 <- update(gca_l2_mod_l1_2, . ~ . + ot3:l1_sum)

l2_l1_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_l1_0, gca_l2_mod_l1_1,
        gca_l2_mod_l1_2, gca_l2_mod_l1_3)
#                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_l2_mod_base   30 183369 183622 -91654   183309                       
# gca_l2_mod_l1_0   31 183368 183629 -91653   183306 2.7975  1    0.09441 .
# gca_l2_mod_l1_1   32 183368 183638 -91652   183304 1.6614  1    0.19742  
# gca_l2_mod_l1_2   33 183369 183647 -91652   183303 0.9824  1    0.32161  
# gca_l2_mod_l1_3   34 183370 183656 -91651   183302 1.8417  1    0.17475 


# add stress effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_stress_0 <- update(gca_l2_mod_base,   . ~ . + condition_sum) 
gca_l2_mod_stress_1 <- update(gca_l2_mod_stress_0, . ~ . + ot1:condition_sum) 
gca_l2_mod_stress_2 <- update(gca_l2_mod_stress_1, . ~ . + ot2:condition_sum) 
gca_l2_mod_stress_3 <- update(gca_l2_mod_stress_2, . ~ . + ot3:condition_sum)

l2_stress_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_stress_0, gca_l2_mod_stress_1,
        gca_l2_mod_stress_2, gca_l2_mod_stress_3)
#                     npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)    
# gca_l2_mod_base       30 183369 183622 -91654   183309                       
# gca_l2_mod_stress_0   31 183369 183630 -91654   183307 1.5751  1    0.20947  
# gca_l2_mod_stress_1   32 183371 183641 -91653   183307 0.3685  1    0.54383  
# gca_l2_mod_stress_2   33 183371 183649 -91653   183305 1.5830  1    0.20832  
# gca_l2_mod_stress_3   34 183370 183657 -91651   183302 3.1447  1    0.07618 .


### ROUTE 1 (proficiency alone)

# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_dele_0 <- update(gca_l2_mod_base,   . ~ . + DELE_z) 
gca_l2_mod_dele_1 <- update(gca_l2_mod_dele_0, . ~ . + ot1:DELE_z) 
gca_l2_mod_dele_2 <- update(gca_l2_mod_dele_1, . ~ . + ot2:DELE_z)
gca_l2_mod_dele_3 <- update(gca_l2_mod_dele_2, . ~ . + ot3:DELE_z)

l2_dele_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_dele_0, gca_l2_mod_dele_1,
        gca_l2_mod_dele_2, gca_l2_mod_dele_3)
#                     npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_l2_mod_base     30 183369 183622 -91654   183309                        
# gca_l2_mod_dele_0   31 183363 183624 -91651   183301 7.7294  1   0.005433 **
# gca_l2_mod_dele_1   32 183363 183632 -91649   183299 2.3905  1   0.122075   
# gca_l2_mod_dele_2   33 183363 183642 -91649   183297 1.2852  1   0.256940   
# gca_l2_mod_dele_3   34 183365 183652 -91649   183297 0.1591  1   0.689971   

# add interaction

gca_l2_mod_dele_int_0 <- update(gca_l2_mod_dele_0,   . ~ . + l1_sum:condition_sum:DELE_z) 
gca_l2_mod_dele_int_1 <- update(gca_l2_mod_dele_int_0, . ~ . + ot1:l1_sum:condition_sum:DELE_z) 
gca_l2_mod_dele_int_2 <- update(gca_l2_mod_dele_int_1, . ~ . + ot2:l1_sum:condition_sum:DELE_z)
gca_l2_mod_dele_int_3 <- update(gca_l2_mod_dele_int_2, . ~ . + ot3:l1_sum:condition_sum:DELE_z)

l2_dele_int_anova <-
  anova(gca_l2_mod_dele_0, gca_l2_mod_dele_int_0, gca_l2_mod_dele_int_1,
        gca_l2_mod_dele_int_2, gca_l2_mod_dele_int_3)
#                       npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_l2_mod_dele_0       31 183363 183624 -91651   183301                        
# gca_l2_mod_dele_int_0   32 183362 183632 -91649   183298 3.1518  1   0.075845 . 
# gca_l2_mod_dele_int_1   33 183364 183642 -91649   183298 0.4593  1   0.497963   
# gca_l2_mod_dele_int_2   34 183357 183644 -91645   183289 8.1576  1   0.004288 **
# gca_l2_mod_dele_int_3   35 183359 183654 -91645   183289 0.0099  1   0.920801

gca_l2_mod_dele_final <- gca_l2_mod_dele_int_2

summary(gca_l2_mod_dele_final)
#                                 Estimate Std. Error t value
# (Intercept)                      1.66874    0.10935  15.260
# ot1                              6.66342    0.32474  20.520
# ot2                             -0.88786    0.23982  -3.702
# ot3                             -1.69238    0.17011  -9.949
# DELE_z                           0.23842    0.08395   2.840
# DELE_z:l1_sum:condition_sum      0.09192    0.04987   1.843
# ot1:DELE_z:l1_sum:condition_sum -0.06753    0.09502  -0.711
# ot2:DELE_z:l1_sum:condition_sum -0.27183    0.09516  -2.856

### ROUTE 2 (L2 use alone)

# add L2 use effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_use_0 <- update(gca_l2_mod_base,  . ~ . + use_z) 
gca_l2_mod_use_1 <- update(gca_l2_mod_use_0, . ~ . + ot1:use_z) 
gca_l2_mod_use_2 <- update(gca_l2_mod_use_1, . ~ . + ot2:use_z)
gca_l2_mod_use_3 <- update(gca_l2_mod_use_2, . ~ . + ot3:use_z)

l2_use_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_use_0, gca_l2_mod_use_1,
        gca_l2_mod_use_2, gca_l2_mod_use_3)
#                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_l2_mod_base    30 183369 183622 -91654   183309                        
# gca_l2_mod_use_0   31 183371 183632 -91654   183309 0.1832  1   0.668622   
# gca_l2_mod_use_1   32 183371 183640 -91653   183307 2.0781  1   0.149423   
# gca_l2_mod_use_2   33 183366 183644 -91650   183300 6.6795  1   0.009753 **
# gca_l2_mod_use_3   34 183365 183651 -91648   183297 2.9724  1   0.084698 . 

# add interaction

gca_l2_mod_use_int_0 <- update(gca_l2_mod_use_2,     . ~ . + l1_sum:condition_sum:use_z) 
gca_l2_mod_use_int_1 <- update(gca_l2_mod_use_int_0, . ~ . + ot1:l1_sum:condition_sum:use_z) 
gca_l2_mod_use_int_2 <- update(gca_l2_mod_use_int_1, . ~ . + ot2:l1_sum:condition_sum:use_z)
gca_l2_mod_use_int_3 <- update(gca_l2_mod_use_int_2, . ~ . + ot3:l1_sum:condition_sum:use_z)

l2_use_int_anova <-
  anova(gca_l2_mod_use_2, gca_l2_mod_use_int_0, gca_l2_mod_use_int_1,
        gca_l2_mod_use_int_2, gca_l2_mod_use_int_3)
#                      npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)   
# gca_l2_mod_use_2       33 183366 183644 -91650   183300                       
# gca_l2_mod_use_int_0   34 183367 183654 -91650   183299 0.6642  1    0.41508  
# gca_l2_mod_use_int_1   35 183366 183661 -91648   183296 3.5456  1    0.05970 .
# gca_l2_mod_use_int_2   36 183364 183668 -91646   183292 3.4130  1    0.06468 .
# gca_l2_mod_use_int_3   37 183365 183677 -91646   183291 1.0641  1    0.30227  

gca_l2_mod_use_final <- gca_l2_mod_use_2

summary(gca_l2_mod_use_final)
#             Estimate Std. Error t value
# (Intercept)  1.60890    0.11266  14.281
# ot1          6.84658    0.33245  20.594
# ot2         -1.05524    0.24614  -4.287
# ot3         -1.69329    0.17003  -9.959
# use_z        0.04044    0.09958   0.406
# ot1:use_z    0.56331    0.26500   2.126
# ot2:use_z   -0.51667    0.19489  -2.651


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
  mget(c(#"mon_cond_anova", "mon_wm_anova", "mon_int_anova",
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





##            NOT MODIFIED


# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_mon <- mon_data %>%
  dplyr::select(time_zero, ot1:ot3, condition_sum) %>% #, ospan
  distinct

dele_dat_l2 <- l2_data %>%
  dplyr::select(l1, time_zero, ot1:ot3, condition_sum) %>%
  distinct %>%
  mutate(l1 = as.character(l1)) %>% 
  expand_grid(., tibble(DELE_z = c(-1, 0, 1)))

use_dat_l2 <- l2_data %>%
  dplyr::select(l1, time_zero, ot1:ot3, condition_sum) %>%
  distinct %>%
  mutate(l1 = as.character(l1)) %>% 
  expand_grid(., tibble(use_z = c(-1, 0, 1)))
  

# write_csv(dele_dat_l2, here::here('dele_dat_l2.csv'))
# dele_dat_l2 <- read_csv(here::here('dele_dat_l2.csv'))

# Get model predictions and SE
fits_all_mon <- predictSE(gca_mod_mon_base, new_dat_mon) %>%  
  as_tibble %>%
  bind_cols(new_dat_mon) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_all_l2_dele <- predictSE(gca_l2_mod_dele_final, dele_dat_l2) %>%
  as_tibble %>%
  bind_cols(dele_dat_l2) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         l1 = fct_recode(l1, EN = "en", MA = "ma"),
         l1 = fct_relevel(l1, "EN", "MA"))

fits_all_l2_use <- predictSE(gca_l2_mod_use_final, use_dat_l2) %>%  
  as_tibble %>%
  bind_cols(use_dat_l2) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         l1 = fct_recode(l1, EN = "en", MA = "ma"),
         l1 = fct_relevel(l1, "EN", "MA"))


# Filter preds at target syllable offset
target_offset_preds_mon <- filter(fits_all_mon, time_zero == 4) %>%
  select(stress = condition_sum, #ospan,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 


target_offset_preds_l2_dele <- filter(fits_all_l2_dele, time_zero == 4) %>% 
  select(l1, stress = condition_sum, DELE = DELE_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(l1)


target_offset_preds_l2_use <- filter(fits_all_l2_use, time_zero == 4) %>%
  select(l1, stress = condition_sum, percent_l2_week = use_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(l1)





# Save models predictions
model_preds <- mget(c(#"fits_all_mon", "fits_all_l2_dele", "fits_all_l2_use", "fits_all_l2_wm",
                      #"target_offset_preds_mon", "target_offset_preds_l2_wm",
                      "target_offset_preds_l2_dele",
                      "target_offset_preds_l2_use"))

save(model_preds,
     file = here("mods", "stress", "gca",
                 "model_preds_l2.Rdata"))

# MON without WM score
saveRDS(target_offset_preds_mon, file = here("mods", "stress", "gca", "model_preds_mon.Rdata"))

# -----------------------------------------------------------------------------


dele_dat_l2 <- l2_data %>%
  dplyr::select(l1, time_zero, ot1:ot3, condition_sum, ospan, DELE_z, use_z) %>%
  distinct

dele_dat_l2 <- dele_dat_l2[ which(dele_dat_l2$l1!='es'), ]
write_csv(dele_dat_l2, here::here('dele_dat_l2.csv'))
dele_dat_l2 <- read_csv(here::here('dele_dat_l2.csv'))

dele_dat_l2 <- expand.grid(
  l1 = c('en', 'ma'), 
  time_zero = 4,
  ot1 = unique(l2_data$ot1),
  ot2 = unique(l2_data$ot2),
  ot3 = unique(l2_data$ot3),
  condition_sum = unique(l2_data$condition_sum),
  DELE_z = c(-1, 0, 1))

fits_all_l2_dele <- predictSE(gca_l2_mod_dele_final, dele_dat_l2) %>%
  as_tibble %>%
  bind_cols(dele_dat_l2) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         l1 = fct_recode(l1, EN = "en", MA = "ma"),
         l1 = fct_relevel(l1, "EN", "MA"))


# Filter preds at target syllable offset
target_offset_preds_dele <- filter(fits_all_l2_dele, time_zero == 4) %>% #
  select(l1, stress = condition_sum, DELE = DELE_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub),
         DELE = c(-1, 0, 1))
  
  
  
  border_1 <- fp_border(width = 1.5)
  border_2 <- fp_border(width = 0.75)
  
  #mean(model_preds$target_offset_preds_l2_dele$DELE) - sd(model_preds$target_offset_preds_l2_dele$DELE)
  
  target_offset_preds_dele %>% 
    mutate(stress = if_else(stress == 1, "Present", "Preterit"),
           l1 = if_else(l1 == 'EN', "English", "Mandarin"),
           l1 = fct_relevel(l1, "English", "Mandarin")) %>% 
    arrange(l1, stress) %>% 
    group_by(l1, stress) %>%
    distinct() %>%
    #ungroup() %>%
    mutate(l1 = blank_same_as_last(as.character(l1)),
           stress = blank_same_as_last(as.character(stress))) %>%
    select(L1 = l1, `Lexical stress` = stress, Proficiency = DELE, Probability = prob,
           LB = prob_lb, UB = prob_ub) %>%
    flextable() %>% 
    width(., j = c(2, 3, 4), width = c(1.1, 1.3, 1.1)) %>%
    font(., fontname = "Times", part = "all") %>%
    fontsize(., size = 11) %>% 
    border_remove(.) %>%  
    border(., part = "header", 
           border.top = border_1,
           border.bottom = border_2) %>% 
    hline_bottom(., part = "body", border = border_1)
  