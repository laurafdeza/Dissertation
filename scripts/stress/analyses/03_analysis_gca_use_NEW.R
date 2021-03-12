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
load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))
# load(paste0(gca_mods_path, "/model_preds_mon.Rdata"))

# Store objects in global env
list2env(gca_mon_mods, globalenv())
list2env(gca_l2_mods, globalenv())
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

wm <- read_csv("./data/clean/ospan_set_z_scores.csv")

stress_50 <- merge(x = stress50, y = wm, by = "participant", all.x=TRUE)

# stress_50 <- na.omit(stress_50)



stress_gc_subset <- stress_50 %>%
  select(., -WM_set) %>%
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

l2_data <- filter(stress_gc_subset, l1 != 'es') %>% 
  mutate(., l1 = fct_relevel(l1, "en", "ma"))

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
  # mod_ot1           9 191650 191726 -95816   191632                          
  # mod_ot2          14 191568 191687 -95770   191540  91.526  5  < 2.2e-16 ***
  # mod_ot3          20 191151 191321 -95556   191111 428.905  6  < 2.2e-16 ***

  
  mod_ot0 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + ot1 | target) +
                       + (1 + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + ot1 + ot2 | target) +
                       + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot3    20 191151 191321 -95556   191111                          
  # mod_ot0    21 190712 190890 -95335   190670 441.128  1  < 2.2e-16 ***
  # mod_ot1a   23 190315 190510 -95135   190269 401.236  2  < 2.2e-16 ***
  # mod_ot2a   26 190168 190388 -95058   190116 152.954  3  < 2.2e-16 ***
  # mod_ot3a   30 190115 190369 -95028   190055  61.043  4  1.751e-12 ***
  
  
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
gca_l2_mod_l1_0 <- update(gca_l2_mod_base,   . ~ . + l1) 
gca_l2_mod_l1_1 <- update(gca_l2_mod_l1_0, . ~ . + ot1:l1) 
gca_l2_mod_l1_2 <- update(gca_l2_mod_l1_1, . ~ . + ot2:l1) 
gca_l2_mod_l1_3 <- update(gca_l2_mod_l1_2, . ~ . + ot3:l1)

l2_l1_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_l1_0, gca_l2_mod_l1_1,
        gca_l2_mod_l1_2, gca_l2_mod_l1_3)
#                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_l2_mod_base   30 188993 189246 -94466   188933                       
# gca_l2_mod_l1_0   31 188992 189254 -94465   188930 2.8402  1    0.09193 .
# gca_l2_mod_l1_1   32 188992 189262 -94464   188928 2.2740  1    0.13156  
# gca_l2_mod_l1_2   33 188992 189271 -94463   188926 1.3905  1    0.23832  
# gca_l2_mod_l1_3   34 188993 189280 -94462   188925 1.5640  1    0.21108   


# add stress effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_stress_0 <- update(gca_l2_mod_base,   . ~ . + condition_sum) 
gca_l2_mod_stress_1 <- update(gca_l2_mod_stress_0, . ~ . + ot1:condition_sum) 
gca_l2_mod_stress_2 <- update(gca_l2_mod_stress_1, . ~ . + ot2:condition_sum) # singular
gca_l2_mod_stress_3 <- update(gca_l2_mod_stress_2, . ~ . + ot3:condition_sum)

l2_stress_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_stress_0, gca_l2_mod_stress_1,
        gca_l2_mod_stress_2, gca_l2_mod_stress_3)
#                     npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)    
# gca_l2_mod_base       30 188993 189246 -94466   188933                       
# gca_l2_mod_stress_0   31 188994 189256 -94466   188932 0.5701  1     0.4502  
# gca_l2_mod_stress_1   32 188996 189267 -94466   188932 0.0232  1     0.8790  
# gca_l2_mod_stress_2   33 188997 189276 -94465   188931 1.1810  1     0.2772  
# gca_l2_mod_stress_3   34 188996 189283 -94464   188928 3.0916  1     0.0787 .


### ROUTE 1 (proficiency alone)

# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_dele_0 <- update(gca_l2_mod_base,   . ~ . + DELE) 
gca_l2_mod_dele_1 <- update(gca_l2_mod_dele_0, . ~ . + ot1:DELE) 
gca_l2_mod_dele_2 <- update(gca_l2_mod_dele_1, . ~ . + ot2:DELE)
gca_l2_mod_dele_3 <- update(gca_l2_mod_dele_2, . ~ . + ot3:DELE)

l2_dele_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_dele_0, gca_l2_mod_dele_1,
        gca_l2_mod_dele_2, gca_l2_mod_dele_3)
#                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_l2_mod_base     30 188993 189246 -94466   188933                        
# gca_l2_mod_dele_0   31 188987 189249 -94463   188925 7.5234  1    0.00609 **
# gca_l2_mod_dele_1   32 188987 189258 -94462   188923 1.9574  1    0.16179   
# gca_l2_mod_dele_2   33 188989 189268 -94461   188923 0.7517  1    0.38595   
# gca_l2_mod_dele_3   34 188990 189278 -94461   188922 0.0976  1    0.75470  


# add interaction

gca_l2_mod_dele_int_0 <- update(gca_l2_mod_dele_0,   . ~ . + l1:condition_sum:DELE) 
gca_l2_mod_dele_int_1 <- update(gca_l2_mod_dele_int_0, . ~ . + ot1:l1:condition_sum:DELE) 
gca_l2_mod_dele_int_2 <- update(gca_l2_mod_dele_int_1, . ~ . + ot2:l1:condition_sum:DELE)
gca_l2_mod_dele_int_3 <- update(gca_l2_mod_dele_int_2, . ~ . + ot3:l1:condition_sum:DELE)

l2_dele_int_anova <-
  anova(gca_l2_mod_dele_0, gca_l2_mod_dele_int_0, gca_l2_mod_dele_int_1,
        gca_l2_mod_dele_int_2, gca_l2_mod_dele_int_3)
#                       npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_l2_mod_dele_0       31 188987 189249 -94463   188925                       
# gca_l2_mod_dele_int_0   33 188990 189269 -94462   188924 1.4052  2    0.49530  
# gca_l2_mod_dele_int_1   35 188991 189287 -94460   188921 3.1764  2    0.20430  
# gca_l2_mod_dele_int_2   37 188992 189305 -94459   188918 2.8605  2    0.23925  
# gca_l2_mod_dele_int_3   39 188989 189319 -94456   188911 6.6570  2    0.03585 *

gca_l2_mod_dele_final <- gca_l2_mod_dele_int_3

### ROUTE 2 (L2 use alone)

# add L2 use effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_mod_use_0 <- update(gca_l2_mod_base,  . ~ . + percent_l2_week) 
gca_l2_mod_use_1 <- update(gca_l2_mod_use_0, . ~ . + ot1:percent_l2_week) 
gca_l2_mod_use_2 <- update(gca_l2_mod_use_1, . ~ . + ot2:percent_l2_week)
gca_l2_mod_use_3 <- update(gca_l2_mod_use_2, . ~ . + ot3:percent_l2_week)

l2_use_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_use_0, gca_l2_mod_use_1,
        gca_l2_mod_use_2, gca_l2_mod_use_3)
#                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_l2_mod_base    30 188993 189246 -94466   188933                       
# gca_l2_mod_use_0   31 188995 189257 -94466   188933 0.1961  1     0.6579  
# gca_l2_mod_use_1   32 188995 189266 -94466   188931 1.1309  1     0.2876  
# gca_l2_mod_use_2   33 188993 189272 -94463   188927 4.7501  1     0.0293 *
# gca_l2_mod_use_3   34 188991 189279 -94462   188923 3.3476  1     0.0673 .

# add interaction

gca_l2_mod_use_int_0 <- update(gca_l2_mod_use_2,     . ~ . + l1:condition_sum:percent_l2_week) 
gca_l2_mod_use_int_1 <- update(gca_l2_mod_use_int_0, . ~ . + ot1:l1:condition_sum:percent_l2_week) 
gca_l2_mod_use_int_2 <- update(gca_l2_mod_use_int_1, . ~ . + ot2:l1:condition_sum:percent_l2_week)
gca_l2_mod_use_int_3 <- update(gca_l2_mod_use_int_2, . ~ . + ot3:l1:condition_sum:percent_l2_week)

l2_use_int_anova <-
  anova(gca_l2_mod_use_0, gca_l2_mod_use_int_0, gca_l2_mod_use_int_1,
        gca_l2_mod_use_int_2, gca_l2_mod_use_int_3)
#                      npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)   
# gca_l2_mod_use_0       31 188995 189257 -94466   188933                         
# gca_l2_mod_use_int_0   35 188996 189292 -94463   188926  6.2237  4   0.183055   
# gca_l2_mod_use_int_1   37 188989 189302 -94457   188915 11.6320  2   0.002979 **
# gca_l2_mod_use_int_2   39 188993 189322 -94457   188915  0.0390  2   0.980703   
# gca_l2_mod_use_int_3   41 188990 189336 -94454   188908  6.9672  2   0.030697 * 

gca_l2_mod_use_final <- gca_l2_mod_use_int_3


### ROUTE 3 (working memory) 
gca_l2_mod_wm_0 <- update(gca_l2_mod_base, . ~ . + ospan) 
gca_l2_mod_wm_1 <- update(gca_l2_mod_wm_0, . ~ . + ot1:ospan) 
gca_l2_mod_wm_2 <- update(gca_l2_mod_wm_1, . ~ . + ot2:ospan) 
gca_l2_mod_wm_3 <- update(gca_l2_mod_wm_2, . ~ . + ot3:ospan)

l2_wm_anova <-
  anova(gca_l2_mod_base, gca_l2_mod_wm_0, gca_l2_mod_wm_1,
        gca_l2_mod_wm_2, gca_l2_mod_wm_3)
#                   npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_full_mod_base   30 188993 189246 -94466   188933                     
# gca_l2_mod_wm_0     31 188994 189257 -94466   188932 0.4262  1     0.5139
# gca_l2_mod_wm_1     32 188996 189267 -94466   188932 0.3709  1     0.5425
# gca_l2_mod_wm_2     33 188998 189277 -94466   188932 0.4146  1     0.5197
# gca_l2_mod_wm_3     34 189000 189287 -94466   188932 0.0196  1     0.8886


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_l2_mod_wm_int_0 <- update(gca_l2_mod_base, . ~ . + l1:condition_sum:ospan)  
gca_l2_mod_wm_int_1 <- update(gca_l2_mod_wm_int_0, . ~ . + ot1:l1:condition_sum:ospan) 
gca_l2_mod_wm_int_2 <- update(gca_l2_mod_wm_int_1, . ~ . + ot2:l1:condition_sum:ospan) 
gca_l2_mod_wm_int_3 <- update(gca_l2_mod_wm_int_2, . ~ . + ot3:l1:condition_sum:ospan) 

l2_wm_int_anova <- anova(gca_l2_mod_base, gca_l2_mod_wm_int_0, gca_l2_mod_wm_int_1, 
                         gca_l2_mod_wm_int_2, gca_l2_mod_wm_int_3)
#                     npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
# gca_l2_mod_base       30 188993 189246 -94466   188933                          
# gca_l2_mod_wm_int_0   32 188993 189263 -94464   188929  4.0251  2  0.1336481    
# gca_l2_mod_wm_int_1   34 188994 189282 -94463   188926  2.4010  2  0.3010364    
# gca_l2_mod_wm_int_2   36 188984 189289 -94456   188912 13.8719  2  0.0009722 ***
# gca_l2_mod_wm_int_3   38 188982 189304 -94453   188906  6.0882  2  0.0476396 * 


gca_l2_mod_wm_final <- gca_l2_mod_wm_int_3


mod_type <- "gca_l2_mod"
mod_spec <- c('_base', 
              "_l1_0", "_l1_1", "_l1_2", "_l1_3", 
              "_stress_0", "_stress_1", "_stress_2", "_stress_3", 
              "_dele_0", "_dele_1", "_dele_2", "_dele_3", 
              "_dele_int_0", "_dele_int_1", "_dele_int_2", "_dele_int_3",
              "_use_0", "_use_1", "_use_2", "_use_3", 
              "_use_int_0", "_use_int_1", "_use_int_2", "_use_int_3",
              "_wm_0", "_wm_1", "_wm_2", "_wm_3",
              "_wm_int_0", "_wm_int_1", "_wm_int_2", "_wm_int_3",
              '_dele_final', '_use_final', '_wm_final')

# Store ind models in list
gca_l2_mods <- mget(c(paste0(mod_type, mod_spec)))

save(gca_l2_mods,
     file = here("mods", "stress", "gca",
                 "gca_l2_mods.Rdata"))


# Save anova model comparisons
nested_model_comparisons <-
  mget(c("mon_cond_anova", "mon_wm_anova", "mon_int_anova",
         'l2_l1_anova', 'l2_stress_anova', 
         'l2_dele_anova', 'l2_dele_int_anova',
         'l2_use_anova', 'l2_use_int_anova',
         'l2_wm_anova', 'l2_wm_int_anova'
  ))

save(nested_model_comparisons,
     file = here("mods", "stress", "gca",
                 "nested_model_comparisons.Rdata"))


}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_mon <- mon_data %>%
  dplyr::select(time_zero, ot1:ot3, condition_sum) %>% #, ospan
  distinct

new_dat_l2 <- l2_data %>%
  dplyr::select(l1, time_zero, ot1:ot3, condition_sum, ospan, DELE, percent_l2_week) %>%
  distinct

write_csv(new_dat_l2, here::here('new_dat_l2.csv'))
new_dat_l2 <- read_csv(here::here('new_dat_l2.csv'))

# Get model predictions and SE
fits_all_mon <- predictSE(gca_mod_mon_base, new_dat_mon) %>%  
  as_tibble %>%
  bind_cols(new_dat_mon) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_all_l2_dele <- predictSE(gca_l2_mod_dele_final, new_dat_l2) %>%
  as_tibble %>%
  bind_cols(new_dat_l2) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(l1, EN = "en", MA = "ma"),
         group = fct_relevel(l1, "EN", "MA"))

fits_all_l2_use <- predictSE(gca_l2_mod_use_final, new_dat_l2) %>%  
  as_tibble %>%
  bind_cols(new_dat_l2) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(l1, EN = "en", MA = "ma"),
         group = fct_relevel(l1, "EN", "MA"))

fits_all_l2_wm <- predictSE(gca_l2_mod_wm_final, new_dat_l2) %>%  
  as_tibble %>%
  bind_cols(new_dat_l2) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(l1, EN = "en", MA = "ma"),
         group = fct_relevel(l1, "EN", "MA"))

# Filter preds at target syllable offset
target_offset_preds_mon <- filter(fits_all_mon, time_zero == 4) %>%
  select(stress = condition_sum, #ospan,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

target_offset_preds_l2_dele <- filter(fits_all_l2_dele, time_zero == 4) %>%
  select(l1, stress = condition_sum, DELE,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(l1)

target_offset_preds_l2_use <- filter(fits_all_l2_use, time_zero == 4) %>%
  select(l1, stress = condition_sum, percent_l2_week,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(l1)

target_offset_preds_l2_wm <- filter(fits_all_l2_wm, time_zero == 4) %>%
  select(l1, stress = condition_sum, wm = ospan,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(l1)

# Save models predictions
model_preds <- mget(c("fits_all_mon", "fits_all_l2_dele", "fits_all_l2_use", "fits_all_l2_wm",
                      "target_offset_preds_mon", "target_offset_preds_l2_dele",
                      "target_offset_preds_l2_use", "target_offset_preds_l2_wm"))

save(model_preds,
     file = here("mods", "stress", "gca",
                 "model_preds.Rdata"))

# MON without WM score
saveRDS(target_offset_preds_mon, file = here("mods", "stress", "gca", "model_preds_mon.Rdata"))

# -----------------------------------------------------------------------------

