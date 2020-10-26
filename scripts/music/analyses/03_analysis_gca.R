#
# Original script by Joseph
# Modified by Cristina to adapt to CrisLau Project
# Last update: 06/12/2019
# Modified by Laura to adapt to Pupurri project
#
# Growth curve analysis ------------------------------------------------------
#
# - Question 1: Does pitch prediction abilities (continuous) influence 
#   SS, IE, AE, IM, and AM (categorical) speakers' abilities to predict verbal tense
#   based on the presence or absence (categorical) of lexical stress?
#     - Pitch and stress are associated in Mandarin speakers
#     - Pitch and stress are not associated in English or Spanish speakers
# - Question 2: Does rhythm prediction abilities (continuous) influence 
#   SS, IE, AE, IM, and AM (categorical) speakers' abilities to predict verbal tense
#   based on the presence or absence (categorical) of lexical stress?
#     - Rhythm and stress are not associated in monolingual speakers (SS)
#     - Rhythm and stress are associated in L2 speakers
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

# Load data
source(here::here("scripts", "02_load_data.R"))

# Get path to saved models
gca_mods_path  <- here("mods", "music", "gca")
 
# Load models as lists
# load(paste0(gca_mods_path, "/ind_mods.Rdata"))
load(paste0(gca_mods_path, "/pitch_mods.Rdata"))
# load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
# load(paste0(gca_mods_path, "/model_preds.Rdata"))
# 
# # Store objects in global env
# list2env(ind_mods, globalenv())
list2env(pitch_mods, globalenv())
# list2env(nested_model_comparisons, globalenv())
# list2env(model_preds, globalenv())

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

auditory <- read_csv("./data/clean/auditory_scores.csv")

music50 <- merge(x = stress50, y = auditory, by = "participant", all.x=TRUE)

music50 <- na.omit(music50)


stress_gc_subset <- music50 %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "mon", "aes", "ams", "ies", "ims"),
            stress_sum = if_else(cond == "1", 1, -1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------







# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + stress_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = stress_gc_subset, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + stress_sum + ot1 | participant) +
             ot2 + (1 + stress_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + stress_sum + ot1 + ot2 | participant) +
             ot3 + (1 + stress_sum + ot1 + ot2 + ot3 | participant))
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3, mod_ot4)
  
# We retain the most complex model: mod_ot4
  ## All participants
  #           Df    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
 #  mod_ot1    9 238591 238669 -119286   238573                          
  # mod_ot2   14 237918 238039 -118945   237890 682.712  5  < 2.2e-16 ***
  # mod_ot3   20 237862 238035 -118911   237822  68.444  6  8.524e-13 ***
  # mod_ot4   21 237593 237775 -118775   237551 271.031  1  < 2.2e-16 ***
  
  
}

# -----------------------------------------------------------------------------





# Individual models -----------------------------------------------------------

#
# only mon
#

gca_mod_mon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +          
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 3e5)),    # 2e4
       REML = F, #na.action = na.exclude,
       data = filter(stress_gc_subset, group == 'mon')) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + stress_sum) # singular
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0, . ~ . + ot1:stress_sum) # singular
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1, . ~ . + ot2:stress_sum) # singular
gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2, . ~ . + ot3:stress_sum) # singular

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2, gca_mod_mon_cond_3) # 
#                      Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_base     30 42823 43032 -21382    42763                     
# gca_mod_mon_cond_0   31 42825 43041 -21382    42763 0.2074  1     0.6488
# gca_mod_mon_cond_1   32 42827 43050 -21381    42763 0.3616  1     0.5476
# gca_mod_mon_cond_2   33 42829 43059 -21381    42763 0.0055  1     0.9411
# gca_mod_mon_cond_3   34 42830 43067 -21381    42762 0.4158  1     0.5190


# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_rhythm_0 <- update(gca_mod_mon_base,     . ~ . + rhythm_dev) # singular
gca_mod_mon_rhythm_1 <- update(gca_mod_mon_rhythm_0, . ~ . + ot1:rhythm_dev) # singular
gca_mod_mon_rhythm_2 <- update(gca_mod_mon_rhythm_1, . ~ . + ot2:rhythm_dev) # singular
gca_mod_mon_rhythm_3 <- update(gca_mod_mon_rhythm_2, . ~ . + ot3:rhythm_dev) # singular

mon_rhythm_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_rhythm_0, gca_mod_mon_rhythm_1,
        gca_mod_mon_rhythm_2, gca_mod_mon_rhythm_3)
#                        Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_base       30 42823 43032 -21382     42763                     
# gca_mod_mon_rhythm_0   31 42825 43041 -21382     42763 0.0354  1     0.8507
# gca_mod_mon_rhythm_1   32 42827 43050 -21382     42763 0.0671  1     0.7956
# gca_mod_mon_rhythm_2   33 42827 43057 -21381     42761 1.8893  1     0.1693
# gca_mod_mon_rhythm_3   34 42829 43066 -21380     42761 0.3489  1     0.5547




# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_pitch_0 <- update(gca_mod_mon_base,    . ~ . + pitch_dev) # singular
gca_mod_mon_pitch_1 <- update(gca_mod_mon_pitch_0, . ~ . + ot1:pitch_dev) # singular
gca_mod_mon_pitch_2 <- update(gca_mod_mon_pitch_1, . ~ . + ot2:pitch_dev) # singular
gca_mod_mon_pitch_3 <- update(gca_mod_mon_pitch_2, . ~ . + ot3:pitch_dev) # singular

mon_pitch_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_pitch_0, gca_mod_mon_pitch_1,
        gca_mod_mon_pitch_2, gca_mod_mon_pitch_3)
#                       Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_base      30 42823 43032 -21382     42763                       
# gca_mod_mon_pitch_0   31 42821 43038 -21380     42759 3.7462  1    0.05293 .
# gca_mod_mon_pitch_1   32 42823 43047 -21380     42759 0.0019  1    0.96568  
# gca_mod_mon_pitch_2   33 42825 43056 -21380     42759 0.0238  1    0.87745  
# gca_mod_mon_pitch_3   34 42821 43058 -21377     42753 6.0752  1    0.01371 *





#######

# only aes

#######

gca_mod_aes_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
       REML = F,
       data = filter(stress_gc_subset, group == "aes")) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_cond_0 <- update(gca_mod_aes_base,   . ~ . + stress_sum) # singular
gca_mod_aes_cond_1 <- update(gca_mod_aes_cond_0, . ~ . + ot1:stress_sum) # singular
gca_mod_aes_cond_2 <- update(gca_mod_aes_cond_1, . ~ . + ot2:stress_sum) # singular
gca_mod_aes_cond_3 <- update(gca_mod_aes_cond_2, . ~ . + ot3:stress_sum) # singular

aes_cond_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_cond_0, gca_mod_aes_cond_1,
        gca_mod_aes_cond_2, gca_mod_aes_cond_3) # 
#                    Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_aes_base   30 47567 47779 -23753    47507                       
# gca_mod_aes_cond_0   31 47568 47787 -23753    47506 0.9101  1    0.34008  
# gca_mod_aes_cond_1   32 47570 47796 -23753    47506 0.0271  1    0.86917  
# gca_mod_aes_cond_2   33 47565 47799 -23750    47499 6.2033  1    0.01275 *
# gca_mod_aes_cond_3   34 47567 47808 -23750    47499 0.1145  1    0.73506  


# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_rhythm_0 <- update(gca_mod_aes_cond_2,   . ~ . + rhythm_dev) # singular
gca_mod_aes_rhythm_1 <- update(gca_mod_aes_rhythm_0, . ~ . + ot1:rhythm_dev) # singular
gca_mod_aes_rhythm_2 <- update(gca_mod_aes_rhythm_1, . ~ . + ot2:rhythm_dev) # singular
gca_mod_aes_rhythm_3 <- update(gca_mod_aes_rhythm_2, . ~ . + ot3:rhythm_dev) # singular

aes_rhythm_anova <-
  anova(gca_mod_aes_cond_2, gca_mod_aes_rhythm_0, gca_mod_aes_rhythm_1,
        gca_mod_aes_rhythm_2, gca_mod_aes_rhythm_3)
#                        Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_aes_cond_2     33 47565 47799 -23750    47499                        
# gca_mod_aes_rhythm_0   34 47560 47800 -23746    47492 7.4644  1   0.006293 **
# gca_mod_aes_rhythm_1   35 47562 47809 -23746    47492 0.0098  1   0.921239   
# gca_mod_aes_rhythm_2   36 47561 47816 -23745    47489 2.4205  1   0.119755   
# gca_mod_aes_rhythm_3   37 47563 47825 -23745    47489 0.0412  1   0.839115   


# add pitch synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_pitch_0 <- update(gca_mod_aes_rhythm_0,    . ~ . + pitch_dev) # singular
gca_mod_aes_pitch_1 <- update(gca_mod_aes_pitch_0, . ~ . + ot1:pitch_dev) # singular
gca_mod_aes_pitch_2 <- update(gca_mod_aes_pitch_1, . ~ . + ot2:pitch_dev) # singular
gca_mod_aes_pitch_3 <- update(gca_mod_aes_pitch_2, . ~ . + ot3:pitch_dev) # singular

aes_pitch_anova <-
  anova(gca_mod_aes_rhythm_0, gca_mod_aes_pitch_0, gca_mod_aes_pitch_1,
        gca_mod_aes_pitch_2, gca_mod_aes_pitch_3)
#                        Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_aes_rhythm_0   34 47560 47800 -23746    47492                     
# gca_mod_aes_pitch_0    35 47562 47809 -23746    47492 0.0000  1     0.9979
# gca_mod_aes_pitch_1    36 47563 47817 -23745    47491 1.3002  1     0.2542
# gca_mod_aes_pitch_2    37 47564 47826 -23745    47490 0.1668  1     0.6830
# gca_mod_aes_pitch_3    38 47566 47835 -23745    47490 0.1129  1     0.7369





#
# only ies
#

gca_mod_ies_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ies")) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_cond_0 <- update(gca_mod_ies_base,     . ~ . + stress_sum) # singular
gca_mod_ies_cond_1 <- update(gca_mod_ies_cond_0,   . ~ . + ot1:stress_sum) # singular
gca_mod_ies_cond_2 <- update(gca_mod_ies_cond_1,   . ~ . + ot2:stress_sum) # singular
gca_mod_ies_cond_3 <- update(gca_mod_ies_cond_2,   . ~ . + ot3:stress_sum) # singular

ies_cond_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_cond_0, gca_mod_ies_cond_1,    
        gca_mod_ies_cond_2, gca_mod_ies_cond_3)
#                      Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ies_base     30 49428 49641 -24684    49368                       
# gca_mod_ies_cond_0   31 49430 49650 -24684    49368 0.1038  1    0.74729  
# gca_mod_ies_cond_1   32 49432 49659 -24684    49368 0.3824  1    0.53632  
# gca_mod_ies_cond_2   33 49434 49668 -24684    49368 0.0885  1    0.76606  
# gca_mod_ies_cond_3   34 49432 49673 -24682    49364 3.4647  1    0.06269 .



# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_rhythm_0 <- update(gca_mod_ies_base,     . ~ . + rhythm_dev) # singular
gca_mod_ies_rhythm_1 <- update(gca_mod_ies_rhythm_0, . ~ . + ot1:rhythm_dev) # singular
gca_mod_ies_rhythm_2 <- update(gca_mod_ies_rhythm_1, . ~ . + ot2:rhythm_dev) # singular
gca_mod_ies_rhythm_3 <- update(gca_mod_ies_rhythm_2, . ~ . + ot3:rhythm_dev) # singular

ies_rhythm_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_rhythm_0, gca_mod_ies_rhythm_1,
        gca_mod_ies_rhythm_2, gca_mod_ies_rhythm_3)
#                        Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ies_base       30 49428 49641  -24684    49368                        
# gca_mod_ies_rhythm_0   31 49429 49649  -24684    49367 0.7195  1   0.396302   
# gca_mod_ies_rhythm_1   32 49431 49658  -24684    49367 0.1098  1   0.740320   
# gca_mod_ies_rhythm_2   33 49426 49660  -24680    49360 7.7533  1   0.005362 **
# gca_mod_ies_rhythm_3   34 49428 49669  -24680    49360 0.0653  1   0.798256

# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_pitch_0 <- update(gca_mod_ies_rhythm_2,    . ~ . + pitch_dev) # singular
gca_mod_ies_pitch_1 <- update(gca_mod_ies_pitch_0, . ~ . + ot1:pitch_dev) # singular
gca_mod_ies_pitch_2 <- update(gca_mod_ies_pitch_1, . ~ . + ot2:pitch_dev) # singular
gca_mod_ies_pitch_3 <- update(gca_mod_ies_pitch_2, . ~ . + ot3:pitch_dev) # singular

ies_pitch_anova <-
  anova(gca_mod_ies_rhythm_2, gca_mod_ies_pitch_0, gca_mod_ies_pitch_1,
        gca_mod_ies_pitch_2, gca_mod_ies_pitch_3)
#                        Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ies_rhythm_2   33 49426 49660  -24680    49360                     
# gca_mod_ies_pitch_0    34 49428 49669  -24680    49360 0.0745  1     0.7848
# gca_mod_ies_pitch_1    35 49428 49676  -24679    49358 1.3457  1     0.2460
# gca_mod_ies_pitch_2    36 49430 49685  -24679    49358 0.0009  1     0.9761
# gca_mod_ies_pitch_3    37 49432 49694  -24679    49358 0.1428  1     0.7055




#
# only ams
#

gca_mod_ams_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ams")) 

# add cond effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_cond_0 <- update(gca_mod_ams_base,     . ~ . + stress_sum) 
gca_mod_ams_cond_1 <- update(gca_mod_ams_cond_0,   . ~ . + ot1:stress_sum) 
gca_mod_ams_cond_2 <- update(gca_mod_ams_cond_1,   . ~ . + ot2:stress_sum) # singular
gca_mod_ams_cond_3 <- update(gca_mod_ams_cond_2,   . ~ . + ot3:stress_sum) # singular

ams_cond_anova <-
  anova(gca_mod_ams_base, gca_mod_ams_cond_0, gca_mod_ams_cond_1,
        gca_mod_ams_cond_2, gca_mod_ams_cond_3)
#                      Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ams_base     30 48201 48413 -24070    48141                       
# gca_mod_ams_cond_0   31 48203 48422 -24070    48141 0.0002  1    0.99002  
# gca_mod_ams_cond_1   32 48205 48431 -24070    48141 0.0573  1    0.81080  
# gca_mod_ams_cond_2   33 48201 48434 -24068    48135 5.7881  1    0.01613 *
# gca_mod_ams_cond_3   34 48203 48443 -24068    48135 0.0116  1    0.91414  



# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_rhythm_0 <- update(gca_mod_ams_cond_2,   . ~ . + rhythm_dev) # singular
gca_mod_ams_rhythm_1 <- update(gca_mod_ams_rhythm_0, . ~ . + ot1:rhythm_dev) # singular
gca_mod_ams_rhythm_2 <- update(gca_mod_ams_rhythm_1, . ~ . + ot2:rhythm_dev) # singular
gca_mod_ams_rhythm_3 <- update(gca_mod_ams_rhythm_2, . ~ . + ot3:rhythm_dev) # singular

ams_rhythm_anova <-
  anova(gca_mod_ams_cond_2, gca_mod_ams_rhythm_0, gca_mod_ams_rhythm_1,
        gca_mod_ams_rhythm_2, gca_mod_ams_rhythm_3)
#                        Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ams_cond_2     33 48201 48434  -24068    48135                       
# gca_mod_ams_rhythm_0   34 48203 48443  -24068    48135 0.0012  1    0.97284  
# gca_mod_ams_rhythm_1   35 48205 48452  -24067    48135 0.3367  1    0.56176  
# gca_mod_ams_rhythm_2   36 48206 48460  -24067    48134 0.9627  1    0.32652  
# gca_mod_ams_rhythm_3   37 48204 48465  -24065    48130 3.8009  1    0.05122 .



# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_pitch_0 <- update(gca_mod_ams_cond_2,  . ~ . + pitch_dev) # singular
gca_mod_ams_pitch_1 <- update(gca_mod_ams_pitch_0, . ~ . + ot1:pitch_dev) # singular
gca_mod_ams_pitch_2 <- update(gca_mod_ams_pitch_1, . ~ . + ot2:pitch_dev) # singular
gca_mod_ams_pitch_3 <- update(gca_mod_ams_pitch_2, . ~ . + ot3:pitch_dev) # singular

ams_pitch_anova <-
  anova(gca_mod_ams_cond_2, gca_mod_ams_pitch_0, gca_mod_ams_pitch_1,
        gca_mod_ams_pitch_2, gca_mod_ams_pitch_3)
#                       Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ams_cond_2    33 48201 48434 -24068    48135                       
# gca_mod_ams_pitch_0   34 48203 48443 -24067    48135 0.3291  1    0.56617  
# gca_mod_ams_pitch_1   35 48201 48448 -24066    48131 3.7538  1    0.05269 .
# gca_mod_ams_pitch_2   36 48202 48457 -24065    48130 0.5339  1    0.46497  
# gca_mod_ams_pitch_3   37 48204 48466 -24065    48130 0.0979  1    0.75439  





#
# only ims
#

gca_mod_ims_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ims")) # singular

# add cond effect to imsercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_cond_0 <- update(gca_mod_ims_base,     . ~ . + stress_sum) # singular
gca_mod_ims_cond_1 <- update(gca_mod_ims_cond_0,   . ~ . + ot1:stress_sum) # singular
gca_mod_ims_cond_2 <- update(gca_mod_ims_cond_1,   . ~ . + ot2:stress_sum) # singular
gca_mod_ims_cond_3 <- update(gca_mod_ims_cond_2,   . ~ . + ot3:stress_sum) # singular

ims_cond_anova <-
  anova(gca_mod_ims_base, gca_mod_ims_cond_0, gca_mod_ims_cond_1, # none singular, and none significant
        gca_mod_ims_cond_2, gca_mod_ims_cond_3)
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_base     30 48559 48771 -24250    48499                     
# gca_mod_ims_cond_0   31 48561 48780 -24250    48499 0.1040  1     0.7471
# gca_mod_ims_cond_1   32 48562 48788 -24249    48498 1.5108  1     0.2190
# gca_mod_ims_cond_2   33 48562 48796 -24248    48496 1.4883  1     0.2225
# gca_mod_ims_cond_3   34 48563 48803 -24247    48495 1.6916  1     0.1934


############## RUN THE NEXT IMS MODELS AGAIN WITH _base AS BASE ########################


# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_rhythm_0 <- update(gca_mod_ims_wm_2,   . ~ . + rhythm_dev) # singular
gca_mod_ims_rhythm_1 <- update(gca_mod_ims_rhythm_0, . ~ . + ot1:rhythm_dev) # singular
gca_mod_ims_rhythm_2 <- update(gca_mod_ims_rhythm_1, . ~ . + ot2:rhythm_dev) # singular
gca_mod_ims_rhythm_3 <- update(gca_mod_ims_rhythm_2, . ~ . + ot3:rhythm_dev) # singular

ims_rhythm_anova <-
  anova(gca_mod_ims_wm_2, gca_mod_ims_rhythm_0, gca_mod_ims_rhythm_1,
        gca_mod_ims_rhythm_2, gca_mod_ims_rhythm_3)
#                        Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ims_wm_2       33 48555 48788  -24244    48489                       
# gca_mod_ims_rhythm_0   34 48557 48797  -24244    48489 0.2341  1    0.62847  
# gca_mod_ims_rhythm_1   35 48554 48802  -24242    48484 4.4538  1    0.03482 *
# gca_mod_ims_rhythm_2   36 48555 48810  -24242    48483 1.3104  1    0.25233  
# gca_mod_ims_rhythm_3   37 48557 48818  -24241    48483 0.3800  1    0.53761  



# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_pitch_0 <- update(gca_mod_ims_rhythm_1,    . ~ . + pitch_dev) # singular
gca_mod_ims_pitch_1 <- update(gca_mod_ims_pitch_0, . ~ . + ot1:pitch_dev) # singular
gca_mod_ims_pitch_2 <- update(gca_mod_ims_pitch_1, . ~ . + ot2:pitch_dev) # singular
gca_mod_ims_pitch_3 <- update(gca_mod_ims_pitch_2, . ~ . + ot3:pitch_dev) # singular

ims_pitch_anova <-
  anova(gca_mod_ims_rhythm_1, gca_mod_ims_pitch_0, gca_mod_ims_pitch_1,
        gca_mod_ims_pitch_2, gca_mod_ims_pitch_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_rhythm_1   35 48554 48802 -24242    48484                        
# gca_mod_ims_pitch_0    36 48554 48809 -24241    48482 2.0275  1   0.154471   
# gca_mod_ims_pitch_1    37 48556 48818 -24241    48482 0.0605  1   0.805631   
# gca_mod_ims_pitch_2    38 48550 48819 -24237    48474 8.3610  1   0.003834 **
# gca_mod_ims_pitch_3    39 48550 48825 -24236    48472 2.1906  1   0.138853   


# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model pitch
gca_full_mod_pitch_base <-
  lmer(eLog ~ 1 + stress_sum * pitch_dev * (ot1 + ot2 + ot3) +          
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 3e5)), 
       data = stress_gc_subset, REML = F) # , na.action = na.exclude



# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_pitch_group_0 <- update(gca_full_mod_pitch_base,    . ~ . + group)
gca_full_mod_pitch_group_1 <- update(gca_full_mod_pitch_group_0, . ~ . + ot1:group)
gca_full_mod_pitch_group_2 <- update(gca_full_mod_pitch_group_1, . ~ . + ot2:group)
gca_full_mod_pitch_group_3 <- update(gca_full_mod_pitch_group_2, . ~ . + ot3:group)

full_pitch_group_anova <-
  anova(gca_full_mod_pitch_base, gca_full_mod_pitch_group_0, gca_full_mod_pitch_group_1,
        gca_full_mod_pitch_group_2, gca_full_mod_pitch_group_3)
#                              Df    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
# gca_full_mod_pitch_base      42 237452 237816 -118684   237368                          
# gca_full_mod_pitch_group_0   46 237435 237833 -118671   237343 25.5573  4  3.886e-05 ***
# gca_full_mod_pitch_group_1   50 237430 237863 -118665   237330 12.5330  4    0.01380 *  
# gca_full_mod_pitch_group_2   54 237410 237877 -118651   237302 28.5113  4  9.824e-06 ***
# gca_full_mod_pitch_group_3   58 237410 237912 -118647   237294  7.9617  4    0.09299 .  



# add 3-way int to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_pitch_int_0 <- update(gca_full_mod_pitch_group_2, . ~ . + stress_sum:pitch_dev:group)
gca_full_mod_pitch_int_1 <- update(gca_full_mod_pitch_int_0,   . ~ . + ot1:stress_sum:pitch_dev:group)
gca_full_mod_pitch_int_2 <- update(gca_full_mod_pitch_int_1,   . ~ . + ot2:stress_sum:pitch_dev:group)
gca_full_mod_pitch_int_3 <- update(gca_full_mod_pitch_int_2,   . ~ . + ot3:stress_sum:pitch_dev:group)

full_pitch_int_anova <-
  anova(gca_full_mod_pitch_group_2, gca_full_mod_pitch_int_0, gca_full_mod_pitch_int_1,
        gca_full_mod_pitch_int_2, gca_full_mod_pitch_int_3)
#                            npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
# gca_full_mod_group_2         54 237410 237877 -118651   237302                        
# gca_full_mod_pitch_int_0     58 237411 237914 -118648   237295  6.4235  4    0.16967  
# gca_full_mod_pitch_int_1     62 237411 237948 -118643   237287  8.3782  4    0.07867 .
# gca_full_mod_pitch_int_2     66 237407 237979 -118638   237275 11.5755  4    0.02080 *
# gca_full_mod_pitch_int_3     70 237414 238021 -118637   237274  0.9770  4    0.91327 



summary(gca_full_mod_pitch_int_2)

# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ams"))
gca_full_mod_pitch_int_2_ams <- update(gca_full_mod_pitch_int_2)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ims"))
gca_full_mod_pitch_int_2_ims <- update(gca_full_mod_pitch_int_2)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ies"))
gca_full_mod_pitch_int_2_ies <- update(gca_full_mod_pitch_int_2)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "aes"))
gca_full_mod_pitch_int_2_aes <- update(gca_full_mod_pitch_int_2)


mod_type <- "gca_full_mod_pitch"
mod_spec <- c("_base", "_group_0", "_group_1", "_group_2", "_group_3",
              "_int_0", "_int_1", "_int_2", "_int_3",
              "_int_2_ams", "_int_2_ims", "_int_2_aes", "_int_2_ies")

# Store ind models in list
pitch_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(pitch_mods,
     file = here("mods", "music", "gca",
                 "pitch_mods.Rdata"))





# Base model rhythm
gca_full_mod_rhythm_base <-
  lmer(eLog ~ 1 + stress_sum * rhythm_dev * (ot1 + ot2 + ot3) +          
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 3e5)), 
       data = stress_gc_subset, REML = F) # , na.action = na.exclude



# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_rhythm_group_0 <- update(gca_full_mod_rhythm_base,    . ~ . + group)
gca_full_mod_rhythm_group_1 <- update(gca_full_mod_rhythm_group_0, . ~ . + ot1:group)
gca_full_mod_rhythm_group_2 <- update(gca_full_mod_rhythm_group_1, . ~ . + ot2:group)
gca_full_mod_rhythm_group_3 <- update(gca_full_mod_rhythm_group_2, . ~ . + ot3:group)

full_rhythm_group_anova <-
  anova(gca_full_mod_rhythm_base, gca_full_mod_rhythm_group_0, gca_full_mod_rhythm_group_1,
        gca_full_mod_rhythm_group_2, gca_full_mod_rhythm_group_3)
#                               Df    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
# gca_full_mod_rhythm_base      42 237443 237807 -118680   237359                          
# gca_full_mod_rhythm_group_0   46 237425 237823 -118666   237333 26.4859  4  2.525e-05 ***
# gca_full_mod_rhythm_group_1   50 237425 237858 -118662   237325  7.8674  4    0.09656 .  
# gca_full_mod_rhythm_group_2   54 237405 237873 -118648   237297 27.8573  4  1.333e-05 ***
# gca_full_mod_rhythm_group_3   58 237404 237907 -118644   237288  8.6875  4    0.06940 . 




# add 3-way int to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_rhythm_int_0 <- update(gca_full_mod_rhythm_group_2, . ~ . + stress_sum:rhythm_dev:group)
gca_full_mod_rhythm_int_1 <- update(gca_full_mod_rhythm_int_0,   . ~ . + ot1:stress_sum:rhythm_dev:group)
gca_full_mod_rhythm_int_2 <- update(gca_full_mod_rhythm_int_1,   . ~ . + ot2:stress_sum:rhythm_dev:group)
gca_full_mod_rhythm_int_3 <- update(gca_full_mod_rhythm_int_2,   . ~ . + ot3:stress_sum:rhythm_dev:group)

full_rhythm_int_anova <-
  anova(gca_full_mod_rhythm_group_2, gca_full_mod_rhythm_int_0, gca_full_mod_rhythm_int_1,
        gca_full_mod_rhythm_int_2, gca_full_mod_rhythm_int_3)
#                      npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
# gca_full_mod_group_2  54 237405 237873 -118648   237297                          
# gca_full_mod_rhythm_int_0     58 237411 237913 -118647   237295  1.9737  4    0.74059    
# gca_full_mod_rhythm_int_1     62 237370 237907 -118623   237246 49.1550  4   5.42e-10 ***
# gca_full_mod_rhythm_int_2     66 237368 237940 -118618   237236  9.8350  4    0.04330 *  
# gca_full_mod_rhythm_int_3     70 237366 237972 -118613   237226 10.3034  4    0.03562 *  


summary(gca_full_mod_rhythm_int_3)


# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ams"))
gca_full_mod_rhythm_int_3_ams <- update(gca_full_mod_rhythm_int_3)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ims"))
gca_full_mod_rhythm_int_3_ims <- update(gca_full_mod_rhythm_int_3)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "ies"))
gca_full_mod_rhythm_int_3_ies <- update(gca_full_mod_rhythm_int_3)

stress_gc_subset %<>% mutate(., group = fct_relevel(group, "aes"))
gca_full_mod_rhythm_int_3_aes <- update(gca_full_mod_rhythm_int_3)


mod_type <- "gca_full_mod_rhythm"
mod_spec <- c("_base", "_group_0", "_group_1", "_group_2", "_group_3",
              "_int_0", "_int_1", "_int_2", "_int_3",
              "_int_3_ams", "_int_3_ims", "_int_3_aes", "_int_3_ies")

# Store ind models in list
rhythm_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(rhythm_mods,
     file = here("mods", "music", "gca",
                 "rhythm_mods.Rdata"))













}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, stress_sum, pitch_dev,
                rhythm_dev) %>%
  distinct

# Get model predictions and SE
fits_all_pitch <- predictSE(gca_full_mod_pitch_int_2, new_dat_all) %>%        #change depending on significance
  as_tibble %>%
  bind_cols(new_dat_all) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, SS = "mon", AE = "aes", IE = "ies", AM = "ams", IM = "ims"))

fits_all_rhythm <- predictSE(gca_full_mod_rhythm_int_3, new_dat_all) %>%        #change depending on significance
  as_tibble %>%
  bind_cols(new_dat_all) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, SS = "mon", AE = "aes", IE = "ies", AM = "ams", IM = "ims"))


# Filter preds at target offset
target_offset_preds_pitch <- filter(fits_all_pitch, time_zero == 4) %>%
  select(group, stress = stress_sum, pitch_dev,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

target_offset_preds_rhythm <- filter(fits_all_rhythm, time_zero == 4) %>%
  select(group, stress = stress_sum, rhythm_dev,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {
# # Build model names programatically
# mod_type <- "gca_mod_"
# mod_spec <- c("_base", "_cond_0", "_cond_1", "_cond_2", "_cond_3", 
#               "_wm_0", "_wm_1", "_wm_2", "_wm_3",
#               "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
#               "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3")
# 
# # Store ind models in list
# ind_mods <- mget(c(#paste0(mod_type, "mon", mod_spec),
#                    paste0(mod_type, "aes", mod_spec),
#                    paste0(mod_type, "ies", mod_spec),
#                    paste0(mod_type, "ams", mod_spec),
#                    paste0(mod_type, "ims", mod_spec)
#                    ))
# 
# save(ind_mods,
#      file = here("mods", "music", "gca",
#                  "ind_mods.Rdata"))
# 
# # Store full (ot1, ot2, ot3, group, single effects) models in list
# full_mods <- mget(c(
#   "gca_full_mod_base",
#   "gca_full_mod_group_0", "gca_full_mod_group_1",
#   "gca_full_mod_group_2", "gca_full_mod_group_3", 
#   "gca_full_mod_int_0", "gca_full_mod_int_1", 
#   "gca_full_mod_int_2", "gca_full_mod_int_3",
#   "gca_full_mod_int_2_ams", "gca_full_mod_int_2_ims",
#   "gca_full_mod_int_2_aes", "gca_full_mod_int_2_ies"))
# 
# save(full_mods,
#      file = here("mods", "music", "gca",
#                  "full_mods.Rdata"))
# 
# # Save anova model comparisons
# nested_model_comparisons <-
#   mget(c("mon_cond_anova", "mon_wm_anova", "mon_rhythm_anova", "mon_pitch_anova",
#          
#          "aes_cond_anova", "aes_wm_anova", "aes_rhythm_anova", "aes_pitch_anova",
#          
#          "ies_cond_anova", "ies_wm_anova", "ies_rhythm_anova", "ies_pitch_anova",
#          
#          "ams_cond_anova", "ams_wm_anova", "ams_rhythm_anova", "ams_pitch_anova",
#          
#          "ims_cond_anova", "ims_wm_anova", "ims_rhythm_anova", "ims_pitch_anova",
#          
#          "full_group_anova", "full_int_anova"))
# 
# save(nested_model_comparisons,
#      file = here("mods", "music", "gca",
#                  "nested_model_comparisons.Rdata"))









# Save models predictions
model_preds <- mget(c("fits_all_pitch", 'fits_all_rhythm', 
                      "target_offset_preds_pitch", 'target_offset_preds_rhythm'))

save(model_preds,
     file = here("mods", "music", "gca",
                 "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------

