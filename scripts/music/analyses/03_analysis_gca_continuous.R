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
gca_mods_path  <- here("mods", "music", "gca", "continuous")
 
# Load models as lists
# load(paste0(gca_mods_path, "/ind_mods.Rdata"))
# load(paste0(gca_mods_path, "/pitch_mods.Rdata"))
# load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
# load(paste0(gca_mods_path, "/model_preds.Rdata"))
# 
# # Store objects in global env
# list2env(ind_mods, globalenv())
# list2env(pitch_mods, globalenv())
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

# _dev = condval
# _sd = condsd


music50 <- music50 %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "mon", "aes", "ams", "ies", "ims"),
            stress_sum = if_else(cond == "1", 1, -1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

mon_music <- filter(music50, l1 == 'es') %>% select(-DELE, -percent_l2_week, 
                                                    -DELE_z, -use_z, -prof, -group)
l2_music <- music50 %>%
  filter(., l1 != 'es') %>% 
  filter(., participant != 'ies04' & participant != 'ies17' & participant != 'ies28' & participant != 'aes32') %>%
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1),
         prof_std = (DELE - mean(DELE))/sd(DELE),
         use_std = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week)) %>%
  select(-DELE, -percent_l2_week, 
         -DELE_z, -use_z, -prof, -group)


# -----------------------------------------------------------------------------







# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + stress_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = mon_music, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + stress_sum + ot1 | participant) +
             ot2 + (1 + stress_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + stress_sum + ot1 + ot2 | participant) +
             ot3 + (1 + stress_sum + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #           Df  AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    9 42875 42938 -21429    42857                         
  # mod_ot2   14 42772 42869 -21372    42744 113.51  5  < 2.2e-16 ***
  # mod_ot3   20 42670 42809 -21315    42630 114.12  6  8.524e-13 *** 
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                         + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                         + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #         npar   AIC   BIC logLik deviance    Chisq Df Pr(>Chisq)    
  # mod_ot3   20 42670 42809 -21315    42630                           
  # mod_ot4   21 42532 42678 -21245    42490 139.6675  1    < 2e-16 ***
  # mod_ot5   23 42339 42500 -21147    42293 196.6082  2    < 2e-16 ***
  # mod_ot6   26 42262 42443 -21105    42210  83.5892  3    < 2e-16 ***
  # mod_ot7   30 42260 42469 -21100    42200   9.9532  4  5.541e-07 ***

  
  
}

# -----------------------------------------------------------------------------





# Individual models -----------------------------------------------------------
if (F) {
#
# only mon
#

gca_mod_mon_base <- mod_ot7
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +          
  #        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 + ot3 | target),
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 3e5)),    # 2e4
  #      REML = F, #na.action = na.exclude,
  #      data = mon_music)


# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + stress_sum) 
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0, . ~ . + ot1:stress_sum)
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1, . ~ . + ot2:stress_sum) 
gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2, . ~ . + ot3:stress_sum) 

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2, gca_mod_mon_cond_3)  
#                      Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_base     30 42260 42469 -21100    42200                     
# gca_mod_mon_cond_0   31 42261 42477 -21100    42199 0.3277  1     0.5670
# gca_mod_mon_cond_1   32 42262 42485 -21099    42198 1.3291  1     0.2490
# gca_mod_mon_cond_2   33 42264 42494 -21099    42198 0.0041  1     0.9490
# gca_mod_mon_cond_3   34 42266 42503 -21099    42198 0.3968  1     0.5288


# BRANCH #1
# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_rhythm_0 <- update(gca_mod_mon_base,     . ~ . + rhythm_dev) 
gca_mod_mon_rhythm_1 <- update(gca_mod_mon_rhythm_0, . ~ . + ot1:rhythm_dev) 
gca_mod_mon_rhythm_2 <- update(gca_mod_mon_rhythm_1, . ~ . + ot2:rhythm_dev) 
gca_mod_mon_rhythm_3 <- update(gca_mod_mon_rhythm_2, . ~ . + ot3:rhythm_dev) 

mon_rhythm_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_rhythm_0, gca_mod_mon_rhythm_1,
        gca_mod_mon_rhythm_2, gca_mod_mon_rhythm_3)
#                        Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_base       30 42260 42469 -21100    42200                     
# gca_mod_mon_rhythm_0   31 42262 42478 -21100    42200 0.1418  1     0.7065
# gca_mod_mon_rhythm_1   32 42264 42487 -21100    42200 0.0013  1     0.9713
# gca_mod_mon_rhythm_2   33 42265 42495 -21100    42199 0.1828  1     0.6690
# gca_mod_mon_rhythm_3   34 42267 42504 -21100    42199 0.3932  1     0.5306


gca_mod_mon_rhythm_int_0 <- update(gca_mod_mon_base,     . ~ . + rhythm_dev:stress_sum) 
gca_mod_mon_rhythm_int_1 <- update(gca_mod_mon_rhythm_int_0, . ~ . + ot1:rhythm_dev:stress_sum) 
gca_mod_mon_rhythm_int_2 <- update(gca_mod_mon_rhythm_int_1, . ~ . + ot2:rhythm_dev:stress_sum) 
gca_mod_mon_rhythm_int_3 <- update(gca_mod_mon_rhythm_int_2, . ~ . + ot3:rhythm_dev:stress_sum) 

mon_rhythm_int_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_rhythm_int_0, gca_mod_mon_rhythm_int_1,
        gca_mod_mon_rhythm_int_2, gca_mod_mon_rhythm_int_3)
# npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
# gca_mod_mon_base           30 42260 42469 -21100    42200                          
# gca_mod_mon_rhythm_int_0   31 42259 42475 -21098    42197  3.2471  1  0.0715518 .  
# gca_mod_mon_rhythm_int_1   32 42246 42469 -21091    42182 14.3812  1  0.0001493 ***
# gca_mod_mon_rhythm_int_2   33 42246 42476 -21090    42180  2.3533  1  0.1250163    
# gca_mod_mon_rhythm_int_3   34 42248 42485 -21090    42180  0.0300  1  0.8624037    


# BRANCH #2
# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_pitch_0 <- update(gca_mod_mon_base,    . ~ . + pitch_dev) 
gca_mod_mon_pitch_1 <- update(gca_mod_mon_pitch_0, . ~ . + ot1:pitch_dev) 
gca_mod_mon_pitch_2 <- update(gca_mod_mon_pitch_1, . ~ . + ot2:pitch_dev) 
gca_mod_mon_pitch_3 <- update(gca_mod_mon_pitch_2, . ~ . + ot3:pitch_dev) 

mon_pitch_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_pitch_0, gca_mod_mon_pitch_1,
        gca_mod_mon_pitch_2, gca_mod_mon_pitch_3)
#                       Df   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_base      30 42260 42469 -21100    42200                        
# gca_mod_mon_pitch_0   31 42253 42469 -21096    42191 8.3130  1   0.003936 **
# gca_mod_mon_pitch_1   32 42255 42478 -21095    42191 0.8621  1   0.353165   
# gca_mod_mon_pitch_2   33 42254 42484 -21094    42188 2.5969  1   0.107071   
# gca_mod_mon_pitch_3   34 42254 42491 -21093    42186 1.9670  1   0.160771 

gca_mod_mon_pitch_int_0 <- update(gca_mod_mon_pitch_0,     . ~ . + pitch_dev:stress_sum) 
gca_mod_mon_pitch_int_1 <- update(gca_mod_mon_pitch_int_0, . ~ . + ot1:pitch_dev:stress_sum) 
gca_mod_mon_pitch_int_2 <- update(gca_mod_mon_pitch_int_1, . ~ . + ot2:pitch_dev:stress_sum) 
gca_mod_mon_pitch_int_3 <- update(gca_mod_mon_pitch_int_2, . ~ . + ot3:pitch_dev:stress_sum) 

mon_pitch_int_anova <-
  anova(gca_mod_mon_pitch_0, gca_mod_mon_pitch_int_0, gca_mod_mon_pitch_int_1,
        gca_mod_mon_pitch_int_2, gca_mod_mon_pitch_int_3)
# npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_mon_pitch_0       31 42253 42469 -21096    42191                     
# gca_mod_mon_pitch_int_0   32 42253 42476 -21095    42189 1.9588  1     0.1616
# gca_mod_mon_pitch_int_1   33 42255 42485 -21094    42189 0.5335  1     0.4651
# gca_mod_mon_pitch_int_2   34 42256 42493 -21094    42188 0.8872  1     0.3462
# gca_mod_mon_pitch_int_3   35 42258 42502 -21094    42188 0.0498  1     0.8233

summary(gca_mod_mon_rhythm_int_1)
#             Estimate Std. Error t value
# (Intercept)   1.6357     0.1417  11.539
# ot1                         4.6060     0.4660   9.883
# ot2                        -1.1433     0.3859  -2.963
# ot3                        -0.9423     0.2666  -3.535
# rhythm_dev:stress_sum       0.9144     0.4035   2.266
# ot1:rhythm_dev:stress_sum   3.5447     0.9279   3.820
summary(gca_mod_mon_pitch_0)
#             Estimate Std. Error t value
# (Intercept)   1.6779     0.1488  11.279
# ot1           4.5813     0.4662   9.828
# ot2          -1.1731     0.3880  -3.024
# ot3          -0.9433     0.2680  -3.520
# pitch_dev     0.8316     0.2330   3.569

mod_type <- "gca_mod_mon"
mod_spec <- c("_base", 
              "_cond_0", "_cond_1", "_cond_2", "_cond_3",
              "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
              "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3",
              "_rhythm_int_0", "_rhythm_int_1", "_rhythm_int_2", "_rhythm_int_3",
              "_pitch_int_0", "_pitch_int_1", "_pitch_int_2", "_pitch_int_3")

# Store ind models in list
mon_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods,
     file = here("mods", "music", "gca", "continuous",
                 "mon_mods.Rdata"))




#
# only EN
#

en_music <- l2_music %>%
  filter(., l1 == 'en')

# Random effects

mod_ot1 <-
  lmer(eLog ~ 1 + ot1 +
         (1 + stress_sum + ot1 | participant),
       control = lmerControl(optimizer = 'bobyqa'),
       data = en_music, weights = 1/wts, REML = F)

mod_ot2 <-
  update(mod_ot1, . ~ . -(1 + stress_sum + ot1 | participant) +
           ot2 + (1 + stress_sum + ot1 + ot2 | participant))

mod_ot3 <-
  update(mod_ot2, . ~ . -(1 + stress_sum + ot1 + ot2 | participant) +
           ot3 + (1 + stress_sum + ot1 + ot2 + ot3 | participant))

anova(mod_ot1, mod_ot2, mod_ot3)
#           Df    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
# mod_ot1    9 90417 90487 -45200    90399                          
# mod_ot2   14 90391 90499 -45182    90363  36.161  5   8.82e-07 ***
# mod_ot3   20 90142 90296 -45051    90102 261.153  6  < 2.2e-16 *** 

mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))

mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))

mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                    + (1 + ot1 + ot2 | target))

mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                    + (1 + ot1 + ot2 + ot3 | target))

anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
#           npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
#   mod_ot3   90142 90296 -45051    90102                          
#   mod_ot4   21 89774 89936 -44866    89732 369.684  1  < 2.2e-16 ***
#   mod_ot5   23 89554 89731 -44754    89508 224.645  2  < 2.2e-16 ***
#   mod_ot6   26 89429 89629 -44688    89377 131.175  3  < 2.2e-16 ***
#   mod_ot7   30 89382 89613 -44661    89322  54.298  4  4.558e-11 ***



# Fixed effects

gca_mod_en_base <- mod_ot7
# lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
#        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
#        (1 + ot1 + ot2 + ot3 | target),
#      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
#      REML = F,
#      data = en_music

# add proficiency synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_prof_0 <- update(gca_mod_en_base,     . ~ . + prof_std) 
gca_mod_en_prof_1 <- update(gca_mod_en_prof_0, . ~ . + ot1:prof_std) 
gca_mod_en_prof_2 <- update(gca_mod_en_prof_1, . ~ . + ot2:prof_std) 
gca_mod_en_prof_3 <- update(gca_mod_en_prof_2, . ~ . + ot3:prof_std) 

en_prof_anova <-
  anova(gca_mod_en_base, gca_mod_en_prof_0, gca_mod_en_prof_1,
        gca_mod_en_prof_2, gca_mod_en_prof_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mod_en_base     30 89382 89613 -44661    89322                       
# gca_mod_en_prof_0   31 89384 89623 -44661    89322 0.3159  1    0.57409  
# gca_mod_en_prof_1   32 89381 89628 -44658    89317 4.9547  1    0.02602 *
# gca_mod_en_prof_2   33 89378 89632 -44656    89312 4.8743  1    0.02726 *
# gca_mod_en_prof_3   34 89379 89641 -44656    89311 0.7857  1    0.37540  


# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_cond_0 <- update(gca_mod_en_prof_2,   . ~ . + stress_sum) 
gca_mod_en_cond_1 <- update(gca_mod_en_cond_0, . ~ . + ot1:stress_sum) 
gca_mod_en_cond_2 <- update(gca_mod_en_cond_1, . ~ . + ot2:stress_sum) 
gca_mod_en_cond_3 <- update(gca_mod_en_cond_2, . ~ . + ot3:stress_sum) 

en_cond_anova <-
  anova(gca_mod_en_prof_2, gca_mod_en_cond_0, gca_mod_en_cond_1,
        gca_mod_en_cond_2, gca_mod_en_cond_3) 
#                     Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_en_prof_2   33 89378 89632 -44656    89312                       
# gca_mod_en_cond_0   34 89380 89642 -44656    89312 0.4687  1    0.49359  
# gca_mod_en_cond_1   35 89381 89651 -44656    89311 0.4470  1    0.50376  
# gca_mod_en_cond_2   36 89377 89654 -44652    89305 6.4153  1    0.01131 *
# gca_mod_en_cond_3   37 89378 89663 -44652    89304 1.1393  1    0.28580  


# BRANCH #1
# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_rhythm_0 <- update(gca_mod_en_cond_2,   . ~ . + rhythm_dev) 
gca_mod_en_rhythm_1 <- update(gca_mod_en_rhythm_0, . ~ . + ot1:rhythm_dev) 
gca_mod_en_rhythm_2 <- update(gca_mod_en_rhythm_1, . ~ . + ot2:rhythm_dev) 
gca_mod_en_rhythm_3 <- update(gca_mod_en_rhythm_2, . ~ . + ot3:rhythm_dev) 

en_rhythm_anova <-
  anova(gca_mod_en_cond_2, gca_mod_en_rhythm_0, gca_mod_en_rhythm_1,
        gca_mod_en_rhythm_2, gca_mod_en_rhythm_3)
#                       Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_en_cond_2     36 89377 89654 -44652    89305                       
# gca_mod_en_rhythm_0   37 89375 89661 -44651    89301 3.2529  1     0.0713 .
# gca_mod_en_rhythm_1   38 89376 89669 -44650    89300 1.5477  1     0.2135  
# gca_mod_en_rhythm_2   39 89378 89678 -44650    89300 0.0968  1     0.7557  
# gca_mod_en_rhythm_3   40 89378 89686 -44649    89298 1.9494  1     0.1627 


gca_mod_en_rhythm_int_0 <- update(gca_mod_en_cond_2,   . ~ . + prof_std:stress_sum:rhythm_dev) 
gca_mod_en_rhythm_int_1 <- update(gca_mod_en_rhythm_int_0, . ~ . + ot1:prof_std:stress_sum:rhythm_dev) 
gca_mod_en_rhythm_int_2 <- update(gca_mod_en_rhythm_int_1, . ~ . + ot2:prof_std:stress_sum:rhythm_dev) 
gca_mod_en_rhythm_int_3 <- update(gca_mod_en_rhythm_int_2, . ~ . + ot3:prof_std:stress_sum:rhythm_dev) 

en_rhythm_int_anova <-
  anova(gca_mod_en_cond_2, gca_mod_en_rhythm_int_0, gca_mod_en_rhythm_int_1,
        gca_mod_en_rhythm_int_2, gca_mod_en_rhythm_int_3)
#                         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_en_cond_2         36 89377 89654 -44652    89305                       
# gca_mod_en_rhythm_int_0   37 89377 89662 -44651    89303 1.8493  1    0.17386  
# gca_mod_en_rhythm_int_1   38 89374 89667 -44649    89298 5.2648  1    0.02176 *
# gca_mod_en_rhythm_int_2   39 89376 89676 -44649    89298 0.0003  1    0.98685  
# gca_mod_en_rhythm_int_3   40 89376 89684 -44648    89296 1.9627  1    0.16122  


# BRANCH #2
# add pitch processing effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_pitch_0 <- update(gca_mod_en_cond_2,   . ~ . + pitch_dev) 
gca_mod_en_pitch_1 <- update(gca_mod_en_pitch_0, . ~ . + ot1:pitch_dev) 
gca_mod_en_pitch_2 <- update(gca_mod_en_pitch_1, . ~ . + ot2:pitch_dev) 
gca_mod_en_pitch_3 <- update(gca_mod_en_pitch_2, . ~ . + ot3:pitch_dev) 

en_pitch_anova <-
  anova(gca_mod_en_cond_2, gca_mod_en_pitch_0, gca_mod_en_pitch_1,
        gca_mod_en_pitch_2, gca_mod_en_pitch_3)
#                       Df    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_en_cond_2    36 89377 89654 -44652    89305                     
# gca_mod_en_pitch_0   37 89378 89664 -44652    89304 0.4518  1     0.5015
# gca_mod_en_pitch_1   38 89380 89673 -44652    89304 0.0307  1     0.8608
# gca_mod_en_pitch_2   39 89382 89683 -44652    89304 0.0451  1     0.8318
# gca_mod_en_pitch_3   40 89383 89692 -44652    89303 0.7578  1     0.3840


gca_mod_en_pitch_int_0 <- update(gca_mod_en_cond_2,   . ~ . + prof_std:stress_sum:pitch_dev) 
gca_mod_en_pitch_int_1 <- update(gca_mod_en_pitch_int_0, . ~ . + ot1:prof_std:stress_sum:pitch_dev) 
gca_mod_en_pitch_int_2 <- update(gca_mod_en_pitch_int_1, . ~ . + ot2:prof_std:stress_sum:pitch_dev) 
gca_mod_en_pitch_int_3 <- update(gca_mod_en_pitch_int_2, . ~ . + ot3:prof_std:stress_sum:pitch_dev) 

en_pitch_int_anova <-
  anova(gca_mod_en_cond_2, gca_mod_en_pitch_int_0, gca_mod_en_pitch_int_1,
        gca_mod_en_pitch_int_2, gca_mod_en_pitch_int_3)
#                        npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_en_cond_2        36 89377 89654 -44652    89305                       
# gca_mod_en_pitch_int_0   37 89375 89660 -44650    89301 3.7721  1    0.05211 .
# gca_mod_en_pitch_int_1   38 89377 89669 -44650    89301 0.4045  1    0.52476  
# gca_mod_en_pitch_int_2   39 89379 89679 -44650    89301 0.0329  1    0.85599  
# gca_mod_en_pitch_int_3   40 89380 89689 -44650    89300 0.1694  1    0.68061  

summary(gca_mod_en_rhythm_int_1)
# Estimate Std. Error t value
# (Intercept)                         1.25888    0.09919  12.691
# ot1                                 5.53020    0.33508  16.504
# ot2                                 0.28946    0.26950   1.074
# ot3                                -1.36985    0.22785  -6.012
# prof_std                            0.09053    0.05853   1.547
# stress_sum                         -0.04794    0.07481  -0.641
# ot1:prof_std                        0.44010    0.20844   2.111
# ot2:prof_std                       -0.30895    0.14329  -2.156
# ot1:stress_sum                      0.16524    0.23457   0.704
# ot2:stress_sum                      0.48299    0.17851   2.706
# prof_std:stress_sum:rhythm_dev      0.64703    0.42345   1.528
# ot1:prof_std:stress_sum:rhythm_dev  2.35704    1.02510   2.299
summary(gca_mod_en_cond_2) # no effect of pitch (maybe gca_mod_en_pitch_int_0 if anything?)
# Estimate Std. Error t value
# (Intercept)     1.25778    0.09927  12.671
# ot1             5.52623    0.33506  16.493
# ot2             0.28886    0.26966   1.071
# ot3            -1.36833    0.22847  -5.989
# prof_std        0.09001    0.05848   1.539
# stress_sum     -0.03643    0.07484  -0.487
# ot1:prof_std    0.44292    0.20953   2.114
# ot2:prof_std   -0.32887    0.14242  -2.309
# ot1:stress_sum  0.20675    0.23327   0.886
# ot2:stress_sum  0.48724    0.17797   2.738


mod_type <- "gca_mod_en"
mod_spec <- c("_base", 
              "_cond_0", "_cond_1", "_cond_2", "_cond_3",
              "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
              "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3",
              "_prof_0", "_prof_1", "_prof_2", "_prof_3",
              "_rhythm_int_0", "_rhythm_int_1", "_rhythm_int_2", "_rhythm_int_3",
              "_pitch_int_0", "_pitch_int_1", "_pitch_int_2", "_pitch_int_3")

# Store ind models in list
en_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(en_mods,
     file = here("mods", "music", "gca", "continuous",
                 "en_mods.Rdata"))




#
# only Ma CH
#

ma_music <- l2_music %>%
  filter(., l1 == 'ma')

# Random effects

mod_ot1 <-
  lmer(eLog ~ 1 + ot1 +
         (1 + stress_sum + ot1 | participant),
       control = lmerControl(optimizer = 'bobyqa'),
       data = ma_music, weights = 1/wts, REML = F)

mod_ot2 <-
  update(mod_ot1, . ~ . -(1 + stress_sum + ot1 | participant) +
           ot2 + (1 + stress_sum + ot1 + ot2 | participant))

mod_ot3 <-
  update(mod_ot2, . ~ . -(1 + stress_sum + ot1 + ot2 | participant) +
           ot3 + (1 + stress_sum + ot1 + ot2 + ot3 | participant))

anova(mod_ot1, mod_ot2, mod_ot3)
#           Df    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
# mod_ot1    9 95518 95588 -47750    95500                          
# mod_ot2   14 95477 95585 -47724    95449  51.274  5  7.599e-10 ***
# mod_ot3   20 95314 95469 -47637    95274 175.060  6  < 2.2e-16 ***

mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))

mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))

mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                    + (1 + ot1 + ot2 | target))

mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                    + (1 + ot1 + ot2 + ot3 | target))

anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
#           npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
#   mod_ot3   20 95314 95469 -47637    95274                           
#   mod_ot4   21 95162 95325 -47560    95120 153.6754  1  < 2.2e-16 ***
#   mod_ot5   23 94970 95149 -47462    94924 195.9337  2  < 2.2e-16 ***
#   mod_ot6   26 94935 95137 -47442    94883  40.7541  3  7.374e-09 ***
#   mod_ot7   30 94935 95167 -47437    94875   8.5841  4    0.07238 .  



# Fixed effects

gca_mod_ma_base <- mod_ot6
# lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
#        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
#        (1 + ot1 + ot2 | target),
#      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
#      REML = F,
#      data = ma_music

# add proficimacy synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_prof_0 <- update(gca_mod_ma_base,     . ~ . + prof_std) 
gca_mod_ma_prof_1 <- update(gca_mod_ma_prof_0, . ~ . + ot1:prof_std) 
gca_mod_ma_prof_2 <- update(gca_mod_ma_prof_1, . ~ . + ot2:prof_std) 
gca_mod_ma_prof_3 <- update(gca_mod_ma_prof_2, . ~ . + ot3:prof_std) 

ma_prof_anova <-
  anova(gca_mod_ma_base, gca_mod_ma_prof_0, gca_mod_ma_prof_1,
        gca_mod_ma_prof_2, gca_mod_ma_prof_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mod_ma_base     26 94935 95137 -47442    94883                       
# gca_mod_ma_prof_0   27 94932 95141 -47439    94878 5.6652  1     0.0173 *
# gca_mod_ma_prof_1   28 94931 95148 -47438    94875 2.5425  1     0.1108  
# gca_mod_ma_prof_2   29 94933 95158 -47438    94875 0.1014  1     0.7502  
# gca_mod_ma_prof_3   30 94934 95167 -47437    94874 0.5995  1     0.4388 


# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_cond_0 <- update(gca_mod_ma_prof_0,   . ~ . + stress_sum) 
gca_mod_ma_cond_1 <- update(gca_mod_ma_cond_0, . ~ . + ot1:stress_sum) 
gca_mod_ma_cond_2 <- update(gca_mod_ma_cond_1, . ~ . + ot2:stress_sum) 
gca_mod_ma_cond_3 <- update(gca_mod_ma_cond_2, . ~ . + ot3:stress_sum) 

ma_cond_anova <-
  anova(gca_mod_ma_prof_0, gca_mod_ma_cond_0, gca_mod_ma_cond_1,
        gca_mod_ma_cond_2, gca_mod_ma_cond_3) 
#                     Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ma_prof_0   27 94932 95141 -47439    94878                       
# gca_mod_ma_cond_0   28 94929 95146 -47436    94873 4.8571  1    0.02753 *
# gca_mod_ma_cond_1   29 94931 95156 -47436    94873 0.0081  1    0.92818  
# gca_mod_ma_cond_2   30 94933 95166 -47436    94873 0.0214  1    0.88368  
# gca_mod_ma_cond_3   31 94934 95175 -47436    94872 0.4953  1    0.48159


# BRANCH #1
# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_rhythm_0 <- update(gca_mod_ma_cond_0,   . ~ . + rhythm_dev) 
gca_mod_ma_rhythm_1 <- update(gca_mod_ma_rhythm_0, . ~ . + ot1:rhythm_dev) 
gca_mod_ma_rhythm_2 <- update(gca_mod_ma_rhythm_1, . ~ . + ot2:rhythm_dev) 
gca_mod_ma_rhythm_3 <- update(gca_mod_ma_rhythm_2, . ~ . + ot3:rhythm_dev) 

ma_rhythm_anova <-
  anova(gca_mod_ma_cond_0, gca_mod_ma_rhythm_0, gca_mod_ma_rhythm_1,
        gca_mod_ma_rhythm_2, gca_mod_ma_rhythm_3)
#                       Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ma_cond_0     28 94929 95146 -47436    94873                          
# gca_mod_ma_rhythm_0   29 94930 95155 -47436    94872  0.4678  1  0.4940202    
# gca_mod_ma_rhythm_1   30 94919 95152 -47430    94859 13.3358  1  0.0002604 ***
# gca_mod_ma_rhythm_2   31 94920 95161 -47429    94858  0.7953  1  0.3725159    
# gca_mod_ma_rhythm_3   32 94922 95170 -47429    94858  0.0247  1  0.8751221   


gca_mod_ma_rhythm_int_0 <- update(gca_mod_ma_rhythm_1,   . ~ . + prof_std:stress_sum:rhythm_dev) 
gca_mod_ma_rhythm_int_1 <- update(gca_mod_ma_rhythm_int_0, . ~ . + ot1:prof_std:stress_sum:rhythm_dev) 
gca_mod_ma_rhythm_int_2 <- update(gca_mod_ma_rhythm_int_1, . ~ . + ot2:prof_std:stress_sum:rhythm_dev) 
gca_mod_ma_rhythm_int_3 <- update(gca_mod_ma_rhythm_int_2, . ~ . + ot3:prof_std:stress_sum:rhythm_dev) 

ma_rhythm_int_anova <-
  anova(gca_mod_ma_rhythm_1, gca_mod_ma_rhythm_int_0, gca_mod_ma_rhythm_int_1,
        gca_mod_ma_rhythm_int_2, gca_mod_ma_rhythm_int_3)
#                         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ma_rhythm_1       30 94919 95152 -47430    94859                       
# gca_mod_ma_rhythm_int_0   31 94920 95160 -47429    94858 1.2100  1    0.27133  
# gca_mod_ma_rhythm_int_1   32 94917 95166 -47427    94853 4.3762  1    0.03644 *
# gca_mod_ma_rhythm_int_2   33 94919 95175 -47427    94853 0.3388  1    0.56055  
# gca_mod_ma_rhythm_int_3   34 94920 95184 -47426    94852 0.7984  1    0.37158  

# BRANCH #2
# add pitch processing effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_pitch_0 <- update(gca_mod_ma_cond_0,   . ~ . + pitch_dev) 
gca_mod_ma_pitch_1 <- update(gca_mod_ma_pitch_0, . ~ . + ot1:pitch_dev) 
gca_mod_ma_pitch_2 <- update(gca_mod_ma_pitch_1, . ~ . + ot2:pitch_dev) 
gca_mod_ma_pitch_3 <- update(gca_mod_ma_pitch_2, . ~ . + ot3:pitch_dev) 

ma_pitch_anova <-
  anova(gca_mod_ma_cond_0, gca_mod_ma_pitch_0, gca_mod_ma_pitch_1,
        gca_mod_ma_pitch_2, gca_mod_ma_pitch_3)
#                       Df    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ma_cond_0    28 94929 95146 -47436    94873                     
# gca_mod_ma_pitch_0   29 94930 95155 -47436    94872 0.8889  1     0.3458
# gca_mod_ma_pitch_1   30 94931 95164 -47436    94871 0.5645  1     0.4525
# gca_mod_ma_pitch_2   31 94933 95174 -47436    94871 0.0851  1     0.7705
# gca_mod_ma_pitch_3   32 94934 95182 -47435    94870 1.5375  1     0.2150


gca_mod_ma_pitch_int_0 <- update(gca_mod_ma_cond_0,   . ~ . + prof_std:stress_sum:pitch_dev) 
gca_mod_ma_pitch_int_1 <- update(gca_mod_ma_pitch_int_0, . ~ . + ot1:prof_std:stress_sum:pitch_dev) 
gca_mod_ma_pitch_int_2 <- update(gca_mod_ma_pitch_int_1, . ~ . + ot2:prof_std:stress_sum:pitch_dev) 
gca_mod_ma_pitch_int_3 <- update(gca_mod_ma_pitch_int_2, . ~ . + ot3:prof_std:stress_sum:pitch_dev) 

ma_pitch_int_anova <-
  anova(gca_mod_ma_cond_0, gca_mod_ma_pitch_int_0, gca_mod_ma_pitch_int_1,
        gca_mod_ma_pitch_int_2, gca_mod_ma_pitch_int_3)
#                        npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ma_cond_0        28 94929 95146 -47436    94873                         
# gca_mod_ma_pitch_int_0   29 94931 95156 -47436    94873  0.1788  1   0.672429   
# gca_mod_ma_pitch_int_1   30 94932 95165 -47436    94872  0.7680  1   0.380852   
# gca_mod_ma_pitch_int_2   31 94923 95164 -47431    94861 10.5754  1   0.001146 **
# gca_mod_ma_pitch_int_3   32 94924 95173 -47430    94860  0.8615  1   0.353325 

summary( gca_mod_ma_rhythm_int_1 )
# Estimate Std. Error t value
# (Intercept)                         0.97103    0.09251  10.497
# ot1                                 4.23407    0.30936  13.687
# ot2                                 0.25581    0.21124   1.211
# ot3                                -0.97697    0.14395  -6.787
# prof_std                            0.14857    0.05826   2.550
# stress_sum                         -0.12161    0.05652  -2.152
# rhythm_dev                          0.91484    0.72378   1.264
# ot1:rhythm_dev                      8.31169    2.14629   3.873
# prof_std:stress_sum:rhythm_dev     -0.56747    0.45613  -1.244
# ot1:prof_std:stress_sum:rhythm_dev -2.24021    1.06654  -2.100
summary(gca_mod_ma_pitch_int_2)
# Estimate Std. Error t value
# (Intercept)                        0.97289    0.09318  10.441
# ot1                                4.25819    0.32287  13.189
# ot2                                0.24707    0.21351   1.157
# ot3                               -0.96121    0.14475  -6.641
# prof_std                           0.14478    0.05743   2.521
# stress_sum                        -0.13778    0.05752  -2.395
# prof_std:stress_sum:pitch_dev     -0.02859    0.15022  -0.190
# ot1:prof_std:stress_sum:pitch_dev  0.37075    0.35455   1.046
# ot2:prof_std:stress_sum:pitch_dev  1.09612    0.33526   3.269


mod_type <- "gca_mod_ma"
mod_spec <- c("_base", 
              "_cond_0", "_cond_1", "_cond_2", "_cond_3",
              "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
              "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3",
              "_prof_0", "_prof_1", "_prof_2", "_prof_3",
              "_rhythm_int_0", "_rhythm_int_1", "_rhythm_int_2", "_rhythm_int_3",
              "_pitch_int_0", "_pitch_int_1", "_pitch_int_2", "_pitch_int_3")

# Store ind models in list
ma_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(ma_mods,
     file = here("mods", "music", "gca", "continuous",
                 "ma_mods.Rdata"))



# correlation pitch & rhythm

library("report") 
library(nortest)

# ASSUMPTION 1. Normality
ggplot(mon_music, aes(x = pitch_dev)) + 
  geom_density()

ggplot(mon_music, aes(x=rhythm_dev)) + 
  geom_density()


ad.test(mon_music$pitch_dev) # Anderson-Darling normality test for large samples
# A = 64.874, p-value < 2.2e-16

ad.test(mon_music$rhythm_dev) # neither normal > Spearman
# A = 336.33, p-value < 2.2e-16

# Spearman's correlation test
report(cor.test(mon_music$rhythm_dev, mon_music$pitch_dev, method = "spearman", exact=FALSE))
# The Spearman's rank correlation rho between mon_music$rhythm_dev and mon_music$pitch_dev 
# is negative, statistically significant, and small (rho = 0.18, S = 6.58e+10, p < .001)

corr_mon <- ggplot(mon_music) +
  aes(x = pitch_dev, y = rhythm_dev) +
  geom_point(size = .8) + #colour = "#0c4c8a"
  xlab("Pitch anticipation") + ylab("Rhythm synchronization") +
  geom_smooth(method=lm) + # for linear
  theme_gray(base_size = 12,
             base_family = "Times") +
  theme(legend.position = "bottom",
        legend.box = "vertical")

# annotate plot
library(grid)
# Create a text
grob <- grobTree(textGrob("rho = 0.18, S = .000, p < .001", x=0.05,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10))) #, fontface="italic"
# Plot
corr_mon_annot <- corr_mon + annotation_custom(grob)


# ASSUMPTION 1. Normality
ggplot(en_music, aes(x = pitch_dev)) + 
  geom_density()

ggplot(en_music, aes(x=rhythm_dev)) + 
  geom_density()


ad.test(en_music$pitch_dev) # Anderson-Darling normality test for large samples
# A = 602.87, p-value < 2.2e-16

ad.test(en_music$rhythm_dev) # neither normal > Spearman
# A = 443.87, p-value < 2.2e-16

# Spearman's correlation test
report(cor.test(en_music$rhythm_dev, en_music$pitch_dev, method = "spearman", exact=FALSE))
# The Spearman's rank correlation rho between en_music$rhythm_dev and en_music$pitch_dev 
# is negative, statistically significant, and small (rho = -0.15, S = 8.58e+11, p < .001)

corr_en <- ggplot(en_music) +
  aes(x = pitch_dev, y = rhythm_dev) +
  geom_point(size = .8) + #colour = "#0c4c8a"
  xlab("Pitch anticipation") + ylab("Rhythm synchronization") +
  geom_smooth(method=lm) + # for linear
  theme_gray(base_size = 12,
             base_family = "Times") +
  theme(legend.position = "bottom",
        legend.box = "vertical")

# annotate plot
# Create a text
grob <- grobTree(textGrob("rho = -0.15, S = .000, p < .001", x=0.45,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10))) #, fontface="italic"
# Plot
corr_en_annot <- corr_en + annotation_custom(grob)




# ASSUMPTION 1. Normality
ggplot(ma_music, aes(x = pitch_dev)) + 
  geom_density()

ggplot(ma_music, aes(x=rhythm_dev)) + 
  geom_density()


ad.test(ma_music$pitch_dev) # Anderson-Darling normality test for large samples
# A = 154.83, p-value < 2.2e-16

ad.test(ma_music$rhythm_dev) # neither normal > Spearman
# A = 617.8, p-value < 2.2e-16

# Spearman's correlation test
report(cor.test(ma_music$rhythm_dev, ma_music$pitch_dev, method = "spearman", exact=FALSE))
# The Spearman's rank correlation rho between ma_music$rhythm_dev and ma_music$pitch_dev 
# is negative, statistically significant, and small (rho = -0.26, S = 1.08e+12, p < .001)

corr_ma <- ggplot(ma_music) +
  aes(x = pitch_dev, y = rhythm_dev) +
  geom_point(size = .8) + #colour = "#0c4c8a"
  xlab("Pitch anticipation") + ylab("Rhythm synchronization") +
  geom_smooth(method=lm) + # for linear
  theme_gray(base_size = 12,
             base_family = "Times") +
  theme(legend.position = "bottom",
        legend.box = "vertical")

# annotate plot
grob <- grobTree(textGrob("rho = -0.26, S = .000, p < .001", x=0.45,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10))) #, fontface="italic"
# Plot
corr_ma_annot <- corr_ma + annotation_custom(grob)


figs_path <- here("figs", "music", "correlation")
ggsave(paste0(figs_path, "/correlation_mon.png"), corr_mon, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/correlation_en.png"), corr_en, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/correlation_ma.png"), corr_ma, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/correlation_mon_annotated.png"), corr_mon_annot, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/correlation_en_annotated.png"), corr_en_annot, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/correlation_ma_annotated.png"), corr_ma_annot, width = 180,
       height = 120, units = "mm", dpi = 600)



# #
# # only L2
# #
# 
# # Random effects
# 
# mod_ot1 <-
#   lmer(eLog ~ 1 + ot1 +
#          (1 + stress_sum + ot1 | participant),
#        control = lmerControl(optimizer = 'bobyqa'),
#        data = l2_music, weights = 1/wts, REML = F)
# 
# mod_ot2 <-
#   update(mod_ot1, . ~ . -(1 + stress_sum + ot1 | participant) +
#            ot2 + (1 + stress_sum + ot1 + ot2 | participant))
# 
# mod_ot3 <-
#   update(mod_ot2, . ~ . -(1 + stress_sum + ot1 + ot2 | participant) +
#            ot3 + (1 + stress_sum + ot1 + ot2 + ot3 | participant))
# 
# anova(mod_ot1, mod_ot2, mod_ot3)
# #           Df    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
# # mod_ot1    9 185950 186026 -92966   185932                          
# # mod_ot2   14 185872 185990 -92922   185844  87.736  5  < 2.2e-16 ***
# # mod_ot3   20 185459 185628 -92710   185419 424.968  6  < 2.2e-16 *** 
# 
# mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
# 
# mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
# 
# mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
#                     + (1 + ot1 + ot2 | target))
# 
# mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
#                     + (1 + ot1 + ot2 + ot3 | target))
# 
# anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
# #           npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
# #   mod_ot3   20 185459 185628 -92710   185419                          
# #   mod_ot4   21 184996 185173 -92477   184954 465.356  1  < 2.2e-16 ***
# #   mod_ot5   23 184615 184809 -92284   184569 385.320  2  < 2.2e-16 ***
# #   mod_ot6   26 184485 184704 -92217   184433 135.484  3  < 2.2e-16 ***
# #   mod_ot7   30 184435 184688 -92188   184375  57.797  4  8.418e-12 ***
# 
# # Fixed effects
# 
# 
# gca_mod_l2_base <- mod_ot7
#   # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
#   #        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
#   #        (1 + ot1 + ot2 + ot3 | target),
#   #      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
#   #      REML = F,
#   #      data = filter(l2_music, group == "aes")) 
# 
# # add proficiency synchronization effect to intercept, linear slope, quadratic, and cubic time terms
# gca_mod_l2_prof_0 <- update(gca_mod_l2_base,     . ~ . + prof_std) 
# gca_mod_l2_prof_1 <- update(gca_mod_l2_prof_0, . ~ . + ot1:prof_std) 
# gca_mod_l2_prof_2 <- update(gca_mod_l2_prof_1, . ~ . + ot2:prof_std) 
# gca_mod_l2_prof_3 <- update(gca_mod_l2_prof_2, . ~ . + ot3:prof_std) 
# 
# l2_prof_anova <-
#   anova(gca_mod_l2_base, gca_mod_l2_prof_0, gca_mod_l2_prof_1,
#         gca_mod_l2_prof_2, gca_mod_l2_prof_3)
# #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# # gca_mod_l2_base     30 184435 184688 -92188   184375                        
# # gca_mod_l2_prof_0   31 184433 184694 -92185   184371 4.8288  1   0.027988 * 
# # gca_mod_l2_prof_1   32 184427 184697 -92182   184363 7.3395  1   0.006746 **
# # gca_mod_l2_prof_2   33 184428 184706 -92181   184362 1.3292  1   0.248949   
# # gca_mod_l2_prof_3   34 184429 184716 -92181   184361 0.8535  1   0.355560
# 
# 
# # add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
# gca_mod_l2_cond_0 <- update(gca_mod_l2_prof_1,   . ~ . + stress_sum) 
# gca_mod_l2_cond_1 <- update(gca_mod_l2_cond_0, . ~ . + ot1:stress_sum) 
# gca_mod_l2_cond_2 <- update(gca_mod_l2_cond_1, . ~ . + ot2:stress_sum) 
# gca_mod_l2_cond_3 <- update(gca_mod_l2_cond_2, . ~ . + ot3:stress_sum) 
# 
# l2_cond_anova <-
#   anova(gca_mod_l2_prof_1, gca_mod_l2_cond_0, gca_mod_l2_cond_1,
#         gca_mod_l2_cond_2, gca_mod_l2_cond_3) 
# #                     Df    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# # gca_mod_aes_prof_1  32 184427 184697 -92182   184363                       
# # gca_mod_l2_cond_0   33 184425 184704 -92180   184359 3.8262  1    0.05046 .
# # gca_mod_l2_cond_1   34 184427 184714 -92180   184359 0.3315  1    0.56476  
# # gca_mod_l2_cond_2   35 184424 184719 -92177   184354 5.3822  1    0.02034 *
# # gca_mod_l2_cond_3   36 184425 184728 -92176   184353 0.8161  1    0.36633
# 
# # add L1 (English, Mandarin Chinese) effect to intercept, linear slope, quadratic, and cubic time terms
# gca_mod_l2_group_0 <- update(gca_mod_l2_cond_2,   . ~ . + l1_sum) 
# gca_mod_l2_group_1 <- update(gca_mod_l2_group_0, . ~ . + ot1:l1_sum) 
# gca_mod_l2_group_2 <- update(gca_mod_l2_group_1, . ~ . + ot2:l1_sum) 
# gca_mod_l2_group_3 <- update(gca_mod_l2_group_2, . ~ . + ot3:l1_sum) 
# 
# l2_group_anova <-
#   anova(gca_mod_l2_cond_2, gca_mod_l2_group_0, gca_mod_l2_group_1,
#         gca_mod_l2_group_2, gca_mod_l2_group_3)
# #                    npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# # gca_mod_l2_cond_2    35 184424 184719 -92177   184354                        
# # gca_mod_l2_group_0   36 184425 184729 -92177   184353 0.4401  1   0.507056   
# # gca_mod_l2_group_1   37 184418 184730 -92172   184344 9.0896  1   0.002571 **
# # gca_mod_l2_group_2   38 184420 184740 -92172   184344 0.3480  1   0.555259   
# # gca_mod_l2_group_3   39 184419 184747 -92170   184341 3.2606  1   0.070965 . 
# 
# 
# # BRANCH #1
# # add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
# gca_mod_l2_rhythm_0 <- update(gca_mod_l2_group_1,   . ~ . + rhythm_dev) 
# gca_mod_l2_rhythm_1 <- update(gca_mod_l2_rhythm_0, . ~ . + ot1:rhythm_dev) 
# gca_mod_l2_rhythm_2 <- update(gca_mod_l2_rhythm_1, . ~ . + ot2:rhythm_dev) 
# gca_mod_l2_rhythm_3 <- update(gca_mod_l2_rhythm_2, . ~ . + ot3:rhythm_dev) 
# 
# l2_rhythm_anova <-
#   anova(gca_mod_l2_group_1, gca_mod_l2_rhythm_0, gca_mod_l2_rhythm_1,
#         gca_mod_l2_rhythm_2, gca_mod_l2_rhythm_3)
# #                       Df    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# # gca_mod_l2_group_1    37 184418 184730 -92172   184344                       
# # gca_mod_l2_rhythm_0   38 184419 184740 -92172   184343 0.7631  1    0.38236  
# # gca_mod_l2_rhythm_1   39 184418 184747 -92170   184340 3.5278  1    0.06035 .
# # gca_mod_l2_rhythm_2   40 184420 184757 -92170   184340 0.0976  1    0.75476  
# # gca_mod_l2_rhythm_3   41 184421 184767 -92170   184339 0.5187  1    0.47142 


# mod_type <- "gca_mod_l2"
# mod_spec <- c("_base", 
#               "_cond_0", "_cond_1", "_cond_2", "_cond_3",
#               "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
#               # "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3",
#               "_prof_0", "_prof_1", "_prof_2", "_prof_3",
#               # "_rhythm_int_0", "_rhythm_int_1", "_rhythm_int_2", "_rhythm_int_3",
#               # "_pitch_int_0", "_pitch_int_1", "_pitch_int_2", "_pitch_int_3"
#               "_group_0", "_group_1", "_group_2", "_group_3"
#               )
# 
# # Store ind models in list
# l2_mods <- mget(c(paste0(mod_type, mod_spec)
# ))
# 
# save(l2_mods,
#      file = here("mods", "music", "gca", "continuous",
#                  "l2_mods.Rdata"))


}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_rhythm <- mon_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, rhythm_dev) %>%
  distinct

mon_pitch <- mon_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, pitch_dev) %>%
  distinct

en_rhythm <- en_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, rhythm_dev) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1)))

en_pitch <- en_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, pitch_dev) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1)))

ma_rhythm <- ma_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, rhythm_dev) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1)))

ma_pitch <- ma_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, pitch_dev) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1)))



# Get model predictions and SE
fits_mon_pitch <- predictSE(gca_mod_mon_pitch_0, mon_pitch) %>%        
  as_tibble %>%
  bind_cols(mon_pitch) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_mon_rhythm <- predictSE(gca_mod_mon_rhythm_int_1, mon_rhythm) %>%        
  as_tibble %>%
  bind_cols(mon_rhythm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_en_pitch <- predictSE(gca_mod_en_cond_2, en_pitch) %>%        
  as_tibble %>%
  bind_cols(en_pitch) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_en_rhythm <- predictSE(gca_mod_en_rhythm_int_1, en_rhythm) %>%        
  as_tibble %>%
  bind_cols(en_rhythm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ma_pitch <- predictSE(gca_mod_ma_pitch_int_2, ma_pitch) %>%        
  as_tibble %>%
  bind_cols(ma_pitch) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ma_rhythm <- predictSE(gca_mod_ma_rhythm_int_1, ma_rhythm) %>%        
  as_tibble %>%
  bind_cols(ma_rhythm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at syllable target onset
syll_onset_preds_mon_pitch <- filter(fits_mon_pitch, time_zero == 4) %>%
  select(stress_sum, pitch_dev,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_mon_rhythm <- filter(fits_mon_rhythm, time_zero == 4) %>%
  select(stress_sum, rhythm_dev,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_en_pitch <- filter(fits_en_pitch, time_zero == 4) %>%
  select(stress_sum, pitch_dev, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_en_rhythm <- filter(fits_en_rhythm, time_zero == 4) %>%
  select(stress_sum, rhythm_dev, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_ma_pitch <- filter(fits_ma_pitch, time_zero == 4) %>%
  select(stress_sum, pitch_dev, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_ma_rhythm <- filter(fits_ma_rhythm, time_zero == 4) %>%
  select(stress_sum, rhythm_dev, prof_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {

# Save anova model comparisons
nested_model_comparisons <-
  mget(c("mon_cond_anova", "mon_rhythm_anova", "mon_pitch_anova",
         "mon_rhythm_int_anova", "mon_pitch_int_anova",
         "en_cond_anova", "en_rhythm_anova", "en_pitch_anova",
         "en_prof_anova", "en_rhythm_int_anova", "en_pitch_int_anova",
         "ma_cond_anova", "ma_rhythm_anova", "ma_pitch_anova",
         "ma_prof_anova", "ma_rhythm_int_anova", "ma_pitch_int_anova"
         ))

save(nested_model_comparisons,
     file = here("mods", "music", "gca", "continuous",
                 "nested_model_comparisons.Rdata"))




# Save models predictions
model_preds <- mget(c("fits_mon_pitch", 'fits_mon_rhythm', 
                      "fits_en_pitch", 'fits_en_rhythm',
                      "fits_ma_pitch", 'fits_ma_rhythm',
                      "syll_onset_preds_mon_pitch", "syll_onset_preds_mon_rhythm",
                      "syll_onset_preds_en_pitch", "syll_onset_preds_en_rhythm",
                      "syll_onset_preds_ma_pitch", "syll_onset_preds_ma_rhythm"))

save(model_preds,
     file = here("mods", "music", "gca", "continuous",
                 "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------

