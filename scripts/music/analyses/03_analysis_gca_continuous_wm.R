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
stress50 <- read_csv(here("data", "clean", "stress_50ms_notstd.csv"))

# Get path to saved models
gca_mods_path  <- here("mods", "music", "gca", "continuous", "wm")
 
# Load models as lists
load(paste0(gca_mods_path, "/es_mods.Rdata"))
load(paste0(gca_mods_path, "/ma_mods.Rdata"))
load(paste0(gca_mods_path, "/en_mods.Rdata"))
# load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
# load(paste0(gca_mods_path, "/model_preds.Rdata"))
# 
# # Store objects in global env
list2env(es_mods, globalenv())
list2env(ma_mods, globalenv())
list2env(en_mods, globalenv())
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

# _dev = condval
# _sd = condsd


music50 <- music50 %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

ss_music <- filter(music50, l1 == 'es') %>% 
  select(-DELE, -percent_l2_week, -prof, -group, -cond) %>%
  mutate( ospan_std = (WM_set - mean(WM_set))/sd(WM_set))

ma_music <- music50 %>%
  filter(., l1 == 'ma') %>% 
  mutate(.,
         ospan_std = (WM_set - mean(WM_set))/sd(WM_set),
         prof_std = (DELE - mean(DELE))/sd(DELE),
         use_std = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week)) %>%
  select(-DELE, -percent_l2_week, -cond, -prof, -group)

en_music <- music50 %>%
  filter(., l1 == 'en') %>% 
  filter(., participant != 'ies04' & participant != 'ies17' & participant != 'ies28' & participant != 'aes32') %>%
  mutate(.,
         ospan_std = (WM_set - mean(WM_set))/sd(WM_set),
         prof_std = (DELE - mean(DELE))/sd(DELE),
         use_std = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week)) %>%
  select(-DELE, -percent_l2_week, -cond, -prof, -group)

en_music <- na.omit(en_music)
ma_music <- na.omit(ma_music)
ss_music <- na.omit(ss_music)

# -----------------------------------------------------------------------------







# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms




#
# only EN
#



# Random effects

if(F) {
  
  # mod_ot0 <-
  #   lmer(eLog ~ 1 + 
  #          (1 | participant),
  #        control = lmerControl(optimizer = 'bobyqa'),
  #        data = en_music, weights = 1/wts, REML = F)
  
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = en_music, weights = 1/wts, REML = F)
  
  # anova(mod_ot0, mod_ot1)
  # #         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # # mod_ot0    3 92772 92795 -46383    92766                         
  # # mod_ot1    6 90491 90537 -45239    90479 2287.1  3  < 2.2e-16 ***

  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + ot1 | participant) +
             ot2 + (1 + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + ot1 + ot2 | participant) +
             ot3 + (1 + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #           Df    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    6 90491 90537 -45239    90479                         
  # mod_ot2   10 90465 90542 -45223    90445  33.31  4  1.032e-06 ***
  # mod_ot3   15 90222 90338 -45096    90192 252.88  5  < 2.2e-16 ***
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                      + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                      + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #           npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  #   mod_ot3   15 90222 90338 -45096    90192                          
  #   mod_ot4   16 89847 89971 -44908    89815 377.149  1  < 2.2e-16 ***
  #   mod_ot5   18 89628 89766 -44796    89592 223.754  2  < 2.2e-16 ***
  #   mod_ot6   21 89500 89662 -44729    89458 133.087  3  < 2.2e-16 ***
  #   mod_ot7   25 89450 89642 -44700    89400  58.682  4  5.487e-12 ***

  
}




if(F) {


# Fixed effects

gca_mod_en_base <- mod_ot7
# lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
#        (1 + ot1 + ot2 + ot3 | participant) +
#        (1 + ot1 + ot2 + ot3 | target),
#      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
#      REML = F,
#      data = en_music

# add ospan effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_ospan_0 <- update(gca_mod_en_base,     . ~ . + ospan_std) 
gca_mod_en_ospan_1 <- update(gca_mod_en_ospan_0, . ~ . + ot1:ospan_std) 
gca_mod_en_ospan_2 <- update(gca_mod_en_ospan_1, . ~ . + ot2:ospan_std) 
gca_mod_en_ospan_3 <- update(gca_mod_en_ospan_2, . ~ . + ot3:ospan_std) 

en_ospan_anova <-
  anova(gca_mod_en_base, gca_mod_en_ospan_0, gca_mod_en_ospan_1,
        gca_mod_en_ospan_2, gca_mod_en_ospan_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mod_en_base    25 89450 89642 -44700    89400                     
# gca_mod_en_ospan_0   26 89451 89651 -44699    89399 1.1813  1     0.2771
# gca_mod_en_ospan_1   27 89452 89661 -44699    89398 0.1918  1     0.6614
# gca_mod_en_ospan_2   28 89452 89668 -44698    89396 2.2987  1     0.1295
# gca_mod_en_ospan_3   29 89453 89676 -44697    89395 1.2951  1     0.2551


# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_prof_0 <- update(gca_mod_en_base,   . ~ . + prof_std) 
gca_mod_en_prof_1 <- update(gca_mod_en_prof_0, . ~ . + ot1:prof_std) 
gca_mod_en_prof_2 <- update(gca_mod_en_prof_1, . ~ . + ot2:prof_std) 
gca_mod_en_prof_3 <- update(gca_mod_en_prof_2, . ~ . + ot3:prof_std) 

en_prof_anova <-
  anova(gca_mod_en_base, gca_mod_en_prof_0, gca_mod_en_prof_1,
        gca_mod_en_prof_2, gca_mod_en_prof_3) 
#                     Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_en_base     25 89450 89642 -44700    89400                       
# gca_mod_en_prof_0   26 89451 89652 -44700    89399 0.3898  1    0.53243  
# gca_mod_en_prof_1   27 89449 89657 -44697    89395 4.8645  1    0.02741 *
# gca_mod_en_prof_2   28 89445 89661 -44694    89389 5.7612  1    0.01638 *
# gca_mod_en_prof_3   29 89446 89669 -44694    89388 1.1488  1    0.28380  


# BRANCH #1
# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_rhythm_0 <- update(gca_mod_en_prof_2,   . ~ . + rhythm_dev) 
gca_mod_en_rhythm_1 <- update(gca_mod_en_rhythm_0, . ~ . + ot1:rhythm_dev) 
gca_mod_en_rhythm_2 <- update(gca_mod_en_rhythm_1, . ~ . + ot2:rhythm_dev) 
gca_mod_en_rhythm_3 <- update(gca_mod_en_rhythm_2, . ~ . + ot3:rhythm_dev) 

en_rhythm_anova <-
  anova(gca_mod_en_prof_2, gca_mod_en_rhythm_0, gca_mod_en_rhythm_1,
        gca_mod_en_rhythm_2, gca_mod_en_rhythm_3)
#                       Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_en_prof_2     28 89445 89661 -44694    89389                     
# gca_mod_en_rhythm_0   29 89445 89668 -44693    89387 2.1067  1     0.1467
# gca_mod_en_rhythm_1   30 89445 89676 -44692    89385 1.9314  1     0.1646
# gca_mod_en_rhythm_2   31 89445 89684 -44692    89383 1.4232  1     0.2329
# gca_mod_en_rhythm_3   32 89447 89693 -44691    89383 0.6623  1     0.4157


gca_mod_en_rhythm_int_0 <- update(gca_mod_en_prof_2,   . ~ . + prof_std:ospan_std:rhythm_dev) 
gca_mod_en_rhythm_int_1 <- update(gca_mod_en_rhythm_int_0, . ~ . + ot1:prof_std:ospan_std:rhythm_dev) 
gca_mod_en_rhythm_int_2 <- update(gca_mod_en_rhythm_int_1, . ~ . + ot2:prof_std:ospan_std:rhythm_dev) 
gca_mod_en_rhythm_int_3 <- update(gca_mod_en_rhythm_int_2, . ~ . + ot3:prof_std:ospan_std:rhythm_dev) 

en_rhythm_int_anova <-
  anova(gca_mod_en_prof_2, gca_mod_en_rhythm_int_0, gca_mod_en_rhythm_int_1,
        gca_mod_en_rhythm_int_2, gca_mod_en_rhythm_int_3)
#                         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_en_prof_2         28 89445 89661 -44694    89389                     
# gca_mod_en_rhythm_int_0   29 89447 89670 -44694    89389 0.0122  1     0.9122
# gca_mod_en_rhythm_int_1   30 89449 89680 -44694    89389 0.0169  1     0.8965
# gca_mod_en_rhythm_int_2   31 89451 89690 -44694    89389 0.0946  1     0.7585
# gca_mod_en_rhythm_int_3   32 89453 89699 -44694    89389 0.0897  1     0.7646


# BRANCH #2
# add pitch processing effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_pitch_0 <- update(gca_mod_en_prof_2,   . ~ . + pitch_dev) 
gca_mod_en_pitch_1 <- update(gca_mod_en_pitch_0, . ~ . + ot1:pitch_dev) 
gca_mod_en_pitch_2 <- update(gca_mod_en_pitch_1, . ~ . + ot2:pitch_dev) 
gca_mod_en_pitch_3 <- update(gca_mod_en_pitch_2, . ~ . + ot3:pitch_dev) 

en_pitch_anova <-
  anova(gca_mod_en_prof_2, gca_mod_en_pitch_0, gca_mod_en_pitch_1,
        gca_mod_en_pitch_2, gca_mod_en_pitch_3)
#                       Df    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_en_prof_2    28 89445 89661 -44694    89389                     
# gca_mod_en_pitch_0   29 89446 89670 -44694    89388 0.6630  1     0.4155
# gca_mod_en_pitch_1   30 89448 89679 -44694    89388 0.0490  1     0.8248
# gca_mod_en_pitch_2   31 89450 89689 -44694    89388 0.0150  1     0.9025
# gca_mod_en_pitch_3   32 89451 89698 -44694    89387 0.9721  1     0.3242


gca_mod_en_pitch_int_0 <- update(gca_mod_en_prof_2,   . ~ . + prof_std:ospan_std:pitch_dev) 
gca_mod_en_pitch_int_1 <- update(gca_mod_en_pitch_int_0, . ~ . + ot1:prof_std:ospan_std:pitch_dev) 
gca_mod_en_pitch_int_2 <- update(gca_mod_en_pitch_int_1, . ~ . + ot2:prof_std:ospan_std:pitch_dev) 
gca_mod_en_pitch_int_3 <- update(gca_mod_en_pitch_int_2, . ~ . + ot3:prof_std:ospan_std:pitch_dev) 

en_pitch_int_anova <-
  anova(gca_mod_en_prof_2, gca_mod_en_pitch_int_0, gca_mod_en_pitch_int_1,
        gca_mod_en_pitch_int_2, gca_mod_en_pitch_int_3)
#                        npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_en_prof_2        28 89445 89661 -44694    89389                     
# gca_mod_en_pitch_int_0   29 89446 89670 -44694    89388 0.5299  1     0.4667
# gca_mod_en_pitch_int_1   30 89448 89679 -44694    89388 0.0734  1     0.7864
# gca_mod_en_pitch_int_2   31 89450 89689 -44694    89388 0.0436  1     0.8345
# gca_mod_en_pitch_int_3   32 89452 89699 -44694    89388 0.1793  1     0.6720

gca_mod_en_final <- gca_mod_en_prof_2
summary(gca_mod_en_final) 
# Estimate Std. Error t value
# (Intercept)     1.24375    0.09893  12.572
# ot1           5.51130    0.33048  16.677
# ot2           0.31310    0.27191   1.151
# ot3          -1.34346    0.23293  -5.768
# prof_std      0.09631    0.05928   1.625
# ot1:prof_std  0.42496    0.21625   1.965
# ot2:prof_std -0.37666    0.15173  -2.482


mod_type <- "gca_mod_en"
mod_spec <- c("_base", 
              "_prof_0", "_prof_1", "_prof_2", "_prof_3",
              "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
              "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3",
              "_ospan_0", "_ospan_1", "_ospan_2", "_ospan_3",
              "_rhythm_int_0", "_rhythm_int_1", "_rhythm_int_2", "_rhythm_int_3",
              "_pitch_int_0", "_pitch_int_1", "_pitch_int_2", "_pitch_int_3",
              "_final")

# Store ind models in list
en_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(en_mods,
     file = here("mods", "music", "gca", "continuous", "wm",
                 "en_mods.Rdata"))




#
# only Ma CH
#



# Random effects

mod_ot1 <-
  lmer(eLog ~ 1 + ot1 +
         (1 + ot1 | participant),
       control = lmerControl(optimizer = 'bobyqa'),
       data = ma_music, weights = 1/wts, REML = F)

mod_ot2 <-
  update(mod_ot1, . ~ . -(1 + ot1 | participant) +
           ot2 + (1 + ot1 + ot2 | participant))

mod_ot3 <-
  update(mod_ot2, . ~ . -(1 + ot1 + ot2 | participant) +
           ot3 + (1 + ot1 + ot2 + ot3 | participant))

anova(mod_ot1, mod_ot2, mod_ot3)
#           Df    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
# mod_ot1    6 95581 95627 -47784    95569                          
# mod_ot2   10 95541 95619 -47761    95521  47.119  4  1.441e-09 ***
# mod_ot3   15 95370 95487 -47670    95340 181.147  5  < 2.2e-16 ***

mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))

mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))

mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                    + (1 + ot1 + ot2 | target))

mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                    + (1 + ot1 + ot2 + ot3 | target))

anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
#           npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
#   mod_ot3   15 95370 95487 -47670    95340                           
#   mod_ot4   16 95213 95337 -47591    95181 159.1892  1  < 2.2e-16 ***
#   mod_ot5   18 95034 95173 -47499    94998 183.4742  2  < 2.2e-16 ***
#   mod_ot6   21 94997 95159 -47477    94955  43.1130  3  2.329e-09 ***
#   mod_ot7   25 94995 95189 -47473    94945   9.3115  4    0.05377 .  



# Fixed effects

gca_mod_ma_base <- mod_ot6
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
  #        (1 + (ot1 + ot2 + ot3) | participant) +
  #        (1 + ot1 + ot2 | target),
  #      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)),
  #      REML = F,
  #      data = ma_music)

# add ospan effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_ospan_0 <- update(gca_mod_ma_base,     . ~ . + ospan_std) 
gca_mod_ma_ospan_1 <- update(gca_mod_ma_ospan_0, . ~ . + ot1:ospan_std) 
gca_mod_ma_ospan_2 <- update(gca_mod_ma_ospan_1, . ~ . + ot2:ospan_std) 
gca_mod_ma_ospan_3 <- update(gca_mod_ma_ospan_2, . ~ . + ot3:ospan_std) 

ma_ospan_anova <-
  anova(gca_mod_ma_base, gca_mod_ma_ospan_0, gca_mod_ma_ospan_1,
        gca_mod_ma_ospan_2, gca_mod_ma_ospan_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mod_ma_base    21 94997 95159 -47477    94955                     
# gca_mod_ma_ospan_0   22 94998 95169 -47477    94954 0.4118  1     0.5210
# gca_mod_ma_ospan_1   23 95000 95178 -47477    94954 0.4574  1     0.4988
# gca_mod_ma_ospan_2   24 95002 95188 -47477    94954 0.0596  1     0.8071
# gca_mod_ma_ospan_3   25 95003 95197 -47477    94953 0.3240  1     0.5692


# add profition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_prof_0 <- update(gca_mod_ma_base,   . ~ . + prof_std) 
gca_mod_ma_prof_1 <- update(gca_mod_ma_prof_0, . ~ . + ot1:prof_std) 
gca_mod_ma_prof_2 <- update(gca_mod_ma_prof_1, . ~ . + ot2:prof_std) 
gca_mod_ma_prof_3 <- update(gca_mod_ma_prof_2, . ~ . + ot3:prof_std) 

ma_prof_anova <-
  anova(gca_mod_ma_base, gca_mod_ma_prof_0, gca_mod_ma_prof_1,
        gca_mod_ma_prof_2, gca_mod_ma_prof_3) 
#                     Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ma_base     21 94997 95159 -47477    94955                       
# gca_mod_ma_prof_0   22 94993 95164 -47474    94949 5.6391  1    0.01756 *
# gca_mod_ma_prof_1   23 94992 95171 -47473    94946 2.5289  1    0.11178  
# gca_mod_ma_prof_2   24 94994 95180 -47473    94946 0.2035  1    0.65191  
# gca_mod_ma_prof_3   25 94995 95189 -47473    94945 0.7321  1    0.39220  

# BRANCH #1
# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_rhythm_0 <- update(gca_mod_ma_prof_0,   . ~ . + rhythm_dev) 
gca_mod_ma_rhythm_1 <- update(gca_mod_ma_rhythm_0, . ~ . + ot1:rhythm_dev) 
gca_mod_ma_rhythm_2 <- update(gca_mod_ma_rhythm_1, . ~ . + ot2:rhythm_dev) 
gca_mod_ma_rhythm_3 <- update(gca_mod_ma_rhythm_2, . ~ . + ot3:rhythm_dev) 

ma_rhythm_anova <-
  anova(gca_mod_ma_prof_0, gca_mod_ma_rhythm_0, gca_mod_ma_rhythm_1,
        gca_mod_ma_rhythm_2, gca_mod_ma_rhythm_3)
#                       Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ma_prof_0     22 94993 95164 -47474    94949                          
# gca_mod_ma_rhythm_0   23 94994 95173 -47474    94948  0.5309  1  0.4662078    
# gca_mod_ma_rhythm_1   24 94982 95169 -47467    94934 13.8844  1  0.0001944 ***
# gca_mod_ma_rhythm_2   25 94984 95178 -47467    94934  0.7787  1  0.3775443    
# gca_mod_ma_rhythm_3   26 94986 95187 -47467    94934  0.0113  1  0.9154395    


gca_mod_ma_rhythm_int_0 <- update(gca_mod_ma_rhythm_1,   . ~ . + ospan_std:prof_std:rhythm_dev) 
gca_mod_ma_rhythm_int_1 <- update(gca_mod_ma_rhythm_int_0, . ~ . + ot1:ospan_std:prof_std:rhythm_dev) 
gca_mod_ma_rhythm_int_2 <- update(gca_mod_ma_rhythm_int_1, . ~ . + ot2:ospan_std:prof_std:rhythm_dev) 
gca_mod_ma_rhythm_int_3 <- update(gca_mod_ma_rhythm_int_2, . ~ . + ot3:ospan_std:prof_std:rhythm_dev) 

ma_rhythm_int_anova <-
  anova(gca_mod_ma_rhythm_1, gca_mod_ma_rhythm_int_0, gca_mod_ma_rhythm_int_1,
        gca_mod_ma_rhythm_int_2, gca_mod_ma_rhythm_int_3)
#                         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ma_rhythm_1       24 94982 95169 -47467    94934                       
# gca_mod_ma_rhythm_int_0   25 94984 95178 -47467    94934 0.0256  1    0.87297  
# gca_mod_ma_rhythm_int_1   26 94986 95188 -47467    94934 0.5844  1    0.44461  
# gca_mod_ma_rhythm_int_2   27 94985 95195 -47466    94931 2.4710  1    0.11596  
# gca_mod_ma_rhythm_int_3   28 94982 95199 -47463    94926 5.2559  1    0.02187 *

# BRANCH #2
# add pitch processing effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_pitch_0 <- update(gca_mod_ma_prof_0,   . ~ . + pitch_dev) 
gca_mod_ma_pitch_1 <- update(gca_mod_ma_pitch_0, . ~ . + ot1:pitch_dev) 
gca_mod_ma_pitch_2 <- update(gca_mod_ma_pitch_1, . ~ . + ot2:pitch_dev) 
gca_mod_ma_pitch_3 <- update(gca_mod_ma_pitch_2, . ~ . + ot3:pitch_dev) 

ma_pitch_anova <-
  anova(gca_mod_ma_prof_0, gca_mod_ma_pitch_0, gca_mod_ma_pitch_1,
        gca_mod_ma_pitch_2, gca_mod_ma_pitch_3)
#                      Df   AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ma_prof_0    22 94993 95164 -47474    94949                     
# gca_mod_ma_pitch_0   23 94994 95173 -47474    94948 0.7631  1     0.3824
# gca_mod_ma_pitch_1   24 94995 95182 -47474    94947 0.6499  1     0.4201
# gca_mod_ma_pitch_2   25 94997 95191 -47474    94947 0.0319  1     0.8582
# gca_mod_ma_pitch_3   26 94998 95200 -47473    94946 1.5766  1     0.2093


gca_mod_ma_pitch_int_0 <- update(gca_mod_ma_prof_0,   . ~ . + ospan_std:prof_std:pitch_dev) 
gca_mod_ma_pitch_int_1 <- update(gca_mod_ma_pitch_int_0, . ~ . + ot1:ospan_std:prof_std:pitch_dev) 
gca_mod_ma_pitch_int_2 <- update(gca_mod_ma_pitch_int_1, . ~ . + ot2:ospan_std:prof_std:pitch_dev) 
gca_mod_ma_pitch_int_3 <- update(gca_mod_ma_pitch_int_2, . ~ . + ot3:ospan_std:prof_std:pitch_dev) 

ma_pitch_int_anova <-
  anova(gca_mod_ma_prof_0, gca_mod_ma_pitch_int_0, gca_mod_ma_pitch_int_1,
        gca_mod_ma_pitch_int_2, gca_mod_ma_pitch_int_3)
#                        npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ma_prof_0        22 94993 95164 -47474    94949                        
# gca_mod_ma_pitch_int_0   23 94987 95166 -47471    94941 7.5525  1   0.005993 **
# gca_mod_ma_pitch_int_1   24 94989 95176 -47471    94941 0.0004  1   0.984839   
# gca_mod_ma_pitch_int_2   25 94989 95183 -47469    94939 2.6607  1   0.102854   
# gca_mod_ma_pitch_int_3   26 94989 95191 -47468    94937 1.7452  1   0.186487   

gca_mod_ma_rhythm_final <- gca_mod_ma_rhythm_1
summary(gca_mod_ma_rhythm_final )
# Estimate Std. Error t value
# (Intercept)    0.96664    0.09144  10.572
# ot1             4.20911    0.29908  14.073
# ot2             0.27295    0.20926   1.304
# ot3            -0.99515    0.14457  -6.884
# prof_std        0.14285    0.05605   2.549
# rhythm_dev      0.91994    0.71949   1.279
# ot1:rhythm_dev  8.36919    2.10481   3.976
gca_mod_ma_pitch_final <- gca_mod_ma_pitch_int_0
summary(gca_mod_ma_pitch_final)
# Estimate Std. Error t value
# (Intercept)   0.95397    0.09061  10.528
# ot1                           4.22625    0.31499  13.417
# ot2                           0.26405    0.20929   1.262
# ot3                          -0.98878    0.14472  -6.832
# prof_std                      0.16820    0.05294   3.177
# prof_std:ospan_std:pitch_dev  0.54971    0.18896   2.909


mod_type <- "gca_mod_ma"
mod_spec <- c("_base", 
              "_ospan_0", "_ospan_1", "_ospan_2", "_ospan_3",
              "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
              "_prof_0", "_prof_1", "_prof_2", "_prof_3",
              "_rhythm_int_0", "_rhythm_int_1", "_rhythm_int_2", "_rhythm_int_3",
              "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3",
              "_pitch_int_0", "_pitch_int_1", "_pitch_int_2", "_pitch_int_3",
              "_pitch_final", "_rhythm_final")

# Store ind models in list
ma_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(ma_mods,
     file = here("mods", "music", "gca", "continuous", "wm",
                 "ma_mods.Rdata"))



}





if(F){

#
# only monolinguals
#



# Random effects

mod_ot1 <-
  lmer(eLog ~ 1 + ot1 +
         (1 + ot1 | participant),
       control = lmerControl(optimizer = 'bobyqa'),
       data = ss_music, weights = 1/wts, REML = F)

mod_ot2 <-
  update(mod_ot1, . ~ . -(1 + ot1 | participant) +
           ot2 + (1 + ot1 + ot2 | participant))

mod_ot3 <-
  update(mod_ot2, . ~ . -(1 + ot1 + ot2 | participant) +
           ot3 + (1 + ot1 + ot2 + ot3 | participant))

anova(mod_ot1, mod_ot2, mod_ot3)
#           Df    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
# mod_ot1    6 42893 42935 -21440    42881                         
# mod_ot2   10 42798 42868 -21389    42778 102.37  4  < 2.2e-16 ***
# mod_ot3   15 42690 42794 -21330    42660 118.67  5  < 2.2e-16 ***

mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))

mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))

mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                    + (1 + ot1 + ot2 | target))

mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                    + (1 + ot1 + ot2 + ot3 | target))

anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
#           npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
#   mod_ot3   15 42690 42794 -21330    42660                           
#   mod_ot4   16 42562 42673 -21265    42530 130.0408  1    < 2e-16 ***
#   mod_ot5   18 42371 42496 -21168    42335 194.7488  2    < 2e-16 ***
#   mod_ot6   21 42295 42442 -21127    42253  81.5610  3    < 2e-16 ***
#   mod_ot7   25 42294 42469 -21122    42244   9.0845  4    0.05902 .  



# Fixed effects

gca_mod_ss_base <- mod_ot6
# lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
#        (1 + (ot1 + ot2 + ot3) | participant) +
#        (1 + ot1 + ot2 | target),
#      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)),
#      REML = F,
#      data = ss_music)

# add ospan effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_ospan_0 <- update(gca_mod_ss_base,     . ~ . + ospan_std) 
gca_mod_ss_ospan_1 <- update(gca_mod_ss_ospan_0, . ~ . + ot1:ospan_std) 
gca_mod_ss_ospan_2 <- update(gca_mod_ss_ospan_1, . ~ . + ot2:ospan_std) 
gca_mod_ss_ospan_3 <- update(gca_mod_ss_ospan_2, . ~ . + ot3:ospan_std) 

ss_ospan_anova <-
  anova(gca_mod_ss_base, gca_mod_ss_ospan_0, gca_mod_ss_ospan_1,
        gca_mod_ss_ospan_2, gca_mod_ss_ospan_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mod_ss_base    21 42295 42442 -21127    42253                     
# gca_mod_ss_ospan_0   22 42297 42451 -21127    42253 0.0391  1     0.8433
# gca_mod_ss_ospan_1   23 42299 42460 -21127    42253 0.1099  1     0.7403
# gca_mod_ss_ospan_2   24 42301 42468 -21127    42253 0.0045  1     0.9463
# gca_mod_ss_ospan_3   25 42303 42477 -21126    42253 0.3138  1     0.5754


# BRANCH #1
# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_rhythm_0 <- update(gca_mod_ss_base,   . ~ . + rhythm_dev) 
gca_mod_ss_rhythm_1 <- update(gca_mod_ss_rhythm_0, . ~ . + ot1:rhythm_dev) 
gca_mod_ss_rhythm_2 <- update(gca_mod_ss_rhythm_1, . ~ . + ot2:rhythm_dev) 
gca_mod_ss_rhythm_3 <- update(gca_mod_ss_rhythm_2, . ~ . + ot3:rhythm_dev) 

ss_rhythm_anova <-
  anova(gca_mod_ss_base, gca_mod_ss_rhythm_0, gca_mod_ss_rhythm_1,
        gca_mod_ss_rhythm_2, gca_mod_ss_rhythm_3)
#                       Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ss_base       21 42295 42442 -21127    42253                     
# gca_mod_ss_rhythm_0   22 42297 42451 -21127    42253 0.2090  1     0.6476
# gca_mod_ss_rhythm_1   23 42299 42459 -21126    42253 0.2064  1     0.6496
# gca_mod_ss_rhythm_2   24 42301 42468 -21126    42253 0.1002  1     0.7516
# gca_mod_ss_rhythm_3   25 42303 42477 -21126    42253 0.2707  1     0.6029


gca_mod_ss_rhythm_int_0 <- update(gca_mod_ss_base,   . ~ . + ospan_std:rhythm_dev) 
gca_mod_ss_rhythm_int_1 <- update(gca_mod_ss_rhythm_int_0, . ~ . + ot1:ospan_std:rhythm_dev) 
gca_mod_ss_rhythm_int_2 <- update(gca_mod_ss_rhythm_int_1, . ~ . + ot2:ospan_std:rhythm_dev) 
gca_mod_ss_rhythm_int_3 <- update(gca_mod_ss_rhythm_int_2, . ~ . + ot3:ospan_std:rhythm_dev) 

ss_rhythm_int_anova <-
  anova(gca_mod_ss_base, gca_mod_ss_rhythm_int_0, gca_mod_ss_rhythm_int_1,
        gca_mod_ss_rhythm_int_2, gca_mod_ss_rhythm_int_3)
#                         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ss_base           21 42295 42442 -21127    42253                     
# gca_mod_ss_rhythm_int_0   22 42297 42451 -21127    42253 0.0374  1     0.8468
# gca_mod_ss_rhythm_int_1   23 42299 42459 -21126    42253 0.7712  1     0.3798
# gca_mod_ss_rhythm_int_2   24 42301 42468 -21126    42253 0.0000  1     0.9945
# gca_mod_ss_rhythm_int_3   25 42302 42477 -21126    42252 0.3022  1     0.5825

# BRANCH #2
# add pitch processing effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_pitch_0 <- update(gca_mod_ss_base,   . ~ . + pitch_dev) 
gca_mod_ss_pitch_1 <- update(gca_mod_ss_pitch_0, . ~ . + ot1:pitch_dev) 
gca_mod_ss_pitch_2 <- update(gca_mod_ss_pitch_1, . ~ . + ot2:pitch_dev) 
gca_mod_ss_pitch_3 <- update(gca_mod_ss_pitch_2, . ~ . + ot3:pitch_dev) 

ss_pitch_anova <-
  anova(gca_mod_ss_base, gca_mod_ss_pitch_0, gca_mod_ss_pitch_1,
        gca_mod_ss_pitch_2, gca_mod_ss_pitch_3)
#                      Df   AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ss_base      21 42295 42442 -21127    42253                        
# gca_mod_ss_pitch_0   22 42289 42443 -21123    42245 8.1963  1   0.004198 **
# gca_mod_ss_pitch_1   23 42290 42450 -21122    42244 1.1559  1   0.282321   
# gca_mod_ss_pitch_2   24 42290 42457 -21121    42242 2.4064  1   0.120842   
# gca_mod_ss_pitch_3   25 42290 42464 -21120    42240 2.1157  1   0.145792   


gca_mod_ss_pitch_int_0 <- update(gca_mod_ss_pitch_0,   . ~ . + ospan_std:pitch_dev) 
gca_mod_ss_pitch_int_1 <- update(gca_mod_ss_pitch_int_0, . ~ . + ot1:ospan_std:pitch_dev) 
gca_mod_ss_pitch_int_2 <- update(gca_mod_ss_pitch_int_1, . ~ . + ot2:ospan_std:pitch_dev) 
gca_mod_ss_pitch_int_3 <- update(gca_mod_ss_pitch_int_2, . ~ . + ot3:ospan_std:pitch_dev) 

ss_pitch_int_anova <-
  anova(gca_mod_ss_pitch_0, gca_mod_ss_pitch_int_0, gca_mod_ss_pitch_int_1,
        gca_mod_ss_pitch_int_2, gca_mod_ss_pitch_int_3)
#                        npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ss_pitch_0       22 42289 42443 -21123    42245                       
# gca_mod_ss_pitch_int_0   23 42289 42450 -21122    42243 1.8190  1    0.17743  
# gca_mod_ss_pitch_int_1   24 42290 42457 -21121    42242 1.2397  1    0.26553  
# gca_mod_ss_pitch_int_2   25 42287 42461 -21118    42237 5.3420  1    0.02082 *
# gca_mod_ss_pitch_int_3   26 42289 42470 -21118    42237 0.0955  1    0.75728  

gca_mod_ss_rhythm_final <- gca_mod_ss_base
summary(gca_mod_ss_rhythm_final )
# Estisste Std. Error t value
# (Intercept)    1.6118     0.1403  11.490
# ot1           4.5187     0.4705   9.604
# ot2          -1.1456     0.3762  -3.045
# ot3          -1.0043     0.2396  -4.191
gca_mod_ss_pitch_final <- gca_mod_ss_pitch_int_2
summary(gca_mod_ss_pitch_final)
# Estisste Std. Error t value
# (Intercept)   1.6590     0.1411  11.760
# ot1                       4.5222     0.4646   9.734
# ot2                      -1.1334     0.3582  -3.164
# ot3                      -0.9984     0.2393  -4.172
# pitch_dev                 0.8855     0.2427   3.649
# pitch_dev:ospan_std       0.8213     0.3090   2.658
# ot1:pitch_dev:ospan_std   1.0369     0.7376   1.406
# ot2:pitch_dev:ospan_std  -2.0813     0.8489  -2.452


mod_type <- "gca_mod_ss"
mod_spec <- c("_base", 
              "_ospan_0", "_ospan_1", "_ospan_2", "_ospan_3",
              "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
              "_rhythm_int_0", "_rhythm_int_1", "_rhythm_int_2", "_rhythm_int_3",
              "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3",
              "_pitch_int_0", "_pitch_int_1", "_pitch_int_2", "_pitch_int_3",
              "_pitch_final", "_rhythm_final")

# Store ind models in list
ss_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(ss_mods,
     file = here("mods", "music", "gca", "continuous", "wm",
                 "ss_mods.Rdata"))



}
# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions

en_rhythm <- en_music %>%
  dplyr::select(time_zero, ot1:ot3, rhythm_dev) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan_std = c(-1, 0, 1)))

en_pitch <- en_music %>%
  dplyr::select(time_zero, ot1:ot3, pitch_dev) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan_std = c(-1, 0, 1)))

ma_rhythm <- ma_music %>%
  dplyr::select(time_zero, ot1:ot3, rhythm_dev) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan_std = c(-1, 0, 1)))

ma_pitch <- ma_music %>%
  dplyr::select(time_zero, ot1:ot3, pitch_dev) %>%
  distinct %>%
  expand_grid(., tibble(prof_std = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan_std = c(-1, 0, 1)))





ss_rhythm <- ss_music %>%
  dplyr::select(time_zero, ot1:ot3, rhythm_dev) %>%
  distinct %>%
  expand_grid(., tibble(
                        ospan_std = c(-1, 0, 1)))

ss_pitch <- ss_music %>%
  dplyr::select(time_zero, ot1:ot3, pitch_dev) %>%
  distinct %>%
  expand_grid(., tibble(
                        ospan_std = c(-1, 0, 1)))



# Get model predictions and SE

fits_en_pitch <- predictSE(gca_mod_en_final, en_pitch) %>%        
  as_tibble %>%
  bind_cols(en_pitch) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_en_rhythm <- predictSE(gca_mod_en_final, en_rhythm) %>%        
  as_tibble %>%
  bind_cols(en_rhythm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


fits_ma_pitch <- predictSE(gca_mod_ma_pitch_final, ma_pitch) %>%        
  as_tibble %>%
  bind_cols(ma_pitch) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ma_rhythm <- predictSE(gca_mod_ma_rhythm_final, ma_rhythm) %>%        
  as_tibble %>%
  bind_cols(ma_rhythm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


fits_ss_pitch <- predictSE(gca_mod_ss_pitch_final, ss_pitch) %>%        
  as_tibble %>%
  bind_cols(ss_pitch) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ss_rhythm <- predictSE(gca_mod_ss_rhythm_final, ss_rhythm) %>%        
  as_tibble %>%
  bind_cols(ss_rhythm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at syllable target onset

syll_onset_preds_en_pitch <- filter(fits_en_pitch, time_zero == 4) %>%
  select(pitch_dev, prof_std, ospan_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_en_rhythm <- filter(fits_en_rhythm, time_zero == 4) %>%
  select(rhythm_dev, prof_std, ospan_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_ma_pitch <- filter(fits_ma_pitch, time_zero == 4) %>%
  select(pitch_dev, prof_std, ospan_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_ma_rhythm <- filter(fits_ma_rhythm, time_zero == 4) %>%
  select(rhythm_dev, prof_std, ospan_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 



syll_onset_preds_ss_pitch <- filter(fits_ss_pitch, time_zero == 4) %>%
  select(pitch_dev, ospan_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_ss_rhythm <- filter(fits_ss_rhythm, time_zero == 4) %>%
  select(rhythm_dev, ospan_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {

# Save anova model comparisons
nested_model_comparisons <-
  mget(c(
         "en_prof_anova", "en_rhythm_anova", "en_pitch_anova",
         "en_ospan_anova", "en_rhythm_int_anova", "en_pitch_int_anova",
         "ma_prof_anova", "ma_rhythm_anova", "ma_pitch_anova",
         "ma_ospan_anova", "ma_rhythm_int_anova", "ma_pitch_int_anova",
         "ss_rhythm_anova", "ss_pitch_anova",
         "ss_ospan_anova", "ss_rhythm_int_anova", "ss_pitch_int_anova"
         ))

save(nested_model_comparisons,
     file = here("mods", "music", "gca", "continuous", "wm",
                 "nested_model_comparisons.Rdata"))




# Save models predictions
model_preds <- mget(c(
                      "fits_en_pitch", 'fits_en_rhythm',
                      "fits_ma_pitch", 'fits_ma_rhythm',
                      "syll_onset_preds_en_pitch", "syll_onset_preds_en_rhythm",
                      "syll_onset_preds_ma_pitch", "syll_onset_preds_ma_rhythm",
                      "fits_ss_pitch", 'fits_ss_rhythm',
                      "syll_onset_preds_ss_pitch", "syll_onset_preds_ss_rhythm"))

save(model_preds,
     file = here("mods", "music", "gca", "continuous", "wm",
                 "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------

