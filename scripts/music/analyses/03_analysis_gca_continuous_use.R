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
gca_mods_path  <- here("mods", "music", "gca", "continuous", "use")
 
# Load models as lists
load(paste0(gca_mods_path, "/ma_mods.Rdata"))
load(paste0(gca_mods_path, "/en_mods.Rdata"))
# load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
# load(paste0(gca_mods_path, "/model_preds.Rdata"))
# 
# # Store objects in global env
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

music50 <- na.omit(music50)

# _dev = condval
# _sd = condsd


music50 <- music50 %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "mon", "aes", "ams", "ies", "ims"),
            stress_sum = if_else(cond == "1", 1, -1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# mon_music <- filter(music50, l1 == 'es') %>% select(-DELE, -percent_l2_week, 
#                                                     -DELE_z, -use_z, -prof, -group)
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




#
# only EN
#

en_music <- l2_music %>%
  filter(., l1 == 'en')

# Random effects

if(F) {
  
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

  
}




if(F) {


# Fixed effects

gca_mod_en_base <- mod_ot7
# lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
#        (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
#        (1 + ot1 + ot2 + ot3 | target),
#      control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
#      REML = F,
#      data = en_music

# add use effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_use_0 <- update(gca_mod_en_base,     . ~ . + use_std) 
gca_mod_en_use_1 <- update(gca_mod_en_use_0, . ~ . + ot1:use_std) 
gca_mod_en_use_2 <- update(gca_mod_en_use_1, . ~ . + ot2:use_std) 
gca_mod_en_use_3 <- update(gca_mod_en_use_2, . ~ . + ot3:use_std) 

en_use_anova <-
  anova(gca_mod_en_base, gca_mod_en_use_0, gca_mod_en_use_1,
        gca_mod_en_use_2, gca_mod_en_use_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mod_en_base    30 89382 89613 -44661    89322                       
# gca_mod_en_use_0   31 89384 89623 -44661    89322 0.5852  1    0.44428  
# gca_mod_en_use_1   32 89382 89629 -44659    89318 3.6959  1    0.05455 .
# gca_mod_en_use_2   33 89384 89638 -44659    89318 0.1008  1    0.75088  
# gca_mod_en_use_3   34 89382 89644 -44657    89314 4.0585  1    0.04395 *


# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_en_cond_0 <- update(gca_mod_en_use_3,   . ~ . + stress_sum) 
gca_mod_en_cond_1 <- update(gca_mod_en_cond_0, . ~ . + ot1:stress_sum) 
gca_mod_en_cond_2 <- update(gca_mod_en_cond_1, . ~ . + ot2:stress_sum) 
gca_mod_en_cond_3 <- update(gca_mod_en_cond_2, . ~ . + ot3:stress_sum) 

en_cond_anova <-
  anova(gca_mod_en_use_3, gca_mod_en_cond_0, gca_mod_en_cond_1,
        gca_mod_en_cond_2, gca_mod_en_cond_3) 
#                     Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_en_use_3    34 89382 89644 -44657    89314                       
# gca_mod_en_cond_0   35 89383 89653 -44657    89313 0.5247  1    0.46885  
# gca_mod_en_cond_1   36 89385 89662 -44656    89313 0.4700  1    0.49299  
# gca_mod_en_cond_2   37 89381 89666 -44653    89307 6.2704  1    0.01228 *
# gca_mod_en_cond_3   38 89381 89674 -44653    89305 1.1365  1    0.28639  


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
# gca_mod_en_cond_2     37 89381 89666 -44653    89307                       
# gca_mod_en_rhythm_0   38 89380 89673 -44652    89304 2.0743  1    0.14980  
# gca_mod_en_rhythm_1   39 89381 89682 -44652    89303 1.3611  1    0.24334  
# gca_mod_en_rhythm_2   40 89383 89691 -44651    89303 0.5565  1    0.45567  
# gca_mod_en_rhythm_3   41 89381 89697 -44650    89299 3.3219  1    0.06836 .


gca_mod_en_rhythm_int_0 <- update(gca_mod_en_cond_2,   . ~ . + use_std:stress_sum:rhythm_dev) 
gca_mod_en_rhythm_int_1 <- update(gca_mod_en_rhythm_int_0, . ~ . + ot1:use_std:stress_sum:rhythm_dev) 
gca_mod_en_rhythm_int_2 <- update(gca_mod_en_rhythm_int_1, . ~ . + ot2:use_std:stress_sum:rhythm_dev) 
gca_mod_en_rhythm_int_3 <- update(gca_mod_en_rhythm_int_2, . ~ . + ot3:use_std:stress_sum:rhythm_dev) 

en_rhythm_int_anova <-
  anova(gca_mod_en_cond_2, gca_mod_en_rhythm_int_0, gca_mod_en_rhythm_int_1,
        gca_mod_en_rhythm_int_2, gca_mod_en_rhythm_int_3)
#                         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_en_cond_2         37 89381 89666 -44653    89307                       
# gca_mod_en_rhythm_int_0   38 89382 89675 -44653    89306 0.1475  1    0.70096  
# gca_mod_en_rhythm_int_1   39 89384 89685 -44653    89306 0.1902  1    0.66273  
# gca_mod_en_rhythm_int_2   40 89383 89692 -44652    89303 2.7680  1    0.09616 .
# gca_mod_en_rhythm_int_3   41 89385 89701 -44652    89303 0.1934  1    0.66009  


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
# gca_mod_en_cond_2    37 89381 89666 -44653    89307                     
# gca_mod_en_pitch_0   38 89382 89675 -44653    89306 0.5639  1     0.4527
# gca_mod_en_pitch_1   39 89384 89685 -44653    89306 0.0042  1     0.9482
# gca_mod_en_pitch_2   40 89386 89694 -44653    89306 0.0323  1     0.8575
# gca_mod_en_pitch_3   41 89387 89703 -44652    89305 0.9840  1     0.3212


gca_mod_en_pitch_int_0 <- update(gca_mod_en_cond_2,   . ~ . + use_std:stress_sum:pitch_dev) 
gca_mod_en_pitch_int_1 <- update(gca_mod_en_pitch_int_0, . ~ . + ot1:use_std:stress_sum:pitch_dev) 
gca_mod_en_pitch_int_2 <- update(gca_mod_en_pitch_int_1, . ~ . + ot2:use_std:stress_sum:pitch_dev) 
gca_mod_en_pitch_int_3 <- update(gca_mod_en_pitch_int_2, . ~ . + ot3:use_std:stress_sum:pitch_dev) 

en_pitch_int_anova <-
  anova(gca_mod_en_cond_2, gca_mod_en_pitch_int_0, gca_mod_en_pitch_int_1,
        gca_mod_en_pitch_int_2, gca_mod_en_pitch_int_3)
#                        npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_en_cond_2        37 89381 89666 -44653    89307                         
# gca_mod_en_pitch_int_0   38 89382 89675 -44653    89306  0.2543  1    0.61409   
# gca_mod_en_pitch_int_1   39 89374 89674 -44648    89296 10.7193  1    0.00106 **
# gca_mod_en_pitch_int_2   40 89373 89682 -44647    89293  2.0741  1    0.14982   
# gca_mod_en_pitch_int_3   41 89374 89690 -44646    89292  1.3067  1    0.25299   

summary(gca_mod_en_cond_2) # no effect of rhythm
# Estimate Std. Error t value
# (Intercept)                        1.26994    0.10119  12.550
# ot1             5.63450    0.33727  16.706
# ot2             0.27851    0.27639   1.008
# ot3            -1.43052    0.22740  -6.291
# use_std         0.06684    0.07696   0.869
# stress_sum     -0.03964    0.07490  -0.529
# ot1:use_std     0.66527    0.26524   2.508
# ot2:use_std    -0.06839    0.17901  -0.382
# ot3:use_std    -0.35861    0.17271  -2.076
# ot1:stress_sum  0.20756    0.23356   0.889
# ot2:stress_sum  0.48345    0.17867   2.706
summary(gca_mod_en_pitch_int_1) 
# Estimate Std. Error t value
# (Intercept)     1.26917    0.10118  12.544
# ot1                            5.63847    0.33693  16.735
# ot2                            0.28063    0.27645   1.015
# ot3                           -1.43402    0.22762  -6.300
# use_std                        0.06670    0.07667   0.870
# stress_sum                    -0.03998    0.07417  -0.539
# ot1:use_std                    0.68100    0.26393   2.580
# ot2:use_std                   -0.06338    0.17944  -0.353
# ot3:use_std                   -0.36205    0.17230  -2.101
# ot1:stress_sum                 0.20354    0.23322   0.873
# ot2:stress_sum                 0.48260    0.17855   2.703
# stress_sum:prof_std:pitch_dev -0.27328    0.12222  -2.236


mod_type <- "gca_mod_en"
mod_spec <- c("_base", 
              "_cond_0", "_cond_1", "_cond_2", "_cond_3",
              "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
              "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3",
              "_use_0", "_use_1", "_use_2", "_use_3",
              "_rhythm_int_0", "_rhythm_int_1", "_rhythm_int_2", "_rhythm_int_3",
              "_pitch_int_0", "_pitch_int_1", "_pitch_int_2", "_pitch_int_3")

# Store ind models in list
en_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(en_mods,
     file = here("mods", "music", "gca", "continuous", "use",
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

gca_mod_ma_base <- 
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + stress_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 | target),
       control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)),
       REML = F,
       data = ma_music)

# add L2 use effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_use_0 <- update(gca_mod_ma_base,     . ~ . + use_std) 
gca_mod_ma_use_1 <- update(gca_mod_ma_use_0, . ~ . + ot1:use_std) 
gca_mod_ma_use_2 <- update(gca_mod_ma_use_1, . ~ . + ot2:use_std) 
gca_mod_ma_use_3 <- update(gca_mod_ma_use_2, . ~ . + ot3:use_std) 

ma_use_anova <-
  anova(gca_mod_ma_base, gca_mod_ma_use_0, gca_mod_ma_use_1,
        gca_mod_ma_use_2, gca_mod_ma_use_3)
#                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mod_ma_base     26 94734 94936 -47341    94682                        
# gca_mod_ma_use_0   27 94736 94946 -47341    94682 0.0115  1   0.914733   
# gca_mod_ma_use_1   28 94736 94953 -47340    94680 2.2930  1   0.129962   
# gca_mod_ma_use_2   29 94730 94955 -47336    94672 7.5881  1   0.005875 **
# gca_mod_ma_use_3   30 94732 94964 -47336    94672 0.5483  1   0.459033   


# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_cond_0 <- update(gca_mod_ma_use_2,   . ~ . + stress_sum) 
gca_mod_ma_cond_1 <- update(gca_mod_ma_cond_0, . ~ . + ot1:stress_sum) 
gca_mod_ma_cond_2 <- update(gca_mod_ma_cond_1, . ~ . + ot2:stress_sum) 
gca_mod_ma_cond_3 <- update(gca_mod_ma_cond_2, . ~ . + ot3:stress_sum) 

ma_cond_anova <-
  anova(gca_mod_ma_use_2, gca_mod_ma_cond_0, gca_mod_ma_cond_1,
        gca_mod_ma_cond_2, gca_mod_ma_cond_3) 
#                     Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ma_use_2    29 94730 94955 -47336    94672                       
# gca_mod_ma_cond_0   30 94731 94964 -47335    94671 1.2654  1    0.26062  
# gca_mod_ma_cond_1   31 94733 94973 -47335    94671 0.0000  1    1.00000  
# gca_mod_ma_cond_2   32 94734 94982 -47335    94670 1.2609  1    0.26148  
# gca_mod_ma_cond_3   33 94732 94988 -47333    94666 3.8038  1    0.05114 .


# BRANCH #1
# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_rhythm_0 <- update(gca_mod_ma_use_2,   . ~ . + rhythm_dev) 
gca_mod_ma_rhythm_1 <- update(gca_mod_ma_rhythm_0, . ~ . + ot1:rhythm_dev) 
gca_mod_ma_rhythm_2 <- update(gca_mod_ma_rhythm_1, . ~ . + ot2:rhythm_dev) 
gca_mod_ma_rhythm_3 <- update(gca_mod_ma_rhythm_2, . ~ . + ot3:rhythm_dev) 

ma_rhythm_anova <-
  anova(gca_mod_ma_use_2, gca_mod_ma_rhythm_0, gca_mod_ma_rhythm_1,
        gca_mod_ma_rhythm_2, gca_mod_ma_rhythm_3)
#                       Df   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ma_use_2      29 94730 94955 -47336    94672                         
# gca_mod_ma_rhythm_0   30 94732 94965 -47336    94672  0.2932  1   0.588146   
# gca_mod_ma_rhythm_1   31 94724 94964 -47331    94662 10.2251  1   0.001385 **
# gca_mod_ma_rhythm_2   32 94726 94974 -47331    94662  0.0239  1   0.877237   
# gca_mod_ma_rhythm_3   33 94728 94983 -47331    94662  0.1599  1   0.689242     


gca_mod_ma_rhythm_int_0 <- update(gca_mod_ma_rhythm_1,   . ~ . + use_std:stress_sum:rhythm_dev) 
gca_mod_ma_rhythm_int_1 <- update(gca_mod_ma_rhythm_int_0, . ~ . + ot1:use_std:stress_sum:rhythm_dev) 
gca_mod_ma_rhythm_int_2 <- update(gca_mod_ma_rhythm_int_1, . ~ . + ot2:use_std:stress_sum:rhythm_dev) 
gca_mod_ma_rhythm_int_3 <- update(gca_mod_ma_rhythm_int_2, . ~ . + ot3:use_std:stress_sum:rhythm_dev) 

ma_rhythm_int_anova <-
  anova(gca_mod_ma_rhythm_1, gca_mod_ma_rhythm_int_0, gca_mod_ma_rhythm_int_1,
        gca_mod_ma_rhythm_int_2, gca_mod_ma_rhythm_int_3)
#                         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ma_rhythm_1       31 94724 94964 -47331    94662                     
# gca_mod_ma_rhythm_int_0   32 94724 94972 -47330    94660 1.5572  1     0.2121
# gca_mod_ma_rhythm_int_1   33 94726 94982 -47330    94660 0.1049  1     0.7460
# gca_mod_ma_rhythm_int_2   34 94728 94991 -47330    94660 0.2970  1     0.5858
# gca_mod_ma_rhythm_int_3   35 94729 95001 -47330    94659 0.3018  1     0.5828

# BRANCH #2
# add pitch processing effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ma_pitch_0 <- update(gca_mod_ma_use_2,   . ~ . + pitch_dev) 
gca_mod_ma_pitch_1 <- update(gca_mod_ma_pitch_0, . ~ . + ot1:pitch_dev) 
gca_mod_ma_pitch_2 <- update(gca_mod_ma_pitch_1, . ~ . + ot2:pitch_dev) 
gca_mod_ma_pitch_3 <- update(gca_mod_ma_pitch_2, . ~ . + ot3:pitch_dev) 

ma_pitch_anova <-
  anova(gca_mod_ma_use_2, gca_mod_ma_pitch_0, gca_mod_ma_pitch_1,
        gca_mod_ma_pitch_2, gca_mod_ma_pitch_3)
#                       Df    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_mod_ma_29 94730 94955 -47336    94672                        
# gca_mod_ma_pitch_0   30 94725 94958 -47332    94665 7.3842  1    0.00658 **
# gca_mod_ma_pitch_1   31 94727 94967 -47332    94665 0.1840  1    0.66797   
# gca_mod_ma_pitch_2   32 94728 94976 -47332    94664 0.6190  1    0.43141   
# gca_mod_ma_pitch_3   33 94730 94986 -47332    94664 0.0111  1    0.91614   


gca_mod_ma_pitch_int_0 <- update(gca_mod_ma_pitch_0,   . ~ . + use_std:stress_sum:pitch_dev) 
gca_mod_ma_pitch_int_1 <- update(gca_mod_ma_pitch_int_0, . ~ . + ot1:use_std:stress_sum:pitch_dev) 
gca_mod_ma_pitch_int_2 <- update(gca_mod_ma_pitch_int_1, . ~ . + ot2:use_std:stress_sum:pitch_dev) 
gca_mod_ma_pitch_int_3 <- update(gca_mod_ma_pitch_int_2, . ~ . + ot3:use_std:stress_sum:pitch_dev) 

ma_pitch_int_anova <-
  anova(gca_mod_ma_pitch_0, gca_mod_ma_pitch_int_0, gca_mod_ma_pitch_int_1,
        gca_mod_ma_pitch_int_2, gca_mod_ma_pitch_int_3)
#                        npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mod_ma_pitch_0       30 94725 94958 -47332    94665                     
# gca_mod_ma_pitch_int_0   31 94726 94966 -47332    94664 1.1783  1     0.2777
# gca_mod_ma_pitch_int_1   32 94728 94976 -47332    94664 0.0008  1     0.9778
# gca_mod_ma_pitch_int_2   33 94728 94984 -47331    94662 1.2916  1     0.2558
# gca_mod_ma_pitch_int_3   34 94730 94994 -47331    94662 0.2563  1     0.6127

summary( gca_mod_ma_rhythm_1 )
# Estimate Std. Error t value
# (Intercept)                         1.44807    0.13215  10.958
# ot1             6.20767    0.36706  16.912
# ot2            -0.66343    0.26259  -2.527
# ot3            -1.55002    0.13135 -11.801
# use_std         0.08338    0.09830   0.848
# rhythm_dev      1.01193    1.18299   0.855
# ot1:use_std     0.40895    0.21812   1.875
# ot2:use_std    -0.50130    0.17569  -2.853
# ot1:rhythm_dev  8.90402    2.67238   3.332
summary(gca_mod_ma_pitch_0)
# Estimate Std. Error t value
# (Intercept)                        1.54285    0.13134  11.747
# ot1          6.21785    0.38226  16.266
# ot2         -0.66297    0.26285  -2.522
# ot3         -1.54927    0.13062 -11.861
# use_std      0.11355    0.09289   1.222
# pitch_dev   -1.09441    0.38743  -2.825
# ot1:use_std  0.43515    0.23560   1.847
# ot2:use_std -0.50147    0.17570  -2.854


mod_type <- "gca_mod_ma"
mod_spec <- c("_base", 
              "_cond_0", "_cond_1", "_cond_2", "_cond_3",
              "_rhythm_0", "_rhythm_1", "_rhythm_2", "_rhythm_3",
              "_use_0", "_use_1", "_use_2", "_use_3",
              "_rhythm_int_0", "_rhythm_int_1", "_rhythm_int_2", "_rhythm_int_3",
              "_pitch_0", "_pitch_1", "_pitch_2", "_pitch_3",
              "_pitch_int_0", "_pitch_int_1", "_pitch_int_2", "_pitch_int_3"
              )

# Store ind models in list
ma_mods <- mget(c(paste0(mod_type, mod_spec)
))

save(ma_mods,
     file = here("mods", "music", "gca", "continuous", "use",
                 "ma_mods.Rdata"))



}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions

en_rhythm <- en_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, rhythm_dev) %>%
  distinct %>%
  expand_grid(., tibble(use_std = c(-1, 0, 1)))

en_pitch <- en_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, pitch_dev) %>%
  distinct %>%
  expand_grid(., tibble(use_std = c(-1, 0, 1)))

ma_rhythm <- ma_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, rhythm_dev) %>%
  distinct %>%
  expand_grid(., tibble(use_std = c(-1, 0, 1)))

ma_pitch <- ma_music %>%
  dplyr::select(time_zero, ot1:ot3, stress_sum, pitch_dev) %>%
  distinct %>%
  expand_grid(., tibble(use_std = c(-1, 0, 1)))



# Get model predictions and SE

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

fits_ma_pitch <- predictSE(gca_mod_ma_pitch_0, ma_pitch) %>%        
  as_tibble %>%
  bind_cols(ma_pitch) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

fits_ma_rhythm <- predictSE(gca_mod_ma_rhythm_1, ma_rhythm) %>%        
  as_tibble %>%
  bind_cols(ma_rhythm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at syllable target onset

syll_onset_preds_en_pitch <- filter(fits_en_pitch, time_zero == 4) %>%
  select(stress_sum, pitch_dev, use_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_en_rhythm <- filter(fits_en_rhythm, time_zero == 4) %>%
  select(stress_sum, rhythm_dev, use_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_ma_pitch <- filter(fits_ma_pitch, time_zero == 4) %>%
  select(stress_sum, pitch_dev, use_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 

syll_onset_preds_ma_rhythm <- filter(fits_ma_rhythm, time_zero == 4) %>%
  select(stress_sum, rhythm_dev, use_std,
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
         "en_cond_anova", "en_rhythm_anova", "en_pitch_anova",
         "en_use_anova", "en_rhythm_int_anova", "en_pitch_int_anova",
         "ma_cond_anova", "ma_rhythm_anova", "ma_pitch_anova",
         "ma_use_anova", "ma_rhythm_int_anova", "ma_pitch_int_anova"
         ))

save(nested_model_comparisons,
     file = here("mods", "music", "gca", "continuous", "use",
                 "nested_model_comparisons.Rdata"))




# Save models predictions
model_preds <- mget(c(
                      "fits_en_pitch", 'fits_en_rhythm',
                      "fits_ma_pitch", 'fits_ma_rhythm',
                      "syll_onset_preds_en_pitch", "syll_onset_preds_en_rhythm",
                      "syll_onset_preds_ma_pitch", "syll_onset_preds_ma_rhythm"))

save(model_preds,
     file = here("mods", "music", "gca", "continuous", "use",
                 "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------

