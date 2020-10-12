#
# Original script by Joseph
# Modified by Cristina to adapt to CrisLau Project
# Last update: 06/12/2019
# Modified by Laura to adapt to Pupurri project
#
# Growth curve analyisis ------------------------------------------------------
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
#     - Rhythm and stress are associated in L2 speakers (the)
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

# Load data
source(here::here("scripts", "02_load_data.R"))

# # Get path to saved models
# gca_mods_path  <- here("reports", "mods", "gca")
# 
# # Load models as lists
# #load(paste0(gca_mods_path, "/ind_mods.Rdata"))
# load(paste0(gca_mods_path, "/full_mods.Rdata"))
# load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
# #load(paste0(gca_mods_path, "/model_preds.Rdata"))
# 
# # Store objects in global env
# list2env(ind_mods, globalenv())
# list2env(full_mods, globalenv())
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

stress50 <- stress50 %>% 
  rename(., linx_stress = cond)

auditory <- read_csv("./data/clean/auditory_scores.csv")

auditory <- auditory %>%
  select(., -X1)

music50 <- merge(x = stress50, y = auditory, by = "participant", all.x=TRUE)
stress50 <- na.omit(stress50)




stress_gc_subset <- music50 %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "mon", "aes", "ams", "ies", "ims"),
            stress_sum = if_else(linx_stress == "1", 1, -1)) %>%           # 1 = present, 2 = preterit        
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
  #           Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
 #  mod_ot1    9 240042 240120 -120012   240024                          
  # mod_ot2   14 239367 239489 -119670   239339 684.486  5  < 2.2e-16 ***
  # mod_ot3   20 239311 239484 -119635   239271  68.534  6  8.167e-13 ***
  # mod_ot4   21 239043 239225 -119501   239001 269.485  1  < 2.2e-16 ***
  
  
}

# -----------------------------------------------------------------------------





# Individual models -----------------------------------------------------------

#
# only mon
#

gca_mod_mon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
       REML = F,
       data = filter(stress_gc_subset, group == "mon")) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + condition_sum) # singular
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0, . ~ . + ot1:condition_sum) # singular
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1, . ~ . + ot2:condition_sum) # singular
gca_mod_mon_cond_3 <- update(gca_mod_mon_cond_2, . ~ . + ot3:condition_sum) # singular

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2, gca_mod_mon_cond_3) # 
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base   30 13994 14171 -6967.1    13934                         
# gca_mod_mon_cond_0 31 13996 14179 -6967.1    13934 0.0590      1     0.8081
# gca_mod_mon_cond_1 32 13996 14185 -6966.2    13932 1.8540      1     0.1733
# gca_mod_mon_cond_2 33 13998 14193 -6966.1    13932 0.1229      1     0.7259
# gca_mod_mon_cond_3 34 14000 14201 -6965.9    13932 0.3998      1     0.5272



 
# Error in anova.merMod(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,  : 
#                         models were not all fitted to the same size of dataset


# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_wm_0 <- update(gca_mod_mon_base,   . ~ . + WM_set) # singular
gca_mod_mon_wm_1 <- update(gca_mod_mon_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_mon_wm_2 <- update(gca_mod_mon_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_mon_wm_3 <- update(gca_mod_mon_wm_2, . ~ . + ot3:WM_set) # singular

mon_wm_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_wm_0, gca_mod_mon_wm_1,
        gca_mod_mon_wm_2, gca_mod_mon_wm_3)
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base 30 13994 14171 -6967.1    13934                           
# gca_mod_mon_wm_0 31 13996 14179 -6967.1    13934 0.0637      1    0.80067  
# gca_mod_mon_wm_1 32 13994 14183 -6965.1    13930 3.8943      1    0.04845 *
# gca_mod_mon_wm_2 33 13996 14191 -6965.1    13930 0.1571      1    0.69180  
# gca_mod_mon_wm_3 34 13998 14198 -6964.9    13930 0.3654      1    0.54551 


# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_spon_0 <- update(gca_mod_mon_wm_1,   . ~ . + spondee) # singular
gca_mod_mon_spon_1 <- update(gca_mod_mon_spon_0, . ~ . + ot1:spondee) # singular
gca_mod_mon_spon_2 <- update(gca_mod_mon_spon_1, . ~ . + ot2:spondee) # singular
gca_mod_mon_spon_3 <- update(gca_mod_mon_spon_2, . ~ . + ot3:spondee) # singular

mon_spon_anova <-
  anova(gca_mod_mon_wm_1, gca_mod_mon_spon_0, gca_mod_mon_spon_1,
        gca_mod_mon_spon_2, gca_mod_mon_spon_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_wm_1   32 13994 14183 -6965.1    13930                         
# gca_mod_mon_spon_0 33 13995 14190 -6964.6    13929 1.1880      1     0.2757
# gca_mod_mon_spon_1 34 13997 14198 -6964.6    13929 0.0010      1     0.9747
# gca_mod_mon_spon_2 35 13999 14205 -6964.3    13929 0.5434      1     0.4610
# gca_mod_mon_spon_3 36 14000 14213 -6964.0    13928 0.5565      1     0.4557

gca_mod_mon_stre_0 <- update(gca_mod_mon_wm_1,   . ~ . + stressed_spondee) # singular
gca_mod_mon_stre_1 <- update(gca_mod_mon_stre_0, . ~ . + ot1:stressed_spondee) # singular
gca_mod_mon_stre_2 <- update(gca_mod_mon_stre_1, . ~ . + ot2:stressed_spondee) # singular
gca_mod_mon_stre_3 <- update(gca_mod_mon_stre_2, . ~ . + ot3:stressed_spondee) # singular

mon_stre_anova <-
  anova(gca_mod_mon_wm_1, gca_mod_mon_stre_0, gca_mod_mon_stre_1,
        gca_mod_mon_stre_2, gca_mod_mon_stre_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_mon_wm_1   32 13994 14183 -6965.1    13930                           
# gca_mod_mon_stre_0 33 13993 14188 -6963.4    13927 3.4473      1    0.06335 .
# gca_mod_mon_stre_1 34 13995 14196 -6963.4    13927 0.0133      1    0.90829  
# gca_mod_mon_stre_2 35 13996 14203 -6963.0    13926 0.7476      1    0.38723  
# gca_mod_mon_stre_3 36 13997 14210 -6962.6    13925 0.8670      1    0.35179 

gca_mod_mon_troc_0 <- update(gca_mod_mon_wm_1,   . ~ . + trochee) # singular
gca_mod_mon_troc_1 <- update(gca_mod_mon_troc_0, . ~ . + ot1:trochee) # singular
gca_mod_mon_troc_2 <- update(gca_mod_mon_troc_1, . ~ . + ot2:trochee) # singular
gca_mod_mon_troc_3 <- update(gca_mod_mon_troc_2, . ~ . + ot3:trochee) # singular

mon_troc_anova <-
  anova(gca_mod_mon_wm_1, gca_mod_mon_troc_0, gca_mod_mon_troc_1,
        gca_mod_mon_troc_2, gca_mod_mon_troc_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_wm_1   32 13994 14183 -6965.1    13930                         
# gca_mod_mon_troc_0 33 13995 14190 -6964.7    13929 0.9383      1     0.3327
# gca_mod_mon_troc_1 34 13997 14198 -6964.4    13929 0.6009      1     0.4382
# gca_mod_mon_troc_2 35 13999 14205 -6964.4    13929 0.0019      1     0.9650
# gca_mod_mon_troc_3 36 14000 14213 -6964.0    13928 0.7239      1     0.3949

gca_mod_mon_unpr_0 <- update(gca_mod_mon_wm_1,   . ~ . + unpredictable) # singular
gca_mod_mon_unpr_1 <- update(gca_mod_mon_unpr_0, . ~ . + ot1:unpredictable) # singular
gca_mod_mon_unpr_2 <- update(gca_mod_mon_unpr_1, . ~ . + ot2:unpredictable) # singular
gca_mod_mon_unpr_3 <- update(gca_mod_mon_unpr_2, . ~ . + ot3:unpredictable) # singular

mon_unpr_anova <-
  anova(gca_mod_mon_wm_1, gca_mod_mon_unpr_0, gca_mod_mon_unpr_1,
        gca_mod_mon_unpr_2, gca_mod_mon_unpr_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_wm_1   32 13994 14183 -6965.1    13930                           
# gca_mod_mon_unpr_0 33 13990 14185 -6962.0    13924 6.2991      1    0.01208 *
# gca_mod_mon_unpr_1 34 13992 14193 -6961.9    13924 0.1002      1    0.75160  
# gca_mod_mon_unpr_2 35 13993 14199 -6961.3    13923 1.2427      1    0.26495  
# gca_mod_mon_unpr_3 36 13994 14207 -6961.1    13922 0.3848      1    0.53505  


# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cdo_0 <- update(gca_mod_mon_unpr_0,   . ~ . + DO_down) # singular
gca_mod_mon_cdo_1 <- update(gca_mod_mon_cdo_0, . ~ . + ot1:DO_down) # singular
gca_mod_mon_cdo_2 <- update(gca_mod_mon_cdo_1, . ~ . + ot2:DO_down) # singular
gca_mod_mon_cdo_3 <- update(gca_mod_mon_cdo_2, . ~ . + ot3:DO_down) # singular

mon_cdo_anova <-
  anova(gca_mod_mon_unpr_0, gca_mod_mon_cdo_0, gca_mod_mon_cdo_1,
        gca_mod_mon_cdo_2, gca_mod_mon_cdo_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_unpr_0 33 13990 14185 -6962.0    13924                           
# gca_mod_mon_cdo_0  34 13989 14190 -6960.5    13921 3.0915      1    0.07870 .
# gca_mod_mon_cdo_1  35 13990 14197 -6960.1    13920 0.7004      1    0.40265  
# gca_mod_mon_cdo_2  36 13988 14200 -6957.7    13916 4.7374      1    0.02951 *
# gca_mod_mon_cdo_3  37 13987 14205 -6956.4    13913 2.6246      1    0.10522  

gca_mod_mon_cup_0 <- update(gca_mod_mon_cdo_2,   . ~ . + DO_up) # singular
gca_mod_mon_cup_1 <- update(gca_mod_mon_cup_0, . ~ . + ot1:DO_up) # singular
gca_mod_mon_cup_2 <- update(gca_mod_mon_cup_1, . ~ . + ot2:DO_up) # singular
gca_mod_mon_cup_3 <- update(gca_mod_mon_cup_2, . ~ . + ot3:DO_up) # singular

mon_cup_anova <-
  anova(gca_mod_mon_cdo_2, gca_mod_mon_cup_0, gca_mod_mon_cup_1,
        gca_mod_mon_cup_2, gca_mod_mon_cup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_mod_mon_cdo_2 36 13988 14200 -6957.7    13916                              
# gca_mod_mon_cup_0 37 13974 14193 -6950.3    13900 14.9243      1  0.0001119 ***
# gca_mod_mon_cup_1 38 13976 14201 -6950.2    13900  0.1558      1  0.6930237    
# gca_mod_mon_cup_2 39 13977 14207 -6949.3    13899  1.7462      1  0.1863506    
# gca_mod_mon_cup_3 40 13979 14215 -6949.3    13899  0.0832      1  0.7729660    

gca_mod_mon_edo_0 <- update(gca_mod_mon_cup_0,   . ~ . + MI_down) # singular
gca_mod_mon_edo_1 <- update(gca_mod_mon_edo_0, . ~ . + ot1:MI_down) # singular
gca_mod_mon_edo_2 <- update(gca_mod_mon_edo_1, . ~ . + ot2:MI_down) # singular
gca_mod_mon_edo_3 <- update(gca_mod_mon_edo_2, . ~ . + ot3:MI_down) # singular

mon_edo_anova <-
  anova(gca_mod_mon_cup_0, gca_mod_mon_edo_0, gca_mod_mon_edo_1,
        gca_mod_mon_edo_2, gca_mod_mon_edo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_mon_cup_0 37 13974 14193 -6950.3    13900                           
# gca_mod_mon_edo_0 38 13981 14205 -6952.5    13905 0.0000      1    1.00000  
# gca_mod_mon_edo_1 39 13978 14209 -6950.3    13900 4.4985      1    0.03393 *
# gca_mod_mon_edo_2 40 13980 14216 -6950.1    13900 0.2770      1    0.59870  
# gca_mod_mon_edo_3 41 13982 14224 -6949.8    13900 0.6804      1    0.40944  

gca_mod_mon_eup_0 <- update(gca_mod_mon_edo_1,   . ~ . + MI_up) # singular
gca_mod_mon_eup_1 <- update(gca_mod_mon_eup_0, . ~ . + ot1:MI_up) # singular
gca_mod_mon_eup_2 <- update(gca_mod_mon_eup_1, . ~ . + ot2:MI_up) # singular
gca_mod_mon_eup_3 <- update(gca_mod_mon_eup_2, . ~ . + ot3:MI_up) # singular

mon_eup_anova <-
  anova(gca_mod_mon_edo_1, gca_mod_mon_eup_0, gca_mod_mon_eup_1,
        gca_mod_mon_eup_2, gca_mod_mon_eup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_mon_edo_1 39 13978 14209 -6950.3    13900                         
# gca_mod_mon_eup_0 40 13979 14215 -6949.4    13899 1.6329      1     0.2013
# gca_mod_mon_eup_1 41 13980 14222 -6949.0    13898 0.9198      1     0.3375
# gca_mod_mon_eup_2 42 13980 14228 -6948.0    13896 2.0311      1     0.1541
# gca_mod_mon_eup_3 43 13981 14235 -6947.7    13895 0.6164      1     0.4324


gca_mod_mon_gdo_0 <- update(gca_mod_mon_edo_1,   . ~ . + SOL_down) # singular
gca_mod_mon_gdo_1 <- update(gca_mod_mon_gdo_0, . ~ . + ot1:SOL_down) # singular
gca_mod_mon_gdo_2 <- update(gca_mod_mon_gdo_1, . ~ . + ot2:SOL_down) # singular
gca_mod_mon_gdo_3 <- update(gca_mod_mon_gdo_2, . ~ . + ot3:SOL_down) # singular

mon_gdo_anova <-
  anova(gca_mod_mon_edo_1, gca_mod_mon_gdo_0, gca_mod_mon_gdo_1,
        gca_mod_mon_gdo_2, gca_mod_mon_gdo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_mon_edo_1 39 13978 14209 -6950.3    13900                         
# gca_mod_mon_gdo_0 40 13981 14217 -6950.6    13901 0.0000      1     1.0000
# gca_mod_mon_gdo_1 41 13981 14223 -6949.7    13899 1.8184      1     0.1775
# gca_mod_mon_gdo_2 42 13983 14231 -6949.7    13899 0.0008      1     0.9779
# gca_mod_mon_gdo_3 43 13983 14237 -6948.6    13897 2.0075      1     0.1565


gca_mod_mon_gup_0 <- update(gca_mod_mon_edo_1,   . ~ . + SOL_up) # singular
gca_mod_mon_gup_1 <- update(gca_mod_mon_gup_0, . ~ . + ot1:SOL_up) # singular
gca_mod_mon_gup_2 <- update(gca_mod_mon_gup_1, . ~ . + ot2:SOL_up) # singular
gca_mod_mon_gup_3 <- update(gca_mod_mon_gup_2, . ~ . + ot3:SOL_up) # singular

mon_gup_anova <-
  anova(gca_mod_mon_edo_1, gca_mod_mon_gup_0, gca_mod_mon_gup_1,
        gca_mod_mon_gup_2, gca_mod_mon_gup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_mon_edo_1 39 13978 14209 -6950.3    13900                           
# gca_mod_mon_gup_0 40 13977 14213 -6948.5    13897 3.5623      1    0.05911 .
# gca_mod_mon_gup_1 41 13979 14221 -6948.6    13897 0.0000      1    1.00000  
# gca_mod_mon_gup_2 42 13981 14229 -6948.3    13897 0.6067      1    0.43602  
# gca_mod_mon_gup_3 43 13981 14235 -6947.4    13895 1.8731      1    0.17112  




#######

# only aes

#######

gca_mod_aes_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 3e5)), 
       REML = F,
       data = filter(stress_gc_subset, group == "aes")) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_cond_0 <- update(gca_mod_aes_base,   . ~ . + condition_sum) # singular
gca_mod_aes_cond_1 <- update(gca_mod_aes_cond_0, . ~ . + ot1:condition_sum) # singular
gca_mod_aes_cond_2 <- update(gca_mod_aes_cond_1, . ~ . + ot2:condition_sum) # singular
gca_mod_aes_cond_3 <- update(gca_mod_aes_cond_2, . ~ . + ot3:condition_sum) # singular

aes_cond_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_cond_0, gca_mod_aes_cond_1,
        gca_mod_aes_cond_2, gca_mod_aes_cond_3) # 
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_aes_base   30 22417 22608 -11178    22357                         
# gca_mod_aes_cond_0 31 22418 22616 -11178    22356 0.1706      1     0.6796
# gca_mod_aes_cond_1 32 22420 22624 -11178    22356 0.5106      1     0.4749
# gca_mod_aes_cond_2 33 22422 22632 -11178    22356 0.1228      1     0.7261
# gca_mod_aes_cond_3 34 22424 22640 -11178    22356 0.0157      1     0.9001


# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_wm_0 <- update(gca_mod_aes_base,   . ~ . + WM_set) # singular
gca_mod_aes_wm_1 <- update(gca_mod_aes_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_aes_wm_2 <- update(gca_mod_aes_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_aes_wm_3 <- update(gca_mod_aes_wm_2, . ~ . + ot3:WM_set) # singular

aes_wm_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_wm_0, gca_mod_aes_wm_1,
        gca_mod_aes_wm_2, gca_mod_aes_wm_3)
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_aes_base 30 22417 22608 -11178    22357                         
# gca_mod_aes_wm_0 31 22419 22616 -11178    22357 0.0031      1     0.9555
# gca_mod_aes_wm_1 32 22420 22624 -11178    22356 0.0408      1     0.8399
# gca_mod_aes_wm_2 33 22422 22633 -11178    22356 0.2917      1     0.5892
# gca_mod_aes_wm_3 34 22422 22639 -11177    22354 1.9778      1     0.1596

# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_spon_0 <- update(gca_mod_aes_base,   . ~ . + spondee) # singular
gca_mod_aes_spon_1 <- update(gca_mod_aes_spon_0, . ~ . + ot1:spondee) # singular
gca_mod_aes_spon_2 <- update(gca_mod_aes_spon_1, . ~ . + ot2:spondee) # singular
gca_mod_aes_spon_3 <- update(gca_mod_aes_spon_2, . ~ . + ot3:spondee) # singular

aes_spon_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_spon_0, gca_mod_aes_spon_1,
        gca_mod_aes_spon_2, gca_mod_aes_spon_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_aes_base   30 22417 22608 -11178    22357                         
# gca_mod_aes_spon_0 31 22417 22615 -11178    22355 1.6053      1     0.2051
# gca_mod_aes_spon_1 32 22417 22621 -11177    22353 1.5435      1     0.2141
# gca_mod_aes_spon_2 33 22419 22630 -11177    22353 0.0223      1     0.8812
# gca_mod_aes_spon_3 34 22420 22637 -11176    22352 1.1325      1     0.2873

gca_mod_aes_stre_0 <- update(gca_mod_aes_base,   . ~ . + stressed_spondee) # singular
gca_mod_aes_stre_1 <- update(gca_mod_aes_stre_0, . ~ . + ot1:stressed_spondee) # singular
gca_mod_aes_stre_2 <- update(gca_mod_aes_stre_1, . ~ . + ot2:stressed_spondee) # singular
gca_mod_aes_stre_3 <- update(gca_mod_aes_stre_2, . ~ . + ot3:stressed_spondee) # singular

aes_stre_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_stre_0, gca_mod_aes_stre_1,
        gca_mod_aes_stre_2, gca_mod_aes_stre_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_aes_base   30 22417 22608 -11178    22357                           
# gca_mod_aes_stre_0 31 22415 22613 -11177    22353 3.3044      1    0.06909 .
# gca_mod_aes_stre_1 32 22416 22620 -11176    22352 1.3310      1    0.24862  
# gca_mod_aes_stre_2 33 22418 22628 -11176    22352 0.0640      1    0.80021  
# gca_mod_aes_stre_3 34 22418 22635 -11175    22350 1.3467      1    0.24585

gca_mod_aes_troc_0 <- update(gca_mod_aes_base,   . ~ . + trochee) # singular
gca_mod_aes_troc_1 <- update(gca_mod_aes_troc_0, . ~ . + ot1:trochee) # singular
gca_mod_aes_troc_2 <- update(gca_mod_aes_troc_1, . ~ . + ot2:trochee) # singular
gca_mod_aes_troc_3 <- update(gca_mod_aes_troc_2, . ~ . + ot3:trochee) # singular

aes_troc_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_troc_0, gca_mod_aes_troc_1,
        gca_mod_aes_troc_2, gca_mod_aes_troc_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_aes_base   30 22417 22608 -11178    22357                         
# gca_mod_aes_troc_0 31 22417 22615 -11178    22355 1.4791      1     0.2239
# gca_mod_aes_troc_1 32 22419 22623 -11178    22355 0.1348      1     0.7135
# gca_mod_aes_troc_2 33 22421 22631 -11178    22355 0.0572      1     0.8110
# gca_mod_aes_troc_3 34 22421 22638 -11176    22353 1.9073      1     0.1673


gca_mod_aes_unpr_0 <- update(gca_mod_aes_base,   . ~ . + unpredictable) # singular
gca_mod_aes_unpr_1 <- update(gca_mod_aes_unpr_0, . ~ . + ot1:unpredictable) # singular
gca_mod_aes_unpr_2 <- update(gca_mod_aes_unpr_1, . ~ . + ot2:unpredictable) # singular
gca_mod_aes_unpr_3 <- update(gca_mod_aes_unpr_2, . ~ . + ot3:unpredictable) # singular

aes_unpr_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_unpr_0, gca_mod_aes_unpr_1,
        gca_mod_aes_unpr_2, gca_mod_aes_unpr_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_aes_base   30 22417 22608 -11178    22357                         
# gca_mod_aes_unpr_0 31 22418 22616 -11178    22356 0.5510      1     0.4579
# gca_mod_aes_unpr_1 32 22420 22624 -11178    22356 0.1015      1     0.7500
# gca_mod_aes_unpr_2 33 22422 22632 -11178    22356 0.2339      1     0.6287
# gca_mod_aes_unpr_3 34 22423 22640 -11178    22355 0.6524      1     0.4192


# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_aes_cdo_0 <- update(gca_mod_aes_base,   . ~ . + DO_down) # singular
gca_mod_aes_cdo_1 <- update(gca_mod_aes_cdo_0, . ~ . + ot1:DO_down) # singular
gca_mod_aes_cdo_2 <- update(gca_mod_aes_cdo_1, . ~ . + ot2:DO_down) # singular
gca_mod_aes_cdo_3 <- update(gca_mod_aes_cdo_2, . ~ . + ot3:DO_down) # singular

aes_cdo_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_cdo_0, gca_mod_aes_cdo_1,
        gca_mod_aes_cdo_2, gca_mod_aes_cdo_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_aes_base  30 22417 22608 -11178    22357                         
# gca_mod_aes_cdo_0 31 22419 22616 -11178    22357 0.0064      1     0.9364
# gca_mod_aes_cdo_1 32 22420 22624 -11178    22356 0.4755      1     0.4905
# gca_mod_aes_cdo_2 33 22421 22632 -11178    22355 0.7104      1     0.3993
# gca_mod_aes_cdo_3 34 22422 22639 -11177    22354 0.8728      1     0.3502

gca_mod_aes_cup_0 <- update(gca_mod_aes_base,   . ~ . + DO_up) # singular
gca_mod_aes_cup_1 <- update(gca_mod_aes_cup_0, . ~ . + ot1:DO_up) # singular
gca_mod_aes_cup_2 <- update(gca_mod_aes_cup_1, . ~ . + ot2:DO_up) # singular
gca_mod_aes_cup_3 <- update(gca_mod_aes_cup_2, . ~ . + ot3:DO_up) # singular

aes_cup_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_cup_0, gca_mod_aes_cup_1,
        gca_mod_aes_cup_2, gca_mod_aes_cup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_mod_aes_base  30 22417 22608 -11178    22357                         
# gca_mod_aes_cup_0 31 22418 22616 -11178    22356 0.0479      1     0.8267
# gca_mod_aes_cup_1 32 22420 22624 -11178    22356 0.6043      1     0.4369
# gca_mod_aes_cup_2 33 22421 22632 -11178    22355 0.8081      1     0.3687
# gca_mod_aes_cup_3 34 22422 22639 -11177    22354 0.6883      1     0.4068

gca_mod_aes_edo_0 <- update(gca_mod_aes_base,   . ~ . + MI_down) # singular
gca_mod_aes_edo_1 <- update(gca_mod_aes_edo_0, . ~ . + ot1:MI_down) # singular
gca_mod_aes_edo_2 <- update(gca_mod_aes_edo_1, . ~ . + ot2:MI_down) # singular
gca_mod_aes_edo_3 <- update(gca_mod_aes_edo_2, . ~ . + ot3:MI_down) # singular

aes_edo_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_edo_0, gca_mod_aes_edo_1,
        gca_mod_aes_edo_2, gca_mod_aes_edo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_aes_base  30 22417 22608 -11178    22357                         
# gca_mod_aes_edo_0 31 22418 22616 -11178    22356 0.1882      1     0.6645
# gca_mod_aes_edo_1 32 22420 22624 -11178    22356 0.4864      1     0.4855
# gca_mod_aes_edo_2 33 22421 22632 -11178    22355 0.6711      1     0.4127
# gca_mod_aes_edo_3 34 22423 22640 -11178    22355 0.0000      1     0.9989

gca_mod_aes_eup_0 <- update(gca_mod_aes_base,   . ~ . + MI_up) # singular
gca_mod_aes_eup_1 <- update(gca_mod_aes_eup_0, . ~ . + ot1:MI_up) # singular
gca_mod_aes_eup_2 <- update(gca_mod_aes_eup_1, . ~ . + ot2:MI_up) # singular
gca_mod_aes_eup_3 <- update(gca_mod_aes_eup_2, . ~ . + ot3:MI_up) # singular

aes_eup_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_eup_0, gca_mod_aes_eup_1,
        gca_mod_aes_eup_2, gca_mod_aes_eup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_aes_base  30 22417 22608 -11178    22357                         
# gca_mod_aes_eup_0 31 22418 22616 -11178    22356 0.2372      1     0.6262
# gca_mod_aes_eup_1 32 22420 22624 -11178    22356 0.0997      1     0.7522
# gca_mod_aes_eup_2 33 22422 22632 -11178    22356 0.3705      1     0.5427
# gca_mod_aes_eup_3 34 22424 22640 -11178    22356 0.1099      1     0.7402

gca_mod_aes_gdo_0 <- update(gca_mod_aes_base,   . ~ . + SOL_down) # singular
gca_mod_aes_gdo_1 <- update(gca_mod_aes_gdo_0, . ~ . + ot1:SOL_down) # singular
gca_mod_aes_gdo_2 <- update(gca_mod_aes_gdo_1, . ~ . + ot2:SOL_down) # singular
gca_mod_aes_gdo_3 <- update(gca_mod_aes_gdo_2, . ~ . + ot3:SOL_down) # singular

aes_gdo_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_gdo_0, gca_mod_aes_gdo_1,
        gca_mod_aes_gdo_2, gca_mod_aes_gdo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_aes_base  30 22417 22608 -11178    22357                         
# gca_mod_aes_gdo_0 31 22419 22616 -11178    22357 0.0150      1     0.9025
# gca_mod_aes_gdo_1 32 22420 22624 -11178    22356 0.4285      1     0.5127
# gca_mod_aes_gdo_2 33 22422 22632 -11178    22356 0.6884      1     0.4067
# gca_mod_aes_gdo_3 34 22423 22640 -11178    22355 0.3218      1     0.5705

gca_mod_aes_gup_0 <- update(gca_mod_aes_base,   . ~ . + SOL_up) # singular
gca_mod_aes_gup_1 <- update(gca_mod_aes_gup_0, . ~ . + ot1:SOL_up) # singular
gca_mod_aes_gup_2 <- update(gca_mod_aes_gup_1, . ~ . + ot2:SOL_up) # singular
gca_mod_aes_gup_3 <- update(gca_mod_aes_gup_2, . ~ . + ot3:SOL_up) # singular

aes_gup_anova <-
  anova(gca_mod_aes_base, gca_mod_aes_gup_0, gca_mod_aes_gup_1,
        gca_mod_aes_gup_2, gca_mod_aes_gup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_aes_base  30 22417 22608 -11178    22357                         
# gca_mod_aes_gup_0 31 22418 22616 -11178    22356 0.1765      1     0.6744
# gca_mod_aes_gup_1 32 22420 22624 -11178    22356 0.0049      1     0.9443
# gca_mod_aes_gup_2 33 22422 22633 -11178    22356 0.0585      1     0.8089
# gca_mod_aes_gup_3 34 22424 22641 -11178    22356 0.0253      1     0.8735




#
# only ies
#

gca_mod_ies_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ies")) # singular

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_cond_0 <- update(gca_mod_ies_base,   . ~ . + condition_sum) # singular
gca_mod_ies_cond_1 <- update(gca_mod_ies_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ies_cond_2 <- update(gca_mod_ies_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ies_cond_3 <- update(gca_mod_ies_cond_2,   . ~ . + ot3:condition_sum) # singular

ies_cond_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_cond_0, gca_mod_ies_cond_1,    
        gca_mod_ies_cond_2, gca_mod_ies_cond_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_base   30 23359 23550 -11650    23299                         
# gca_mod_ies_cond_0 31 23361 23558 -11650    23299 0.0365      1     0.8485
# gca_mod_ies_cond_1 32 23362 23566 -11649    23298 1.3143      1     0.2516
# gca_mod_ies_cond_2 33 23364 23574 -11649    23298 0.0138      1     0.9065
# gca_mod_ies_cond_3 34 23364 23581 -11648    23296 1.4464      1     0.2291


# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_wm_0 <- update(gca_mod_ies_base,   . ~ . + WM_set) # singular
gca_mod_ies_wm_1 <- update(gca_mod_ies_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_ies_wm_2 <- update(gca_mod_ies_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_ies_wm_3 <- update(gca_mod_ies_wm_2, . ~ . + ot3:WM_set) # singular

ies_wm_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_wm_0, gca_mod_ies_wm_1,
        gca_mod_ies_wm_2, gca_mod_ies_wm_3)
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_base 30 23359 23550 -11650    23299                         
# gca_mod_ies_wm_0 31 23359 23557 -11649    23297 1.6632      1     0.1972
# gca_mod_ies_wm_1 32 23361 23565 -11649    23297 0.3165      1     0.5737
# gca_mod_ies_wm_2 33 23360 23571 -11647    23294 2.5894      1     0.1076
# gca_mod_ies_wm_3 34 23362 23578 -11647    23294 0.4919      1     0.4831

# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_spon_0 <- update(gca_mod_ies_base,   . ~ . + spondee) # singular
gca_mod_ies_spon_1 <- update(gca_mod_ies_spon_0, . ~ . + ot1:spondee) # singular
gca_mod_ies_spon_2 <- update(gca_mod_ies_spon_1, . ~ . + ot2:spondee) # singular
gca_mod_ies_spon_3 <- update(gca_mod_ies_spon_2, . ~ . + ot3:spondee) # singular

ies_spon_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_spon_0, gca_mod_ies_spon_1,
        gca_mod_ies_spon_2, gca_mod_ies_spon_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_base   30 23359 23550 -11650    23299                         
# gca_mod_ies_spon_0 31 23360 23557 -11649    23298 1.1810      1     0.2772
# gca_mod_ies_spon_1 32 23362 23566 -11649    23298 0.0112      1     0.9156
# gca_mod_ies_spon_2 33 23364 23574 -11649    23298 0.1178      1     0.7314
# gca_mod_ies_spon_3 34 23365 23581 -11648    23297 0.9223      1     0.3369

gca_mod_ies_stre_0 <- update(gca_mod_ies_base,   . ~ . + stressed_spondee) # singular
gca_mod_ies_stre_1 <- update(gca_mod_ies_stre_0, . ~ . + ot1:stressed_spondee) # singular
gca_mod_ies_stre_2 <- update(gca_mod_ies_stre_1, . ~ . + ot2:stressed_spondee) # singular
gca_mod_ies_stre_3 <- update(gca_mod_ies_stre_2, . ~ . + ot3:stressed_spondee) # singular

ies_stre_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_stre_0, gca_mod_ies_stre_1,
        gca_mod_ies_stre_2, gca_mod_ies_stre_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ies_base   30 23359 23550 -11650    23299                         
# gca_mod_ies_stre_0 31 23359 23556 -11648    23297 2.0416      1     0.1530
# gca_mod_ies_stre_1 32 23361 23565 -11648    23297 0.1615      1     0.6878
# gca_mod_ies_stre_2 33 23363 23573 -11648    23297 0.0654      1     0.7982
# gca_mod_ies_stre_3 34 23364 23580 -11648    23296 1.1963      1     0.2741

gca_mod_ies_troc_0 <- update(gca_mod_ies_base,   . ~ . + trochee) # singular
gca_mod_ies_troc_1 <- update(gca_mod_ies_troc_0, . ~ . + ot1:trochee) # singular
gca_mod_ies_troc_2 <- update(gca_mod_ies_troc_1, . ~ . + ot2:trochee) # singular
gca_mod_ies_troc_3 <- update(gca_mod_ies_troc_2, . ~ . + ot3:trochee) # singular

ies_troc_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_troc_0, gca_mod_ies_troc_1,
        gca_mod_ies_troc_2, gca_mod_ies_troc_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_base   30 23359 23550 -11650    23299                         
# gca_mod_ies_troc_0 31 23361 23558 -11649    23299 0.4856      1     0.4859
# gca_mod_ies_troc_1 32 23363 23566 -11649    23299 0.0029      1     0.9570
# gca_mod_ies_troc_2 33 23364 23574 -11649    23298 0.2192      1     0.6396
# gca_mod_ies_troc_3 34 23365 23582 -11648    23297 1.3470      1     0.2458

gca_mod_ies_unpr_0 <- update(gca_mod_ies_base,   . ~ . + unpredictable) # singular
gca_mod_ies_unpr_1 <- update(gca_mod_ies_unpr_0, . ~ . + ot1:unpredictable) # singular
gca_mod_ies_unpr_2 <- update(gca_mod_ies_unpr_1, . ~ . + ot2:unpredictable) # singular
gca_mod_ies_unpr_3 <- update(gca_mod_ies_unpr_2, . ~ . + ot3:unpredictable) # singular

ies_unpr_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_unpr_0, gca_mod_ies_unpr_1,
        gca_mod_ies_unpr_2, gca_mod_ies_unpr_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_base   30 23359 23550 -11650    23299                           
# gca_mod_ies_unpr_0 31 23361 23558 -11650    23299 0.0528      1    0.81824  
# gca_mod_ies_unpr_1 32 23363 23566 -11649    23299 0.2888      1    0.59098  
# gca_mod_ies_unpr_2 33 23364 23574 -11649    23298 0.3056      1    0.58037  
# gca_mod_ies_unpr_3 34 23363 23580 -11648    23295 3.2673      1    0.07067 .

# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ies_cdo_0 <- update(gca_mod_ies_base,   . ~ . + DO_down) # singular
gca_mod_ies_cdo_1 <- update(gca_mod_ies_cdo_0, . ~ . + ot1:DO_down) # singular
gca_mod_ies_cdo_2 <- update(gca_mod_ies_cdo_1, . ~ . + ot2:DO_down) # singular
gca_mod_ies_cdo_3 <- update(gca_mod_ies_cdo_2, . ~ . + ot3:DO_down) # singular

ies_cdo_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_cdo_0, gca_mod_ies_cdo_1,
        gca_mod_ies_cdo_2, gca_mod_ies_cdo_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ies_base  30 23359 23550 -11650    23299                         
# gca_mod_ies_cdo_0 31 23359 23556 -11648    23297 2.5312      1     0.1116
# gca_mod_ies_cdo_1 32 23360 23564 -11648    23296 0.0174      1     0.8950
# gca_mod_ies_cdo_2 33 23362 23573 -11648    23296 0.0050      1     0.9439
# gca_mod_ies_cdo_3 34 23364 23581 -11648    23296 0.4044      1     0.5248

gca_mod_ies_cup_0 <- update(gca_mod_ies_base,   . ~ . + DO_up) # singular
gca_mod_ies_cup_1 <- update(gca_mod_ies_cup_0, . ~ . + ot1:DO_up) # singular
gca_mod_ies_cup_2 <- update(gca_mod_ies_cup_1, . ~ . + ot2:DO_up) # singular
gca_mod_ies_cup_3 <- update(gca_mod_ies_cup_2, . ~ . + ot3:DO_up) # singular

ies_cup_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_cup_0, gca_mod_ies_cup_1,
        gca_mod_ies_cup_2, gca_mod_ies_cup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_mod_ies_base  30 23359 23550 -11650    23299                         
# gca_mod_ies_cup_0 31 23361 23558 -11650    23299 0.1568      1     0.6921
# gca_mod_ies_cup_1 32 23362 23566 -11649    23298 0.8798      1     0.3482
# gca_mod_ies_cup_2 33 23362 23573 -11648    23296 1.5485      1     0.2134
# gca_mod_ies_cup_3 34 23364 23581 -11648    23296 0.0271      1     0.8693

gca_mod_ies_edo_0 <- update(gca_mod_ies_base,   . ~ . + MI_down) # singular
gca_mod_ies_edo_1 <- update(gca_mod_ies_edo_0, . ~ . + ot1:MI_down) # singular
gca_mod_ies_edo_2 <- update(gca_mod_ies_edo_1, . ~ . + ot2:MI_down) # singular
gca_mod_ies_edo_3 <- update(gca_mod_ies_edo_2, . ~ . + ot3:MI_down) # singular

ies_edo_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_edo_0, gca_mod_ies_edo_1,
        gca_mod_ies_edo_2, gca_mod_ies_edo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ies_base  30 23359 23550 -11650    23299                         
# gca_mod_ies_edo_0 31 23361 23558 -11650    23299 0.1453      1     0.7031
# gca_mod_ies_edo_1 32 23362 23566 -11649    23298 1.0287      1     0.3105
# gca_mod_ies_edo_2 33 23363 23573 -11649    23297 0.7971      1     0.3720
# gca_mod_ies_edo_3 34 23365 23581 -11648    23297 0.1767      1     0.6743

gca_mod_ies_eup_0 <- update(gca_mod_ies_base,   . ~ . + MI_up) # singular
gca_mod_ies_eup_1 <- update(gca_mod_ies_eup_0, . ~ . + ot1:MI_up) # singular
gca_mod_ies_eup_2 <- update(gca_mod_ies_eup_1, . ~ . + ot2:MI_up) # singular
gca_mod_ies_eup_3 <- update(gca_mod_ies_eup_2, . ~ . + ot3:MI_up) # singular

ies_eup_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_eup_0, gca_mod_ies_eup_1,
        gca_mod_ies_eup_2, gca_mod_ies_eup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ies_base  30 23359 23550 -11650    23299                         
# gca_mod_ies_eup_0 31 23361 23558 -11650    23299 0.0129      1     0.9096
# gca_mod_ies_eup_1 32 23361 23564 -11648    23297 2.4404      1     0.1182
# gca_mod_ies_eup_2 33 23362 23573 -11648    23296 0.1226      1     0.7263
# gca_mod_ies_eup_3 34 23364 23581 -11648    23296 0.0809      1     0.7761

gca_mod_ies_gdo_0 <- update(gca_mod_ies_base,   . ~ . + SOL_down) # singular
gca_mod_ies_gdo_1 <- update(gca_mod_ies_gdo_0, . ~ . + ot1:SOL_down) # singular
gca_mod_ies_gdo_2 <- update(gca_mod_ies_gdo_1, . ~ . + ot2:SOL_down) # singular
gca_mod_ies_gdo_3 <- update(gca_mod_ies_gdo_2, . ~ . + ot3:SOL_down) # singular

ies_gdo_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_gdo_0, gca_mod_ies_gdo_1,
        gca_mod_ies_gdo_2, gca_mod_ies_gdo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ies_base  30 23359 23550 -11650    23299                         
# gca_mod_ies_gdo_0 31 23361 23558 -11650    23299 0.0636      1     0.8009
# gca_mod_ies_gdo_1 32 23363 23566 -11649    23299 0.3692      1     0.5434
# gca_mod_ies_gdo_2 33 23363 23573 -11648    23297 1.7727      1     0.1831
# gca_mod_ies_gdo_3 34 23365 23581 -11648    23297 0.0000      1     0.9946

gca_mod_ies_gup_0 <- update(gca_mod_ies_base,   . ~ . + SOL_up) # singular
gca_mod_ies_gup_1 <- update(gca_mod_ies_gup_0, . ~ . + ot1:SOL_up) # singular
gca_mod_ies_gup_2 <- update(gca_mod_ies_gup_1, . ~ . + ot2:SOL_up) # singular
gca_mod_ies_gup_3 <- update(gca_mod_ies_gup_2, . ~ . + ot3:SOL_up) # singular

ies_gup_anova <-
  anova(gca_mod_ies_base, gca_mod_ies_gup_0, gca_mod_ies_gup_1,
        gca_mod_ies_gup_2, gca_mod_ies_gup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ies_base  30 23359 23550 -11650    23299                         
# gca_mod_ies_gup_0 31 23360 23557 -11649    23298 1.2053      1     0.2723
# gca_mod_ies_gup_1 32 23362 23566 -11649    23298 0.1373      1     0.7110
# gca_mod_ies_gup_2 33 23364 23574 -11649    23298 0.1552      1     0.6936
# gca_mod_ies_gup_3 34 23364 23581 -11648    23296 1.2303      1     0.2674




#
# only ams
#

gca_mod_ams_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ams")) # singular

# add cond effect to amsercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_cond_0 <- update(gca_mod_ams_base,   . ~ . + condition_sum) # singular
gca_mod_ams_cond_1 <- update(gca_mod_ams_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ams_cond_2 <- update(gca_mod_ams_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ams_cond_3 <- update(gca_mod_ams_cond_2,   . ~ . + ot3:condition_sum) # singular

ams_cond_anova <-
  anova(gca_mod_ams_base, gca_mod_ams_cond_0, gca_mod_ams_cond_1,
        gca_mod_ams_cond_2, gca_mod_ams_cond_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_base   30 26129 26324 -13035    26069                           
# gca_mod_ams_cond_0 31 26130 26331 -13034    26068 1.5074      1    0.21954  
# gca_mod_ams_cond_1 32 26129 26336 -13032    26065 3.1551      1    0.07569 .
# gca_mod_ams_cond_2 33 26130 26344 -13032    26064 0.3882      1    0.53322  
# gca_mod_ams_cond_3 34 26129 26349 -13030    26061 3.4212      1    0.06437 .


# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_wm_0 <- update(gca_mod_ams_base,   . ~ . + WM_set) # singular
gca_mod_ams_wm_1 <- update(gca_mod_ams_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_ams_wm_2 <- update(gca_mod_ams_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_ams_wm_3 <- update(gca_mod_ams_wm_2, . ~ . + ot3:WM_set) # singular

ams_wm_anova <-
  anova(gca_mod_ams_base, gca_mod_ams_wm_0, gca_mod_ams_wm_1,
        gca_mod_ams_wm_2, gca_mod_ams_wm_3)
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_base 30 26129 26324 -13035    26069                           
# gca_mod_ams_wm_0 31 26131 26332 -13035    26069 0.0476      1    0.82730  
# gca_mod_ams_wm_1 32 26129 26336 -13032    26065 4.2586      1    0.03905 *
# gca_mod_ams_wm_2 33 26130 26344 -13032    26064 0.9086      1    0.34048  
# gca_mod_ams_wm_3 34 26131 26352 -13032    26063 0.5803      1    0.44620  

# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_spon_0 <- update(gca_mod_ams_wm_1,   . ~ . + spondee) # singular
gca_mod_ams_spon_1 <- update(gca_mod_ams_spon_0, . ~ . + ot1:spondee) # singular
gca_mod_ams_spon_2 <- update(gca_mod_ams_spon_1, . ~ . + ot2:spondee) # singular
gca_mod_ams_spon_3 <- update(gca_mod_ams_spon_2, . ~ . + ot3:spondee) # singular

ams_spon_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_spon_0, gca_mod_ams_spon_1,
        gca_mod_ams_spon_2, gca_mod_ams_spon_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_wm_1   32 26129 26336 -13032    26065                         
# gca_mod_ams_spon_0 33 26131 26345 -13032    26065 0.3026      1     0.5822
# gca_mod_ams_spon_1 34 26131 26352 -13032    26063 1.2445      1     0.2646
# gca_mod_ams_spon_2 35 26133 26360 -13031    26063 0.7851      1     0.3756
# gca_mod_ams_spon_3 36 26134 26368 -13031    26062 0.0513      1     0.8208

gca_mod_ams_stre_0 <- update(gca_mod_ams_wm_1,   . ~ . + stressed_spondee) # singular
gca_mod_ams_stre_1 <- update(gca_mod_ams_stre_0, . ~ . + ot1:stressed_spondee) # singular
gca_mod_ams_stre_2 <- update(gca_mod_ams_stre_1, . ~ . + ot2:stressed_spondee) # singular
gca_mod_ams_stre_3 <- update(gca_mod_ams_stre_2, . ~ . + ot3:stressed_spondee) # singular

ams_stre_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_stre_0, gca_mod_ams_stre_1,
        gca_mod_ams_stre_2, gca_mod_ams_stre_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ams_wm_1   32 26129 26336 -13032    26065                         
# gca_mod_ams_stre_0 33 26130 26344 -13032    26064 0.3830      1     0.5360
# gca_mod_ams_stre_1 34 26131 26352 -13032    26063 1.4756      1     0.2245
# gca_mod_ams_stre_2 35 26133 26360 -13032    26063 0.0871      1     0.7679
# gca_mod_ams_stre_3 36 26135 26368 -13032    26063 0.0122      1     0.9119

gca_mod_ams_troc_0 <- update(gca_mod_ams_wm_1,   . ~ . + trochee) # singular
gca_mod_ams_troc_1 <- update(gca_mod_ams_troc_0, . ~ . + ot1:trochee) # singular
gca_mod_ams_troc_2 <- update(gca_mod_ams_troc_1, . ~ . + ot2:trochee) # singular
gca_mod_ams_troc_3 <- update(gca_mod_ams_troc_2, . ~ . + ot3:trochee) # singular

ams_troc_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_troc_0, gca_mod_ams_troc_1,
        gca_mod_ams_troc_2, gca_mod_ams_troc_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_wm_1   32 26129 26336 -13032    26065                         
# gca_mod_ams_troc_0 33 26130 26344 -13032    26064 0.8886      1     0.3459
# gca_mod_ams_troc_1 34 26132 26352 -13032    26064 0.4995      1     0.4797
# gca_mod_ams_troc_2 35 26133 26360 -13032    26063 0.1708      1     0.6794
# gca_mod_ams_troc_3 36 26134 26368 -13031    26062 1.0715      1     0.3006

gca_mod_ams_unpr_0 <- update(gca_mod_ams_wm_1,   . ~ . + unpredictable) # singular
gca_mod_ams_unpr_1 <- update(gca_mod_ams_unpr_0, . ~ . + ot1:unpredictable) # singular
gca_mod_ams_unpr_2 <- update(gca_mod_ams_unpr_1, . ~ . + ot2:unpredictable) # singular
gca_mod_ams_unpr_3 <- update(gca_mod_ams_unpr_2, . ~ . + ot3:unpredictable) # singular

ams_unpr_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_unpr_0, gca_mod_ams_unpr_1,
        gca_mod_ams_unpr_2, gca_mod_ams_unpr_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_wm_1   32 26129 26336 -13032    26065                         
# gca_mod_ams_unpr_0 33 26129 26343 -13032    26063 1.6712      1     0.1961
# gca_mod_ams_unpr_1 34 26130 26350 -13031    26062 1.3694      1     0.2419
# gca_mod_ams_unpr_2 35 26131 26358 -13031    26061 0.7374      1     0.3905
# gca_mod_ams_unpr_3 36 26131 26365 -13030    26059 1.8456      1     0.1743

# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ams_cdo_0 <- update(gca_mod_ams_wm_1,   . ~ . + DO_down) # singular
gca_mod_ams_cdo_1 <- update(gca_mod_ams_cdo_0, . ~ . + ot1:DO_down) # singular
gca_mod_ams_cdo_2 <- update(gca_mod_ams_cdo_1, . ~ . + ot2:DO_down) # singular
gca_mod_ams_cdo_3 <- update(gca_mod_ams_cdo_2, . ~ . + ot3:DO_down) # singular

ams_cdo_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_cdo_0, gca_mod_ams_cdo_1,
        gca_mod_ams_cdo_2, gca_mod_ams_cdo_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ams_wm_1  32 26129 26336 -13032    26065                         
# gca_mod_ams_cdo_0 33 26131 26345 -13032    26065 0.0028      1     0.9581
# gca_mod_ams_cdo_1 34 26133 26353 -13032    26065 0.0519      1     0.8198
# gca_mod_ams_cdo_2 35 26135 26362 -13032    26065 0.0009      1     0.9760
# gca_mod_ams_cdo_3 36 26136 26370 -13032    26064 0.6283      1     0.4280

gca_mod_ams_cup_0 <- update(gca_mod_ams_wm_1,   . ~ . + DO_up) # singular
gca_mod_ams_cup_1 <- update(gca_mod_ams_cup_0, . ~ . + ot1:DO_up) # singular
gca_mod_ams_cup_2 <- update(gca_mod_ams_cup_1, . ~ . + ot2:DO_up) # singular
gca_mod_ams_cup_3 <- update(gca_mod_ams_cup_2, . ~ . + ot3:DO_up) # singular

ams_cup_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_cup_0, gca_mod_ams_cup_1,
        gca_mod_ams_cup_2, gca_mod_ams_cup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_mod_ams_wm_1  32 26129 26336 -13032    26065                         
# gca_mod_ams_cup_0 33 26131 26345 -13032    26065 0.1274      1     0.7212
# gca_mod_ams_cup_1 34 26133 26353 -13032    26065 0.1633      1     0.6861
# gca_mod_ams_cup_2 35 26135 26362 -13032    26065 0.0018      1     0.9657
# gca_mod_ams_cup_3 36 26136 26369 -13032    26064 0.8292      1     0.3625

gca_mod_ams_edo_0 <- update(gca_mod_ams_wm_1,   . ~ . + MI_down) # singular
gca_mod_ams_edo_1 <- update(gca_mod_ams_edo_0, . ~ . + ot1:MI_down) # singular
gca_mod_ams_edo_2 <- update(gca_mod_ams_edo_1, . ~ . + ot2:MI_down) # singular
gca_mod_ams_edo_3 <- update(gca_mod_ams_edo_2, . ~ . + ot3:MI_down) # singular

ams_edo_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_edo_0, gca_mod_ams_edo_1,
        gca_mod_ams_edo_2, gca_mod_ams_edo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ams_wm_1  32 26129 26336 -13032    26065                         
# gca_mod_ams_edo_0 33 26131 26345 -13032    26065 0.2809      1     0.5961
# gca_mod_ams_edo_1 34 26132 26353 -13032    26064 0.3617      1     0.5476
# gca_mod_ams_edo_2 35 26134 26361 -13032    26064 0.2903      1     0.5900
# gca_mod_ams_edo_3 36 26135 26368 -13031    26063 1.2224      1     0.2689

gca_mod_ams_eup_0 <- update(gca_mod_ams_wm_1,   . ~ . + MI_up) # singular
gca_mod_ams_eup_1 <- update(gca_mod_ams_eup_0, . ~ . + ot1:MI_up) # singular
gca_mod_ams_eup_2 <- update(gca_mod_ams_eup_1, . ~ . + ot2:MI_up) # singular
gca_mod_ams_eup_3 <- update(gca_mod_ams_eup_2, . ~ . + ot3:MI_up) # singular

ams_eup_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_eup_0, gca_mod_ams_eup_1,
        gca_mod_ams_eup_2, gca_mod_ams_eup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ams_wm_1  32 26129 26336 -13032    26065                         
# gca_mod_ams_eup_0 33 26130 26344 -13032    26064 0.9963      1     0.3182
# gca_mod_ams_eup_1 34 26131 26352 -13032    26063 0.8077      1     0.3688
# gca_mod_ams_eup_2 35 26133 26360 -13031    26063 0.2957      1     0.5866
# gca_mod_ams_eup_3 36 26134 26368 -13031    26062 0.3444      1     0.5573

gca_mod_ams_gdo_0 <- update(gca_mod_ams_wm_1,   . ~ . + SOL_down) # singular
gca_mod_ams_gdo_1 <- update(gca_mod_ams_gdo_0, . ~ . + ot1:SOL_down) # singular
gca_mod_ams_gdo_2 <- update(gca_mod_ams_gdo_1, . ~ . + ot2:SOL_down) # singular
gca_mod_ams_gdo_3 <- update(gca_mod_ams_gdo_2, . ~ . + ot3:SOL_down) # singular

ams_gdo_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_gdo_0, gca_mod_ams_gdo_1,
        gca_mod_ams_gdo_2, gca_mod_ams_gdo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ams_wm_1  32 26129 26336 -13032    26065                         
# gca_mod_ams_gdo_0 33 26131 26345 -13032    26065 0.1592      1     0.6899
# gca_mod_ams_gdo_1 34 26133 26353 -13032    26065 0.0014      1     0.9705
# gca_mod_ams_gdo_2 35 26134 26361 -13032    26064 0.3078      1     0.5790
# gca_mod_ams_gdo_3 36 26136 26370 -13032    26064 0.4237      1     0.5151

gca_mod_ams_gup_0 <- update(gca_mod_ams_wm_1,   . ~ . + SOL_up) # singular
gca_mod_ams_gup_1 <- update(gca_mod_ams_gup_0, . ~ . + ot1:SOL_up) # singular
gca_mod_ams_gup_2 <- update(gca_mod_ams_gup_1, . ~ . + ot2:SOL_up) # singular
gca_mod_ams_gup_3 <- update(gca_mod_ams_gup_2, . ~ . + ot3:SOL_up) # singular

ams_gup_anova <-
  anova(gca_mod_ams_wm_1, gca_mod_ams_gup_0, gca_mod_ams_gup_1,
        gca_mod_ams_gup_2, gca_mod_ams_gup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ams_wm_1  32 26129 26336 -13032    26065                         
# gca_mod_ams_gup_0 33 26131 26345 -13032    26065 0.2018      1     0.6533
# gca_mod_ams_gup_1 34 26132 26353 -13032    26064 0.3811      1     0.5370
# gca_mod_ams_gup_2 35 26134 26361 -13032    26064 0.4636      1     0.4959
# gca_mod_ams_gup_3 36 26135 26369 -13032    26063 0.4820      1     0.4875



#
# only ims
#

gca_mod_ims_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ims")) # singular

# add cond effect to imsercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_cond_0 <- update(gca_mod_ims_base,   . ~ . + condition_sum) # singular
gca_mod_ims_cond_1 <- update(gca_mod_ims_cond_0,   . ~ . + ot1:condition_sum) # singular
gca_mod_ims_cond_2 <- update(gca_mod_ims_cond_1,   . ~ . + ot2:condition_sum) # singular
gca_mod_ims_cond_3 <- update(gca_mod_ims_cond_2,   . ~ . + ot3:condition_sum) # singular

ims_cond_anova <-
  anova(gca_mod_ims_base, gca_mod_ims_cond_0, gca_mod_ims_cond_1, # none singular, and none significant
        gca_mod_ims_cond_2, gca_mod_ims_cond_3)
#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_base   30 26777 26971 -13358    26717                           
# gca_mod_ims_cond_0 31 26777 26978 -13357    26715 2.0613      1    0.15108  
# gca_mod_ims_cond_1 32 26779 26986 -13357    26715 0.0909      1    0.76310  
# gca_mod_ims_cond_2 33 26778 26992 -13356    26712 2.9194      1    0.08752 .
# gca_mod_ims_cond_3 34 26779 27000 -13356    26711 0.3822      1    0.53642  


# add verbal WM effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_wm_0 <- update(gca_mod_ims_base,   . ~ . + WM_set) # singular
gca_mod_ims_wm_1 <- update(gca_mod_ims_wm_0, . ~ . + ot1:WM_set) # singular
gca_mod_ims_wm_2 <- update(gca_mod_ims_wm_1, . ~ . + ot2:WM_set) # singular
gca_mod_ims_wm_3 <- update(gca_mod_ims_wm_2, . ~ . + ot3:WM_set) # singular

ims_wm_anova <-
  anova(gca_mod_ims_base, gca_mod_ims_wm_0, gca_mod_ims_wm_1,
        gca_mod_ims_wm_2, gca_mod_ims_wm_3)
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_base 30 26777 26971 -13358    26717                           
# gca_mod_ims_wm_0 31 26777 26978 -13358    26715 1.7702      1    0.18336  
# gca_mod_ims_wm_1 32 26778 26986 -13357    26714 0.4471      1    0.50373  
# gca_mod_ims_wm_2 33 26780 26994 -13357    26714 0.7681      1    0.38082  
# gca_mod_ims_wm_3 34 26776 26996 -13354    26708 6.1125      1    0.01342 *

# add rhythm synchronization effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_spon_0 <- update(gca_mod_ims_wm_3,   . ~ . + spondee) # singular
gca_mod_ims_spon_1 <- update(gca_mod_ims_spon_0, . ~ . + ot1:spondee) # singular
gca_mod_ims_spon_2 <- update(gca_mod_ims_spon_1, . ~ . + ot2:spondee) # singular
gca_mod_ims_spon_3 <- update(gca_mod_ims_spon_2, . ~ . + ot3:spondee) # singular

ims_spon_anova <-
  anova(gca_mod_ims_wm_3, gca_mod_ims_spon_0, gca_mod_ims_spon_1,
        gca_mod_ims_spon_2, gca_mod_ims_spon_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_wm_3   34 26776 26996 -13354    26708                              
# gca_mod_ims_spon_0 35 26775 27002 -13352    26705  2.9666      1  0.0850021 .  
# gca_mod_ims_spon_1 36 26774 27008 -13351    26702  2.2845      1  0.1306705    
# gca_mod_ims_spon_2 37 26764 27004 -13345    26690 11.9674      1  0.0005414 ***
# gca_mod_ims_spon_3 38 26766 27012 -13345    26690  0.4673      1  0.4942376

gca_mod_ims_stre_0 <- update(gca_mod_ims_spon_2,   . ~ . + stressed_spondee) # singular
gca_mod_ims_stre_1 <- update(gca_mod_ims_stre_0, . ~ . + ot1:stressed_spondee) # singular
gca_mod_ims_stre_2 <- update(gca_mod_ims_stre_1, . ~ . + ot2:stressed_spondee) # singular
gca_mod_ims_stre_3 <- update(gca_mod_ims_stre_2, . ~ . + ot3:stressed_spondee) # singular

ims_stre_anova <-
  anova(gca_mod_ims_spon_2, gca_mod_ims_stre_0, gca_mod_ims_stre_1,
        gca_mod_ims_stre_2, gca_mod_ims_stre_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# gca_mod_ims_spon_2 37 26764 27004 -13345    26690                         
# gca_mod_ims_stre_0 38 26765 27011 -13344    26689 1.5936      1     0.2068
# gca_mod_ims_stre_1 39 26767 27020 -13344    26689 0.1390      1     0.7093
# gca_mod_ims_stre_2 40 26768 27028 -13344    26688 0.3065      1     0.5799
# gca_mod_ims_stre_3 41 26769 27035 -13344    26687 1.0726      1     0.3004

gca_mod_ims_troc_0 <- update(gca_mod_ims_spon_2,   . ~ . + trochee) # singular
gca_mod_ims_troc_1 <- update(gca_mod_ims_troc_0, . ~ . + ot1:trochee) # singular
gca_mod_ims_troc_2 <- update(gca_mod_ims_troc_1, . ~ . + ot2:trochee) # singular
gca_mod_ims_troc_3 <- update(gca_mod_ims_troc_2, . ~ . + ot3:trochee) # singular

ims_troc_anova <-
  anova(gca_mod_ims_spon_2, gca_mod_ims_troc_0, gca_mod_ims_troc_1,
        gca_mod_ims_troc_2, gca_mod_ims_troc_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_spon_2 37 26764 27004 -13345    26690                         
# gca_mod_ims_troc_0 38 26765 27012 -13345    26689 1.2426      1     0.2650
# gca_mod_ims_troc_1 39 26767 27020 -13345    26689 0.0051      1     0.9431
# gca_mod_ims_troc_2 40 26769 27028 -13344    26689 0.0788      1     0.7790
# gca_mod_ims_troc_3 41 26770 27036 -13344    26688 0.9450      1     0.3310

gca_mod_ims_unpr_0 <- update(gca_mod_ims_spon_2,   . ~ . + unpredictable) # singular
gca_mod_ims_unpr_1 <- update(gca_mod_ims_unpr_0, . ~ . + ot1:unpredictable) # singular
gca_mod_ims_unpr_2 <- update(gca_mod_ims_unpr_1, . ~ . + ot2:unpredictable) # singular
gca_mod_ims_unpr_3 <- update(gca_mod_ims_unpr_2, . ~ . + ot3:unpredictable) # singular

ims_unpr_anova <-
  anova(gca_mod_ims_spon_2, gca_mod_ims_unpr_0, gca_mod_ims_unpr_1,
        gca_mod_ims_unpr_2, gca_mod_ims_unpr_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_spon_2 37 26764 27004 -13345    26690                         
# gca_mod_ims_unpr_0 38 26765 27012 -13345    26689 1.2734      1     0.2591
# gca_mod_ims_unpr_1 39 26766 27018 -13344    26688 1.5809      1     0.2086
# gca_mod_ims_unpr_2 40 26766 27026 -13343    26686 1.0543      1     0.3045
# gca_mod_ims_unpr_3 41 26768 27034 -13343    26686 0.8641      1     0.3526

# add pitch anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ims_cdo_0 <- update(gca_mod_ims_spon_2,   . ~ . + DO_down) # singular
gca_mod_ims_cdo_1 <- update(gca_mod_ims_cdo_0, . ~ . + ot1:DO_down) # singular
gca_mod_ims_cdo_2 <- update(gca_mod_ims_cdo_1, . ~ . + ot2:DO_down) # singular
gca_mod_ims_cdo_3 <- update(gca_mod_ims_cdo_2, . ~ . + ot3:DO_down) # singular

ims_cdo_anova <-
  anova(gca_mod_ims_spon_2, gca_mod_ims_cdo_0, gca_mod_ims_cdo_1,
        gca_mod_ims_cdo_2, gca_mod_ims_cdo_3)
#                    Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ims_spon_2 37 26764 27004 -13345    26690                           
# gca_mod_ims_cdo_0  38 26765 27012 -13344    26689 1.4018      1    0.23642  
# gca_mod_ims_cdo_1  39 26767 27020 -13344    26689 0.0033      1    0.95419  
# gca_mod_ims_cdo_2  40 26768 27028 -13344    26688 0.5443      1    0.46067  
# gca_mod_ims_cdo_3  41 26766 27032 -13342    26684 4.2548      1    0.03914 *

gca_mod_ims_cup_0 <- update(gca_mod_ims_cdo_3,   . ~ . + DO_up) # singular
gca_mod_ims_cup_1 <- update(gca_mod_ims_cup_0, . ~ . + ot1:DO_up) # singular
gca_mod_ims_cup_2 <- update(gca_mod_ims_cup_1, . ~ . + ot2:DO_up) # singular
gca_mod_ims_cup_3 <- update(gca_mod_ims_cup_2, . ~ . + ot3:DO_up) # singular

ims_cup_anova <-
  anova(gca_mod_ims_cdo_3, gca_mod_ims_cup_0, gca_mod_ims_cup_1,
        gca_mod_ims_cup_2, gca_mod_ims_cup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_mod_ims_cdo_3 41 26766 27032 -13342    26684                           
# gca_mod_ims_cup_0 42 26767 27040 -13342    26683 0.9544      1    0.32860  
# gca_mod_ims_cup_1 43 26765 27044 -13340    26679 3.9103      1    0.04799 *
# gca_mod_ims_cup_2 44 26765 27050 -13338    26677 2.3516      1    0.12516  
# gca_mod_ims_cup_3 45 26767 27059 -13338    26677 0.0108      1    0.91708  

gca_mod_ims_edo_0 <- update(gca_mod_ims_cup_1,   . ~ . + MI_down) # singular
gca_mod_ims_edo_1 <- update(gca_mod_ims_edo_0, . ~ . + ot1:MI_down) # singular
gca_mod_ims_edo_2 <- update(gca_mod_ims_edo_1, . ~ . + ot2:MI_down) # singular
gca_mod_ims_edo_3 <- update(gca_mod_ims_edo_2, . ~ . + ot3:MI_down) # singular

ims_edo_anova <-
  anova(gca_mod_ims_cup_1, gca_mod_ims_edo_0, gca_mod_ims_edo_1,
        gca_mod_ims_edo_2, gca_mod_ims_edo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ims_cup_1 43 26765 27044 -13340    26679                         
# gca_mod_ims_edo_0 44 26767 27052 -13339    26679 0.4579      1     0.4986
# gca_mod_ims_edo_1 45 26768 27059 -13339    26678 1.3284      1     0.2491
# gca_mod_ims_edo_2 46 26769 27068 -13339    26677 0.1884      1     0.6643
# gca_mod_ims_edo_3 47 26771 27076 -13339    26677 0.0576      1     0.8103

gca_mod_ims_eup_0 <- update(gca_mod_ims_cup_1,   . ~ . + MI_up) # singular
gca_mod_ims_eup_1 <- update(gca_mod_ims_eup_0, . ~ . + ot1:MI_up) # singular
gca_mod_ims_eup_2 <- update(gca_mod_ims_eup_1, . ~ . + ot2:MI_up) # singular
gca_mod_ims_eup_3 <- update(gca_mod_ims_eup_2, . ~ . + ot3:MI_up) # singular

ims_eup_anova <-
  anova(gca_mod_ims_cup_1, gca_mod_ims_eup_0, gca_mod_ims_eup_1,
        gca_mod_ims_eup_2, gca_mod_ims_eup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ims_cup_1 43 26765 27044 -13340    26679                         
# gca_mod_ims_eup_0 44 26767 27053 -13340    26679 0.0768      1     0.7817
# gca_mod_ims_eup_1 45 26769 27061 -13340    26679 0.0158      1     0.8999
# gca_mod_ims_eup_2 46 26771 27070 -13340    26679 0.0201      1     0.8872
# gca_mod_ims_eup_3 47 26773 27078 -13340    26679 0.0186      1     0.8914


gca_mod_ims_gdo_0 <- update(gca_mod_ims_cup_1,   . ~ . + SOL_down) # singular
gca_mod_ims_gdo_1 <- update(gca_mod_ims_gdo_0, . ~ . + ot1:SOL_down) # singular
gca_mod_ims_gdo_2 <- update(gca_mod_ims_gdo_1, . ~ . + ot2:SOL_down) # singular
gca_mod_ims_gdo_3 <- update(gca_mod_ims_gdo_2, . ~ . + ot3:SOL_down) # singular

ims_gdo_anova <-
  anova(gca_mod_ims_cup_1, gca_mod_ims_gdo_0, gca_mod_ims_gdo_1,
        gca_mod_ims_gdo_2, gca_mod_ims_gdo_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ims_cup_1 43 26765 27044 -13340    26679                         
# gca_mod_ims_gdo_0 44 26767 27052 -13339    26679 0.6279      1     0.4281
# gca_mod_ims_gdo_1 45 26769 27061 -13339    26679 0.0009      1     0.9764
# gca_mod_ims_gdo_2 46 26771 27069 -13339    26679 0.0729      1     0.7871
# gca_mod_ims_gdo_3 47 26772 27077 -13339    26678 0.2664      1     0.6058

gca_mod_ims_gup_0 <- update(gca_mod_ims_cup_1,   . ~ . + SOL_up) # singular
gca_mod_ims_gup_1 <- update(gca_mod_ims_gup_0, . ~ . + ot1:SOL_up) # singular
gca_mod_ims_gup_2 <- update(gca_mod_ims_gup_1, . ~ . + ot2:SOL_up) # singular
gca_mod_ims_gup_3 <- update(gca_mod_ims_gup_2, . ~ . + ot3:SOL_up) # singular

ims_gup_anova <-
  anova(gca_mod_ims_cup_1, gca_mod_ims_gup_0, gca_mod_ims_gup_1,
        gca_mod_ims_gup_2, gca_mod_ims_gup_3)
#                   Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)  
# gca_mod_ims_cup_1 43 26765 27044 -13340    26679                         
# gca_mod_ims_gup_0 44 26767 27052 -13340    26679 0.3662      1     0.5451
# gca_mod_ims_gup_1 45 26768 27060 -13339    26678 1.1236      1     0.2891
# gca_mod_ims_gup_2 46 26770 27068 -13339    26678 0.0045      1     0.9463
# gca_mod_ims_gup_3 47 26771 27076 -13339    26677 0.6734      1     0.4119

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_full_mod_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +          
         (1 + stress_sum * (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 2e4)), 
       data = stress_gc_subset, REML = F, na.action = na.exclude)
# Warning message:
#   In commonArgs(par, fn, control, environment()) :
#   maxfun < 10 * length(par)^2 is not recommended.

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_group_0 <- update(gca_full_mod_base,    . ~ . + group)
gca_full_mod_group_1 <- update(gca_full_mod_group_0, . ~ . + ot1:group)
gca_full_mod_group_2 <- update(gca_full_mod_group_1, . ~ . + ot2:group)
gca_full_mod_group_3 <- update(gca_full_mod_group_2, . ~ . + ot3:group)

full_group_anova <-
  anova(gca_full_mod_base, gca_full_mod_group_0, gca_full_mod_group_1,
        gca_full_mod_group_2, gca_full_mod_group_3)
#                        Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# gca_full_mod_base      51 238356 238799 -119127   238254                          
# gca_full_mod_group_0   55 238335 238812 -119113   238225 29.3362  4  6.680e-06 ***
# gca_full_mod_group_1   59 238333 238845 -119108   238215  9.7470  4    0.04491 *  
# gca_full_mod_group_2   63 238307 238853 -119090   238181 34.5150  4  5.843e-07 ***
# gca_full_mod_group_3   67 238306 238887 -119086   238172  9.2752  4    0.05458 .

mod_type <- "gca_full_mod"
mod_spec <- c("_base", "_group_0", "_group_1", "_group_2", "_group_3")

group_mods <- mget(c(paste0(mod_type, mod_spec)))

save(group_mods,
     file = here("mods", "music", "gca",
                 "group_mods.Rdata"))

gca_full_mod_pitch <- update(gca_full_mod_group_2, . ~ . + pitch_dev)
gca_full_mod_pitch_1 <- update(gca_full_mod_pitch, . ~ . + ot1:pitch_dev)
gca_full_mod_pitch_2 <- update(gca_full_mod_pitch_1, . ~ . + ot2:pitch_dev)
gca_full_mod_pitch_3 <- update(gca_full_mod_pitch_2, . ~ . + ot3:pitch_dev)

full_pitch_anova <-
  anova(gca_full_mod_group_2, gca_full_mod_pitch, gca_full_mod_pitch_1, 
        gca_full_mod_pitch_2, gca_full_mod_pitch_3)





 * rhythm_dev * WM_set * 




# add 3-way int to intercept, linear slope, quadratic, and cubic time terms
gca_full_mod_int_0 <- update(gca_full_mod_group_3, . ~ . + stress_sum:pitch_dev:rhythm_dev:WM_set:group)
gca_full_mod_int_1 <- update(gca_full_mod_int_0,   . ~ . + ot1:coda_sum:stress_sum:pitch_dev:rhythm_dev:WM_set:group)
gca_full_mod_int_2 <- update(gca_full_mod_int_1,   . ~ . + ot2:coda_sum:stress_sum:pitch_dev:rhythm_dev:WM_set:group)
gca_full_mod_int_3 <- update(gca_full_mod_int_2,   . ~ . + ot3:coda_sum:stress_sum:pitch_dev:rhythm_dev:WM_set:group)

full_int_anova <-
  anova(gca_full_mod_group_3, gca_full_mod_int_0, gca_full_mod_int_1,
        gca_full_mod_int_2, gca_full_mod_int_3)

# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "int"))
gca_full_mod_int_relevel <- update(gca_full_mod_int_3)




}

# -----------------------------------------------------------------------------





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, condition_sum) %>%
  distinct

# Get model predictions and SE
fits_all <- predictSE(gca_full_mod_group_3, new_dat_all) %>%        #change depending on significance
  as_tibble %>%
  bind_cols(new_dat_all) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, SS = "mon", AE = "aes", IE = "ies", AM = "ams", IM = "ims"))

# Filter preds at target offset
target_offset_preds <- filter(fits_all, time_zero == 4) %>%
  select(group, cond = condition_sum,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {
# Build model names programatically
mod_type <- "gca_mod_"
mod_spec <- c("_base", "_cond_0", "_cond_1", "_cond_2", "_cond_3", 
              "_wm_0", "_wm_1", "_wm_2", "_wm_3",
              "_spon_0", "_spon_1", "_spon_2", "_spon_3",
              "_stre_0", "_stre_1", "_stre_2", "_stre_3",
              "_troc_0", "_troc_1", "_troc_2", "_troc_3",
              "_unpr_0", "_unpr_1", "_unpr_2", "_unpr_3",
              "_cdo_0", "_cdo_1", "_cdo_2", "_cdo_3",
              "_cup_0", "_cup_1", "_cup_2", "_cup_3",
              "_edo_0", "_edo_1", "_edo_2", "_edo_3",
              "_eup_0", "_eup_1", "_eup_2", "_eup_3",
              "_gdo_0", "_gdo_1", "_gdo_2", "_gdo_3",
              "_gup_0", "_gup_1", "_gup_2", "_gup_3")

# Store ind models in list
ind_mods <- mget(c(#paste0(mod_type, "mon", mod_spec),
                   paste0(mod_type, "aes", mod_spec),
                   paste0(mod_type, "ies", mod_spec),
                   paste0(mod_type, "ams", mod_spec),
                   paste0(mod_type, "ims", mod_spec)
                   ))

save(ind_mods,
     file = here("mods", "music", "gca",
                 "ind_mods.Rdata"))

# Store full (ot1, ot2, ot3, group, single effects) models in list
full_mods <- mget(c(
  "gca_full_mod_base", "gca_full_mod_group_0", "gca_full_mod_group_1",
  "gca_full_mod_group_2", "gca_full_mod_group_3", "gca_full_mod_unpr",
  "gca_full_mod_cup_0", "gca_full_mod_wm_1", "gca_full_mod_edo_1",
  "gca_full_mod_cup_1", "gca_full_mod_cdo_2", "gca_full_mod_spon_2",
  "gca_full_mod_wm_3", "gca_full_mod_cdo_3"))

save(full_mods,
     file = here("mods", "music", "gca",
                 "full_mods.Rdata"))

# Save anova model comparisons
nested_model_comparisons <-
  mget(c(#"mon_cond_anova", "mon_wm_anova", "mon_spon_anova", "mon_stre_anova",
         # "mon_troc_anova", "mon_unpr_anova", "mon_cdo_anova", "mon_cup_anova",
         # "mon_edo_anova", "mon_eup_anova", "mon_gdo_anova", "mon_gup_anova",
         
         "aes_cond_anova", "aes_wm_anova", "aes_spon_anova", "aes_stre_anova",
         "aes_troc_anova", "aes_unpr_anova", "aes_cdo_anova", "aes_cup_anova",
         "aes_edo_anova", "aes_eup_anova", "aes_gdo_anova", "aes_gup_anova",
         
         "ies_cond_anova", "ies_wm_anova", "ies_spon_anova", "ies_stre_anova",
         "ies_troc_anova", "ies_unpr_anova", "ies_cdo_anova", "ies_cup_anova",
         "ies_edo_anova", "ies_eup_anova", "ies_gdo_anova", "ies_gup_anova",
         
         "ams_cond_anova", "ams_wm_anova", "ams_spon_anova", "ams_stre_anova",
         "ams_troc_anova", "ams_unpr_anova", "ams_cdo_anova", "ams_cup_anova",
         "ams_edo_anova", "ams_eup_anova", "ams_gdo_anova", "ams_gup_anova",
         
         "ims_cond_anova", "ims_wm_anova", "ims_spon_anova", "ims_stre_anova",
         "ims_troc_anova", "ims_unpr_anova", "ims_cdo_anova", "ims_cup_anova",
         "ims_edo_anova", "ims_eup_anova", "ims_gdo_anova", "ims_gup_anova",
         
         "full_group_anova", "full_sineff_1_anova", "full_sineff_2_anova"))

save(nested_model_comparisons,
     file = here("mods", "music", "gca",
                 "nested_model_comparisons.Rdata"))









# Save models predictions
model_preds <- mget(c("fits_all", "target_offset_preds"))

save(model_preds,
     file = here("reports", "mods", "gca",
                 "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------

