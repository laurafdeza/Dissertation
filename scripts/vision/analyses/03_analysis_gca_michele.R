#
#
# Growth curve analysis ------------------------------------------------------
#
# - Question 1: Do visuospatial prediction abilities (continuous) influence 
#   Spanish speakers' abilities to predict verbal tense
#   based on the presence or absence (categorical) of lexical stress?
# - Question 2: Do verbal and visuospatial processing speed influence linguistic prediction?
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

# Load data
source(here::here("scripts", "00_load_libs.R"))
# source(here::here("scripts", "02_load_data.R"))
stress50 <- read_csv(here("data", "clean", "stress_50ms_final_onsetc3updated.csv"))


# Get path to saved models
gca_mods_path  <- here("mods", "vision", "gca", "cont_speed_verb")

# Load models as lists
load(paste0(gca_mods_path, "/mon_mods_zero0.Rdata")) # gca_mon_ospan_int_1, gca_mon_corsirt_int_2
load(paste0(gca_mods_path, "/mon_mods_zerosc.Rdata")) # gca_mon_wm_in_1, gca_mon_carsc_1

load(paste0(gca_mods_path, "/model_preds_zero0.Rdata"))
# load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))

# Store objects in global env
#list2env(mon_mods, globalenv())
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

wm <- read_csv("./data/clean/wm_processing_speed.csv")
# vision <- read_csv("./data/clean/vision_scores_nooutliers.csv") # pred car
# corsi <- read_csv("./data/clean/corsi_z_scores.csv")
wm_score <- read_csv("./data/clean/ospan_set_z_scores.csv")

vision50 <- left_join(x = stress50, y = wm, by = "participant", all.x=TRUE)
# vision50 <- left_join(x = vision50, y = vision, by = "participant", all.x=TRUE)
# vision50 <- left_join(x = vision50, y = corsi, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = wm_score, by = "participant", all.x=TRUE)


vision50 <- vision50 %>%
  filter(., time_zero >= -4 & time_zero <= 4) %>%
  mutate(., #l1 = fct_relevel(l1, "es", "en", "ma"),
            stress_sum = if_else(cond == "1", -1, 1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 2, prefix = "ot")


# -----------------------------------------------------------------------------





# MON

mon_vision <- filter(vision50, l1 == 'es') %>% select(-DELE, -percent_l2_week, 
                                                      -prof, -group)

mon_vision <- na.omit(mon_vision)


# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 + ot2 +
           (1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = mon_vision, weights = 1/wts, REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . + (1 | target))
  
  mod_ot3 <- update(mod_ot2, . ~ . -(1 | target) +
                      + (1 + ospan | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # mod_ot1    5 23204 23236 -11597    23194                         
  # mod_ot2    6 23126 23164 -11557    23114 79.958  1  < 2.2e-16 ***
  # mod_ot3    8 23101 23152 -11543    23085 28.588  2  6.198e-07 ***
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_mon_base <- mod_ot3
  # lmer(eLog ~ 1 + ot1 + ot2 +    
  #        (1 | participant) +
  #        (1 + ospan | target),
  #      control = lmerControl(optimizer = 'bobyqa'), #, optCtrl = list(maxfun = 3e5)
  #      data = stress_gc_subset, REML = F)    # , na.action = na.exclude 

# add stress effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_stress_0 <- update(gca_mon_base,    . ~ . + stress_sum)
gca_mon_stress_1 <- update(gca_mon_stress_0, . ~ . + ot1:stress_sum)
gca_mon_stress_2 <- update(gca_mon_stress_1, . ~ . + ot2:stress_sum)

mon_stress_anova <-
  anova(gca_mon_base, gca_mon_stress_0, gca_mon_stress_1,
        gca_mon_stress_2)
                 # npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base        8 23101 23152 -11543    23085                     
# gca_mon_stress_0    9 23103 23160 -11543    23085 0.0854  1     0.7701
# gca_mon_stress_1   10 23103 23167 -11542    23083 1.8184  1     0.1775
# gca_mon_stress_2   11 23105 23175 -11542    23083 0.1812  1     0.6704


# add verbal wm effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_wm_0 <- update(gca_mon_base,    . ~ . + ospan)
gca_mon_wm_1 <- update(gca_mon_wm_0, . ~ . + ot1:ospan)
gca_mon_wm_2 <- update(gca_mon_wm_1, . ~ . + ot2:ospan)

mon_wm_anova <-
  anova(gca_mon_base, gca_mon_wm_0, gca_mon_wm_1,
        gca_mon_wm_2)
#                       npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base    8 23101 23152 -11543    23085                       
# gca_mon_wm_0    9 23103 23160 -11543    23085 0.1282  1    0.72035  
# gca_mon_wm_1   10 23102 23166 -11541    23082 3.0255  1    0.08197 .
# gca_mon_wm_2   11 23102 23172 -11540    23080 1.8075  1    0.17881


# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_wm_int_0 <- update(gca_mon_base, . ~ . + stress_sum:ospan)
gca_mon_wm_int_1 <- update(gca_mon_wm_int_0,   . ~ . + ot1:stress_sum:ospan)
gca_mon_wm_int_2 <- update(gca_mon_wm_int_1,   . ~ . + ot2:stress_sum:ospan)

mon_wm_int_anova <-
  anova(gca_mon_base, gca_mon_wm_int_0, gca_mon_wm_int_1,
        gca_mon_wm_int_2)
#                  npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mon_base        8 23101 23152 -11543    23085                       
# gca_mon_wm_int_0    9 23102 23159 -11542    23084 1.2179  1    0.26978  
# gca_mon_wm_int_1   10 23100 23163 -11540    23080 4.1818  1    0.04086 *
# gca_mon_wm_int_2   11 23101 23171 -11540    23079 0.4481  1    0.50325  


summary(gca_mon_wm_int_1)
#                      Estimate Std. Error t value
# (Intercept)           0.66153    0.11672   5.668
# ot1                   2.37318    0.13968  16.990
# ot2                   0.56377    0.13590   4.148
# stress_sum:ospan     -0.08615    0.07102  -1.213
# ot1:stress_sum:ospan -0.28411    0.13868  -2.049

car::vif(gca_mon_wm_int_1)
# ot1                     ot2           stress_sum:ospan 
# 1.016263             1.017342             1.003383 
# ot1:stress_sum:ospan 
# 1.007339



# save models  
mod_type <- "gca_mon"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", 
              "_wm_0", "_wm_1", "_wm_2", 
              "_wm_int_0", "_wm_int_1", "_wm_int_2" 
              )

# Store ind models in list
mon_mods_michele <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods_michele,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "mon_mods_michele.Rdata"))

}

# ---------------------------------------------





# EN

en_vision <- filter(vision50, l1 == 'en') %>% select(-percent_l2_week, 
                                                      -prof, -group) %>%
  mutate(DELE_z = (DELE - mean(DELE))/sd(DELE))

en_vision <- na.omit(en_vision)


# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 + ot2 +
           (1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = en_vision, weights = 1/wts, REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . + (1 | target))
  
  mod_ot3 <- update(mod_ot2, . ~ . -(1 | target) +
                      + (1 + ospan | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)
  # mod_ot1    5 49274 49309 -24632    49264                          
  # mod_ot2    6 49167 49210 -24578    49155 108.537  1  < 2.2e-16 ***
  # mod_ot3    8 49114 49171 -24549    49098  57.323  2  3.569e-13 ***
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_en_base <- mod_ot3
  # lmer(eLog ~ 1 + ot1 + ot2 +    
  #        (1 | participant) +
  #        (1 + ospan | target),
  #      control = lmerControl(optimizer = 'bobyqa'), #, optCtrl = list(maxfun = 3e5)
  #      data = stress_gc_subset, REML = F)    # , na.action = na.exclude 
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_stress_0 <- update(gca_en_base,    . ~ . + stress_sum)
  gca_en_stress_1 <- update(gca_en_stress_0, . ~ . + ot1:stress_sum)
  gca_en_stress_2 <- update(gca_en_stress_1, . ~ . + ot2:stress_sum)
  
  en_stress_anova <-
    anova(gca_en_base, gca_en_stress_0, gca_en_stress_1,
          gca_en_stress_2)
  #                 npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_base        8 49114 49171  -24549    49098                        
  # gca_en_stress_0    9 49114 49177  -24548    49096 2.5294  1   0.111743   
  # gca_en_stress_1   10 49108 49178  -24544    49088 7.9330  1   0.004854 **
  # gca_en_stress_2   11 49109 49186  -24543    49087 1.0503  1   0.305432
  
  
  # add verbal wm effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_wm_0 <- update(gca_en_stress_1,    . ~ . + ospan)
  gca_en_wm_1 <- update(gca_en_wm_0, . ~ . + ot1:ospan)
  gca_en_wm_2 <- update(gca_en_wm_1, . ~ . + ot2:ospan)
  
  en_wm_anova <-
    anova(gca_en_stress_1, gca_en_wm_0, gca_en_wm_1,
          gca_en_wm_2)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_en_stress_1   10 49108 49178 -24544    49088                     
  # gca_en_wm_0       11 49110 49187 -24544    49088 0.0620  1     0.8033
  # gca_en_wm_1       12 49110 49195 -24543    49086 1.4364  1     0.2307
  # gca_en_wm_2       13 49112 49204 -24543    49086 0.0941  1     0.7590
  
  
  # add prof effect to intercept, linear slope, quadratic, and cubic time terms
  gca_en_prof_0 <- update(gca_en_stress_1,    . ~ . + DELE_z)
  gca_en_prof_1 <- update(gca_en_prof_0, . ~ . + ot1:DELE_z)
  gca_en_prof_2 <- update(gca_en_prof_1, . ~ . + ot2:DELE_z)
  
  en_prof_anova <-
    anova(gca_en_stress_1, gca_en_prof_0, gca_en_prof_1,
          gca_en_prof_2)
  #                 npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)
  # gca_en_stress_1   10 49108 49178 -24544    49088                          
  # gca_en_prof_0     11 49109 49187 -24544    49087  0.5719  1     0.4495    
  # gca_en_prof_1     12 49094 49179 -24535    49070 17.1027  1  3.541e-05 ***
  # gca_en_prof_2     13 49096 49188 -24535    49070  0.1307  1     0.7177
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_wm_i_0 <- update(gca_en_prof_1, . ~ . + stress_sum:ospan)
  gca_en_wm_i_1 <- update(gca_en_wm_i_0,   . ~ . + ot1:stress_sum:ospan)
  gca_en_wm_i_2 <- update(gca_en_wm_i_1,   . ~ . + ot2:stress_sum:ospan)
  
  en_wm_i_anova <-
    anova(gca_en_prof_1, gca_en_wm_i_0, gca_en_wm_i_1,
          gca_en_wm_i_2)
  #               npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_prof_1   12 49094 49179  -24535    49070                     
  # gca_en_wm_i_0   13 49096 49188  -24535    49070 0.2123  1     0.6450
  # gca_en_wm_i_1   14 49096 49195  -24534    49068 1.4891  1     0.2224
  # gca_en_wm_i_2   15 49098 49204  -24534    49068 0.6762  1     0.4109

  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_wm_in_0 <- update(gca_en_prof_1, . ~ . + stress_sum:DELE_z)
  gca_en_wm_in_1 <- update(gca_en_wm_in_0,   . ~ . + ot1:stress_sum:DELE_z)
  gca_en_wm_in_2 <- update(gca_en_wm_in_1,   . ~ . + ot2:stress_sum:DELE_z)
  
  en_wm_in_anova <-
    anova(gca_en_prof_1, gca_en_wm_in_0, gca_en_wm_in_1,
          gca_en_wm_in_2)
  #                npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_prof_1    12 49094 49179 -24535    49070                       
  # gca_en_wm_in_0   13 49092 49184 -24533    49066 4.3166  1    0.03774 *
  # gca_en_wm_in_1   14 49091 49190 -24532    49063 2.1696  1    0.14076  
  # gca_en_wm_in_2   15 49093 49199 -24532    49063 0.2091  1    0.64744
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_wm_nt_0 <- update(gca_en_wm_in_0, . ~ . + ospan:DELE_z)
  gca_en_wm_nt_1 <- update(gca_en_wm_nt_0,   . ~ . + ot1:ospan:DELE_z)
  gca_en_wm_nt_2 <- update(gca_en_wm_nt_1,   . ~ . + ot2:ospan:DELE_z)
  
  en_wm_nt_anova <-
    anova(gca_en_wm_in_0, gca_en_wm_nt_0, gca_en_wm_nt_1,
          gca_en_wm_nt_2)
  #                npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_en_wm_in_0   13 49092 49184 -24533    49066                     
  # gca_en_wm_nt_0   14 49093 49192 -24532    49065 0.7281  1     0.3935
  # gca_en_wm_nt_1   15 49094 49200 -24532    49064 0.7250  1     0.3945
  # gca_en_wm_nt_2   16 49096 49209 -24532    49064 0.0052  1     0.9424
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_en_wm_int_0 <- update(gca_en_wm_in_0, . ~ . + stress_sum:ospan:DELE_z)
  gca_en_wm_int_1 <- update(gca_en_wm_int_0,   . ~ . + ot1:stress_sum:ospan:DELE_z)
  gca_en_wm_int_2 <- update(gca_en_wm_int_1,   . ~ . + ot2:stress_sum:ospan:DELE_z)
  
  en_wm_int_anova <-
    anova(gca_en_wm_in_0, gca_en_wm_int_0, gca_en_wm_int_1,
          gca_en_wm_int_2)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_en_wm_in_0    13 49092 49184 -24533    49066                       
  # gca_en_wm_int_0   14 49090 49189 -24531    49062 3.4319  1    0.06395 .
  # gca_en_wm_int_1   15 49088 49194 -24529    49058 4.2491  1    0.03927 *
  # gca_en_wm_int_2   16 49088 49202 -24528    49056 1.5914  1    0.20712
  
  
  summary(gca_en_wm_int_1)
  #                             Estimate Std. Error t value
  # (Intercept)                  0.04624    0.09548   0.484
  # ot1                          1.36392    0.09537  14.302
  # ot2                          0.72638    0.09645   7.532
  # stress_sum                   0.11952    0.07590   1.575
  # DELE_z                       0.02905    0.05714   0.508
  # ot1:stress_sum               0.32326    0.09712   3.328
  # ot1:DELE_z                   0.40361    0.09789   4.123
  # stress_sum:DELE_z            0.05170    0.03285   1.574
  # stress_sum:DELE_z:ospan     -0.04526    0.02856  -1.585
  # ot1:stress_sum:DELE_z:ospan -0.17599    0.08531  -2.063
  
  car::vif(gca_en_wm_int_1)
  # ot1                         ot2 
  # 1.004665                    1.008748 
  # stress_sum                      DELE_z 
  # 1.005863                    1.008943 
  # ot1:stress_sum                  ot1:DELE_z 
  # 1.043842                    1.010844 
  # stress_sum:DELE_z     stress_sum:DELE_z:ospan 
  # 1.058220                    1.080031 
  # ot1:stress_sum:DELE_z:ospan 
  # 1.061725 
  
  
  
  # save models  
  mod_type <- "gca_en"
  mod_spec <- c("_base", 
                "_stress_0", "_stress_1", "_stress_2", 
                "_prof_0", "_prof_1", "_prof_2",
                "_wm_0", "_wm_1", "_wm_2", 
                "_wm_i_0", "_wm_i_1", "_wm_i_2",
                "_wm_in_0", "_wm_in_1", "_wm_in_2",
                "_wm_nt_0", "_wm_nt_1", "_wm_nt_2",
                "_wm_int_0", "_wm_int_1", "_wm_int_2" 
  )
  
  # Store ind models in list
  en_mods_michele <- mget(c(paste0(mod_type, mod_spec)
  ))
  
  save(en_mods_michele,
       file = here("mods", "vision", "gca", "cont_speed_verb",
                   "en_mods_michele.Rdata"))
  
}






# Ma CH

ma_vision <- filter(vision50, l1 == 'ma') %>% select(-percent_l2_week, 
                                                     -prof, -group) %>%
  mutate(DELE_z = (DELE - mean(DELE))/sd(DELE))

ma_vision <- na.omit(ma_vision)


# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 + ot2 +
           (1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = ma_vision, weights = 1/wts, REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . + (1 | target))
  
  anova(mod_ot1, mod_ot2)
  #         npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)
  # mod_ot1    5 51616 51651 -25803    51606                         
  # mod_ot2    6 51587 51630 -25788    51575 30.711  1  2.994e-08 ***
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_ma_base <- mod_ot2
  # lmer(eLog ~ 1 + ot1 + ot2 +    
  #        (1 | participant) +
  #        (1 | target),
  #      control = lmerControl(optimizer = 'bobyqa'), #, optCtrl = list(maxfun = 3e5)
  #      data = stress_gc_subset, REML = F)    # , na.action = na.exclude 
  
  # add stress effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_stress_0 <- update(gca_ma_base,    . ~ . + stress_sum)
  gca_ma_stress_1 <- update(gca_ma_stress_0, . ~ . + ot1:stress_sum)
  gca_ma_stress_2 <- update(gca_ma_stress_1, . ~ . + ot2:stress_sum)
  
  ma_stress_anova <-
    anova(gca_ma_base, gca_ma_stress_0, gca_ma_stress_1,
          gca_ma_stress_2)
  #                 npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_base        6 51587 51630 -25788    51575                       
  # gca_ma_stress_0    7 51583 51633 -25785    51569 5.7097  1    0.01687 *
  # gca_ma_stress_1    8 51585 51642 -25784    51569 0.3658  1    0.54530  
  # gca_ma_stress_2    9 51587 51651 -25784    51569 0.0544  1    0.81565 
  
  
  # add verbal wm effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_wm_0 <- update(gca_ma_stress_0,    . ~ . + ospan)
  gca_ma_wm_1 <- update(gca_ma_wm_0, . ~ . + ot1:ospan)
  gca_ma_wm_2 <- update(gca_ma_wm_1, . ~ . + ot2:ospan)
  
  ma_wm_anova <-
    anova(gca_ma_stress_0, gca_ma_wm_0, gca_ma_wm_1,
          gca_ma_wm_2)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_ma_stress_0    7 51583 51633 -25785    51569                     
  # gca_ma_wm_0        8 51585 51642 -25785    51569 0.1777  1     0.6733
  # gca_ma_wm_1        9 51587 51651 -25785    51569 0.0274  1     0.8685
  # gca_ma_wm_2       10 51589 51660 -25784    51569 0.1474  1     0.7011
  
  
  # add prof effect to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_prof_0 <- update(gca_ma_stress_0,    . ~ . + DELE_z)
  gca_ma_prof_1 <- update(gca_ma_prof_0, . ~ . + ot1:DELE_z)
  gca_ma_prof_2 <- update(gca_ma_prof_1, . ~ . + ot2:DELE_z)
  
  ma_prof_anova <-
    anova(gca_ma_stress_0, gca_ma_prof_0, gca_ma_prof_1,
          gca_ma_prof_2)
  #                 npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)
  # gca_ma_stress_0    7 51583 51633 -25785    51569                       
  # gca_ma_prof_0      8 51580 51637 -25782    51564  5.0739  1    0.02429 *
  # gca_ma_prof_1      9 51581 51645 -25781    51563  1.4728  1    0.22490  
  # gca_ma_prof_2     10 51583 51654 -25781    51563  0.2190  1    0.63979 
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_wm_i_0 <- update(gca_ma_prof_0, . ~ . + stress_sum:ospan)
  gca_ma_wm_i_1 <- update(gca_ma_wm_i_0,   . ~ . + ot1:stress_sum:ospan)
  gca_ma_wm_i_2 <- update(gca_ma_wm_i_1,   . ~ . + ot2:stress_sum:ospan)
  
  ma_wm_i_anova <-
    anova(gca_ma_prof_0, gca_ma_wm_i_0, gca_ma_wm_i_1,
          gca_ma_wm_i_2)
  #               npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_ma_prof_0    8 51580 51637 -25782    51564                     
  # gca_ma_wm_i_0    9 51580 51644 -25781    51562 2.5305  1     0.1117
  # gca_ma_wm_i_1   10 51581 51652 -25780    51561 0.8718  1     0.3505
  # gca_ma_wm_i_2   11 51581 51659 -25779    51559 2.2578  1     0.1329
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_wm_in_0 <- update(gca_ma_prof_0, . ~ . + stress_sum:DELE_z)
  gca_ma_wm_in_1 <- update(gca_ma_wm_in_0,   . ~ . + ot1:stress_sum:DELE_z)
  gca_ma_wm_in_2 <- update(gca_ma_wm_in_1,   . ~ . + ot2:stress_sum:DELE_z)
  
  ma_wm_in_anova <-
    anova(gca_ma_prof_0, gca_ma_wm_in_0, gca_ma_wm_in_1,
          gca_ma_wm_in_2)
  #                npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_ma_prof_0     8 51580 51637 -25782    51564                     
  # gca_ma_wm_in_0    9 51582 51646 -25782    51564 0.6210  1     0.4307
  # gca_ma_wm_in_1   10 51584 51655 -25782    51564 0.0213  1     0.8841
  # gca_ma_wm_in_2   11 51586 51664 -25782    51564 0.0746  1     0.7847
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_wm_nt_0 <- update(gca_ma_prof_0, . ~ . + ospan:DELE_z)
  gca_ma_wm_nt_1 <- update(gca_ma_wm_nt_0,   . ~ . + ot1:ospan:DELE_z)
  gca_ma_wm_nt_2 <- update(gca_ma_wm_nt_1,   . ~ . + ot2:ospan:DELE_z)
  
  ma_wm_nt_anova <-
    anova(gca_ma_prof_0, gca_ma_wm_nt_0, gca_ma_wm_nt_1,
          gca_ma_wm_nt_2)
  #                npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_ma_prof_0     8 51580 51637 -25782    51564                     
  # gca_ma_wm_nt_0    9 51582 51646 -25782    51564 0.2727  1     0.6015
  # gca_ma_wm_nt_1   10 51583 51654 -25782    51563 0.9374  1     0.3330
  # gca_ma_wm_nt_2   11 51585 51663 -25781    51563 0.2695  1     0.6037
  
  
  # add 2-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_ma_wm_int_0 <- update(gca_ma_prof_0, . ~ . + stress_sum:ospan:DELE_z)
  gca_ma_wm_int_1 <- update(gca_ma_wm_int_0,   . ~ . + ot1:stress_sum:ospan:DELE_z)
  gca_ma_wm_int_2 <- update(gca_ma_wm_int_1,   . ~ . + ot2:stress_sum:ospan:DELE_z)
  
  ma_wm_int_anova <-
    anova(gca_ma_prof_0, gca_ma_wm_int_0, gca_ma_wm_int_1,
          gca_ma_wm_int_2)
  #                 npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_ma_prof_0      8 51580 51637 -25782    51564                     
  # gca_ma_wm_int_0    9 51581 51645 -25782    51563 1.3452  1     0.2461
  # gca_ma_wm_int_1   10 51582 51653 -25781    51562 1.3193  1     0.2507
  # gca_ma_wm_int_2   11 51583 51661 -25780    51561 1.0616  1     0.3029
  
  
  car::vif(gca_ma_prof_0)
  #  ot1        ot2 stress_sum     DELE_z 
  # 1.000577   1.000703   1.001258   1.000017 
  
  
  
  # save models  
  mod_type <- "gca_ma"
  mod_spec <- c("_base", 
                "_stress_0", "_stress_1", "_stress_2", 
                "_prof_0", "_prof_1", "_prof_2",
                "_wm_0", "_wm_1", "_wm_2", 
                "_wm_i_0", "_wm_i_1", "_wm_i_2",
                "_wm_in_0", "_wm_in_1", "_wm_in_2",
                "_wm_nt_0", "_wm_nt_1", "_wm_nt_2",
                "_wm_int_0", "_wm_int_1", "_wm_int_2" 
  )
  
  # Store ind models in list
  ma_mods_michele <- mget(c(paste0(mod_type, mod_spec)
  ))
  
  save(ma_mods_michele,
       file = here("mods", "vision", "gca", "cont_speed_verb",
                   "ma_mods_michele.Rdata"))
  
}












  
# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_wm <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot2, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(ospan = c(-1, 0, 1)))

# Get model predictions and SE
fits_mon_wm <- predictSE(gca_mon_wm_int_1, mon_wm) %>%        
  as_tibble %>%
  bind_cols(mon_wm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_mon_wm<- filter(fits_mon_wm, time_zero == 4) %>%
  select(stress = stress_sum, ospan, 
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 




# Create design dataframe for predictions
en_wm <- en_vision %>%
  dplyr::select(time_zero, ot1:ot2, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(DELE_z = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan = c(-1, 0, 1)))

# Get model predictions and SE
fits_en_wm <- predictSE(gca_en_wm_int_1, en_wm) %>%        
  as_tibble %>%
  bind_cols(en_wm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_en_wm <- filter(fits_en_wm, time_zero == 4) %>%
  select(stress = stress_sum, ospan, DELE_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 




# Create design dataframe for predictions
ma_wm <- ma_vision %>%
  dplyr::select(time_zero, ot1:ot2, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(DELE_z = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan = c(-1, 0, 1)))

# Get model predictions and SE
fits_ma_wm <- predictSE(gca_ma_prof_0, ma_wm) %>%        
  as_tibble %>%
  bind_cols(ma_wm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_ma_wm <- filter(fits_ma_wm, time_zero == 4) %>%
  select(stress = stress_sum, ospan, DELE_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 





# Save models predictions
model_preds_michele <- mget(c('fits_mon_wm', 'fits_en_wm',
                            'fits_ma_wm', 
                      "preds_mon_wm", "preds_en_wm",
                      "preds_ma_wm"
                      ))

save(model_preds_michele,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "model_preds_michele.Rdata"))



# -----------------------------------------------------------------------------






# Save models -----------------------------------------------------------------

if(F) {
  
  # Save anova model comparisons
  # nested_model_comparisons <-
  #   mget(c("mon_stress_anova", "mon_car_anova", "mon_corsi_anova",
  #          "mon_car_int_anova", "mon_corsi_int_anova"
  #   ))
  # 
  # save(nested_model_comparisons,
  #      file = here("mods", "vision", "gca", "continuous",
  #                  "nested_model_comparisons.Rdata"))
  
  
  
  
  
  
}

# -----------------------------------------------------------------------------

