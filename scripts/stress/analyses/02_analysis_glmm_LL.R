#
# Original script by Joseph
# Modified by Cristina to adapt to CrisLau Project
# Then by Laura 
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
#source(here::here("scripts", "02_load_data.R"))
stress50 <- read_csv(here("data", "clean", "stress_50ms_final_onsetc3updated.csv"))

# Get path to saved models
gca_mods_path  <- here("mods", "stress", "glmm", "LL_changes")

# Load models as lists
load(paste0(gca_mods_path, "/gca_mon_mods.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_l1ot2.Rdata")) #gca_l2_mods only dele parameter significant
#load(paste0(gca_mods_path, "/nested_model_l2.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata")) #mon
load(paste0(gca_mods_path, "/model_preds_l1ot2.Rdata")) #L2

# # Store objects in global env
# list2env(gca_mon_mods, globalenv())
# list2env(gca_l2_mods, globalenv())
# list2env(gca_l2_mods_sum, globalenv())
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
#    - We can select the appropriate time course subset by selecting the
#      target syllable onset, bin 4 (200ms / 50 = 4), and keeping an
#      equal number of bins on each side:
#                     8 7 6 5 4 3 2 1 X 1 2 3 4 5 6 7 8
#                                     ^
#                     center of time course (bin 4)
#
#
# Number of bins:     1  2  3  4 5 6 7 8 9 10 11 12 13 14 15 16 17
# Actual bin number: -4 -3 -2 -1 0 1 2 3 4  5  6  7  8  9 10 11 12



# stress_50 <- na.omit(stress50)


stress_glmm_subset <- stress50 %>%
  # select(., -WM_set) %>%
  filter(., time_zero == 4) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma"),
         condition_sum = if_else(cond == "1", -1, 1)
         )       # 1 = present (now -1), 2 = past (now 1)



# -----------------------------------------------------------------------------







#################### MONOLINGUAL SPEAKERS ########################################

# Build up random effects to test time terms
if(F){
  
  mon_data <- filter(stress_glmm_subset, l1 == 'es') %>% select(-DELE, -percent_l2_week)
  
  mod_ot1 <-
    lmer(eLog ~ 1 + 
           (1 | participant),
         data = mon_data, weights = 1/wts, REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . + (1 | target))
  
  mod_ot3 <- update(mod_ot2, . ~ . -(1 | participant) +
                      (1 + condition_sum | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    3 2574.6 2587.2 -1284.3   2568.6                          
  # mod_ot2    4 2562.7 2579.3 -1277.3   2554.7 13.9740  1  0.0001854 ***
  # mod_ot3    6 2561.0 2586.1 -1274.5   2549.0  5.6302  2  0.0598983 .

  
}



# Individual model MON -----------------------------------------------------------

glmm_mon_base <- mod_ot2
  # lmer(eLog ~ 1 +
  #        (1 + condition_sum | participant) +
  #        (1 | target),
  #      REML = F,
  #      data = filter(mon_data)) 

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
glmm_mon_cond <- update(glmm_mon_base,   . ~ . + condition_sum)

mon_cond_anova <-
  anova(glmm_mon_base, glmm_mon_cond)
#                 Df    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# glmm_mon_base    4 2562.7 2579.3 -1277.3   2554.7                     
# glmm_mon_cond    5 2564.4 2585.3 -1277.2   2554.4 0.2066  1     0.6494



mod_type <- "glmm_mon"
mod_spec <- c('_base', "_cond")
              

# Store ind models in list
glmm_mon_mods <- mget(c(paste0(mod_type, mod_spec)))

save(glmm_mon_mods,
     file = here("mods", "stress", "glmm", "LL_changes",
                 "glmm_mon_mods.Rdata"))

 

#################### L2 SPEAKERS ########################################

l2_data <- stress_glmm_subset%>%
  filter(., l1 != 'es') %>% 
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1))

##### PROFICIENCY

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){

  mod_ot1 <-
    lmer(eLog ~ 1 +
           (1 | participant),
         data = l2_data, weights = 1/wts, 
         control=lmerControl(optimizer="bobyqa",
                             optCtrl=list(maxfun=2e5)),
         REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . + (1 | target))
  
  mod_ot3 <- update(mod_ot2, . ~ . -(1 | participant) +
                    + (1 + condition_sum | participant))
                     
  mod_ot4 <- update(mod_ot3, . ~ . -(1 | target) +
                    + (1 + DELE_z + l1_sum | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3, mod_ot4)
  #          npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    3 11061 11078 -5527.6    11055                         
  # mod_ot2    4 11040 11062 -5516.0    11032 23.0987  1  1.539e-06 ***
  # mod_ot3    6 11027 11061 -5507.7    11015 16.7236  2  0.0002336 ***
  # mod_ot4   11 11029 11090 -5503.4    11007  8.4565  5  0.1328086  
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
# Base model
glmm_l2_dele_base <- mod_ot3
   #lmer(eLog ~ 1 +
   #       (1 + condition_sum | participant) +
   #       (1 | target), 
   #       control=lmerControl(optimizer="bobyqa",
   #                           optCtrl=list(maxfun=2e5)),
   #     data = l2_data, REML = F)


# add condition effect to intercept, linear slope, quadratic, and cubic time terms
glmm_l2_dele_condition <- update(glmm_l2_dele_base, . ~ . + condition_sum) 

# add group effect to intercept, linear slope, quadratic, and cubic time terms
glmm_l2_dele_l1 <- update(glmm_l2_dele_condition, . ~ . + l1_sum) 

# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
glmm_l2_dele_dele <- update(glmm_l2_dele_l1,   . ~ . + DELE_z) 

anova(glmm_l2_dele_base, glmm_l2_dele_condition, glmm_l2_dele_l1,
        glmm_l2_dele_dele)
#                        npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# glmm_l2_dele_base         6 11027 11061 -5507.7    11015                       
# glmm_l2_dele_condition    7 11026 11065 -5506.1    11012 3.2317  1    0.07223 .
# glmm_l2_dele_l1           8 11026 11070 -5504.8    11010 2.5229  1    0.11221  
# glmm_l2_dele_dele         9 11027 11078 -5504.6    11009 0.4898  1    0.48403

# Add 2-way interactions
glmm_l2_dele_int_1 <- update(glmm_l2_dele_base,     . ~ . + l1_sum:DELE_z) 
glmm_l2_dele_int_2 <- update(glmm_l2_dele_int_1, . ~ . + l1_sum:condition_sum) 
glmm_l2_dele_int_3 <- update(glmm_l2_dele_int_2, . ~ . + DELE_z:condition_sum)

anova(glmm_l2_dele_int_1, glmm_l2_dele_int_2, glmm_l2_dele_int_3)
#                    npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# glmm_l2_dele_int_1    7 11029 11068 -5507.5    11015                     
# glmm_l2_dele_int_2    8 11031 11076 -5507.4    11015 0.1237  1     0.7251
# glmm_l2_dele_int_3    9 11032 11082 -5507.0    11014 0.9035  1     0.3418

# add 3-way interaction 
glmm_l2_dele_int_4 <- update(glmm_l2_dele_int_3,     . ~ . + l1_sum:DELE_z:condition_sum) 

anova(glmm_l2_dele_int_3, glmm_l2_dele_int_4)
#                    npar   AIC   BIC  logLik deviance Chisq Df Pr(>Chisq)   
# glmm_l2_dele_int_3    9 11032 11082 -5507.0    11014                    
# glmm_l2_dele_int_4   10 11032 11088 -5505.9    11012  2.21  1     0.1371

glmm_l2_dele_final <- glmm_l2_dele_base

summary(glmm_l2_dele_final)

mod_type <- "glmm_l2_dele"
mod_spec <- c('_base', "_l1", "_dele", "_condition", 
              "_int_1", "_int_2", "_int_3", "_int_4", 
              '_final') 

# Store ind models in list
glmm_l2_mods_dele <- mget(c(paste0(mod_type, mod_spec)))

save(glmm_l2_mods_dele,
     file = here("mods", "stress", "glmm", "LL_changes",
                 "glmm_l2_mods_dele.Rdata"))



}

# -----------------------------------------------------------------------------







##### USE

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  
  mod_ot1 <-
    lmer(eLog ~ 1 +
           (1 | participant),
         data = l2_data, weights = 1/wts, 
         control=lmerControl(optimizer="bobyqa",
                             optCtrl=list(maxfun=2e5)),
         REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . + (1 | target))
  
  mod_ot3 <- update(mod_ot2, . ~ . -(1 | participant) +
                      + (1 + condition_sum | participant))
  
  mod_ot4 <- update(mod_ot2, . ~ . -(1 | target) +
                      + (1 + use_z + l1_sum | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3, mod_ot4)
  #          npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    3 11061 11078 -5527.6    11055                         
  # mod_ot2    4 11040 11062 -5516.0    11032 23.099  1  1.539e-06 ***
  # mod_ot3    6 11027 11061 -5507.7    11015 16.724  2  0.0002336 ***
  # mod_ot4    9 11044 11094 -5513.0    11026  0.000  3  1.0000000    
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
  # Base model
  glmm_l2_use_base <- mod_ot3
  #lmer(eLog ~ 1 +
  #       (1 + condition_sum | participant) +
  #       (1 | target), 
  #       control=lmerControl(optimizer="bobyqa",
  #                           optCtrl=list(maxfun=2e5)),
  #     data = l2_data, REML = F)
  
  
  # add condition effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_use_condition <- update(glmm_l2_use_base, . ~ . + condition_sum) 
  
  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_use_l1 <- update(glmm_l2_use_condition, . ~ . + l1_sum) 
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_use_use <- update(glmm_l2_use_l1,   . ~ . + use_z) 
  
  anova(glmm_l2_use_base, glmm_l2_use_condition, glmm_l2_use_l1,
        glmm_l2_use_use)
  #                       npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # glmm_l2_use_base         6 11027 11061 -5507.7    11015                       
  # glmm_l2_use_condition    7 11026 11065 -5506.1    11012 3.2317  1    0.07223 .
  # glmm_l2_use_l1           8 11026 11070 -5504.8    11010 2.5229  1    0.11221  
  # glmm_l2_use_use          9 11022 11073 -5502.3    11004 5.0452  1    0.02469 *
  
  # Add 2-way interactions
  glmm_l2_use_int_1 <- update(glmm_l2_use_use,     . ~ . + l1_sum:use_z) 
  glmm_l2_use_int_2 <- update(glmm_l2_use_int_1, . ~ . + l1_sum:condition_sum) 
  glmm_l2_use_int_3 <- update(glmm_l2_use_int_2, . ~ . + use_z:condition_sum)
  
  anova(glmm_l2_use_use, glmm_l2_use_int_1, glmm_l2_use_int_2, glmm_l2_use_int_3)
  #                   npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # glmm_l2_use_use      9 11022 11073 -5502.3    11004                     
  # glmm_l2_use_int_1   10 11024 11080 -5502.3    11004 0.0134  1     0.9079
  # glmm_l2_use_int_2   11 11026 11088 -5502.3    11004 0.0065  1     0.9356
  # glmm_l2_use_int_3   12 11036 11104 -5506.2    11012 0.0000  1     1.0000
  
  # add 3-way interaction 
  glmm_l2_use_int_4 <- update(glmm_l2_use_use,     . ~ . + l1_sum:use_z:condition_sum) 
  
  anova(glmm_l2_use_use, glmm_l2_use_int_4)
  #                   npar   AIC   BIC  logLik deviance Chisq Df Pr(>Chisq)   
  # glmm_l2_use_use      9 11022 11073 -5502.3    11004                     
  # glmm_l2_use_int_4   10 11024 11080 -5502.1    11004 0.3296  1     0.5659
  
  glmm_l2_use_final <- glmm_l2_use_use
  
  summary(glmm_l2_use_final)
  
  mod_type <- "glmm_l2_use"
  mod_spec <- c('_base', "_l1", "_use", "_condition", 
                "_int_1", "_int_2", "_int_3", "_int_4", 
                '_final') 
  
  # Store ind models in list
  glmm_l2_mods_use <- mget(c(paste0(mod_type, mod_spec)))
  
  save(glmm_l2_mods_use,
       file = here("mods", "stress", "glmm", "LL_changes",
                   "glmm_l2_mods_use.Rdata"))
  
  
  
}

# -----------------------------------------------------------------------------







##### L2 EXPERIENCE

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  
  mod_ot1 <-
    lmer(eLog ~ 1 +
           (1 | participant),
         data = l2_data, weights = 1/wts, 
         control=lmerControl(optimizer="bobyqa",
                             optCtrl=list(maxfun=2e5)),
         REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . + (1 | target))
  
  mod_ot3 <- update(mod_ot2, . ~ . -(1 | target) +
                      + (1 + use_z + DELE_z + l1_sum | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar    AIC   BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    3  11061 11078 -5527.6    11055                         
  # mod_ot2    4  11040 11062 -5516.0    11032  23.0987  1  1.539e-06 ***
  # mod_ot3   13  11050 11122 -5511.9    11024   8.2891  9     0.5053 
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
  # Base model
  glmm_l2_expl1_base <- mod_ot2
  #lmer(eLog ~ 1 +
  #       (1 | participant) +
  #       (1 | target), 
  #       control=lmerControl(optimizer="bobyqa",
  #                           optCtrl=list(maxfun=2e5)),
  #     data = l2_data, REML = F)
  
  # add use effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_expl1_use <- update(glmm_l2_expl1_base,   . ~ . + use_z) 
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_expl1_dele <- update(glmm_l2_expl1_use,   . ~ . + DELE_z) 
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_expl1_dele <- update(glmm_l2_expl1_use,   . ~ . + DELE_z) 
  
  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_expl1_l1 <- update(glmm_l2_expl1_dele, . ~ . + l1_sum) 
  
  anova(glmm_l2_expl1_base, glmm_l2_expl1_use, glmm_l2_expl1_dele, glmm_l2_expl1_l1)
  #                    npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # glmm_l2_expl1_base    4 11040 11062 -5516.0    11032                       
  # glmm_l2_expl1_use     5 11040 11068 -5514.8    11030 2.4431  1    0.11804  
  # glmm_l2_expl1_dele    6 11037 11071 -5512.5    11025 4.5704  1    0.03253 *
  # glmm_l2_expl1_l1      7 11036 11075 -5511.1    11022 2.8627  1    0.09065 .
  
  # Add 2-way interactions
  glmm_l2_expl1_int_1 <- update(glmm_l2_expl1_dele,  . ~ . + l1_sum:use_z) 
  glmm_l2_expl1_int_2 <- update(glmm_l2_expl1_int_1, . ~ . + l1_sum:DELE_z) 
  glmm_l2_expl1_int_3 <- update(glmm_l2_expl1_int_2, . ~ . + DELE_z:use_z)
  
  anova(glmm_l2_expl1_dele, glmm_l2_expl1_int_1, glmm_l2_expl1_int_2, glmm_l2_expl1_int_3)
  #                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # glmm_l2_expl1_dele     6 11037 11071 -5512.5    11025                     
  # glmm_l2_expl1_int_1    7 11038 11077 -5512.0    11024 0.9953  1     0.3185
  # glmm_l2_expl1_int_2    8 11040 11085 -5512.0    11024 0.1229  1     0.7259
  # glmm_l2_expl1_int_3    9 11042 11092 -5511.8    11024 0.3921  1     0.5312
  
  # add 3-way interaction 
  glmm_l2_expl1_int_4 <- update(glmm_l2_expl1_dele,     . ~ . + l1_sum:use_z:DELE_z) 
  
  anova(glmm_l2_expl1_dele, glmm_l2_expl1_int_4)
  #                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # glmm_l2_expl1_dele     6 11037 11071 -5512.5    11025                     
  # glmm_l2_expl1_int_4    7 11039 11078 -5512.4    11025 0.2173  1     0.6411
  
  glmm_l2_expl1_final <- glmm_l2_expl1_dele
  
  summary(glmm_l2_expl1_final)
  
  mod_type <- "glmm_l2_expl1"
  mod_spec <- c('_base', "_dele", "_use", "_l1", 
                "_int_1", "_int_2", "_int_3", "_int_4", 
                '_final') 
  
  # Store ind models in list
  glmm_l2_mods_expl1 <- mget(c(paste0(mod_type, mod_spec)))
  
  save(glmm_l2_mods_expl1,
       file = here("mods", "stress", "glmm", "LL_changes",
                   "glmm_l2_mods_expl1.Rdata"))
  
  
  
}






## with condition instead of L1

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  
  mod_ot1 <-
    lmer(eLog ~ 1 +
           (1 | participant),
         data = l2_data, weights = 1/wts, 
         control=lmerControl(optimizer="bobyqa",
                             optCtrl=list(maxfun=2e5)),
         REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . + (1 | target))
  
  mod_ot3 <- update(mod_ot2, . ~ . - (1 | participant) +
                      + (1 + condition_sum | participant))
  
  mod_ot4 <- update(mod_ot3, . ~ . -(1 | target) +
                      + (1 + use_z + DELE_z | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3, mod_ot4)
  #         npar    AIC   BIC  logLik deviance    Chisq Df Pr(>Chisq)    
  # mod_ot1    3  11061 11078 -5527.6    11055                         
  # mod_ot2    4  11040 11062 -5516.0    11032  23.0987  1  1.539e-06 ***
  # mod_ot3    6  11027 11061 -5507.7    11015  16.7236  2  0.0002336 ***
  # mod_ot4   11  11036 11097 -5506.7    11014   1.8573  5  0.8685112
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
  # Base model
  glmm_l2_expcond_base <- mod_ot3
  #lmer(eLog ~ 1 +
  #       (1 + condition_sum | participant) +
  #       (1 | target), 
  #       control=lmerControl(optimizer="bobyqa",
  #                           optCtrl=list(maxfun=2e5)),
  #     data = l2_data, REML = F)
  
  # add use effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_expcond_use <- update(glmm_l2_expcond_base,   . ~ . + use_z) 
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_expcond_dele <- update(glmm_l2_expcond_use,   . ~ . + DELE_z) 
  
  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  glmm_l2_expcond_condition <- update(glmm_l2_expcond_dele, . ~ . + condition_sum) 
  
  anova(glmm_l2_expcond_base, glmm_l2_expcond_use, 
        glmm_l2_expcond_dele, glmm_l2_expcond_condition)
  #                           npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # glmm_l2_expcond_base         6 11027 11061 -5507.7    11015                         
  # glmm_l2_expcond_use          7 11033 11072 -5509.7    11019  0.000  1  1.0000000    
  # glmm_l2_expcond_dele         8 11023 11068 -5503.7    11007 11.954  1  0.0005454 ***
  # glmm_l2_expcond_condition    9 11028 11079 -5505.1    11010  0.000  1  1.0000000 
  
  # Add 2-way interactions
  glmm_l2_expcond_int_1 <- update(glmm_l2_expcond_dele,  . ~ . + condition_sum:use_z) 
  glmm_l2_expcond_int_2 <- update(glmm_l2_expcond_int_1, . ~ . + condition_sum:DELE_z) 
  glmm_l2_expcond_int_3 <- update(glmm_l2_expcond_int_2, . ~ . + DELE_z:use_z)
  
  anova(glmm_l2_expcond_dele, glmm_l2_expcond_int_1, glmm_l2_expcond_int_2, glmm_l2_expcond_int_3)
  #                       npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
  # glmm_l2_expcond_dele     8 11023 11068 -5503.7    11007                     
  # glmm_l2_expcond_int_1    9 11024 11075 -5503.2    11006 0.9969  1     0.3181
  # glmm_l2_expcond_int_2   10 11024 11080 -5502.0    11004 2.3114  1     0.1284
  # glmm_l2_expcond_int_3   11 11026 11087 -5501.9    11004 0.2059  1     0.6500
  
  # add 3-way interaction 
  glmm_l2_expcond_int_4 <- update(glmm_l2_expcond_dele,     . ~ . + condition_sum:use_z:DELE_z) 
  
  anova(glmm_l2_expcond_dele, glmm_l2_expcond_int_4)
  #                       npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # glmm_l2_expcond_dele     8 11023 11068 -5503.7    11007                     
  # glmm_l2_expcond_int_4    9 11025 11076 -5503.7    11007 0.0085  1     0.9268
  
  glmm_l2_expcond_final <- glmm_l2_expcond_dele
  
  summary(glmm_l2_expcond_final)
  
  mod_type <- "glmm_l2_expcond"
  mod_spec <- c('_base', "_dele", "_use", "_l1", 
                "_int_1", "_int_2", "_int_3", "_int_4", 
                '_final') 
  
  # Store ind models in list
  glmm_l2_mods_expcond <- mget(c(paste0(mod_type, mod_spec)))
  
  save(glmm_l2_mods_expcond,
       file = here("mods", "stress", "glmm", "LL_changes",
                   "glmm_l2_mods_expcond.Rdata"))
  
  
  
}

# -----------------------------------------------------------------------------







# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_mon <- mon_data %>%
  dplyr::select(time_zero, ot1:ot2) %>% 
  distinct
  
# Get model predictions and SE
fits_all_mon <- predictSE(gca_mon_mods$gca_mod_mon_base, new_dat_mon) %>%  
  as_tibble %>%
  bind_cols(new_dat_mon) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at first syllable offset
target_offset_preds_mon <- filter(fits_all_mon, time_zero == 4) %>%
  select(elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))


# -----------------------------------------------------------------------------


new_l2_dat <- l2_data %>%
  dplyr::select(l1_sum, time_zero, ot1:ot2) %>%
  distinct %>%
  # mutate(l1_sum = as.character(l1_sum)) %>% 
  expand_grid(., tibble(DELE_z = c(-1, 0, 1))) %>%
  expand_grid(., tibble(use_z = c(-1, 0, 1)))

fits_all_l2 <- predictSE(gca_l2_mod_final, new_l2_dat) %>%
  as_tibble %>%
  bind_cols(new_l2_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at first syllable offset
target_offset_preds_l2 <- filter(fits_all_l2, time_zero == 4) %>% #
  select(l1 = l1_sum, DELE = DELE_z, `L2 use` = use_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))
 


model_preds <- mget(c("fits_all_mon", 
                      # "fits_all_l2", 
  "target_offset_preds_mon" 
  # "target_offset_preds_l2"
  ))

save(model_preds,
     file = here("mods", "stress", "gca", "LL_changes",
                 "model_preds.Rdata"))



fits_all_l2_ints <- predictSE(gca_l2_mod_full, new_l2_dat) %>%
  as_tibble %>%
  bind_cols(new_l2_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

target_onset_preds_l2_ints <- filter(fits_all_l2_ints, time_zero == 4) %>% #
  select(L1 = l1_sum, DELE = DELE_z, `L2 use` = use_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))

model_preds <- mget(c("fits_all_l2_ints", 
                      "target_onset_preds_l2_ints"))

save(model_preds,
     file = here("mods", "stress", "gca", "LL_changes",
                 "model_preds_ints.Rdata"))





# -------------------------------------------------------------------------------------------------------
gca_l2_base <- mod_ot5
#lmer(eLog ~ 1 + (ot1 + ot2) +         
#       (1 + ot1 + ot2 | participant) +
#       (1 + DELE_z + use_z + l1_sum + ot1 | target), 
#     control = lmerControl(optimizer = 'bobyqa',
#                           optCtrl = list(maxfun = 2e4)),
#     data = l2_data, REML = F)

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_l2_l1_0 <- update(gca_l2_base, . ~ . + l1_sum) 
gca_l2_l1_1 <- update(gca_l2_l1_0, . ~ . + ot1:l1_sum) 
gca_l2_l1_2 <- update(gca_l2_l1_1, . ~ . + ot2:l1_sum) 

l2_l1_anova_all <-
  anova(gca_l2_base, gca_l2_l1_0, gca_l2_l1_1,
        gca_l2_l1_2)
#                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_l2_base   25 100675 100870 -50312   100625                    
# gca_l2_l1_0   26 100677 100880 -50313   100625 0.000  1     1.0000
# gca_l2_l1_1   27 100676 100887 -50311   100622 2.690  1     0.1010
# gca_l2_l1_2   28 100678 100896 -50311   100622 0.618  1     0.4318


# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_dele_0 <- update(gca_l2_base,   . ~ . + DELE_z) 
gca_l2_dele_1 <- update(gca_l2_dele_0, . ~ . + ot1:DELE_z) 
gca_l2_dele_2 <- update(gca_l2_dele_1, . ~ . + ot2:DELE_z)

l2_dele_anova_all <-
  anova(gca_l2_base, gca_l2_dele_0, gca_l2_dele_1,
        gca_l2_dele_2)
#               npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_l2_base     25 100675 100870 -50312   100625                       
# gca_l2_dele_0   26 100674 100876 -50311   100622 3.0880  1    0.07887 .
# gca_l2_dele_1   27 100669 100880 -50308   100615 6.4988  1    0.01079 *
# gca_l2_dele_2   28 100671 100889 -50308   100615 0.0801  1    0.77722 



# add L2 use effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_use_0 <- update(gca_l2_dele_1,  . ~ . + use_z) 
gca_l2_use_1 <- update(gca_l2_use_0, . ~ . + ot1:use_z) 
gca_l2_use_2 <- update(gca_l2_use_1, . ~ . + ot2:use_z)

l2_use_anova_all <-
  anova(gca_l2_dele_1, gca_l2_use_0, gca_l2_use_1,
        gca_l2_use_2)
#               npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
# gca_l2_dele_1   27 100669 100880 -50308   100615                         
# gca_l2_use_0    28 100682 100900 -50313   100626  0.000  1  1.0000000    
# gca_l2_use_1    29 100669 100895 -50306   100611 14.961  1  0.0001098 ***
# gca_l2_use_2    30 100673 100906 -50306   100613  0.000  1  1.0000000




# Add interactions one by one
gca_l2_i_0 <- update(gca_l2_use_1,  . ~ . + l1_sum:DELE_z) 
gca_l2_i_1 <- update(gca_l2_i_0, . ~ . + ot1:l1_sum:DELE_z) 
gca_l2_i_2 <- update(gca_l2_i_1, . ~ . + ot2:l1_sum:DELE_z)

anova(gca_l2_use_1, gca_l2_i_0, gca_l2_i_1,
      gca_l2_i_2)
#                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_l2_use_1   29 100669 100895 -50306   100611                     
# gca_l2_i_0     30 100674 100908 -50307   100614 0.0000  1     1.0000
# gca_l2_i_1     31 100675 100917 -50307   100613 0.5501  1     0.4583
# gca_l2_i_2     32 100677 100926 -50306   100613 0.3818  1     0.5366


gca_l2_in_0 <- update(gca_l2_use_1,   . ~ . + l1_sum:use_z) 
gca_l2_in_1 <- update(gca_l2_in_0, . ~ . + ot1:l1_sum:use_z) 
gca_l2_in_2 <- update(gca_l2_in_1, . ~ . + ot2:l1_sum:use_z)

anova(gca_l2_use_1, gca_l2_in_0, gca_l2_in_1,
      gca_l2_in_2)
#              npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_l2_use_1   29 100669 100895 -50306   100611                         
# gca_l2_in_0    30 100675 100909 -50307   100615  0.000  1  1.0000000    
# gca_l2_in_1    31 100681 100923 -50310   100619  0.000  1  1.0000000    
# gca_l2_in_2    32 100672 100921 -50304   100608 11.418  1  0.0007275 ***



gca_l2_it_0 <- update(gca_l2_in_2,     . ~ . + DELE_z:use_z) 
gca_l2_it_1 <- update(gca_l2_it_0, . ~ . + ot1:DELE_z:use_z) 
gca_l2_it_2 <- update(gca_l2_it_1, . ~ . + ot2:DELE_z:use_z)

anova(gca_l2_in_2, gca_l2_it_0, gca_l2_it_1,
      gca_l2_it_2)
#                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_l2_in_2   32 100672 100921 -50304   100608                          
# gca_l2_it_0   33 100685 100942 -50309   100619  0.0000  1  1.0000000    
# gca_l2_it_1   34 100673 100938 -50303   100605 13.7579  1  0.0002079 ***
# gca_l2_it_2   35 100675 100948 -50302   100605  0.4534  1  0.5007372 


gca_l2_nt_0 <- update(gca_l2_it_1,     . ~ . + l1_sum:DELE_z:use_z) 
gca_l2_nt_1 <- update(gca_l2_nt_0, . ~ . + ot1:l1_sum:DELE_z:use_z) 
gca_l2_nt_2 <- update(gca_l2_nt_1, . ~ . + ot2:l1_sum:DELE_z:use_z)

anova(gca_l2_it_1, gca_l2_nt_0, gca_l2_nt_1,
      gca_l2_nt_2)
#             npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
# gca_l2_it_1   34 100673 100938 -50303   100605                         
# gca_l2_nt_0   35 100675 100948 -50303   100605  0.180  1  0.6713390    
# gca_l2_nt_1   36 100684 100964 -50306   100612  0.000  1  1.0000000    
# gca_l2_nt_2   37 100674 100963 -50300   100600 11.485  1  0.0007018 ***

mod_type <- "gca_l2"
mod_spec <- c('_base', 
              "_l1_0", "_l1_1", "_l1_2", 
              "_dele_0", "_dele_1", "_dele_2", 
              "_use_0", "_use_1", "_use_2", 
              "_i_0", "_i_1", "_i_2", # Dele x l1
              "_in_0", "_in_1", "_in_2", #use x l1
              "_it_0", "_it_1", "_it_2", # dele x use
              "_nt_0", "_nt_1", "_nt_2")

# Store ind models in list
gca_l2 <- mget(c(paste0(mod_type, mod_spec)))

save(gca_l2,
     file = here("mods", "stress", "gca", "LL_changes",
                 "gca_l2_l1ot2.Rdata"))

fits_all_l2_l1ot2 <- predictSE(gca_l2_nt_2, new_l2_dat) %>%
  as_tibble %>%
  bind_cols(new_l2_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

target_onset_preds_l2_l1ot2 <- filter(fits_all_l2_l1ot2, time_zero == 4) %>% #
  select(L1 = l1_sum, `Spanish proficiency` = DELE_z, `Spanish use` = use_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))

model_preds_l1ot2 <- mget(c("fits_all_l2_l1ot2", 
                      "target_onset_preds_l2_l1ot2"))

save(model_preds_l1ot2,
     file = here("mods", "stress", "gca", "LL_changes",
                 "model_preds_l1ot2.Rdata"))


car::vif(gca_l2_nt_2)
# ot1                     ot2                  DELE_z                   use_z 
# 1.145707                1.025393                1.178469                1.300949 
# ot1:DELE_z               ot1:use_z            use_z:l1_sum            DELE_z:use_z 
# 1.345421                1.566706                1.562949                1.417314 
# ot1:use_z:l1_sum        ot2:use_z:l1_sum        ot1:DELE_z:use_z     DELE_z:use_z:l1_sum 
# 1.615377                1.368916                1.673724                1.473754 
# ot1:DELE_z:use_z:l1_sum ot2:DELE_z:use_z:l1_sum 
# 1.525169                1.354092 





