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
gca_mods_path  <- here("mods", "stress", "gca", "LL_changes")

# Load models as lists
load(paste0(gca_mods_path, "/gca_mon_simple.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_dele_simple.Rdata")) 
load(paste0(gca_mods_path, "/gca_l2_use_simple.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_exp_simple.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata")) #mon    



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

stress50 %>%
  filter(time_zero > -5 & time_zero < 5) %>%
  ggplot(., aes(x = time_zero, y = target_prop, color = l1, fill = l1)) +
  geom_vline(xintercept = 4, lty = 3) +
  geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
               stroke = 0.5, pch = 21, show.legend = FALSE) +
  scale_color_discrete(name="L1",
                       breaks = c('es', 'en', 'ma'),
                       labels = c('Spanish', "English", 'Chinese')) +
  scale_x_continuous(breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8),
                     labels = c("-400", "-300", "-200", "-100", "0", "100", "200", "300", "400")) +
  # ggtitle("Time course per verbal tense") +
  # xlab("Time in 50 ms bins (0 = marker time before accounting for 200 ms processing)") +
  labs(x = "Time relative to final syllable onset (ms)",
       y = "Proportion of fixations on target",
       caption = "Mean +/- 95% CI") +
  # annotate("text", x = 3.65, y = 0.53, label = '200ms',
  #                       angle = 90, size = 3, hjust = 0) +
  theme_grey(base_size = 10, base_family = "Times") +
  theme(legend.position = 'bottom')



stress_gc_subset <- stress50 %>%
  # select(., -WM_set) %>%
  filter(., time_zero >= -4 & time_zero <= 4) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma"),
         condition_sum = if_else(cond == "1", -1, 1)) %>%       # 1 = present (now -1), 2 = past (now 1)
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")



# -----------------------------------------------------------------------------







#################### MONOLINGUAL SPEAKERS ########################################

# Build up random effects to test time terms
if(F){
  
  mon_data <- filter(stress_gc_subset, l1 == 'es') %>% select(-DELE, -percent_l2_week)
  
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),   # , optCtrl=list(maxfun=2e5)
         data = mon_data, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . + ot2)
  
  mod_ot3 <- update(mod_ot2, . ~ . + (1 | target))
  
  anova(mod_ot1, mod_ot2, mod_ot3) 
  #         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    4 24019 24045 -12006    24011                         
  # mod_ot2    5 24003 24035 -11996    23993 18.626  1   1.59e-05 ***
  # mod_ot3    6 23918 23957 -11953    23906 86.305  1  < 2.2e-16 ***
  
}



# Fixed effects -----------------------------------------------------------

gca_mon_base <- mod_ot3
# lmer(eLog ~ 1 + ot1 + ot2 +
#        (1 | participant) +
#        (1 | target),
#      control = lmerControl(optimizer = 'bobyqa'), 
#      REML = F,
#      data = filter(mon_data)) 

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_cond_0 <- update(gca_mon_base,   . ~ . + condition_sum)
gca_mon_cond_1 <- update(gca_mon_cond_0,   . ~ . + ot1:condition_sum)
gca_mon_cond_2 <- update(gca_mon_cond_1,   . ~ . + ot2:condition_sum)

mon_cond_anova <-
  anova(gca_mon_base, gca_mon_cond_0, gca_mon_cond_1,
        gca_mon_cond_2)
#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mon_base      6 23918 23957 -11953    23906                     
# gca_mon_cond_0    7 23920 23965 -11953    23906 0.0005  1     0.9826
# gca_mon_cond_1    8 23920 23971 -11952    23904 2.2549  1     0.1332
# gca_mon_cond_2    9 23922 23979 -11952    23904 0.0914  1     0.7625

gca_mon_final <- gca_mon_base


summary(gca_mon_final)
confint(gca_mon_final)
#                           2.5 %     97.5 %
# .sig01                0.4029908 0.7433551
# .sig02      0.1894206 0.4524980
# .sigma      2.9857553 3.1154799
# (Intercept) 0.6479812 0.6601830
# ot1         2.3821537 2.3972766
# ot2         0.5431761 0.5582188

car::vif(gca_mon_final)
# ot1                  ot2  
# 1.014012        1.014012 
performance::r2_nakagawa(gca_mon_final, by_group = FALSE, tolerance = 1e-05)
# R2 for Mixed Models

# Conditional R2: 0.102
# Marginal R2: 0.064

rsq::rsq(gca_mon_final)
# $model
# [1] 0.1127688

# $fixed
# [1] 0.05992289
# 
# $random
# [1] 0.05284591
rsq::rsq(gca_mon_final, adj = TRUE)
# $model
# [1] 0.1123569

# $fixed
# [1] 0.05948646
# 
# $random
# [1] 0.05287044


mod_type <- "gca_mon"
mod_spec <- c('_base',
              '_cond_0', '_cond_1', '_cond_2',
              '_final')


# Store ind models in list
gca_mon_simple <- mget(c(paste0(mod_type, mod_spec)))

save(gca_mon_simple,
     file = here("mods", "stress", "gca", "LL_changes",
                 "gca_mon_simple.Rdata"))



#################### L2 SPEAKERS ########################################

l2_data <- stress_gc_subset%>%
  filter(., l1 != 'es') %>% 
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1),
         use_z = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week),
         DELE_z = (DELE - mean(DELE))/sd(DELE) #,
         #ospan = (WM_set - mean(WM_set))/sd(WM_set)
         ) 


### PROFICIENCY

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 | participant),
         control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=2e5)),
         data = l2_data, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . + ot2 )
  
  anova(mod_ot1, mod_ot2) 
  
  #         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    4 100955 100986 -50474   100947                         
  # mod_ot2    5 100883 100922 -50437   100873 73.681  1  < 2.2e-16 ***
  
  mod_ot3 <- update(mod_ot2, . ~ . + (1 | target))
  
  anova(mod_ot2, mod_ot3) 
  #         npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot2    5 100883 100922 -50437   100873                         
  # mod_ot3    6 100766 100813 -50377   100754 119.44  1  < 2.2e-16 ***

  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_l2_dele_base <- mod_ot3
    # lmer(eLog ~ 1 + (ot1 + ot2) +
    #        (1 | participant) +
    #        (1 | target),
    #      control = lmerControl(optimizer = 'bobyqa',
    #                            optCtrl = list(maxfun = 2e5)),
    #      data = l2_data, REML = F)
  
  # add condition effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_dele_cond_0 <- update(gca_l2_dele_base, . ~ . + condition_sum) 
  gca_l2_dele_cond_1 <- update(gca_l2_dele_cond_0, . ~ . + ot1:condition_sum) 
  gca_l2_dele_cond_2 <- update(gca_l2_dele_cond_1, . ~ . + ot2:condition_sum) 
  
  l2_dele_cond_anova <-
    anova(gca_l2_dele_base, gca_l2_dele_cond_0, gca_l2_dele_cond_1,
          gca_l2_dele_cond_2)
  #                    npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_dele_base      6 100766 100813 -50377   100754                       
  # gca_l2_dele_cond_0    7 100763 100818 -50375   100749 4.6014  1    0.03195 *
  # gca_l2_dele_cond_1    8 100760 100823 -50372   100744 5.0641  1    0.02443 *
  # gca_l2_dele_cond_2    9 100762 100832 -50372   100744 0.0445  1    0.83293
  
  
  
  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_dele_l1_0 <- update(gca_l2_dele_cond_1, . ~ . + l1_sum) 
  gca_l2_dele_l1_1 <- update(gca_l2_dele_l1_0, . ~ . + ot1:l1_sum) 
  gca_l2_dele_l1_2 <- update(gca_l2_dele_l1_1, . ~ . + ot2:l1_sum) 
  
  l2_dele_l1_anova <-
    anova(gca_l2_dele_cond_1, gca_l2_dele_l1_0, gca_l2_dele_l1_1,
          gca_l2_dele_l1_2)
  #                    npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_dele_cond_1    8 100760 100823 -50372   100744                       
  # gca_l2_dele_l1_0      9 100762 100832 -50372   100744 0.0180  1    0.89316  
  # gca_l2_dele_l1_1     10 100761 100839 -50371   100741 2.8480  1    0.09149 .
  # gca_l2_dele_l1_2     11 100763 100849 -50371   100741 0.3314  1    0.56486
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_dele_dele_0 <- update(gca_l2_dele_cond_1,   . ~ . + DELE_z) 
  gca_l2_dele_dele_1 <- update(gca_l2_dele_dele_0, . ~ . + ot1:DELE_z) 
  gca_l2_dele_dele_2 <- update(gca_l2_dele_dele_1, . ~ . + ot2:DELE_z)
  
  l2_dele_dele_anova <-
    anova(gca_l2_dele_cond_1, gca_l2_dele_dele_0, gca_l2_dele_dele_1,
          gca_l2_dele_dele_2)
  #                    npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_dele_cond_1    8 100760 100823 -50372   100744                          
  # gca_l2_dele_dele_0    9 100758 100828 -50370   100740  3.9628  1  0.0465165 *  
  # gca_l2_dele_dele_1   10 100748 100826 -50364   100728 12.2044  1  0.0004768 ***
  # gca_l2_dele_dele_2   11 100750 100836 -50364   100728  0.1284  1  0.7201023 
  
  
  # Add interactions one by one
  gca_l2_dele_l1dele_0 <- update(gca_l2_dele_dele_1,     . ~ . + l1_sum:DELE_z) 
  gca_l2_dele_l1dele_1 <- update(gca_l2_dele_l1dele_0, . ~ . + ot1:l1_sum:DELE_z) 
  gca_l2_dele_l1dele_2 <- update(gca_l2_dele_l1dele_1, . ~ . + ot2:l1_sum:DELE_z)
  
  anova(gca_l2_dele_dele_1, gca_l2_dele_l1dele_0, gca_l2_dele_l1dele_1,
        gca_l2_dele_l1dele_2)
  #                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_dele_dele_1     10 100748 100826 -50364   100728                       
  # gca_l2_dele_l1dele_0   11 100748 100833 -50363   100726 2.3043  1    0.12901  
  # gca_l2_dele_l1dele_1   12 100746 100840 -50361   100722 3.5009  1    0.06134 .
  # gca_l2_dele_l1dele_2   13 100748 100849 -50361   100722 0.0914  1    0.76244 
  
  
  gca_l2_dele_l1cond_0 <- update(gca_l2_dele_dele_1,     . ~ . + l1_sum:condition_sum) 
  gca_l2_dele_l1cond_1 <- update(gca_l2_dele_l1cond_0, . ~ . + ot1:l1_sum:condition_sum) 
  gca_l2_dele_l1cond_2 <- update(gca_l2_dele_l1cond_1, . ~ . + ot2:l1_sum:condition_sum)
  
  anova(gca_l2_dele_dele_1, gca_l2_dele_l1cond_0, gca_l2_dele_l1cond_1,
        gca_l2_dele_l1cond_2)
  #                      npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_l2_dele_dele_1     10 100748 100826 -50364   100728                     
  # gca_l2_dele_l1cond_0   11 100748 100834 -50363   100726 2.1913  1     0.1388
  # gca_l2_dele_l1cond_1   12 100747 100841 -50362   100723 2.6351  1     0.1045
  # gca_l2_dele_l1cond_2   13 100749 100850 -50361   100723 0.5506  1     0.4581
  
  
  
  gca_l2_dele_delecond_0 <- update(gca_l2_dele_dele_1,     . ~ . + DELE_z:condition_sum) 
  gca_l2_dele_delecond_1 <- update(gca_l2_dele_delecond_0, . ~ . + ot1:DELE_z:condition_sum) 
  gca_l2_dele_delecond_2 <- update(gca_l2_dele_delecond_1, . ~ . + ot2:DELE_z:condition_sum)
  
  anova(gca_l2_dele_dele_1, gca_l2_dele_delecond_0, gca_l2_dele_delecond_1,
        gca_l2_dele_delecond_2)
  #                        npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_dele_dele_1       10 100748 100826 -50364   100728                       
  # gca_l2_dele_delecond_0   11 100744 100830 -50361   100722 5.6304  1    0.01765 *
  # gca_l2_dele_delecond_1   12 100744 100838 -50360   100720 2.1640  1    0.14127  
  # gca_l2_dele_delecond_2   13 100746 100848 -50360   100720 0.0002  1    0.99015  
  
  
  
  gca_l2_dele_int_0 <- update(gca_l2_dele_delecond_0, . ~ . + l1_sum:DELE_z:condition_sum) 
  gca_l2_dele_int_1 <- update(gca_l2_dele_int_0, . ~ . + ot1:l1_sum:DELE_z:condition_sum) 
  gca_l2_dele_int_2 <- update(gca_l2_dele_int_1, . ~ . + ot2:l1_sum:DELE_z:condition_sum)
  
  anova(gca_l2_dele_delecond_0, gca_l2_dele_int_0, gca_l2_dele_int_1,
        gca_l2_dele_int_2)
  #                        npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_dele_delecond_0   11 100744 100830 -50361   100722                     
  # gca_l2_dele_int_0        12 100744 100838 -50360   100720 2.0115  1     0.1561
  # gca_l2_dele_int_1        13 100746 100847 -50360   100720 0.6762  1     0.4109
  # gca_l2_dele_int_2        14 100747 100856 -50360   100719 0.5234  1     0.4694
  
  
  gca_l2_dele_final <- gca_l2_dele_delecond_0
  
  summary(gca_l2_dele_final)
  confint(gca_l2_dele_final) #default "profile": computing a likelihood profile and finding the appropriate cutoffs based on the likelihood ratio test
  #                           2.5 %     97.5 %
  # .sig01               0.28452345 0.41807962
  # .sig02               0.22706334 0.41424777
  # .sigma               3.15076191 3.21697749
  # (Intercept)          0.04531402 0.05106354
  # ot1                  1.18130835 1.19027634
  # ot2                  0.60815022 0.61410468
  # condition_sum        0.10815292 0.11262853
  # DELE_z               0.07175546 0.07506920
  # ot1:condition_sum    0.15303488 0.15877913
  # ot1:DELE_z           0.23000472 0.23578460
  # condition_sum:DELE_z 0.05117113 0.05305324
  confint(gca_l2_dele_final, method = 'Wald') # approximating the confidence intervals (of fixed-effect parameters only; all variance-covariance parameters CIs will be returned as NA) based on the estimated local curvature of the likelihood surface
  #                             2.5 %     97.5 %
  # .sig01                         NA         NA
  # .sig02                         NA         NA
  # .sigma                         NA         NA
  # (Intercept)          -0.083457436 0.17455785
  # ot1                   1.052928723 1.31079063
  # ot2                   0.478966130 0.73783732
  # condition_sum         0.006786088 0.20988390
  # DELE_z               -0.002495342 0.14627839
  # ot1:condition_sum     0.024337246 0.28220326
  # ot1:DELE_z            0.100244024 0.36023850
  # condition_sum:DELE_z  0.008921262 0.09357506
  
  car::vif(gca_l2_dele_final)
  # ot1                  ot2        condition_sum 
  # 1.002369             1.002530             1.001593 
  # DELE_z    ot1:condition_sum           ot1:DELE_z 
  # 1.003313             1.003838             1.004037 
  # condition_sum:DELE_z 
  # 1.000949
  performance::r2_nakagawa(gca_l2_dele_final, by_group = FALSE, tolerance = 1e-05)
  # R2 for Mixed Models
  
  # Conditional R2: 0.041
  # Marginal R2: 0.021
  
  rsq::rsq(gca_l2_dele_final)
  # $model
  # [1] 0.06484025
  # 
  # $fixed
  # [1] 0.03632173
  # 
  # $random
  # [1] 0.02851852
  rsq::rsq(gca_l2_dele_final, adj = TRUE)
  # $model
  # [1] 0.06447496
  # 
  # $fixed
  # [1] 0.03594529
  # 
  # $random
  # [1] 0.02852966
  
  
  mod_type <- "gca_l2_dele"
  mod_spec <- c('_base',
                '_cond_0', '_cond_1', '_cond_2',
                '_l1_0', '_l1_1', '_l1_2',
                '_dele_0', '_dele_1', '_dele_2',
                '_cond_0', '_cond_1', '_cond_2',
                '_l1dele_0', '_l1dele_1', '_l1dele_2',
                '_l1cond_0', '_l1cond_1', '_l1cond_2',
                '_delecond_0', '_delecond_1', '_delecond_2',
                '_int_0', '_int_1', '_int_2',
                '_final') 
  
  # Store ind models in list
  gca_l2_dele_simple <- mget(c(paste0(mod_type, mod_spec)))
  
  save(gca_l2_dele_simple,
       file = here("mods", "stress", "gca", "LL_changes",
                   "gca_l2_dele_simple.Rdata"))
  
  
  
  
}

### USE

# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_l2_use_base <- mod_ot3
  # lmer(eLog ~ 1 + (ot1 + ot2) +
  #        (1 | participant) +
  #        (1 | target),
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 2e5)),
  #      data = l2_data, REML = F)
  
  # add condition effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_use_cond_0 <- update(gca_l2_use_base, . ~ . + condition_sum) 
  gca_l2_use_cond_1 <- update(gca_l2_use_cond_0, . ~ . + ot1:condition_sum) 
  gca_l2_use_cond_2 <- update(gca_l2_use_cond_1, . ~ . + ot2:condition_sum) 
  
  l2_use_cond_anova <-
    anova(gca_l2_use_base, gca_l2_use_cond_0, gca_l2_use_cond_1,
          gca_l2_use_cond_2)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_use_base      6 100766 100813 -50377   100754                       
  # gca_l2_use_cond_0    7 100763 100818 -50375   100749 4.6014  1    0.03195 *
  # gca_l2_use_cond_1    8 100760 100823 -50372   100744 5.0641  1    0.02443 *
  # gca_l2_use_cond_2    9 100762 100832 -50372   100744 0.0445  1    0.83293
  
  
  
  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_use_l1_0 <- update(gca_l2_use_cond_1, . ~ . + l1_sum) 
  gca_l2_use_l1_1 <- update(gca_l2_use_l1_0, . ~ . + ot1:l1_sum) 
  gca_l2_use_l1_2 <- update(gca_l2_use_l1_1, . ~ . + ot2:l1_sum) 
  
  l2_use_l1_anova <-
    anova(gca_l2_use_cond_1, gca_l2_use_l1_0, gca_l2_use_l1_1,
          gca_l2_use_l1_2)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_use_cond_1    8 100760 100823 -50372   100744                       
  # gca_l2_use_l1_0      9 100762 100832 -50372   100744 0.0180  1    0.89316  
  # gca_l2_use_l1_1     10 100761 100839 -50371   100741 2.8480  1    0.09149 .
  # gca_l2_use_l1_2     11 100763 100849 -50371   100741 0.3314  1    0.56486 
  
  
  
  # add use effect to intercept, linear slope, quadratic, and cubic time terms
  
  gca_l2_use_use_0 <- update(gca_l2_use_cond_1,   . ~ . + use_z) 
  gca_l2_use_use_1 <- update(gca_l2_use_use_0, . ~ . + ot1:use_z) 
  gca_l2_use_use_2 <- update(gca_l2_use_use_1, . ~ . + ot2:use_z)
  
  l2_use_use_anova <-
    anova(gca_l2_use_cond_1, gca_l2_use_use_0, gca_l2_use_use_1,
          gca_l2_use_use_2)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_use_cond_1    8 100760 100823 -50372   100744                        
  # gca_l2_use_use_0     9 100762 100832 -50372   100744 0.2890  1   0.590849   
  # gca_l2_use_use_1    10 100756 100834 -50368   100736 7.8599  1   0.005054 **
  # gca_l2_use_use_2    11 100755 100841 -50367   100733 3.0390  1   0.081288 . 
  
  # add interactions
  gca_l2_use_l1use_0 <- update(gca_l2_use_use_1,    . ~ . + l1_sum:use_z) 
  gca_l2_use_l1use_1 <- update(gca_l2_use_l1use_0, . ~ . + ot1:l1_sum:use_z) 
  gca_l2_use_l1use_2 <- update(gca_l2_use_l1use_1, . ~ . + ot2:l1_sum:use_z)
  
  anova(gca_l2_use_use_1, gca_l2_use_l1use_0, gca_l2_use_l1use_1,
        gca_l2_use_l1use_2)
  #                    npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_use_use_1     10 100756 100834 -50368   100736                     
  # gca_l2_use_l1use_0   11 100757 100843 -50368   100735 1.0165  1     0.3133
  # gca_l2_use_l1use_1   12 100759 100852 -50367   100735 0.5223  1     0.4699
  # gca_l2_use_l1use_2   13 100759 100860 -50367   100733 1.4066  1     0.2356
  
  gca_l2_use_l1cond_0 <- update(gca_l2_use_use_1,    . ~ . + l1_sum:condition_sum) 
  gca_l2_use_l1cond_1 <- update(gca_l2_use_l1cond_0, . ~ . + ot1:l1_sum:condition_sum) 
  gca_l2_use_l1cond_2 <- update(gca_l2_use_l1cond_1, . ~ . + ot2:l1_sum:condition_sum)
  
  anova(gca_l2_use_use_1, gca_l2_use_l1cond_0, gca_l2_use_l1cond_1,
        gca_l2_use_l1cond_2)
  #                     npar    AIC    BIC logLik deviance    Chisq Df Pr(>Chisq)    
  # gca_l2_use_use_1      10 100756 100834 -50368   100736                     
  # gca_l2_use_l1cond_0   11 100756 100842 -50367   100734 2.1974  1     0.1382
  # gca_l2_use_l1cond_1   12 100755 100849 -50366   100731 2.5375  1     0.1112
  # gca_l2_use_l1cond_2   13 100757 100858 -50365   100731 0.4800  1     0.4884
  
  gca_l2_use_usecond_0 <- update(gca_l2_use_use_1,     . ~ . + condition_sum:use_z) 
  gca_l2_use_usecond_1 <- update(gca_l2_use_usecond_0, . ~ . + ot1:condition_sum:use_z) 
  gca_l2_use_usecond_2 <- update(gca_l2_use_usecond_1, . ~ . + ot2:condition_sum:use_z)
  
  anova(gca_l2_use_use_1, gca_l2_use_usecond_0, gca_l2_use_usecond_1,
        gca_l2_use_usecond_2)
  #                      npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_use_use_1       10 100756 100834 -50368   100736                       
  # gca_l2_use_usecond_0   11 100752 100838 -50365   100730 5.8751  1    0.01536 *
  # gca_l2_use_usecond_1   12 100750 100843 -50363   100726 4.6936  1    0.03028 *
  # gca_l2_use_usecond_2   13 100748 100850 -50361   100722 3.0720  1    0.07965 .
  
  gca_l2_use_int_0 <- update(gca_l2_use_usecond_1, . ~ . + l1_sum:use_z:condition_sum) 
  gca_l2_use_int_1 <- update(gca_l2_use_int_0, . ~ . + ot1:l1_sum:use_z:condition_sum) 
  gca_l2_use_int_2 <- update(gca_l2_use_int_1, . ~ . + ot2:l1_sum:use_z:condition_sum)
  
  anova(gca_l2_use_usecond_1, gca_l2_use_int_0, gca_l2_use_int_1,
        gca_l2_use_int_2)
  #                      npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_l2_use_usecond_1   12 100750 100843 -50363   100726                       
  # gca_l2_use_int_0       13 100746 100847 -50360   100720 5.7869  1    0.01615 *
  # gca_l2_use_int_1       14 100746 100855 -50359   100718 2.0530  1    0.15191  
  # gca_l2_use_int_2       15 100744 100861 -50357   100714 3.3269  1    0.06815 .
  
  gca_l2_use_final <- gca_l2_use_int_0
  
  summary(gca_l2_use_final)
  confint(gca_l2_use_final)
  #                                  2.5 %      97.5 %
  # .sig01                      0.29538355  0.43078216
  # .sig02                      0.22420813  0.40977395
  # .sigma                      3.15008888  3.21629118
  # (Intercept)                 0.04219306  0.04813174
  # ot1                         1.17453889  1.18384393
  # ot2                         0.60347240  0.60962774
  # condition_sum               0.09921669  0.10381378
  # use_z                      -0.02698325 -0.02350678
  # ot1:condition_sum           0.15209783  0.15803306
  # ot1:use_z                   0.18158837  0.18754462
  # condition_sum:use_z         0.04086556  0.04288544
  # ot1:condition_sum:use_z    -0.13355761 -0.12771281
  # condition_sum:use_z:l1_sum  0.05422822  0.05626015
  confint(gca_l2_use_final, method = 'Wald') 
  #                             2.5 %     97.5 %
  # .sig01                         NA         NA
  # .sig02                         NA         NA
  # .sigma                         NA         NA
  # (Intercept)                -0.086568806  0.171457583
  # ot1                         1.046011391  1.304247757
  # ot2                         0.474328972  0.733152073
  # condition_sum              -0.001557811  0.200376611
  # use_z                      -0.102672704  0.048999090
  # ot1:condition_sum           0.023306381  0.281390199
  # ot1:use_z                   0.054508196  0.309179620
  # condition_sum:use_z        -0.003016547  0.084918114
  # ot1:condition_sum:use_z    -0.260366402 -0.006255299
  # condition_sum:use_z:l1_sum  0.010065228  0.098562611
  car::vif(gca_l2_use_final)
  # ot1                        ot2              condition_sum 
  # 1.005329                   1.002154                   1.005103 
  # use_z          ot1:condition_sum                  ot1:use_z 
  # 1.003876                   1.005591                   1.009940 
  # condition_sum:use_z    ot1:condition_sum:use_z condition_sum:use_z:l1_sum 
  # 1.105098                   1.017300                   1.095712
  performance::r2_nakagawa(gca_l2_use_final, by_group = FALSE, tolerance = 1e-05)
  # R2 for Mixed Models
  
  # Conditional R2: 0.042
  # Marginal R2: 0.021
  
  rsq::rsq(gca_l2_use_final)
  # $model
  # [1] 0.06516387
  
  # $fixed
  # [1] 0.03582148
  # 
  # $random
  # [1] 0.02934239
  rsq::rsq(gca_l2_use_final, adj = TRUE)
  # $model
  # [1] 0.06469431
  # 
  # $fixed
  # [1] 0.03533718
  # 
  # $random
  # [1] 0.02935713
  
  mod_type <- "gca_l2_use"
  mod_spec <- c('_base',
                '_cond_0', '_cond_1', '_cond_2',
                '_l1_0', '_l1_1', '_l1_2',
                '_use_0', '_use_1', '_use_2',
                '_cond_0', '_cond_1', '_cond_2',
                '_l1use_0', '_l1use_1', '_l1use_2',
                '_l1cond_0', '_l1cond_1', '_l1cond_2',
                '_usecond_0', '_usecond_1', '_usecond_2',
                '_int_0', '_int_1', '_int_2',
                '_final') 
  
  # Store ind models in list
  gca_l2_use_simple <- mget(c(paste0(mod_type, mod_spec)))
  
  save(gca_l2_use_simple,
       file = here("mods", "stress", "gca", "LL_changes",
                   "gca_l2_use_simple.Rdata"))
  
  
  
  
  
  
  
  # # Save anova model comparisons
  # nested_model_comparisons <-
  #   mget(c('l2_l1_anova', 'l2_dele_anova', 
  #          'l2_use_anova', 'l2_int_anova'
  #   ))
  # 
  # save(nested_model_comparisons,
  #      file = here("mods", "stress", "gca", "LL_changes",
  #                  "nested_model_comparisons_l2.Rdata"))
  
  
}

# -----------------------------------------------------------------------------






### L2 EXP

# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_l2_exp_base <- mod_ot3
  # lmer(eLog ~ 1 + (ot1 + ot2) +
  #        (1 | participant) +
  #        (1 | target),
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 2e5)),
  #      data = l2_data, REML = F)
  
  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_exp_l1_0 <- update(gca_l2_exp_base, . ~ . + l1_sum) 
  gca_l2_exp_l1_1 <- update(gca_l2_exp_l1_0, . ~ . + ot1:l1_sum) 
  gca_l2_exp_l1_2 <- update(gca_l2_exp_l1_1, . ~ . + ot2:l1_sum) 
  
  l2_exp_l1_anova <-
    anova(gca_l2_exp_base, gca_l2_exp_l1_0, gca_l2_exp_l1_1,
          gca_l2_exp_l1_2)
  #                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_exp_base    6 100766 100813 -50377   100754                       
  # gca_l2_exp_l1_0    7 100768 100822 -50377   100754 0.0147  1    0.90340  
  # gca_l2_exp_l1_1    8 100767 100829 -50375   100751 2.9692  1    0.08486 .
  # gca_l2_exp_l1_2    9 100769 100839 -50375   100751 0.3591  1    0.54899 
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_exp_dele_0 <- update(gca_l2_exp_base, . ~ . + DELE_z) 
  gca_l2_exp_dele_1 <- update(gca_l2_exp_dele_0, . ~ . + ot1:DELE_z) 
  gca_l2_exp_dele_2 <- update(gca_l2_exp_dele_1, . ~ . + ot2:DELE_z) 
  
  l2_exp_dele_anova <-
    anova(gca_l2_exp_base, gca_l2_exp_dele_0, gca_l2_exp_dele_1,
          gca_l2_exp_dele_2)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_use_base      6 100766 100813 -50377   100754                       
  # gca_l2_use_dele_0    7 100764 100819 -50375   100750  3.8740  1  0.0490391 *  
  # gca_l2_exp_dele_1    8 100754 100816 -50369   100738 11.9497  1  0.0005466 ***
  # gca_l2_exp_dele_2    9 100756 100826 -50369   100738  0.0902  1  0.7639364 
  
  # add use effect to intercept, linear slope, quadratic, and cubic time terms
  
  gca_l2_exp_use_0 <- update(gca_l2_exp_dele_1,   . ~ . + use_z) 
  gca_l2_exp_use_1 <- update(gca_l2_exp_use_0, . ~ . + ot1:use_z) 
  gca_l2_exp_use_2 <- update(gca_l2_exp_use_1, . ~ . + ot2:use_z)
  
  l2_exp_use_anova <-
    anova(gca_l2_exp_dele_1, gca_l2_exp_use_0, gca_l2_exp_use_1,
          gca_l2_exp_use_2)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_exp_dele_1    8 100754 100816 -50369   100738                       
  # gca_l2_exp_use_0     9 100755 100825 -50368   100737 1.2932  1    0.25546  
  # gca_l2_exp_use_1    10 100753 100831 -50366   100733 4.1106  1    0.04261 *
  # gca_l2_exp_use_2    11 100752 100838 -50365   100730 2.8931  1    0.08896 . 
  
  # add interactions
  gca_l2_exp_l1use_0 <- update(gca_l2_exp_use_1,    . ~ . + l1_sum:use_z) 
  gca_l2_exp_l1use_1 <- update(gca_l2_exp_l1use_0, . ~ . + ot1:l1_sum:use_z) 
  gca_l2_exp_l1use_2 <- update(gca_l2_exp_l1use_1, . ~ . + ot2:l1_sum:use_z)
  
  anova(gca_l2_exp_use_1, gca_l2_exp_l1use_0, gca_l2_exp_l1use_1,
        gca_l2_exp_l1use_2)
  #                    npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_exp_use_1     10 100753 100831 -50366   100733                     
  # gca_l2_exp_l1use_0   11 100753 100839 -50366   100731 1.2847  1     0.2570
  # gca_l2_exp_l1use_1   12 100755 100848 -50365   100731 0.6416  1     0.4231
  # gca_l2_exp_l1use_2   13 100756 100857 -50365   100730 1.0876  1     0.2970
  
  gca_l2_exp_l1dele_0 <- update(gca_l2_exp_use_1,    . ~ . + l1_sum:DELE_z) 
  gca_l2_exp_l1dele_1 <- update(gca_l2_exp_l1dele_0, . ~ . + ot1:l1_sum:DELE_z) 
  gca_l2_exp_l1dele_2 <- update(gca_l2_exp_l1dele_1, . ~ . + ot2:l1_sum:DELE_z)
  
  anova(gca_l2_exp_use_1, gca_l2_exp_l1dele_0, gca_l2_exp_l1dele_1,
        gca_l2_exp_l1dele_2)
  #                     npar    AIC    BIC logLik deviance    Chisq Df Pr(>Chisq)    
  # gca_l2_exp_use_1      10 100753 100831 -50366   100733                       
  # gca_l2_exp_l1dele_0   11 100752 100838 -50365   100730 2.5560  1    0.10988  
  # gca_l2_exp_l1dele_1   12 100750 100844 -50363   100726 3.9303  1    0.04742 *
  # gca_l2_exp_l1dele_2   13 100752 100853 -50363   100726 0.1396  1    0.70865
  
  gca_l2_exp_usedele_0 <- update(gca_l2_exp_l1dele_1,     . ~ . + DELE_z:use_z) 
  gca_l2_exp_usedele_1 <- update(gca_l2_exp_usedele_0, . ~ . + ot1:DELE_z:use_z) 
  gca_l2_exp_usedele_2 <- update(gca_l2_exp_usedele_1, . ~ . + ot2:DELE_z:use_z)
  
  anova(gca_l2_exp_l1dele_1, gca_l2_exp_usedele_0, gca_l2_exp_usedele_1,
        gca_l2_exp_usedele_2)
  #                      npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_exp_l1dele_1    12 100750 100844 -50363   100726                       
  # gca_l2_exp_usedele_0   13 100752 100853 -50363   100726 0.0063  1    0.93686  
  # gca_l2_exp_usedele_1   14 100749 100858 -50361   100721 5.0525  1    0.02459 *
  # gca_l2_exp_usedele_2   15 100751 100868 -50361   100721 0.0204  1    0.88635  
  
  gca_l2_exp_int_0 <- update(gca_l2_exp_usedele_1, . ~ . + l1_sum:use_z:DELE_z) 
  gca_l2_exp_int_1 <- update(gca_l2_exp_int_0, . ~ . + ot1:l1_sum:use_z:DELE_z) 
  gca_l2_exp_int_2 <- update(gca_l2_exp_int_1, . ~ . + ot2:l1_sum:use_z:DELE_z)
  
  anova(gca_l2_exp_usedele_1, gca_l2_exp_int_0, gca_l2_exp_int_1,
        gca_l2_exp_int_2)
  #                      npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_l2_exp_usedele_1   14 100749 100858 -50361   100721                     
  # gca_l2_exp_int_0       15 100751 100868 -50360   100721 0.5044  1     0.4776
  # gca_l2_exp_int_1       16 100752 100876 -50360   100720 1.0247  1     0.3114
  # gca_l2_exp_int_2       17 100753 100886 -50360   100719 0.3247  1     0.5688
  
  gca_l2_exp_final <- gca_l2_exp_usedele_1
  
  summary(gca_l2_exp_final)
  confint(gca_l2_exp_final)
  #                          2.5 %       97.5 %
  # .sig01             0.276859817  0.409793848
  # .sig02             0.239311245  0.433198655
  # .sigma             3.150748982  3.216964384
  # (Intercept)        0.047407702  0.053604714
  # ot1                1.129403208  1.138163549
  # ot2                0.601854648  0.608035500
  # DELE_z             0.079131297  0.082653240
  # use_z             -0.051481957 -0.047997464
  # ot1:DELE_z         0.198372848  0.204563718
  # ot1:use_z          0.146968525  0.153136797
  # DELE_z:l1_sum      0.063234397  0.066641585
  # DELE_z:use_z      -0.006129523 -0.002737731
  # ot1:DELE_z:l1_sum -0.143004958 -0.136923452
  # ot1:DELE_z:use_z   0.146778498  0.152693412
  confint(gca_l2_exp_final, method = 'Wald') 
  #                             2.5 %     97.5 %
  # .sig01                         NA         NA
  # .sig02                         NA         NA
  # .sigma                         NA         NA
  # (Intercept)                -0.086215208  0.18155772
  # ot1                0.997233285  1.26260155
  # ot2                0.472726182  0.73152364
  # DELE_z             0.003526831  0.15503661
  # use_z             -0.126872799  0.02420414
  # ot1:DELE_z         0.064670295  0.33260026
  # ot1:use_z          0.016142425  0.27832692
  # DELE_z:l1_sum     -0.010179316  0.13693820
  # DELE_z:use_z      -0.079413766  0.06744241
  # ot1:DELE_z:l1_sum -0.273231585 -0.01225882
  # ot1:DELE_z:use_z   0.018837960  0.27522069
  car::vif(gca_l2_exp_final)
  # ot1                        ot2          DELE_z 
  # 1.061785          1.002174          1.072901 
  # use_z        ot1:DELE_z         ot1:use_z 
  # 1.074554          1.066631          1.070900 
  # DELE_z:l1_sum      DELE_z:use_z ot1:DELE_z:l1_sum 
  # 1.011646          1.007131          1.012945 
  # ot1:DELE_z:use_z 
  # 1.066214
  performance::r2_nakagawa(gca_l2_exp_final, by_group = FALSE, tolerance = 1e-05)
  # R2 for Mixed Models
  
  # Conditional R2: 0.041
  # Marginal R2: 0.020
 
  rsq::rsq(gca_l2_exp_final)
  # $model
  # [1] 0.06562273
  
  # $fixed
  # [1] 0.03661266
  # 
  # $random
  # [1] 0.02901007
  
  rsq::rsq(gca_l2_exp_final, adj = TRUE)
  # $model
  # [1] 0.06510123
  
  # $fixed
  # [1] 0.03607496
  # 
  # $random
  # [1] 0.02902626
  
  
  mod_type <- "gca_l2_exp"
  mod_spec <- c('_base',
                '_dele_0', '_dele_1', '_dele_2',
                '_l1_0', '_l1_1', '_l1_2',
                '_use_0', '_use_1', '_use_2',
                '_l1use_0', '_l1use_1', '_l1use_2',
                '_l1dele_0', '_l1dele_1', '_l1dele_2',
                '_usedele_0', '_usedele_1', '_usedele_2',
                '_int_0', '_int_1', '_int_2',
                '_final') 
  
  # Store ind models in list
  gca_l2_exp_simple <- mget(c(paste0(mod_type, mod_spec)))
  
  save(gca_l2_exp_simple,
       file = here("mods", "stress", "gca", "LL_changes",
                   "gca_l2_exp_simple.Rdata"))
}





# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_mon <- mon_data %>%
  dplyr::select(time_zero, ot1:ot2, condition_sum) %>% 
  distinct

# Get model predictions and SE
fits_all_mon <- predictSE(gca_mon_final, new_dat_mon) %>%  
  as_tibble %>%
  bind_cols(new_dat_mon) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at first syllable offset
target_onset_preds_mon <- filter(fits_all_mon, time_zero == 4) %>%
  select(elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))


# -----------------------------------------------------------------------------


new_dele_dat <- l2_data %>%
  dplyr::select(l1_sum, condition_sum, time_zero, ot1:ot2) %>%
  distinct %>%
  # mutate(l1_sum = as.character(l1_sum)) %>% 
  expand_grid(., tibble(DELE_z = c(-1, 0, 1))) 

fits_all_dele <- predictSE(gca_l2_dele_final, new_dele_dat) %>%
  as_tibble %>%
  bind_cols(new_dele_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at first syllable offset
target_onset_preds_dele <- filter(fits_all_dele, time_zero == 4) %>% #
  select(`L1 experience` = l1_sum, `L2 proficiency` = DELE_z, `Stress pattern` = condition_sum,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))





new_use_dat <- l2_data %>%
  dplyr::select(l1_sum, condition_sum, time_zero, ot1:ot3) %>%
  distinct %>%
  # mutate(l1_sum = as.character(l1_sum)) %>% 
  expand_grid(., tibble(use_z = c(-1, 0, 1))) 

fits_all_use <- predictSE(gca_l2_use_final, new_use_dat) %>%
  as_tibble %>%
  bind_cols(new_use_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at first syllable offset
target_onset_preds_use <- filter(fits_all_use, time_zero == 4) %>% #
  select(`L1 experience` = l1_sum, `L2 use` = use_z, `Stress pattern` = condition_sum,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))




new_exp_dat <- l2_data %>%
  dplyr::select(l1_sum, time_zero, ot1:ot3) %>%
  distinct %>%
  # mutate(l1_sum = as.character(l1_sum)) %>% 
  expand_grid(., tibble(DELE_z = c(-1, 0, 1))) %>%
  expand_grid(., tibble(use_z = c(-1, 0, 1))) 

fits_all_exp <- predictSE(gca_l2_exp_final, new_exp_dat) %>%
  as_tibble %>%
  bind_cols(new_exp_dat) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)


# Filter preds at first syllable offset
target_onset_preds_exp <- filter(fits_all_exp, time_zero == 4) %>% #
  select(`L1 experience` = l1_sum, `L2 use` = use_z, `L2 proficiency` = DELE_z,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub))



model_preds_simple <- mget(c("fits_all_mon", "fits_all_dele", "fits_all_use", "fits_all_exp",
                             "target_onset_preds_mon", "target_onset_preds_dele", 
                             "target_onset_preds_use", "target_onset_preds_exp"))

save(model_preds_simple,
     file = here("mods", "stress", "gca", "LL_changes",
                 "model_preds_simple.Rdata"))









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

