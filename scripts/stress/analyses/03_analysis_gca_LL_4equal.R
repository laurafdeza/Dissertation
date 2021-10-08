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
#source(here::here("scripts", "02_load_data.R"))
stress50 <- read_csv(here("data", "clean", "stress_50ms_final_onsetc3updated.csv"))

# Get path to saved models
gca_mods_path  <- here("mods", "stress", "gca", "LL_changes")

# Load models as lists
load(paste0(gca_mods_path, "/gca_mon_mods_4equal.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_dele_4equal.Rdata")) 
load(paste0(gca_mods_path, "/gca_l2_use_4equal.Rdata"))
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
  #         npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    9 24000 24057 -11991    23982                          
  # mod_ot2   14 23981 24070 -11976    23953 29.3437  5  1.985e-05 ***
  # mod_ot3   20 23989 24116 -11974    23949  3.6882  6     0.7188   
  
  
  mod_ot0 <- update(mod_ot2, . ~ . + (1 | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + ot1 | target) +
                       + (1 + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + ot1 + ot2 | target) +
                       + ot3 + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot2, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot2    14 23981 24070 -11976    23953                          
  # mod_ot0    15 23886 23981 -11928    23856 97.0551  1  < 2.2e-16 ***
  # mod_ot1a   17 23868 23976 -11917    23834 21.8694  2  1.783e-05 ***
  # mod_ot2a   20 23869 23996 -11914    23829  5.1894  3     0.1584    
  # mod_ot3a   25 23876 24035 -11913    23826  2.4741  5     0.7804    
  
}



# Fixed effects -----------------------------------------------------------

gca_mod_mon_base <- mod_ot1a
# lmer(eLog ~ 1 + (ot1 + ot2) +
#        (1 + ot1 + ot2 | participant) +
#        (1 + ot1 | target),
#      control = lmerControl(optimizer = 'bobyqa'), 
#      REML = F,
#      data = filter(mon_data)) 

# add condition (paroxytone, oxytone) effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_mon_cond_0 <- update(gca_mod_mon_base,   . ~ . + condition_sum)
gca_mod_mon_cond_1 <- update(gca_mod_mon_cond_0,   . ~ . + ot1:condition_sum)
gca_mod_mon_cond_2 <- update(gca_mod_mon_cond_1,   . ~ . + ot2:condition_sum)

mon_cond_anova <-
  anova(gca_mod_mon_base, gca_mod_mon_cond_0, gca_mod_mon_cond_1,
        gca_mod_mon_cond_2)
#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_mon_base     17 23868 23976 -11917    23834                     
# gca_mod_mon_cond_0   18 23870 23984 -11917    23834 0.0617  1     0.8039
# gca_mod_mon_cond_1   19 23870 23991 -11916    23832 1.7096  1     0.1910
# gca_mod_mon_cond_2   20 23872 23999 -11916    23832 0.3252  1     0.5685

gca_mod_mon_final <- gca_mod_mon_base


mod_type <- "gca_mod_mon"
mod_spec <- c('_base',
              '_cond_0', '_cond_1', '_cond_2',
              '_final')


# Store ind models in list
gca_mon_mods_4equal <- mget(c(paste0(mod_type, mod_spec)))

save(gca_mon_mods_4equal,
     file = here("mods", "stress", "gca", "LL_changes",
                 "gca_mon_mods_4equal.Rdata"))



#################### L2 SPEAKERS ########################################

l2_data <- stress_gc_subset%>%
  filter(., l1 != 'es') %>% 
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1))


### PROFICIENCY

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + condition_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=2e5)),   
         data = l2_data, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + condition_sum + ot1 | participant) +
             ot2 + (1 + condition_sum + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + condition_sum + ot1 + ot2 | participant) +
             ot3 + (1 + condition_sum + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  
  #         npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    9 100923 100993 -50452   100905                          
  # mod_ot2   14 100881 100990 -50426   100853 52.2873  5   4.71e-10 **
  # mod_ot3   20 100888 101044 -50424   100848  4.0451  6     0.6706  
  
  
  
  mod_ot0 <- update(mod_ot2, . ~ . + (1 + DELE_z + l1_sum | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 + DELE_z + l1_sum | target) + 
                       + (1 + DELE_z + l1_sum + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + DELE_z + l1_sum + ot1 | target) +
                       + (1 + DELE_z + l1_sum + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + DELE_z + l1_sum + ot1 + ot2 | target) +
                       + ot3 + (1 + DELE_z + l1_sum + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot2, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot2    14 100881 100990 -50426   100853                           
  # mod_ot0    20 100740 100896 -50350   100700 152.4527  6     <2e-16 ***
  # mod_ot1a   24 100569 100756 -50260   100521 179.1601  4     <2e-16 ***
  # mod_ot2a   29 100576 100802 -50259   100518   2.8668  5     0.7205    
  # mod_ot3a   36 100586 100866 -50257   100514   4.4918  7     0.7217   
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_l2_dele_base <- #mod_ot1a
    lmer(eLog ~ 1 + (ot1 + ot2) +
           (1 + condition_sum + ot1 + ot2 | participant) +
           (1 + DELE_z + l1_sum + ot1 | target),
         control = lmerControl(optimizer = 'bobyqa',
                               optCtrl = list(maxfun = 2e5)),
         data = l2_data, REML = F)
  
  
  
  # add condition effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_dele_cond_0 <- update(gca_l2_dele_base, . ~ . + condition_sum) 
  gca_l2_dele_cond_1 <- update(gca_l2_dele_cond_0, . ~ . + ot1:condition_sum) 
  gca_l2_dele_cond_2 <- update(gca_l2_dele_cond_1, . ~ . + ot2:condition_sum) 
  
  l2_dele_cond_anova <-
    anova(gca_l2_dele_base, gca_l2_dele_cond_0, gca_l2_dele_cond_1,
          gca_l2_dele_cond_2)
  #                    npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_dele_base     24 102124 102312 -51038   102076                       
  # gca_l2_dele_cond_0   25 102122 102317 -51036   102072 4.2863  1    0.03842 *
  # gca_l2_dele_cond_1   26 102124 102327 -51036   102072 0.2875  1    0.59185  
  # gca_l2_dele_cond_2   27 102124 102334 -51035   102070 2.3444  1    0.12573 
  
  
  
  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_dele_l1_0 <- update(gca_l2_dele_cond_0, . ~ . + l1_sum) 
  gca_l2_dele_l1_1 <- update(gca_l2_dele_l1_0, . ~ . + ot1:l1_sum) 
  gca_l2_dele_l1_2 <- update(gca_l2_dele_l1_1, . ~ . + ot2:l1_sum) 
  
  l2_dele_l1_anova <-
    anova(gca_l2_dele_cond_0, gca_l2_dele_l1_0, gca_l2_dele_l1_1,
          gca_l2_dele_l1_2)
  #                    npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_dele_cond_0   25 102122 102317 -51036   102072                     
  # gca_l2_dele_l1_0     26 102124 102326 -51036   102072 0.3476  1     0.5554
  # gca_l2_dele_l1_1     27 102125 102335 -51035   102071 1.2222  1     0.2689
  # gca_l2_dele_l1_2     28 102125 102344 -51035   102069 1.1610  1     0.2813
  
  
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  
  gca_l2_dele_dele_0 <- update(gca_l2_dele_cond_0,   . ~ . + DELE_z) 
  gca_l2_dele_dele_1 <- update(gca_l2_dele_dele_0, . ~ . + ot1:DELE_z) 
  gca_l2_dele_dele_2 <- update(gca_l2_dele_dele_1, . ~ . + ot2:DELE_z)
  
  l2_dele_dele_anova <-
    anova(gca_l2_dele_cond_0, gca_l2_dele_dele_0, gca_l2_dele_dele_1,
          gca_l2_dele_dele_2)
  #                     npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_dele_cond_0   25 102122 102317 -51036   102072                       
  # gca_l2_dele_dele_0   26 102120 102323 -51034   102068 4.1316  1    0.04209 *
  # gca_l2_dele_dele_1   27 102118 102329 -51032   102064 3.7016  1    0.05436 .
  # gca_l2_dele_dele_2   28 102120 102338 -51032   102064 0.7679  1    0.38086 
  
  
  
  # Add interactions one by one
  gca_l2_dele_l1dele_0 <- update(gca_l2_dele_dele_0,     . ~ . + l1_sum:DELE_z) 
  gca_l2_dele_l1dele_1 <- update(gca_l2_dele_l1dele_0, . ~ . + ot1:l1_sum:DELE_z) 
  gca_l2_dele_l1dele_2 <- update(gca_l2_dele_l1dele_1, . ~ . + ot2:l1_sum:DELE_z)
  
  anova(gca_l2_dele_dele_0, gca_l2_dele_l1dele_0, gca_l2_dele_l1dele_1,
        gca_l2_dele_l1dele_2)
  #                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_dele_dele_0     26 102120 102323 -51034   102068                       
  # gca_l2_dele_l1dele_0   27 102119 102329 -51032   102065 3.5363  1    0.06004 .
  # gca_l2_dele_l1dele_1   28 102120 102338 -51032   102064 1.0058  1    0.31590  
  # gca_l2_dele_l1dele_2   29 102121 102348 -51032   102063 0.0297  1    0.86326  
  
  
  gca_l2_dele_l1cond_0 <- update(gca_l2_dele_dele_0,     . ~ . + l1_sum:condition_sum) 
  gca_l2_dele_l1cond_1 <- update(gca_l2_dele_l1cond_0, . ~ . + ot1:l1_sum:condition_sum) 
  gca_l2_dele_l1cond_2 <- update(gca_l2_dele_l1cond_1, . ~ . + ot2:l1_sum:condition_sum)
  
  anova(gca_l2_dele_dele_0, gca_l2_dele_l1cond_0, gca_l2_dele_l1cond_1,
        gca_l2_dele_l1cond_2)
  #                      npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_l2_dele_dele_0     26 102120 102323 -51034   102068                     
  # gca_l2_dele_l1cond_0   27 102122 102332 -51034   102068 0.1735  1     0.6770
  # gca_l2_dele_l1cond_1   28 102124 102342 -51034   102068 0.0205  1     0.8862
  # gca_l2_dele_l1cond_2   29 102126 102352 -51034   102068 0.0121  1     0.9126
  
  
  
  gca_l2_dele_delecond_0 <- update(gca_l2_dele_dele_0,     . ~ . + DELE_z:condition_sum) 
  gca_l2_dele_delecond_1 <- update(gca_l2_dele_delecond_0, . ~ . + ot1:DELE_z:condition_sum) 
  gca_l2_dele_delecond_2 <- update(gca_l2_dele_delecond_1, . ~ . + ot2:DELE_z:condition_sum)
  
  anova(gca_l2_dele_dele_0, gca_l2_dele_delecond_0, gca_l2_dele_delecond_1,
        gca_l2_dele_delecond_2)
  #                        npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_dele_dele_0       26 102120 102323 -51034   102068                       
  # gca_l2_dele_delecond_0   27 102122 102332 -51034   102068 0.3131  1    0.57578  
  # gca_l2_dele_delecond_1   28 102120 102338 -51032   102064 4.0628  1    0.04384 *
  # gca_l2_dele_delecond_2   29 102122 102348 -51032   102064 0.0145  1    0.90415 
  
  
  
  gca_l2_dele_int_0 <- update(gca_l2_dele_delecond_1, . ~ . + l1_sum:DELE_z:condition_sum) 
  gca_l2_dele_int_1 <- update(gca_l2_dele_int_0, . ~ . + ot1:l1_sum:DELE_z:condition_sum) 
  gca_l2_dele_int_2 <- update(gca_l2_dele_int_1, . ~ . + ot2:l1_sum:DELE_z:condition_sum)
  
  anova(gca_l2_dele_delecond_1, gca_l2_dele_int_0, gca_l2_dele_int_1,
        gca_l2_dele_int_2)
  #                        npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_dele_delecond_1   28 102120 102338 -51032   102064                       
  # gca_l2_dele_int_0        29 102120 102346 -51031   102062 1.2587  1    0.26190  
  # gca_l2_dele_int_1        30 102119 102353 -51030   102059 3.2390  1    0.07191 .
  # gca_l2_dele_int_2        31 102121 102362 -51029   102059 0.4916  1    0.48321 
  
  
  gca_l2_dele_final <- gca_l2_dele_delecond_1
  
  
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
  gca_l2_dele_4equal <- mget(c(paste0(mod_type, mod_spec)))
  
  save(gca_l2_dele_4equal,
       file = here("mods", "stress", "gca", "LL_changes",
                   "gca_l2_dele_4equal.Rdata"))
  
  
  
  
}

### USE

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  
  
  mod_ot0 <- update(mod_ot2, . ~ . + (1 + use_z + l1_sum | target))
  
  mod_ot1a <- update(mod_ot0, . ~ . -(1 + use_z + l1_sum | target) + 
                       + (1 + use_z + l1_sum + ot1 | target))
  
  mod_ot2a <- update(mod_ot1a, . ~ . -(1 + use_z + l1_sum + ot1 | target) +
                       + (1 + use_z + l1_sum + ot1 + ot2 | target))
  
  mod_ot3a <- update(mod_ot2a, . ~ . -(1 + use_z + l1_sum + ot1 + ot2 | target) +
                       + ot3 + (1 + use_z + l1_sum + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot2, mod_ot0, mod_ot1a, mod_ot2a, mod_ot3a)
  #          npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot2    14 100881 100990 -50426   100853                           
  # mod_ot0    20 100698 100853 -50329   100658 194.96  6     <2e-16 ***
  # mod_ot1a   24 100561 100748 -50257   100513 144.45  4     <2e-16 ***
  # mod_ot2a   29 100708 100934 -50325   100650   0.00  5          1    
  # mod_ot3a   36 100578 100858 -50253   100506 143.64  7     <2e-16 ***
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_l2_use_base <- mod_ot3a
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
  #        (1 + condition_sum + ot1 + ot2 | participant) +
  #        (1 + use_z + l1_sum + ot1 + ot2 + ot3 | target),
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 2e5)),
  #      data = l2_data, REML = F)
  
  
  
  # add condition effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_use_cond_0 <- update(gca_l2_use_base, . ~ . + condition_sum) 
  gca_l2_use_cond_1 <- update(gca_l2_use_cond_0, . ~ . + ot1:condition_sum) 
  gca_l2_use_cond_2 <- update(gca_l2_use_cond_1, . ~ . + ot2:condition_sum) 
  gca_l2_use_cond_3 <- update(gca_l2_use_cond_2, . ~ . + ot3:condition_sum) 
  
  l2_use_cond_anova <-
    anova(gca_l2_use_base, gca_l2_use_cond_0, gca_l2_use_cond_1,
          gca_l2_use_cond_2,  gca_l2_use_cond_3)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_use_base     36 100578 100858 -50253   100506                       
  # gca_l2_use_cond_0   37 100575 100863 -50250   100501 5.3797  1    0.02037 *
  # gca_l2_use_cond_1   38 100575 100871 -50249   100499 1.9063  1    0.16738  
  # gca_l2_use_cond_2   39 100577 100881 -50249   100499 0.0326  1    0.85669
  # gca_l2_use_cond_3   40 100579 100890 -50249   100499 0.0364  1    0.84871  
  
  
  
  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_use_l1_0 <- update(gca_l2_use_cond_0, . ~ . + l1_sum) 
  gca_l2_use_l1_1 <- update(gca_l2_use_l1_0, . ~ . + ot1:l1_sum) 
  gca_l2_use_l1_2 <- update(gca_l2_use_l1_1, . ~ . + ot2:l1_sum) 
  gca_l2_use_l1_3 <- update(gca_l2_use_l1_2, . ~ . + ot3:l1_sum) 
  
  l2_use_l1_anova <-
    anova(gca_l2_use_cond_0, gca_l2_use_l1_0, gca_l2_use_l1_1,
          gca_l2_use_l1_2, gca_l2_use_l1_3)
  #                    npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # gca_l2_use_cond_0   37 100575 100863 -50250   100501                     
  # gca_l2_use_l1_0     38 100577 100873 -50250   100501 0.0035  1     0.9527
  # gca_l2_use_l1_1     39 100576 100880 -50249   100498 2.3084  1     0.1287
  # gca_l2_use_l1_2     40 100577 100889 -50249   100497 0.9289  1     0.3352
  # gca_l2_use_l1_3     41 100578 100897 -50248   100496 1.5500  1     0.2131
  
  
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  
  gca_l2_use_use_0 <- update(gca_l2_use_cond_0,   . ~ . + use_z) 
  gca_l2_use_use_1 <- update(gca_l2_use_use_0, . ~ . + ot1:use_z) 
  gca_l2_use_use_2 <- update(gca_l2_use_use_1, . ~ . + ot2:use_z)
  gca_l2_use_use_3 <- update(gca_l2_use_use_2, . ~ . + ot3:use_z)
  
  l2_use_use_anova <-
    anova(gca_l2_use_cond_0, gca_l2_use_use_0, gca_l2_use_use_1,
          gca_l2_use_use_2, gca_l2_use_use_3)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_use_cond_0   37 100575 100863 -50250   100501                        
  # gca_l2_use_use_0    38 100576 100872 -50250   100500 0.3630  1   0.546833   
  # gca_l2_use_use_1    39 100571 100875 -50247   100493 6.7517  1   0.009366 **
  # gca_l2_use_use_2    40 100571 100882 -50245   100491 2.7243  1   0.098834 . 
  # gca_l2_use_use_3    41 100570 100890 -50244   100488 2.5690  1   0.108976 
  
  
  
  
  gca_l2_use_l1use_0 <- update(gca_l2_use_use_1,    . ~ . + l1_sum:use_z) 
  gca_l2_use_l1use_1 <- update(gca_l2_use_l1use_0, . ~ . + ot1:l1_sum:use_z) 
  gca_l2_use_l1use_2 <- update(gca_l2_use_l1use_1, . ~ . + ot2:l1_sum:use_z)
  gca_l2_use_l1use_3 <- update(gca_l2_use_l1use_2, . ~ . + ot3:l1_sum:use_z)
  
  anova(gca_l2_use_use_1, gca_l2_use_l1use_0, gca_l2_use_l1use_1,
        gca_l2_use_l1use_2, gca_l2_use_l1use_3)
  #                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_use_use_1     39 100571 100875 -50247   100493                       
  # gca_l2_use_l1use_0   40 100572 100884 -50246   100492 1.2692  1    0.25991  
  # gca_l2_use_l1use_1   41 100572 100891 -50245   100490 2.5474  1    0.11047  
  # gca_l2_use_l1use_2   42 100576 100903 -50246   100492 0.0000  1    1.00000  
  # gca_l2_use_l1use_3   43 100574 100909 -50244   100488 3.8223  1    0.05058 .
  
  
  
  
  gca_l2_use_l1cond_0 <- update(gca_l2_use_use_1,    . ~ . + l1_sum:condition_sum) 
  gca_l2_use_l1cond_1 <- update(gca_l2_use_l1cond_0, . ~ . + ot1:l1_sum:condition_sum) 
  gca_l2_use_l1cond_2 <- update(gca_l2_use_l1cond_1, . ~ . + ot2:l1_sum:condition_sum)
  gca_l2_use_l1cond_3 <- update(gca_l2_use_l1cond_2, . ~ . + ot3:l1_sum:condition_sum)
  
  anova(gca_l2_use_use_1, gca_l2_use_l1cond_0, gca_l2_use_l1cond_1,
        gca_l2_use_l1cond_2, gca_l2_use_l1cond_3)
  #                     npar    AIC    BIC logLik deviance    Chisq Df Pr(>Chisq)    
  # gca_l2_use_use_1      39 100571 100875 -50247   100493                           
  # gca_l2_use_l1cond_0   40 100573 100885 -50247   100493   0.1397  1     0.7085    
  # gca_l2_use_l1cond_1   41 100709 101029 -50314   100627   0.0000  1     1.0000    
  # gca_l2_use_l1cond_2   42 100574 100901 -50245   100490 137.3043  1     <2e-16 ***
  # gca_l2_use_l1cond_3   43 100576 100911 -50245   100490   0.1385  1     0.7098  
  
  
  
  
  gca_l2_use_usecond_0 <- update(gca_l2_use_l1cond_2,     . ~ . + condition_sum:use_z) 
  gca_l2_use_usecond_1 <- update(gca_l2_use_usecond_0, . ~ . + ot1:condition_sum:use_z) 
  gca_l2_use_usecond_2 <- update(gca_l2_use_usecond_1, . ~ . + ot2:condition_sum:use_z)
  gca_l2_use_usecond_3 <- update(gca_l2_use_usecond_2, . ~ . + ot3:condition_sum:use_z)
  
  anova(gca_l2_use_l1cond_2, gca_l2_use_usecond_0, gca_l2_use_usecond_1,
        gca_l2_use_usecond_2, gca_l2_use_usecond_3)
  #                      npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_use_l1cond_2    42 100574 100901 -50245   100490                       
  # gca_l2_use_usecond_0   43 100575 100910 -50245   100489 0.7256  1     0.3943  
  # gca_l2_use_usecond_1   44 100572 100915 -50242   100484 4.9363  1     0.0263 *
  # gca_l2_use_usecond_2   45 100572 100923 -50241   100482 1.9808  1     0.1593  
  # gca_l2_use_usecond_3   46 100574 100933 -50241   100482 0.1021  1     0.7493 
  
  
  
  gca_l2_use_int_0 <- update(gca_l2_use_usecond_1, . ~ . + l1_sum:use_z:condition_sum) 
  gca_l2_use_int_1 <- update(gca_l2_use_int_0, . ~ . + ot1:l1_sum:use_z:condition_sum) 
  gca_l2_use_int_2 <- update(gca_l2_use_int_1, . ~ . + ot2:l1_sum:use_z:condition_sum)
  gca_l2_use_int_3 <- update(gca_l2_use_int_2, . ~ . + ot3:l1_sum:use_z:condition_sum)
  
  
  anova(gca_l2_use_usecond_1, gca_l2_use_int_0, gca_l2_use_int_1,
        gca_l2_use_int_2, gca_l2_use_int_3)
  #                      npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_l2_use_usecond_1   44 100572 100915 -50242   100484                     
  # gca_l2_use_int_0       45 100573 100924 -50241   100483 1.3538  1     0.2446
  # gca_l2_use_int_1       46 100574 100932 -50241   100482 1.0561  1     0.3041
  # gca_l2_use_int_2       47 100574 100941 -50240   100480 1.7048  1     0.1917
  # gca_l2_use_int_3       48 100574 100948 -50239   100478 2.4597  1     0.1168
  
  gca_l2_use_final <- gca_l2_use_usecond_1
  
  
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
  gca_l2_use_4equal <- mget(c(paste0(mod_type, mod_spec)))
  
  save(gca_l2_use_4equal,
       file = here("mods", "stress", "gca", "LL_changes",
                   "gca_l2_use_4equal.Rdata"))
  
  
  
  
  
  
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







# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
new_dat_mon <- mon_data %>%
  dplyr::select(time_zero, ot1:ot2, condition_sum) %>% 
  distinct

# Get model predictions and SE
fits_all_mon <- predictSE(gca_mod_mon_final, new_dat_mon) %>%  
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




model_preds_4equal <- mget(c("fits_all_mon", "fits_all_dele", "fits_all_use",
                             "target_onset_preds_mon", "target_onset_preds_dele", "target_onset_preds_use"))

save(model_preds_4equal,
     file = here("mods", "stress", "gca", "LL_changes",
                 "model_preds_4equal.Rdata"))









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

