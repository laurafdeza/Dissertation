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
load(paste0(gca_mods_path, "/gca_mon_mods.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_l1ot2.Rdata")) #gca_l2_mods only dele parameter significant
#load(paste0(gca_mods_path, "/nested_model_l2.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata")) #mon
load(paste0(gca_mods_path, "/model_preds_l1ot2.Rdata")) #L2



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


stress_gc_subset <- stress50 %>%
  # select(., -WM_set) %>%
  filter(., time_zero >= -8 & time_zero <= 4) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma")#,
         #   condition_sum = if_else(cond == "1", -1, 1)
         ) %>%       # 1 = present (now -1), 2 = past (now 1)
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")



# -----------------------------------------------------------------------------


 

#################### L2 SPEAKERS ########################################

l2_data <- stress_gc_subset%>%
  filter(., l1 != 'es') %>% 
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1))

# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){

  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=2e5)),   
         data = l2_data, weights = 1/wts, REML = F)
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + ot1 | participant) +
             ot2 + (1 + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + ot1 + ot2 | participant) +
             ot3 + (1 + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  
  #         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod_ot1    6  46473 146522 -73230   146461                          
  # mod_ot2   10 146355 146436 -73167   146335 125.883  4  < 2.2e-16 ***
  # mod_ot3   15 146332 146455 -73151   146302  32.695  5  4.326e-06 ***
  

          mod_ot4 <- update(mod_ot3, . ~ . + (1 + DELE_z + use_z + l1_sum | target))
          
          mod_ot5 <- update(mod_ot4, . ~ . -(1 + DELE_z + use_z + l1_sum | target) + 
                               (1 + DELE_z + use_z + l1_sum + ot1 | target))
          
          mod_ot6 <- update(mod_ot5, . ~ . -(1 + DELE_z + use_z + l1_sum + ot1 | target) +
                               + (1 + DELE_z + use_z + l1_sum + ot1 + ot2 | target))
          
          mod_ot7 <- update(mod_ot6, . ~ . -(1 + DELE_z + use_z + l1_sum + ot1 + ot2 | target) +
                               + (1 + DELE_z + use_z + l1_sum + ot1 + ot2 + ot3 | target)) 
                    
          
          anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
          #         npar    AIC    BIC logLik deviance    Chisq Df Pr(>Chisq)    
          # mod_ot3   15 146332 146455 -73151   146302                           
          # mod_ot4   25 146121 146325 -73035   146071 231.4748 10  < 2.2e-16 ***
          # mod_ot5   30 146079 146324 -73010   146019  51.1566  5  8.033e-10 ***
          # mod_ot6   36 146074 146368 -73001   146002  17.3872  6   0.007961 ** 
          # mod_ot7   43 146084 146435 -72999   145998   4.2813  7   0.746864  
  
  
  
}

# -----------------------------------------------------------------------------





# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_l2_400_base <- mod_ot6
   #lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +         
   #       (1 + ot1 + ot2 + ot3 | participant) +
   #       (1 + DELE_z + use_z + l1_sum + ot1 + ot2 | target), 
   #     control = lmerControl(optimizer = 'bobyqa',
   #                           optCtrl = list(maxfun = 2e5)),
   #     data = l2_data, REML = F)

# add group effect to intercept, linear slope, quadratic, and cubic time terms
gca_l2_400_l1_0 <- update(gca_l2_400_base, . ~ . + l1_sum) 
gca_l2_400_l1_1 <- update(gca_l2_400_l1_0, . ~ . + ot1:l1_sum) 
gca_l2_400_l1_2 <- update(gca_l2_400_l1_1, . ~ . + ot2:l1_sum) 
gca_l2_400_l1_3 <- update(gca_l2_400_l1_2, . ~ . + ot3:l1_sum) 

l2_l1_anova <-
  anova(gca_l2_400_base, gca_l2_400_l1_0, gca_l2_400_l1_1,
        gca_l2_400_l1_2, gca_l2_400_l1_3)
#                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
# gca_l2_400_base   36 146074 146368 -73001   146002                     
# gca_l2_400_l1_0   37 146076 146378 -73001   146002 0.1269  1     0.7217
# gca_l2_400_l1_1   38 146078 146388 -73001   146002 0.4609  1     0.4972
# gca_l2_400_l1_2   39 146077 146396 -73000   145999 2.2917  1     0.1301
# gca_l2_400_l1_3   40 146078 146405 -72999   145998 0.8264  1     0.3633

mod_type <- "gca_l2_400"
mod_spec <- c('_base', 
              "_l1_0", "_l1_1", "_l1_2", "_l1_3")

# , 
#               "_dele_0", "_dele_1", "_dele_2", 
#               "_use_0", "_use_1", "_use_2", 
#               "_int_0", "_int_1", "_int_2",
#               '_final') 

# Store ind models in list
gca_l2_400msbeforeonset <- mget(c(paste0(mod_type, mod_spec)))

save(gca_l2_400msbeforeonset,
     file = here("mods", "stress", "gca", "LL_changes",
                 "gca_l2_400msbeforeonset.Rdata"))



# add proficiency effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_400_dele_0 <- update(gca_l2_400_base,   . ~ . + DELE_z) 
gca_l2_400_dele_1 <- update(gca_l2_400_dele_0, . ~ . + ot1:DELE_z) 
gca_l2_400_dele_2 <- update(gca_l2_400_dele_1, . ~ . + ot2:DELE_z)

l2_dele_anova <-
  anova(gca_l2_400_base, gca_l2_400_dele_0, gca_l2_400_dele_1,
        gca_l2_400_dele_2)
#                     npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
# gca_l2_400_base     25 100704 100899 -50327   100654                          
# gca_l2_400_dele_0   26 100712 100915 -50330   100660  0.0000  1  1.0000000    
# gca_l2_400_dele_1   27 100699 100910 -50323   100645 15.0735  1  0.0001034 ***
# gca_l2_400_dele_2   28 100701 100919 -50322   100645  0.7087  1  0.3998799 



# add L2 use effect to intercept, linear slope, quadratic, and cubic time terms

gca_l2_400_use_0 <- update(gca_l2_400_dele_1,  . ~ . + use_z) 
gca_l2_400_use_1 <- update(gca_l2_400_use_0, . ~ . + ot1:use_z) 
gca_l2_400_use_2 <- update(gca_l2_400_use_1, . ~ . + ot2:use_z)

l2_use_anova <-
  anova(gca_l2_400_dele_1, gca_l2_400_use_0, gca_l2_400_use_1,
        gca_l2_400_use_2)
#                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
# gca_l2_mod_dele_1   27 100699 100910 -50323   100645                       
# gca_l2_mod_use_0    28 100700 100918 -50322   100644 1.3387  1    0.24726  
# gca_l2_mod_use_1    29 100700 100926 -50321   100642 2.2794  1    0.13111  
# gca_l2_mod_use_2    30 100698 100932 -50319   100638 3.9104  1    0.04799 *


        # Add interactions one by one
        gca_l2_mod_i_0 <- update(gca_l2_mod_use_2,     . ~ . + l1_sum:DELE_z) 
        gca_l2_mod_i_1 <- update(gca_l2_mod_i_0, . ~ . + ot1:l1_sum:DELE_z) 
        gca_l2_mod_i_2 <- update(gca_l2_mod_i_1, . ~ . + ot2:l1_sum:DELE_z)
        
        anova(gca_l2_mod_use_2, gca_l2_mod_i_0, gca_l2_mod_i_1,
              gca_l2_mod_i_2)
        #                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
        # gca_l2_mod_use_2   30 100698 100932 -50319   100638                        
        # gca_l2_mod_i_0     31 100707 100949 -50322   100645 0.0000  1   1.000000   
        # gca_l2_mod_i_1     32 100700 100949 -50318   100636 8.9327  1   0.002801 **
        # gca_l2_mod_i_2     33 100701 100958 -50317   100635 1.1081  1   0.292489   
        
        
        gca_l2_mod_in_0 <- update(gca_l2_mod_i_1,     . ~ . + l1_sum:use_z) 
        gca_l2_mod_in_1 <- update(gca_l2_mod_in_0, . ~ . + ot1:l1_sum:use_z) 
        gca_l2_mod_in_2 <- update(gca_l2_mod_in_1, . ~ . + ot2:l1_sum:use_z)
        
        anova(gca_l2_mod_i_1, gca_l2_mod_in_0, gca_l2_mod_in_1,
              gca_l2_mod_in_2)
        #                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
        # gca_l2_mod_i_1    32 100700 100949 -50318   100636                       
        # gca_l2_mod_in_0   33 100702 100959 -50318   100636 0.2277  1     0.6333  
        # gca_l2_mod_in_1   34 100699 100964 -50316   100631 4.7269  1     0.0297 *
        # gca_l2_mod_in_2   35 100701 100973 -50315   100631 0.5246  1     0.4689  
        
        
        
        gca_l2_mod_it_0 <- update(gca_l2_mod_in_1,     . ~ . + DELE_z:use_z) 
        gca_l2_mod_it_1 <- update(gca_l2_mod_it_0, . ~ . + ot1:DELE_z:use_z) 
        gca_l2_mod_it_2 <- update(gca_l2_mod_it_1, . ~ . + ot2:DELE_z:use_z)
        
        anova(gca_l2_mod_in_1, gca_l2_mod_it_0, gca_l2_mod_it_1,
              gca_l2_mod_it_2)
        #                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
        # gca_l2_mod_in_1   34 100699 100964 -50316   100631                        
        # gca_l2_mod_it_0   35 100710 100983 -50320   100640 0.0000  1   1.000000   
        # gca_l2_mod_it_1   36 100710 100990 -50319   100638 2.6782  1   0.101733   
        # gca_l2_mod_it_2   37 100702 100991 -50314   100628 9.5322  1   0.002019 **


        
        gca_l2_mod_nt_0 <- update(gca_l2_mod_it_2,     . ~ . + l1_sum:DELE_z:use_z) 
        gca_l2_mod_nt_1 <- update(gca_l2_mod_nt_0, . ~ . + ot1:l1_sum:DELE_z:use_z) 
        gca_l2_mod_nt_2 <- update(gca_l2_mod_nt_1, . ~ . + ot2:l1_sum:DELE_z:use_z)
        
        anova(gca_l2_mod_it_2, gca_l2_mod_nt_0, gca_l2_mod_nt_1,
              gca_l2_mod_nt_2)
        #                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
        # gca_l2_mod_it_2   37 100702 100991 -50314   100628                       
        # gca_l2_mod_nt_0   38 100704 101000 -50314   100628 0.4083  1    0.52284  
        # gca_l2_mod_nt_1   39 100702 101006 -50312   100624 3.3326  1    0.06792 .
        # gca_l2_mod_nt_2   40 100704 101016 -50312   100624 0.0496  1    0.82371 
        
        gca_l2_mod_full <- gca_l2_mod_it_2
        
        summary(gca_l2_mod_full)
        #                   Estimate Std. Error t value
        # (Intercept)        0.06278    0.07005   0.896
        # ot1                1.31089    0.14902   8.797
        # ot2                0.67789    0.07742   8.756
        # DELE_z             0.10075    0.06055   1.664
        # use_z             -0.07046    0.07013  -1.005
        # ot1:DELE_z         0.32036    0.11890   2.694
        # ot1:use_z          0.25774    0.13327   1.934
        # ot2:use_z          0.18267    0.10097   1.809
        # DELE_z:l1_sum      0.05794    0.04937   1.174
        # use_z:l1_sum       0.02026    0.05562   0.364
        # DELE_z:use_z      -0.01097    0.06955  -0.158
        # ot1:DELE_z:l1_sum -0.15724    0.10594  -1.484
        # ot1:use_z:l1_sum   0.19421    0.11956   1.624
        # ot1:DELE_z:use_z   0.25436    0.14943   1.702
        # ot2:DELE_z:use_z   0.02563    0.10697   0.240
        
        
# add 3-way interaction directly

gca_l2_mod_int_0 <- update(gca_l2_mod_use_2,     . ~ . + l1_sum:DELE_z:use_z) 
gca_l2_mod_int_1 <- update(gca_l2_mod_int_0, . ~ . + ot1:l1_sum:DELE_z:use_z) 
gca_l2_mod_int_2 <- update(gca_l2_mod_int_1, . ~ . + ot2:l1_sum:DELE_z:use_z)

l2_int_anova <-
  anova(gca_l2_mod_use_2, gca_l2_mod_int_0, gca_l2_mod_int_1,
        gca_l2_mod_int_2)
#                      npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)   
# gca_l2_mod_use_2   30 100698 100932 -50319   100638                     
# gca_l2_mod_int_0   31 100699 100941 -50319   100637 0.2564  1     0.6126
# gca_l2_mod_int_1   32 100700 100949 -50318   100636 1.5717  1     0.2100
# gca_l2_mod_int_2   33 100711 100968 -50322   100645 0.0000  1     1.0000

gca_l2_mod_final <- gca_l2_mod_use_2

summary(gca_l2_mod_final)
#             Estimate Std. Error t value
# (Intercept)  0.06372    0.06967   0.915
# ot1          1.35457    0.14908   9.086
# ot2          0.67896    0.07619   8.912
# DELE_z       0.10416    0.05693   1.830
# use_z       -0.06765    0.06533  -1.036
# ot1:DELE_z   0.23245    0.11182   2.079
# ot1:use_z    0.21189    0.12212   1.735
# ot2:use_z    0.16292    0.09139   1.783




      # models with interactions one by one
      mod_spec <- c("_i_0", "_i_1", "_i_2", # Dele x l1
                    "_in_0", "_in_1", "_in_2", #use x l1
                    "_it_0", "_it_1", "_it_2", # dele x use
                    "_nt_0", "_nt_1", "_nt_2", # all three
                    '_full') 
      
      # Store ind models in list
      gca_l2_mods_ints <- mget(c(paste0(mod_type, mod_spec)))
      
      save(gca_l2_mods,
           file = here("mods", "stress", "gca", "LL_changes",
                       "gca_l2_mods_ints.Rdata"))


# Save anova model comparisons
nested_model_comparisons <-
  mget(c('l2_l1_anova', 'l2_dele_anova', 
         'l2_use_anova', 'l2_int_anova'
  ))

save(nested_model_comparisons,
     file = here("mods", "stress", "gca", "LL_changes",
                 "nested_model_comparisons_l2.Rdata"))


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





