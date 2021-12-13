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
load(paste0(gca_mods_path, "/mon_mods_onlypred.Rdata")) # gca_mon_ospan_int_1, gca_mon_corsirt_int_2

load(paste0(gca_mods_path, "/model_preds_onlypred.Rdata"))
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
vision <- read_csv("./data/clean/vision_scores_nooutliers.csv") # pred car
corsi <- read_csv("./data/clean/corsi_z_scores.csv")
wm_score <- read_csv("./data/clean/ospan_set_z_scores.csv")

vision50 <- left_join(x = stress50, y = wm, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = vision, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = corsi, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = wm_score, by = "participant", all.x=TRUE)

mon_vision <- filter(vision50, l1 == 'es') %>% select(-DELE, -percent_l2_week, 
                                                      -prof, -group)


mon_vision <- na.omit(mon_vision)




mon_vision <- mon_vision %>%
  filter(., time_zero >= -4 & time_zero <= 4) %>%
  mutate(., #l1 = fct_relevel(l1, "es", "en", "ma"),
            stress_sum = if_else(cond == "1", -1, 1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 2, prefix = "ot")




# -----------------------------------------------------------------------------







# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){
  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 + ot2 +
           (1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = mon_vision, weights = 1/wts, REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . + (1 | target))
  
  anova(mod_ot1, mod_ot2)
  #         npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # mod_ot1    5 14454 14483 -7221.8    14444                         
  # mod_ot2    6 14363 14398 -7175.3    14351 92.906  1  < 2.2e-16 ***
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_mon_base <- mod_ot2
  # lmer(eLog ~ 1 + ot1 + ot2 +    
  #        (1 | participant) +
  #        (1 | target),
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
# gca_mon_base        6 14363 14398 -7175.3    14351                     
# gca_mon_stress_0    7 14364 14405 -7175.1    14350 0.4526  1     0.5011
# gca_mon_stress_1    8 14365 14412 -7174.5    14349 1.2131  1     0.2707
# gca_mon_stress_2    9 14367 14420 -7174.5    14349 0.0239  1     0.8772


# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_0 <- update(gca_mon_base,    . ~ . + car_dev)
gca_mon_car_1 <- update(gca_mon_car_0, . ~ . + ot1:car_dev)
gca_mon_car_2 <- update(gca_mon_car_1, . ~ . + ot2:car_dev)

mon_car_anova <-
  anova(gca_mon_base, gca_mon_car_0, gca_mon_car_1,
        gca_mon_car_2)
#               npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base     6 14363 14398 -7175.3    14351                       
# gca_mon_car_0    7 14364 14404 -7174.8    14350 1.0896  1    0.29656  
# gca_mon_car_1    8 14361 14408 -7172.5    14345 4.5690  1    0.03256 *
# gca_mon_car_2    9 14363 14415 -7172.3    14345 0.2672  1    0.60524  

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_int_0 <- update(gca_mon_car_1,    . ~ . + car_dev:stress_sum)
gca_mon_car_int_1 <- update(gca_mon_car_int_0, . ~ . + ot1:car_dev:stress_sum)
gca_mon_car_int_2 <- update(gca_mon_car_int_1, . ~ . + ot2:car_dev:stress_sum)

mon_car_int_anova <-
  anova(gca_mon_car_1, gca_mon_car_int_0, gca_mon_car_int_1,
        gca_mon_car_int_2)
#                   npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_car_1        8 14361 14408 -7172.5    14345                     
# gca_mon_car_int_0    9 14363 14416 -7172.5    14345 0.0002  1     0.9898
# gca_mon_car_int_1   10 14365 14424 -7172.4    14345 0.0588  1     0.8084
# gca_mon_car_int_2   11 14367 14431 -7172.4    14345 0.0665  1     0.7965



car::vif(gca_mon_car_1)
# ot1          ot2         car_dev      ot1:car_dev 
# 1.024691    1.013401    1.003478    1.014514 





}

# -----------------------------------------------------------------------------


# save models  
mod_type <- "gca_mon"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", 
              "_car_0", "_car_1", "_car_2",
              "_car_int_0", "_car_int_1", "_car_int_2"
              )

# Store ind models in list
mon_mods_onlypred <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods_onlypred,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "mon_mods_onlypred.Rdata"))

  
#########    CORRELATION
library(nortest)
library("report") 

# ASSUMPTION 1. Normality
ggplot(mon_vision, aes(x=eLog)) + 
  geom_density()
# if error for graphics, run dev.off() and then run code again

ggplot(mon_vision, aes(x=car_dev)) + 
  geom_density()


ad.test(mon_vision$eLog) # Anderson-Darling normality test for large samples
# A = 430.26, p-value < 2.2e-16

ad.test(mon_vision$car_dev) # neither normal > Spearman
# A = 133.58, p-value < 2.2e-16


# Spearman correlation between 2 variables
cor(mon_vision$eLog, mon_vision$car_dev, method = "spearman")
# 0.04358656

corr_all <- ggplot(mon_vision) +
  aes(x = eLog, y = car_dev) +
  geom_point(size = .8) + #colour = "#0c4c8a"
  xlab("Language predictioin") + ylab("Visuospatial predictioin") +
  geom_smooth(method=lm) + # for linear
  theme_gray(base_size = 12,
             base_family = "Times") +
  theme(legend.position = "bottom",
        legend.box = "vertical")

figs_path <- here("figs", "vision", "gca", "cont_speed_verb")
ggsave(paste0(figs_path, "/correlation_prediction.png"), corr_all, width = 180,
       height = 120, units = "mm", dpi = 600)


test <- cor.test(mon_vision$eLog, mon_vision$car_dev, method = "spearman", exact=FALSE)
test
# S = 2775871427, p-value = 0.02648
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.04358656 
report(test)














# add visuospatial proc speed effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_0 <- update(gca_mon_base,    . ~ . + corsi)
gca_mon_corsi_1 <- update(gca_mon_corsi_0, . ~ . + ot1:corsi)
gca_mon_corsi_2 <- update(gca_mon_corsi_1, . ~ . + ot2:corsi)

mon_corsi_anova <-
  anova(gca_mon_base, gca_mon_corsi_0, gca_mon_corsi_1,
        gca_mon_corsi_2)
#                 npar   AIC   BIC  logLik deviance   Chisq Df Pr(>Chisq)
# gca_mon_base       6 14363 14398 -7175.3    14351                        
# gca_mon_corsi_0    7 14362 14403 -7173.9    14348  2.7958  1  0.0945112 .  
# gca_mon_corsi_1    8 14355 14402 -7169.5    14339  8.7102  1  0.0031644 ** 
# gca_mon_corsi_2    9 14346 14399 -7164.0    14328 11.1517  1  0.0008396 ***

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_corsi_int_0 <- update(gca_mon_corsi_2, . ~ . + stress_sum:corsi)
gca_mon_corsi_int_1 <- update(gca_mon_corsi_int_0,   . ~ . + ot1:stress_sum:corsi)
gca_mon_corsi_int_2 <- update(gca_mon_corsi_int_1,   . ~ . + ot2:stress_sum:corsi)

mon_corsi_int_anova <-
  anova(gca_mon_corsi_2, gca_mon_corsi_int_0, gca_mon_corsi_int_1,
        gca_mon_corsi_int_2)
#                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_corsi_2        9 14346 14399 -7164.0    14328                       
# gca_mon_corsi_int_0   10 14348 14406 -7163.9    14328 0.0675  1    0.79497  
# gca_mon_corsi_int_1   11 14349 14413 -7163.3    14327 1.3167  1    0.25118  
# gca_mon_corsi_int_2   12 14346 14416 -7160.9    14322 4.8563  1    0.02755 *


# BRANCH #1
# add verbal proc time effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_wm_0 <- update(gca_mon_corsi_int_2,    . ~ . + ospan)
gca_mon_wm_1 <- update(gca_mon_wm_0, . ~ . + ot1:ospan)
gca_mon_wm_2 <- update(gca_mon_wm_1, . ~ . + ot2:ospan)

mon_wm_anova <-
  anova(gca_mon_corsi_int_2, gca_mon_wm_0, gca_mon_wm_1,
        gca_mon_wm_2)
#                       npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_corsi_int_2     12 14346 14416 -7160.9    14322                       
# gca_mon_wm_0            13 14347 14423 -7160.5    14321 0.6621  1     0.4158  
# gca_mon_wm_1            14 14346 14428 -7158.9    14318 3.2690  1     0.0706 .
# gca_mon_wm_2            15 14347 14435 -7158.7    14317 0.3513  1     0.5534

gca_mon_wm_i_0 <- update(gca_mon_corsi_int_2, . ~ . + stress_sum:ospan)
gca_mon_wm_i_1 <- update(gca_mon_wm_i_0,   . ~ . + ot1:stress_sum:ospan)
gca_mon_wm_i_2 <- update(gca_mon_wm_i_1,   . ~ . + ot2:stress_sum:ospan)

mon_wm_i_anova <-
  anova(gca_mon_corsi_int_2, gca_mon_wm_i_0, gca_mon_wm_i_1,
        gca_mon_wm_i_2)
#                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_corsi_int_2   12 14346 14416 -7160.9    14322                       
# gca_mon_wm_i_0        13 14348 14424 -7160.8    14322 0.0861  1    0.76919  
# gca_mon_wm_i_1        14 14347 14429 -7159.3    14319 2.9744  1    0.08459 .
# gca_mon_wm_i_2        15 14348 14436 -7158.8    14318 1.0345  1    0.30909  


gca_mon_wm_in_0 <- update(gca_mon_corsi_int_2,    . ~ . + ospan:corsi)
gca_mon_wm_in_1 <- update(gca_mon_wm_in_0,   . ~ . + ot1:ospan:corsi)
gca_mon_wm_in_2 <- update(gca_mon_wm_in_1,   . ~ . + ot2:ospan:corsi)

mon_wm_in_anova <-
  anova(gca_mon_corsi_int_2, gca_mon_wm_in_0, gca_mon_wm_in_1,
        gca_mon_wm_in_2)
#                     npar   AIC   BIC  logLik deviance    Chisq Df Pr(>Chisq)  
# gca_mon_corsi_int_2   12 14346 14416 -7160.9    14322                       
# gca_mon_wm_in_0       13 14347 14423 -7160.4    14321 0.8050  1    0.36962  
# gca_mon_wm_in_1       14 14342 14424 -7157.1    14314 6.6257  1    0.01005 *
# gca_mon_wm_in_2       15 14344 14432 -7157.1    14314 0.0455  1    0.83103 


# add 3-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_wm_int_0 <- update(gca_mon_wm_in_1, . ~ . + stress_sum:ospan:corsi)
gca_mon_wm_int_1 <- update(gca_mon_wm_int_0,   . ~ . + ot1:stress_sum:ospan:corsi)
gca_mon_wm_int_2 <- update(gca_mon_wm_int_1,   . ~ . + ot2:stress_sum:ospan:corsi)

mon_wm_int_anova <-
  anova(gca_mon_wm_in_1, gca_mon_wm_int_0, gca_mon_wm_int_1,
        gca_mon_wm_int_2)
#                  npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# gca_mon_wm_in_1    14 14342 14424 -7157.1    14314                     
# gca_mon_wm_int_0   15 14344 14432 -7157.1    14314 0.1609  1     0.6883
# gca_mon_wm_int_1   16 14345 14439 -7156.6    14313 0.9804  1     0.3221
# gca_mon_wm_int_2   17 14347 14447 -7156.5    14313 0.1332  1     0.7151


# BRANCH #2
# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_carsc_0 <- update(gca_mon_corsi_int_2,    . ~ . + car_dev)
gca_mon_carsc_1 <- update(gca_mon_carsc_0, . ~ . + ot1:car_dev)
gca_mon_carsc_2 <- update(gca_mon_carsc_1, . ~ . + ot2:car_dev)

mon_carsc_anova <-
  anova(gca_mon_corsi_int_2, gca_mon_carsc_0, gca_mon_carsc_1,
        gca_mon_carsc_2)
#                     npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_corsi_int_2   12 14346 14416 -7160.9    14322                       
# gca_mon_carsc_0       13 14346 14422 -7160.0    14320 1.6584  1    0.19782  
# gca_mon_carsc_1       14 14342 14424 -7157.0    14314 6.0023  1    0.01429 *
# gca_mon_carsc_2       15 14343 14431 -7156.5    14313 1.0416  1    0.30744

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_carsc_i_0 <- update(gca_mon_carsc_1,    . ~ . + car_dev:stress_sum)
gca_mon_carsc_i_1 <- update(gca_mon_carsc_i_0, . ~ . + ot1:car_dev:stress_sum)
gca_mon_carsc_i_2 <- update(gca_mon_carsc_i_1, . ~ . + ot2:car_dev:stress_sum)

mon_carsc_i_anova <-
  anova(gca_mon_carsc_1, gca_mon_carsc_i_0, gca_mon_carsc_i_1,
        gca_mon_carsc_i_2)
#                   npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_carsc_1     14 14342 14424 -7157.0    14314                     
# gca_mon_carsc_i_0   15 14344 14432 -7157.0    14314 0.0008  1     0.9776
# gca_mon_carsc_i_1   16 14346 14440 -7156.9    14314 0.1589  1     0.6902
# gca_mon_carsc_i_2   17 14348 14448 -7156.9    14314 0.0001  1     0.9918


gca_mon_carsc_in_0 <- update(gca_mon_carsc_1,     . ~ . + car_dev:corsi)
gca_mon_carsc_in_1 <- update(gca_mon_carsc_in_0, . ~ . + ot1:car_dev:corsi)
gca_mon_carsc_in_2 <- update(gca_mon_carsc_in_1, . ~ . + ot2:car_dev:corsi)

mon_carsc_in_anova <-
  anova(gca_mon_carsc_1, gca_mon_carsc_in_0, gca_mon_carsc_in_1,
        gca_mon_carsc_in_2)
#                    npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_carsc_1      14 14342 14424 -7157.0    14314                       
# gca_mon_carsc_in_0   15 14341 14428 -7155.3    14311 3.4579  1    0.06295 .
# gca_mon_carsc_in_1   16 14342 14436 -7155.3    14310 0.0331  1    0.85560  
# gca_mon_carsc_in_2   17 14344 14444 -7155.3    14310 0.0171  1    0.89602


gca_mon_carsc_int_0 <- update(gca_mon_carsc_1, . ~ . + stress_sum:car_dev:corsi)
gca_mon_carsc_int_1 <- update(gca_mon_carsc_int_0,   . ~ . + ot1:stress_sum:car_dev:corsi)
gca_mon_carsc_int_2 <- update(gca_mon_carsc_int_1,   . ~ . + ot2:stress_sum:car_dev:corsi)

mon_carsc_int_anova <-
  anova(gca_mon_carsc_1, gca_mon_carsc_int_0, gca_mon_carsc_int_1,
        gca_mon_carsc_int_2)
#                     npar   AIC   BIC  logLik deviance   Chisq Df Pr(>Chisq)  
# gca_mon_carsc_1       14 14342 14424 -7157.0    14314                     
# gca_mon_carsc_int_0   15 14344 14432 -7157.0    14314 0.0433  1     0.8352
# gca_mon_carsc_int_1   16 14344 14438 -7156.1    14312 1.7692  1     0.1835
# gca_mon_carsc_int_2   17 14346 14446 -7156.1    14312 0.0061  1     0.9378


car::vif(gca_mon_wm_in_1)
# ot1                  ot2                corsi 
# 1.058848             1.067298             1.866491 
# ot1:corsi            ot2:corsi     corsi:stress_sum 
# 2.003956             1.083209             1.041388 
# corsi:ospan ot1:corsi:stress_sum ot2:corsi:stress_sum 
# 1.827475             1.051210             1.083087 
# ot1:corsi:ospan 
# 1.983169 

car::vif(gca_mon_carsc_1)
# ot1                     ot2                     corsi 
# 1.068870             1.065002             1.027411 
# car_dev            ot1:corsi            ot2:corsi 
# 1.007993             1.070206             1.069940 
# corsi:stress_sum          ot1:car_dev ot1:corsi:stress_sum 
# 1.043949             1.028380             1.044958 
# ot2:corsi:stress_sum 
# 1.084214 



# save models  
mod_type <- "gca_mon"
mod_spec <- c("_base", 
              "_stress_0", "_stress_1", "_stress_2", 
              "_corsi_0", "_corsi_1", "_corsi_2", 
              "_corsi_int_0", "_corsi_int_1", "_corsi_int_2", 
              "_wm_0", "_wm_1", "_wm_2", 
              "_wm_i_0", "_wm_i_1", "_wm_i_2", 
              "_wm_in_0", "_wm_in_1", "_wm_in_2", 
              "_wm_int_0", "_wm_int_1", "_wm_int_2", 
              "_carsc_0", "_carsc_1", "_carsc_2", 
              "_carsc_i_0", "_carsc_i_1", "_carsc_i_2",
              "_carsc_in_0", "_carsc_in_1", "_carsc_in_2",
              "_carsc_int_0", "_carsc_int_1", "_carsc_int_2"
)

# Store ind models in list
mon_mods_zerosc <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods_zerosc,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "mon_mods_zerosc.Rdata"))


# ---------------------------------------------



  
# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_verbalps <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot2, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(corsi_rt = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan_rt = c(-1, 0, 1)))

# Get model predictions and SE
fits_mon_verbalps <- predictSE(gca_mon_ospan_int_1, mon_verbalps) %>%        
  as_tibble %>%
  bind_cols(mon_verbalps) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_mon_verbalps <- filter(fits_mon_verbalps, time_zero == 4) %>%
  select(stress = stress_sum, ospan_rt, corsi_rt,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 





# Create design dataframe for predictions
mon_car <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot2, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(corsi_rt = c(-1, 0, 1))) %>%
  expand_grid(., tibble(car_dev = c(-1, 0, 1)))

# Get model predictions and SE
fits_mon_car <- predictSE(gca_mon_corsirt_int_2, mon_car) %>%        
  as_tibble %>%
  bind_cols(mon_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_mon_car <- filter(fits_mon_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, corsi_rt,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 





# Create design dataframe for predictions
mon_verbalwm <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot2, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(corsi = c(-1, 0, 1))) %>%
  expand_grid(., tibble(ospan = c(-1, 0, 1)))

# Get model predictions and SE
fits_mon_verbalwm <- predictSE(gca_mon_wm_in_1, mon_verbalwm) %>%        
  as_tibble %>%
  bind_cols(mon_verbalwm) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_mon_verbalwm <- filter(fits_mon_verbalwm, time_zero == 4) %>%
  select(stress = stress_sum, ospan, corsi,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 




# Create design dataframe for predictions
mon_carsc <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot2, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(corsi = c(-1, 0, 1))) %>%
  expand_grid(., tibble(car_dev = c(-1, 0, 1)))

# Get model predictions and SE
fits_mon_carsc <- predictSE(gca_mon_carsc_1, mon_carsc) %>%        
  as_tibble %>%
  bind_cols(mon_carsc) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_mon_carsc <- filter(fits_mon_carsc, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, corsi,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 


# Save models predictions
model_preds_zero0 <- mget(c('fits_mon_verbalps', 'fits_mon_car',
                            'fits_mon_verbalwm', 'fits_mon_carsc',
                      "preds_mon_verbalps", "preds_mon_car",
                      "preds_mon_verbalwm", "preds_mon_carsc"
                      ))

save(model_preds_zero0,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "model_preds_zero0.Rdata"))



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








# l2_data <- stress_gc_subset%>%
#   filter(., l1 != 'es') %>% 
#   mutate(., l1_sum = if_else(l1 == 'en', -1, 1),
#          use_z = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week),
#          DELE_z = (DELE - mean(DELE))/sd(DELE)
#   ) 

