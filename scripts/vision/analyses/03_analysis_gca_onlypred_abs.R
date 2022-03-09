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


vision <- read_csv("./data/clean/vision_scores_nooutliers.csv") # pred car

vision$car_dev <- abs(vision$car_dev)

vision50 <- left_join(x = stress50, y = vision, by = "participant", all.x=TRUE)

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
    lmer(eLog ~ 1 + 
           (1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = mon_vision, weights = 1/wts, REML = F)
  
  mod_ot2 <- update(mod_ot1, . ~ . -(1 | participant) + ot1 + (1 + ot1 | participant))
  
  mod_ot3 <- update(mod_ot2, . ~ . -(1 + ot1 | participant) + ot2 + (1 + ot1 | participant))
  
  anova(mod_ot1, mod_ot2, mod_ot3)
  #         npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot1    3 22630 22649 -11312    22624                          
  # mod_ot2    6 22394 22432 -11191    22382 241.921  3  < 2.2e-16 ***
  # mod_ot3    7 22376 22420 -11181    22362  20.034  1  7.608e-06 ***
  
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  #mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) + (1 + ot1 + ot2 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5) #, mod_ot6)
  #         npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot3    7 22376 22420 -11181    22362                         
  # mod_ot4    8 22290 22340 -11137    22274 88.716  1  < 2.2e-16 ***
  # mod_ot5   10 22278 22341 -11129    22258 15.672  2  0.0003952 ***
  
  
}

# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_mon_base <- mod_ot5
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
# gca_mon_base       10 22278 22341 -11129    22258                     
# gca_mon_stress_0   11 22280 22349 -11129    22258 0.1275  1     0.7210
# gca_mon_stress_1   12 22280 22356 -11128    22256 1.4732  1     0.2248
# gca_mon_stress_2   13 22282 22364 -11128    22256 0.4400  1     0.5071


# add visuospatial anticipation effect to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_0 <- update(gca_mon_base,    . ~ . + car_dev)
gca_mon_car_1 <- update(gca_mon_car_0, . ~ . + ot1:car_dev)
gca_mon_car_2 <- update(gca_mon_car_1, . ~ . + ot2:car_dev)

mon_car_anova <-
  anova(gca_mon_base, gca_mon_car_0, gca_mon_car_1,
        gca_mon_car_2)
#               npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)
# gca_mon_base    10 22278 22341 -11129    22258                       
# gca_mon_car_0   11 22276 22345 -11127    22254 3.9637  1    0.04649 *
# gca_mon_car_1   12 22278 22354 -11127    22254 0.0657  1    0.79772  
# gca_mon_car_2   13 22280 22362 -11127    22254 0.0252  1    0.87375  

# add 2-way int to intercept, linear slope, quadratic, and cubic time terms
gca_mon_car_int_0 <- update(gca_mon_car_0,    . ~ . + car_dev:stress_sum)
gca_mon_car_int_1 <- update(gca_mon_car_int_0, . ~ . + ot1:car_dev:stress_sum)
gca_mon_car_int_2 <- update(gca_mon_car_int_1, . ~ . + ot2:car_dev:stress_sum)

mon_car_int_anova <-
  anova(gca_mon_car_0, gca_mon_car_int_0, gca_mon_car_int_1,
        gca_mon_car_int_2)
#                   npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# gca_mon_car_0       11 22276 22345 -11127    22254                     
# gca_mon_car_int_0   12 22278 22354 -11127    22254 0.1613  1     0.6879
# gca_mon_car_int_1   13 22280 22362 -11127    22254 0.0249  1     0.8746
# gca_mon_car_int_2   14 22282 22370 -11127    22254 0.2376  1     0.6260


car::vif(gca_mon_car_0)
#      ot1      ot2  car_dev 
# 1.003850 1.004470 1.000694 

summary(gca_mon_car_0)
# Estimate Std. Error t value
# (Intercept)  0.8758     0.1580   5.544
# ot1           2.4976     0.3081   8.107
# ot2           0.6452     0.1380   4.676
# car_dev      -1.1492     0.5520  -2.082

confint(gca_mon_car_0)
#                    2.5 %     97.5 %
# (Intercept)  0.874852321  0.8927308
# ot1          2.495704244  2.5316941
# ot2          0.644167382  0.6619923
# car_dev     -1.152942371 -1.0856081

performance::r2_nakagawa(gca_mon_car_0, by_group = FALSE, tolerance = 1e-05)
# Conditional R2: 0.137
# Marginal R2: 0.072

rsq::rsq(gca_mon_car_0)
# $model
# [1] 0.1466243
# 
# $fixed
# [1] 0.06620766
# 
# $random
# [1] 0.08041666

rsq::rsq(gca_mon_car_0, adj = TRUE)
# $model
# [1] 0.1459873
# 
# $fixed
# [1] 0.06551062
# 
# $random
# [1] 0.08047668

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
mon_mods_onlypred_abs <- mget(c(paste0(mod_type, mod_spec)
))

save(mon_mods_onlypred_abs,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "mon_mods_onlypred_abs.Rdata"))

  
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




  
# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_car <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot2, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(car_dev = c(-1, 0, 1)))
  

# Get model predictions and SE
fits_mon_car <- predictSE(gca_mon_car_0, mon_car) %>%        
  as_tibble %>%
  bind_cols(mon_car) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

# Filter preds at target syllable offset
preds_mon_car <- filter(fits_mon_car, time_zero == 4) %>%
  select(stress = stress_sum, car_dev, 
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) 


# Save models predictions
model_preds_onlypred_abs <- mget(c('fits_mon_car',
                      "preds_mon_car"
                      ))

save(model_preds_onlypred_abs,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "model_preds_onlypred_abs.Rdata"))



# -----------------------------------------------------------------------------
