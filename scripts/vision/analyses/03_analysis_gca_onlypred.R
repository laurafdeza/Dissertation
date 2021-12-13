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




  
# Model predictions for plotting ---------------------------------------------

# Create design dataframe for predictions
mon_car <- mon_vision %>%
  dplyr::select(time_zero, ot1:ot2, stress_sum) %>%
  distinct %>%
  expand_grid(., tibble(car_dev = c(-1, 0, 1)))
  

# Get model predictions and SE
fits_mon_car <- predictSE(gca_mon_car_1, mon_car) %>%        
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
model_preds_onlypred <- mget(c('fits_mon_car',
                      "preds_mon_car"
                      ))

save(model_preds_onlypred,
     file = here("mods", "vision", "gca", "cont_speed_verb",
                 "model_preds_onlypred.Rdata"))



# -----------------------------------------------------------------------------


mon_vision %>%
  summarise(., mean(car_dev),
               sd(car_dev),
               min(car_dev),
               max(car_dev),
               median(car_dev))
# `mean(car_dev)` `sd(car_dev)` `min(car_dev)` `max(car_dev)` `median(car_dev)`
# <dbl>         <dbl>          <dbl>          <dbl>             <dbl>
#   1         -0.0144         0.194         -0.451          0.212            0.0713              





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

