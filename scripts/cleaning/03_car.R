# Clean car prediction scores
# select only complete cases
# remove outliers
# calculate descriptives
# remove trials where button press happened after 400 ms of car's supposed reappearance
# calculate random effects per participant

# Last update: 2021-04-07

source(here::here("scripts", "00_load_libs.R"))
source(here::here("scripts", "01_helpers.R"))

car <- read_csv(here("data", 'clean', 'car_sec_long.csv'))

# How much data we lose by selecting only good trials (!= NA & <400ms) in proportion
mean(complete.cases(car$mean_dev_sc))
# 0.7536688

car <- car %>%
  select(., -X1) %>%
  filter(., complete.cases(car$mean_dev_sc))

d <- density(car$mean_dev_sc)
plot(d, main="Density of visuospatial prediction\ntime deviations from ms 0,\nwhen car should reappear")

# eliminate outliers
Q <- quantile(car$mean_dev_sc, probs=c(.25, .75), na.rm = FALSE)

plot(density((Q)))

iqr <- IQR(car$mean_dev_sc)

car <- car %>%
  filter(., car$mean_dev_sc > (Q[1] - 1.5*iqr) & 
           car$mean_dev_sc < (Q[2]+1.5*iqr)) %>%
  filter(., mean_dev_sc < 0.4) %>%
  mutate(direction_sum = ifelse(direction == 'right', -1, 1),
         speed = fct_relevel(speed, "slow", "medium", "fast"))

plot(density(car$mean_dev_sc))


car %>%
  group_by(., direction, speed) %>%
  summarize(., mean = mean(mean_dev_sc),
            sd = sd(mean_dev_sc))
#   direction speed     mean    sd
# 1 left      fast    0.0686 0.181
# 2 left      medium  0.0459 0.184
# 3 left      slow   -0.0190 0.271
# 4 right     fast    0.0395 0.166
# 5 right     medium  0.0517 0.188
# 6 right     slow   -0.0357 0.242

car %>%
  group_by(., direction) %>%
  summarize(., mean = mean(mean_dev_sc),
            sd = sd(mean_dev_sc))
#   direction   mean    sd
# 1 left      0.0302 0.221
# 2 right     0.0159 0.207

car %>%
  group_by(., speed) %>%
  summarize(., mean = mean(mean_dev_sc),
            sd = sd(mean_dev_sc))
#   speed     mean    sd
# 1 fast    0.0546 0.174
# 2 medium  0.0488 0.185
# 3 slow   -0.0275 0.257


car %>%
  summarize(., mean = mean(mean_dev_sc),
            sd = sd(mean_dev_sc))
#     mean    sd
# 1 0.0231 0.214


car_glm <-  lmer(mean_dev_sc ~ speed + direction_sum + (1 | subject_id),
                 data = car)

car <- car %>%
  mutate(speed = fct_relevel(speed, 'fast', 'medium', 'slow'))

car_fast <- update(car_glm)


summary(car_glm)
#                  Estimate Std. Error t value
# (Intercept)     -0.010219   0.016781  -0.609
# speedmedium      0.104417   0.013419   7.781
# speedfast        0.097063   0.013755   7.056
# direction_sum    0.011337   0.005533   2.049

confint(car_glm)
#                       2.5 %     97.5 %
# .sig01         0.1461431327 0.19244006
# .sigma         0.1260163329 0.14303929
# (Intercept)   -0.0430682299 0.02284137
# speedmedium    0.0780435825 0.13068439
# speedfast      0.0701092996 0.12397897
# direction_sum  0.0004937954 0.02216338


summary(car_fast)
#                  Estimate Std. Error t value
# (Intercept)      0.086844   0.017214   5.045
# speedmedium      0.007355   0.013985   0.526
# speedslow       -0.097063   0.013755  -7.056
# direction_sum    0.011337   0.005533   2.049

confint(car_fast)
#                       2.5 %      97.5 %
# .sig01         0.1461431327  0.19244006
# .sigma         0.1260163329  0.14303929
# (Intercept)    0.0531608229  0.12079936
# speedmedium   -0.0200762461  0.03471913
# speedslow     -0.1239789705 -0.07010930
# direction_sum  0.0004937954  0.02216338



pretty_fixed_effects <- car_glm %>%  
  tidy_lme4() %>% 
  mutate(p = format_pval(p), 
         Parameter = fix_param_names(Parameter)) %>% 
  mutate_each(funs(format_fixef_num), Estimate:t) %>% 
  rename(`_t_` = t, `_p_` = p) 

pretty_fixed_effects %>% 
  select(-effect) %>%
  knitr::kable(format = "pandoc", align = str_tokenize("lrrrr")) 

# Parameter              Estimate      SE            _t_      _p_
# ----------------  -------------  ------  -------------  -------
# Intercept          &minus;0.010   0.017   &minus;0.609     .543
# speedmedium               0.104   0.013          7.781   < .001
# speedfast                 0.097   0.014          7.056   < .001
# direction_sum             0.011   0.006          2.049     .040


pretty_fixed_effects <- car_fast %>%  
  tidy_lme4() %>% 
  mutate(p = format_pval(p), 
         Parameter = fix_param_names(Parameter)) %>% 
  mutate_each(funs(format_fixef_num), Estimate:t) %>% 
  rename(`_t_` = t, `_p_` = p) 

pretty_fixed_effects %>% 
  select(-effect) %>%
  knitr::kable(format = "pandoc", align = str_tokenize("lrrrr")) 
# Parameter              Estimate      SE            _t_      _p_
# ----------------  -------------  ------  -------------  -------
# Intercept                 0.087   0.017          5.045   < .001
# speedmedium               0.007   0.014          0.526     .599
# speedslow          &minus;0.097   0.014   &minus;7.056   < .001
# direction_sum             0.011   0.006          2.049     .040





car_ranef <- ranef(car_glm) %>% as_tibble()  

car_sel <- car_ranef %>%
  select(., grp, condval, condsd) %>%
  rename(., participant = grp,
         car_dev = condval,
         car_sd = condsd)

car_sel$participant <- str_replace(car_sel$participant, "ae", "aes")
car_sel$participant <- str_replace(car_sel$participant, "ie", "ies")
car_sel$participant <- str_replace(car_sel$participant, "am", "ams")
car_sel$participant <- str_replace(car_sel$participant, "im", "ims")
car_sel$participant <- str_replace(car_sel$participant, "mo", "mon")

write.csv(car_sel,'./data/clean/vision_scores_nooutliers-400.csv', row.names = F)

