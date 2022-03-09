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
  # select(., -X1) %>%
  filter(., complete.cases(car$mean_dev_sc))

d <- density(car$mean_dev_sc)
plot(d, main="Density of visuospatial prediction\ntime deviations from ms 0,\nwhen car should reappear")

# eliminate outliers
mean(car$mean_dev_sc)
sd(car$mean_dev_sc)

sd2 <- sd(car$mean_dev_sc) * 2

sd2below <- mean(car$mean_dev_sc) - sd2
sd2above <- mean(car$mean_dev_sc) + sd2

car_tidy <- car %>%
  filter(., car$mean_dev_sc > sd2below & 
           car$mean_dev_sc < sd2above) %>%
  mutate(direction_sum = ifelse(direction == 'right', -1, 1),
         speed = fct_relevel(speed, "slow", "medium", "fast"))

plot(density(car_tidy$mean_dev_sc))

# calculate amount of data lost
(count(car) - count(car_tidy)) * 100 / count(car)
# 5.146036

# descriptives of tidy data
car_tidy %>%
  group_by(., direction, speed) %>%
  summarize(., mean = mean(mean_dev_sc),
            sd = sd(mean_dev_sc))
#   direction speed    mean    sd
# 1 left      fast   0.0662 0.313
# 2 left      medium 0.0887 0.218
# 3 left      fast   0.0720 0.184
# 4 right     slow   0.0460 0.271
# 5 right     medium 0.0878 0.219
# 6 right     fast   0.0520 0.178

car_tidy %>%
  group_by(., direction) %>%
  summarize(., mean = mean(mean_dev_sc),
            sd = sd(mean_dev_sc))
#   direction   mean    sd
# 1 left      0.0751 0.250
# 2 right     0.0616 0.230

car_tidy %>%
  group_by(., speed) %>%
  summarize(., mean = mean(mean_dev_sc),
            sd = sd(mean_dev_sc))
#   speed     mean    sd
# 1 fast    0.0563 0.293
# 2 medium  0.0882 0.218
# 3 fast    0.0622 0.181


car_tidy %>%
  summarize(., mean = mean(mean_dev_sc),
            sd = sd(mean_dev_sc))
#     mean    sd
# 1 0.0685 0.240


car_glm <-  lmer(mean_dev_sc ~ speed + direction_sum + (1 | subject_id),
                 data = car_tidy)

car_tidy <- car_tidy %>%
  mutate(speed = fct_relevel(speed, 'fast', 'medium', 'slow'))

car_fast <- update(car_glm)


summary(car_glm)
#                 Estimate Std. Error t value
# (Intercept)     0.065524   0.018797   3.486
# speedmedium     0.072036   0.013881   5.190
# speedfast       0.042031   0.014412   2.916
# direction_sum   0.014678   0.005782   2.539

confint(car_glm)
#                      2.5 %     97.5 %
# .sig01         0.176648484 0.22938161
# .sigma         0.138002613 0.15581943
# (Intercept)    0.028682794 0.10248720
# speedmedium    0.044730178 0.09921844
# speedfast      0.013723234 0.07024148
# direction_sum  0.003339145 0.02599405


summary(car_fast)
#                Estimate Std. Error t value
# (Intercept)    0.107555   0.019616   5.483
# speedmedium    0.030005   0.014766   2.032
# speedslow     -0.042031   0.014412  -2.916
# direction_sum  0.014678   0.005782   2.539

confint(car_fast)
#                      2.5 %      97.5 %
# .sig01         0.176648484  0.22938161
# .sigma         0.138002613  0.15581943
# (Intercept)    0.069165350  0.14625742
# speedmedium    0.001080125  0.05890611
# speedslow     -0.070241477 -0.01372323
# direction_sum  0.003339145  0.02599405



pretty_fixed_effects <- car_glm %>%  
  tidy_lme4() %>% 
  mutate(p = format_pval(p), 
         Parameter = fix_param_names(Parameter)) %>% 
  mutate_each(funs(format_fixef_num), Estimate:t) %>% 
  rename(`_t_` = t, `_p_` = p) 

pretty_fixed_effects %>% 
  select(-effect) %>%
  knitr::kable(format = "pandoc", align = str_tokenize("lrrrr")) 

# Parameter        Estimate      SE     _t_      _p_
# --------------  ---------  ------  ------  -------
# Intercept           0.066   0.019   3.486   < .001
# speedmedium         0.072   0.014   5.190   < .001
# speedfast           0.042   0.014   2.916     .004
# direction_sum       0.015   0.006   2.539     .011


pretty_fixed_effects <- car_fast %>%  
  tidy_lme4() %>% 
  mutate(p = format_pval(p), 
         Parameter = fix_param_names(Parameter)) %>% 
  mutate_each(funs(format_fixef_num), Estimate:t) %>% 
  rename(`_t_` = t, `_p_` = p) 

pretty_fixed_effects %>% 
  select(-effect) %>%
  knitr::kable(format = "pandoc", align = str_tokenize("lrrrr")) 
# Parameter            Estimate      SE            _t_      _p_
# --------------  -------------  ------  -------------  -------
# Intercept               0.108   0.020          5.483   < .001
# speedmedium             0.030   0.015          2.032     .042
# speedslow        &minus;0.042   0.014   &minus;2.916     .004
# direction_sum           0.015   0.006          2.539     .011





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

write.csv(car_sel,'./data/clean/vision_scores_nooutliers.csv', row.names = F)

car_sel <- read.csv('./data/clean/vision_scores_nooutliers.csv')
# distribution plot for visuospatial prediction deviation times
p <- ggplot(car_sel, aes(x=car_dev)) + 
  geom_density() +
  geom_vline(aes(xintercept=mean(car_dev)),
             color="black", linetype="dashed", size=.5) +
  geom_vline(aes(xintercept=sd(car_dev)),
             color="black", linetype="dotted", size=.5) +
  geom_vline(aes(xintercept=-sd(car_dev)),
             color="black", linetype="dotted", size=.5) +
  scale_x_continuous(breaks=round(seq(-.40, .40, .1), 2)) +
  labs(y = 'Density',
       x = 'Deviation time from reappearance millisecond (0.00)') +
  theme_grey(base_size = 12, base_family = "Times")
  
ggsave('devtime_density.png',
       plot = p, dpi = 600, device = "png",
       path = "./figs/vision/",
       height = 3.5, width = 4.5, units = 'in')

ggsave('devtime_density_long.png',
       plot = p, dpi = 600, device = "png",
       path = "./figs/vision/",
       height = 2.5, width = 4.5, units = 'in')

