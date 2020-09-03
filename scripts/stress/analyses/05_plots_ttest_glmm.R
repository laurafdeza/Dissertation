
source(here::here("scripts", "00_load_libs.R"))



# T-test plots

stress_ttest <- readRDS('./mods/stress/stress_ttest.rds')
stress_ttest$cond <- as.factor(stress_ttest$cond)

# We will plot the models
# This will be almost exactly the same as
# the previous plot, but it will use the
# confidence interval from the test we
# actually conducted

(stress_ttest %>%
  ungroup(.) %>%
  mutate(., group = factor(group, levels = c("mon", 'aes', 'ies', 'ams', 'ims')),
            condition = factor(cond, levels = c("present", "preterit"))) %>%
  ggplot(., aes(x = group, y = estimate, color = cond,
                group = interaction(group, cond), dodge = cond)) +
    geom_hline(yintercept = 0.5, lty = 3) +
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                       position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    geom_point(position = position_dodge(width = 0.75), size = 2.75, color = 'grey90') +
    ylim(0.3, 1.0) +
    scale_x_discrete(labels = c("ESS", 'Adv EN', 'Int EN', 'Adv MA', 'Int MA')) +
    scale_color_discrete(name = "Condition", labels = c("Present", "Preterit")) +
    labs(title = 'Mean fixations at target syllable offset',
         y = 'Target fixations', x = '', caption = 'Mean and lower-bound 99% CI') +
    theme_minimal(base_size = 12, base_family = 'Times') +
    theme(plot.title = element_text(hjust = 0.65)) -> stressTargetFixMOD1)

ggsave('stressTargetFixMOD1.png', plot = stressTargetFixMOD1, dpi = 600, device = "png",          # already commented by C.
        path = "./figs/stress/",
        height = 4, width = 7, unit = 'in')










# load GLMM model

prop_0_mod_final <- readRDS('./mods/stress/glmm/12_prop_0_mod_final.rds')



# GLMM plots
# I changed the slices but I'm not sure the rows I selected are the ones that should be selected

et_ci <- confint(prop_0_mod_final, method = "Wald", level = 0.99) %>%
  as.data.frame(.) %>%
  slice(., 5:10) %>%   # slice changed here
  rename(., ciLow = `0.5 %`, ciHi = `99.5 %`)

stressFixModP0 <- broom.mixed::tidy(prop_0_mod_final) %>% slice(1:6) %>%   # and slice changed here too
  cbind(., et_ci) %>%
  mutate(., term = recode(term, `(Intercept)` = '(Intercept)',  
                          groupae = 'aes',
                          groupie = 'ies',
                          groupam = 'ams',
                          groupim = 'ims'),
         term = factor(term, levels = c('aes',
                                        'ies',
                                        'ams',
                                        'ims',
                                        '(Intercept)'))) %>%
  ggplot(., aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_errorbarh(aes(xmin = ciLow, xmax = ciHi), height = 0.2, size = 0.65) +
  geom_point(size = 3) +
  geom_point(size = 2, color = 'lightgrey') +
  labs(y = 'Term', x = 'Estimate +/- 95% CI') +
  theme_bw(base_size = 15, base_family = 'Times') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave('stressP0.png', plot = stressFixModP0, dpi = 600, device = "png",                   
        path = "./figs/stress/",
        height = 2.5, width = 7, unit = 'in')






# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
prop_0_mod_final %>%
#  group_by(., group, condition, participant) %>%
#  summarise(., meanFix = mean(target_prop)) %>%
  ggplot(., aes(x = group, y = mean(target_prop), 
                dodge = condition_sum, color = condition_sum, 
                group = interaction(group, condition_sum))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Group') + 
#  scale_x_discrete(labels = c("mon", "aes", "ies", "ams", "ims")) +         # reorder?
  ggtitle('Mean target fixations as a function of group\nand target type') +    
#  scale_color_brewer(palette = "Set1", name = '', labels = c('Present', 'Preterit')) +       #reorder
  theme_bw(base_size = 16, base_family = 'Times') -> stress_rel_target_fix





source(here::here("scripts", "02_load_data.R"))

df_stress <- stress10 %>%
  filter(., time_zero == 20) %>%
  mutate(., condition_sum = if_else(cond == "1", 1, -1))

df_stress$condition_sum <- as.factor(as.numeric(df_stress$condition_sum))

# If cond does not apply, do not run

stressFixModsHLS <- df_stress %>%
  mutate(., group = factor(group, levels = c('mon', 'aes', 'ies', 'ams', 'ims'))) %>%
  group_by(., group, participant) %>%
#  summarise(., meanFix = mean(target_prop)) %>%
  ggplot(., aes(x = group, y = mean(target_prop),
                group = interaction(group, condition_sum), dodge = condition_sum, color = condition_sum)) +           # change?
  geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
               position = position_dodge(width = 0.5), color = 'black',
               size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.95)) +
  stat_summary(fun.y = mean, geom = 'point', size = 2.75,
               position = position_dodge(width = 0.5), show.legend = TRUE) +
  labs(y = 'Target fixations', x = '', caption = 'Mean +/- 99% bootstrap CI.',
       title = "Mean fixations at target syllable offset") +
  coord_cartesian(ylim = c(0.4, 1)) +
#  scale_x_discrete(labels = c('ESS', 'Adv EN', 'Int EN', 'Adv MA', 'Int EN')) +
  scale_color_brewer(name = '', palette = 'Set1', labels = c("Present", "Preterit")) +
  theme_minimal(base_size = 12, base_family = 'Times')

ggsave('stressFixModsHLS.png', plot = stressFixModsHLS, dpi = 600, device = "png",                  # already commented by C.
        path = "./figs/stress/",
        height = 4, width = 7, unit = 'in')











# Calculate mean target fixation as a function of group, condition, WM
# for each participant. We will plot the mean and calculate the
# bootstrapped 95% confidence interval and plot it all.

stressFixModsP1 <- df_stress %>%
  mutate(., group = factor(group, levels = c("mon", 'aes', 'ies', 'ams', 'ims'))) %>%
  group_by(., group, condition_sum, participant) %>%
  summarise(., meanFix = mean(target_prop)) %>%
  ggplot(., aes(x = group, y = meanFix, shape = condition_sum,
                group = interaction(group, condition_sum), dodge = condition_sum, color = group)) +
  geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
               position = position_dodge(width = 0.5), color = 'black',
               size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.99)) +
  stat_summary(fun.y = mean, geom = 'point', size = 2.75,
               position = position_dodge(width = 0.5)) +
  ylim(0, 1) +
  labs(y = '% Correct', x = NULL, caption = '') +
  scale_x_discrete(labels = c('mon', 'aes', 'ies', 'ams', 'ims')) +
  scale_color_manual(name = '', values = c('grey90', 'grey75', 'grey55', 'white', 'black'), guide = FALSE) +
  scale_shape_manual(name = '', values = c(16, 17), labels = c('Present', 'Preterit')) +
  guides(shape = guide_legend(override.aes = list(shape = c(1, 2), color = 'black'))) +
  theme_bw(base_size = 15, base_family = 'Times') +
  theme(legend.position = c(0.26, 0.14),
        legend.box.just = "left",
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.75, 'lines'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())



stressFixModsP2 <- df_stress %>%
  filter(., cond == 1) %>%
  mutate(., group = factor(group, levels = c("mon", 'aes', 'ies', 'ams', 'ims'))) %>%
  group_by(., group, WM_set, participant) %>%
  summarise(., meanFix = mean(target_prop)) %>%
  ggplot(., aes(x = group, y = meanFix, color = WM_set,
                group = interaction(group, WM_set), dodge = WM_set)) +
  geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
               position = position_dodge(width = 0.5), color = 'black',
               size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.99)) +
  stat_summary(fun.y = mean, geom = 'point', size = 2.75,
               position = position_dodge(width = 0.5)) +
  ylim(0, 1) +
  labs(y = '% Correct', x = NULL, caption = '') +
  scale_x_discrete(labels = c('mon', 'aes', 'ies', 'ams', 'ims')) +
  scale_color_continuous(name = "WM score") +
  guides(shape = guide_legend(override.aes = list(shape = c(1, 2), color = 'black'))) +
  theme_bw(base_size = 15, base_family = 'Times') +
  theme(legend.position = 'none',
        legend.box.just = "left",
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.75, 'lines'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Present tense")




stressFixModsP3 <- df_stress %>%
  filter(., cond != 1) %>%
  mutate(., group = factor(group, levels = c("mon", 'aes', 'ies', 'ams', 'ims'))) %>%
  group_by(., group, WM_set, participant) %>%
  summarise(., meanFix = mean(target_prop)) %>%
  ggplot(., aes(x = group, y = meanFix, color = WM_set,
                group = interaction(group, WM_set), dodge = WM_set)) +
  geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
               position = position_dodge(width = 0.5), color = 'black',
               size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.99)) +
  stat_summary(fun.y = mean, geom = 'point', size = 2.75,
               position = position_dodge(width = 0.5)) +
  ylim(0, 1) +
  labs(y = '% Correct', x = NULL, caption = '') +
  scale_x_discrete(labels = c('mon', 'aes', 'ies', 'ams', 'ims')) +
  scale_color_continuous(name = "WM score") +
  guides(shape = guide_legend(override.aes = list(shape = c(1, 2), color = 'black'))) +
  theme_bw(base_size = 15, base_family = 'Times') +
  theme(#legend.position = c(0.26, 0.14),
    legend.box.just = "left",
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(0.75, 'lines'),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Preterit tense")


# arrange WM plots together
stress_WM_eyet_plot <- stressFixModsP2 + stressFixModsP3  


# arrange all plots together
#stress_eyet_plot <- stressFixModsP1 + stressFixModsP2   

ggsave('stressFixModsP1.png',
       plot = stressFixModsP1, dpi = 600, device = "png",
       path = "./figs/stress/",
       height = 3.5, width = 8.5, units = 'in')


ggsave('stress_WM_eyet_plot.png',
       plot = stress_WM_eyet_plot, dpi = 600, device = "png",
      path = "./figs/stress/",
       height = 3.5, width = 8.5, units = 'in')

