











# T-test plots

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
    scale_x_discrete(labels = c("mon", 'aes', 'ies', 'ams', 'ims')) +
    labs(title = 'Mean fixations at target syllable offset',
         y = 'Target fixations', x = '', caption = 'Mean and lower-bound 99% CI') +
    theme_minimal(base_size = 12, base_family = 'Times') -> stressTargetFixMOD1)

# ggsave('stressTargetFixMOD1.png', plot = stressTargetFixMOD1, dpi = 600, device = "png",          # already commented by C.
#         path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#         height = 4, width = 7, unit = 'in')
















# GLMM plots


et_ci <- confint(prop_0_mod_final, method = "Wald", level = 0.99) %>%
  as.data.frame(.) %>%
  slice(., 13:16) %>%
  rename(., ciLow = `0.5 %`, ciHi = `99.5 %`)

stressFixModP0 <- broom::tidy(prop_0_mod_final) %>% slice(1:4) %>%
  cbind(., et_ci) %>%
  mutate(., term = recode(term, `(Intercept)` = '(Intercept)',         # not sure if something else needs to be changed for terms
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

# ggsave('stressP0.png', plot = stressFixModP0, dpi = 600, device = "png",                   # already commented by C.
#         path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#         height = 2.5, width = 7, unit = 'in')



# If cond does not apply, do not run

stressFixModsHLS <- df_stress %>%
  mutate(., group = factor(group, levels = c('mon', 'aes', 'ies', 'ams', 'ims')) %>%
  group_by(., group, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix,
                group = interaction(group, coda), dodge = coda, color = coda)) +           # change?
  geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
               position = position_dodge(width = 0.5), color = 'black',
               size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.95)) +
  stat_summary(fun.y = mean, geom = 'point', size = 2.75,
               position = position_dodge(width = 0.5), show.legend = TRUE) +
  labs(y = 'Target fixations', x = '', caption = 'Mean +/- 99% bootstrap CI.',
       title = "Mean fixations at target syllable offset") +
  coord_cartesian(ylim = c(0.4, 1)) +
  scale_x_discrete(labels = c('mon', 'aes', 'ies', 'ams', 'ims')) +
  scale_color_brewer(name = '', palette = 'Set1', labels = c("Present", "Preterit")) +
  theme_minimal(base_size = 12, base_family = 'Times')

# ggsave('stressFixModsHLS.png', plot = stressFixModsHLS, dpi = 600, device = "png",                  # already commented by C.
#         path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#         height = 4, width = 7, unit = 'in')











# Calculate mean target fixation as a function of group, condition,
# for each participant. We will plot the mean and calculate the
# bootstrapped 95% confidence interval and plot it all.

stressFixModsP1 <- df_stress %>%
  mutate(., group = factor(group, levels = c("mon", 'aes', 'ies', 'ams', 'ims'))) %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix, shape = condition,
                group = interaction(group, condition), dodge = condition, color = group)) +
  geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
               position = position_dodge(width = 0.5), color = 'black',
               size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.99)) +
  stat_summary(fun.y = mean, geom = 'point', size = 2.75,
               position = position_dodge(width = 0.5)) +
  ylim(0, 1) +
  labs(y = '% Correct', x = NULL, caption = '') +
  scale_x_discrete(labels = c('mon', 'aes', 'ies', 'ams', 'ims')) +
  scale_color_manual(name = '', values = c('grey90', 'grey75', 'grey55'), guide = FALSE) +
  scale_shape_manual(name = '', values = c(16, 17), labels = c('Present', 'Preterit')) +
  guides(shape = guide_legend(override.aes = list(shape = c(1, 2), color = 'black'))) +
  theme_bw(base_size = 15, base_family = 'Times') +
  theme(legend.position = c(0.26, 0.14),
        legend.box.just = "left",
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.75, 'lines'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# arrange all plots together
#stress_eyet_plot <- stressFixModsP1 + stressFixModsP2               # not necessary since I only have one plot

# ggsave('stress_eyet_plot.png',                                                    # already commented by C.
#        plot = stress_eyet_plot, dpi = 600, device = "png",
#       path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#        height = 3.5, width = 8.5, units = 'in')

