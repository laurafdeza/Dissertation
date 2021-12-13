# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
source(here::here("scripts", "01_helpers.R"))
#stress50 <- read_csv(here("data", "clean", "stress_50ms_allparticipants.csv"))


# Get path to saved models
gca_mods_path  <- here("mods", "music", "gca", "continuous", "wm")

# Load models as lists
load(paste0(gca_mods_path, "/ss_mods.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))


# Set path for saving figs
figs_path <- here("figs", "music", "gca")

# -----------------------------------------------------------------------------






# # Plot raw data ---------------------------------------------------------------
# 
# stress50$cond <- as.factor(as.character(stress50$cond))
# 
# condition_names <- c(
#   `1` = 'Present',
#   `2` = 'Preterit'
# )
# 
# 
# stress_p1 <- stress50 %>%
#     #na.omit(.) %>%
#     filter(., time_zero >= -10, time_zero <= 12, l1 != 'ma') %>%
#     mutate(., l1 = fct_relevel(l1, "es", "en")) %>%
#     ggplot(., aes(x = time_zero, y = target_prop, fill = l1)) +
#     facet_grid(cond ~ ., labeller = as_labeller(condition_names)) +
#     geom_hline(yintercept = 0.5, color = 'grey40', lty = 3) +
#     geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
#     geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
#     stat_summary(fun.y = "mean", geom = "line", size = 1) +
#     stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
#                  stroke = 0.5, pch = 21) +
#     scale_fill_brewer(palette = 'Set1', name = "L1",
#                        labels = c("Spanish", "English")) +
#     scale_x_continuous(breaks = c(-10, 0, 10),
#                        labels = c("-500", "0", "500")) +
#     labs(y = 'Proportion of target fixations',
#          x = 'Time relative to target syllable offset (ms)',
#          caption = "Mean +/- 95% CI") +
#     annotate("text", x = 3.3, y = 0.02, label = '200ms',
#              angle = 90, size = 3, hjust = 0, family = 'Times') +
#     theme_grey(base_size = 12, base_family = "Times") +
#     theme(legend.position = 'bottom')
# 
# l1_names <- c(
#   `es` = 'Spanish speakers',
#   `en` = 'English speakers'
# )
# 
# stress_p2 <- stress50 %>%
#   #na.omit(.) %>%
#   filter(., time_zero >= -10, time_zero <= 12, l1 != 'ma') %>%
#   mutate(., l1 = fct_relevel(l1, "es", "en")) %>%
#   ggplot(., aes(x = time_zero, y = target_prop, fill = cond)) +
#   facet_grid(l1 ~ ., labeller = as_labeller(l1_names)) +
#   geom_hline(yintercept = 0.5, color = 'grey40', lty = 3) +
#   geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
#   geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
#   stat_summary(fun.y = "mean", geom = "line", size = 1) +
#   stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
#                stroke = 0.5, pch = 21) +
#   scale_fill_brewer(palette = 'Set1', name = "Tense",
#                     labels = c("Present", "Preterite")) +
#   scale_x_continuous(breaks = c(-10, 0, 10),
#                      labels = c("-500", "0", "500")) +
#   labs(y = 'Proportion of target fixations',
#        x = 'Time relative to target syllable offset (ms)',
#        caption = "Mean +/- 95% CI") +
#   annotate("text", x = 3.3, y = 0.02, label = '200ms',
#            angle = 90, size = 3, hjust = 0, family = 'Times') +
#   theme_grey(base_size = 12, base_family = "Times") +
#   theme(legend.position = 'bottom')
# 
# 
# ggsave('timecourse_l1.png',
#        plot = stress_p1, dpi = 600, device = "png",
#        path = figs_path,
#        height = 3.5, width = 4.5, units = 'in')
# ggsave('timecourse_tense.png',
#        plot = stress_p2, dpi = 600, device = "png",
#        path = figs_path,
#        height = 3.5, width = 4.5, units = 'in')



# -----------------------------------------------------------------------------




# Plot GCA --------------------------------------------------------------------

# mon

model_preds$fits_ss_pitch %>%
  #mutate(ospan_std = as.factor(ospan_std)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin, color = pitch_dev)) +
  facet_grid(. ~ ospan_std) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  # stat_summary(fun.y = "mean", geom = "line", size = 1) +
  # stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
  #              alpha = 0.2)
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       lty = 'Stress condition') +
  ggtitle('Spanish use') +
  theme_grey(base_size = 10, base_family = "Times") +
  theme(
        legend.position = 'bottom',
        plot.title = element_text(size = 9, hjust = 0.5),
        plot.margin = margin(t = 5, l = 5, r = 24))





# # L2 speakers
# 
# library(grid)
# library(gridExtra)
# 
# # all variables 
# plot <- fits_all_en %>%
#   mutate(Proficiency = as.factor(DELE_z),
#          `Spanish use` = as.factor(use_z),
#          Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
#          Stress = fct_relevel(Stress, 'Paroxytone')) %>%
#   ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
#                 lty = Stress)) +
#   facet_grid(Proficiency ~ `Spanish use`) +
#   geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
#   geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
#   geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
#   geom_line(size = 0.35) +
#   scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
#                      labels = c("-200", "0", "200", "400", "600")) +
#   labs(x = "Time (ms) relative to target syllable offset",
#        y = "Empirical logit of looks to target",
#        lty = 'Stress condition') +
#   ggtitle('Spanish use') +
#   theme_grey(base_size = 10, base_family = "Times") + 
#   theme(
#         legend.position = 'bottom',
#         plot.title = element_text(size = 9, hjust = 0.5),
#         plot.margin = margin(t = 5, l = 5, r = 24))
# 
# plot2 <- grid.arrange(plot,     #en_gca_plot <- 
#              bottom = textGrob('Proficiency', rot = 270,
#                                x = .97, y = 4.2, gp = gpar(fontsize = 9,
#                                                            fontfamily = 'Times')))
# 
# 
# # Or Save using the 'export' button instead
# ggsave(paste0(figs_path, "/en_gca_all.png"), plot2, width = 180,
#        height = 120, units = "mm", dpi = 600)
# 
# 
# # stress
# stress_en <- fits_all_en %>%
#   mutate(Proficiency = as.factor(DELE_z),
#          `Spanish use` = as.factor(use_z),
#          Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
#          Stress = fct_relevel(Stress, 'Paroxytone')) %>%
#   ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
#                 color = Stress)) +
#   # facet_grid(Proficiency ~ `Spanish use`) +
#   geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
#   geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
#   # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
#   stat_summary(fun.y = "mean", geom = "line", size = 1) +  
#   stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
#                alpha = 0.2) +
#   scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
#                      labels = c("-200", "0", "200", "400", "600")) +
#   labs(x = "Time (ms) relative to target syllable offset",
#        y = "Empirical logit of looks to target",
#        color = 'Stress condition') +
#   theme_grey(base_size = 10, base_family = "Times") + 
#   theme(
#     legend.position = 'bottom')
# 
# ggsave(paste0(figs_path, "/en_gca_stress.png"), stress_en, width = 180,
#        height = 100, units = "mm", dpi = 600)
# 
# # Proficiency
# prof_en <- fits_all_en %>%
#   mutate(Proficiency = as.factor(DELE_z),
#          `Spanish use` = as.factor(use_z),
#          Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
#          Stress = fct_relevel(Stress, 'Paroxytone')) %>%
#   ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
#                 lty = Proficiency)) +
#   geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
#   geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
#   # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
#   stat_summary(fun.y = "mean", geom = "line", size = .5) +  
#   stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
#                alpha = 0.2) +
#   scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
#                      labels = c("-200", "0", "200", "400", "600")) +
#   labs(x = "Time (ms) relative to target syllable offset",
#        y = "Empirical logit of looks to target") +
#   theme_grey(base_size = 10, base_family = "Times") + 
#   theme(
#     legend.position = 'bottom')
# 
# ggsave(paste0(figs_path, "/en_gca_prof.png"), prof_en, width = 180,
#        height = 100, units = "mm", dpi = 600)
# 
# # Spanish use
# use_en <- fits_all_en %>%
#   mutate(Proficiency = as.factor(DELE_z),
#          `Spanish use` = as.factor(use_z),
#          Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
#          Stress = fct_relevel(Stress, 'Paroxytone')) %>%
#   ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
#                 lty = `Spanish use`)) +
#   geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
#   geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
#   # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
#   stat_summary(fun.y = "mean", geom = "line", size = .5) +  
#   stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
#                alpha = 0.2) +
#   scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
#                      labels = c("-200", "0", "200", "400", "600")) +
#   labs(x = "Time (ms) relative to target syllable offset",
#        y = "Empirical logit of looks to target") +
#   theme_grey(base_size = 10, base_family = "Times") + 
#   theme(
#     legend.position = 'bottom')
# 
# ggsave(paste0(figs_path, "/en_gca_use.png"), use_en, width = 180,
#        height = 100, units = "mm", dpi = 600)

