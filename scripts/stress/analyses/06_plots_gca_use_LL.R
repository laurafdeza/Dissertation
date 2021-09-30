# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
source(here::here("scripts", "01_helpers.R"))
#source(here::here("scripts", "02_load_data.R"))
stress50 <- read_csv(here("data", "clean", "stress_50ms_final_onsetc3updated.csv"))

# Get path to saved models
gca_mods_path  <- here("mods", "stress", "gca", "LL_changes")

# Load models as lists and store in global environment
load(paste0(gca_mods_path, "/gca_mon_mods.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_mods.Rdata"))
load(paste0(gca_mods_path, "/nested_model_l2.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_mods_ints.Rdata"))
# load(paste0(gca_mods_path, "/model_preds_ints.Rdata"))
load(paste0(gca_mods_path, "/model_preds_l1ot2.Rdata"))

list2env(gca_mon_mods, globalenv())
list2env(gca_l2_mods, globalenv())
list2env(model_preds, globalenv())
list2env(gca_l2_mods_ints, globalenv())
list2env(model_preds_ints, globalenv())

# Set path for saving figs
figs_path <- here("figs", "stress", "gca", "LL_changes")

# -----------------------------------------------------------------------------






# Plot raw data ---------------------------------------------------------------

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
#     filter(., time_zero >= -10, time_zero <= 20) %>%
#     mutate(., l1 = fct_relevel(l1, "es", "en", "ma")) %>%
#     ggplot(., aes(x = time_zero, y = target_prop, fill = l1)) +
#     facet_grid(. ~ cond, labeller = as_labeller(condition_names)) +
#     geom_hline(yintercept = 0.5, color = 'white', size = 3) +
#     geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
#     geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
#     stat_summary(fun.y = "mean", geom = "line", size = 1) +  
#     stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
#                  stroke = 0.5, pch = 21) +
#     scale_fill_brewer(palette = 'Set1', name = "L1",
#                        labels = c("ES", "EN", "MA")) +
#     scale_x_continuous(breaks = c(-10, 0, 10, 20),
#                        labels = c("-500", "0", "500", "1000")) +
#     labs(y = 'Proportion of target fixations',
#          x = 'Time relative to target syllable offset (ms)',
#          caption = "Mean +/- 95% CI") +
#     annotate("text", x = 3.3, y = 0.02, label = '200ms',
#              angle = 90, size = 3, hjust = 0) +
#     theme_grey(base_size = 12, base_family = "Times")
# 
# l1_names <- c(
#   `es` = 'Spanish speakers',
#   `en` = 'English speakers',
#   `ma` = 'Mandarin speakers'
# )
# 
# stress_p2 <- stress50 %>%        
#   #na.omit(.) %>%
#   filter(., time_zero >= -10, time_zero <= 20) %>%
#   mutate(., l1 = fct_relevel(l1, "es", "en", "ma")) %>%
#   ggplot(., aes(x = time_zero, y = target_prop, fill = cond)) +
#   facet_grid(. ~ l1, labeller = as_labeller(l1_names)) +
#   geom_hline(yintercept = 0.5, color = 'white', size = 3) +
#   geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
#   geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
#   stat_summary(fun.y = "mean", geom = "line", size = 1) +  
#   stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
#                stroke = 0.5, pch = 21) +
#   scale_fill_brewer(palette = 'Set1', name = "Tense",
#                     labels = c("Present", "Preterit")) +
#   scale_x_continuous(breaks = c(-10, 0, 10, 20),
#                      labels = c("-500", "0", "500", "1000")) +
#   labs(y = 'Proportion of target fixations',
#        x = 'Time relative to target syllable offset (ms)',
#        caption = "Mean +/- 95% CI") +
#   annotate("text", x = 3.3, y = 0.02, label = '200ms',
#            angle = 90, size = 3, hjust = 0) +
#   theme_grey(base_size = 12, base_family = "Times")
# 
# 
# ggsave('stress_l1.png',
#        plot = stress_p1, dpi = 600, device = "png",
#        path = figs_path,
#        height = 3.5, width = 8.5, units = 'in')
# ggsave('stress_tense.png',
#        plot = stress_p1, dpi = 600, device = "png",
#        path = figs_path,
#        height = 3.5, width = 8.5, units = 'in')
# 
# 
# 
# # -----------------------------------------------------------------------------




# Plot GCA nested-interactions --------------------------------------------------------------------


# L2 speakers

# All three variables 
base_l2_ints <- fits_all_l2_ints %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         # Stress = if_else(condition_sum == 1, "Present", "Preterite"),
         # Stress = fct_relevel(Stress, 'Present'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish proficiency`)) + #fill = Stress, color = Stress, 
  facet_grid(`Spanish use` ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  # scale_color_brewer(palette = "Set1", name = "Condition") +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points")) #cm

library(grid)
library(gridExtra)

stress_l2_ints <- grid.arrange(base_l2_ints,     #en_gca_plot <- 
                          bottom = textGrob('Spanish use', rot = 270,
                                            x = 0.98, y = 4.3, gp = gpar(fontsize = 9,
                                                                      fontfamily = 'Times')))

ggsave(paste0(figs_path, "/stress_l2_ints.png"), stress_l2_ints, width = 180,
       height = 120, units = "mm", dpi = 600)


# Proficiency and use
prof_plot <- fits_all_l2_ints %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         # Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
         # Stress = fct_relevel(Stress, 'Paroxytone')
         ) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish proficiency`)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom')

use_plot <- fits_all_l2_ints %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         # Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
         # Stress = fct_relevel(Stress, 'Paroxytone')
  ) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish use`)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom')

ggsave(paste0(figs_path, "/prof_ints.png"), prof_plot, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/use_ints.png"), use_plot, width = 180,
       height = 120, units = "mm", dpi = 600)


# 2-way interactions
profxl1 <- fits_all_l2_ints %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish proficiency`)) + 
  facet_grid(. ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"))



usexl1 <- fits_all_l2_ints %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish use`)) + 
  facet_grid(. ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"))


delexuse <- fits_all_l2_ints %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish proficiency`)) + 
  facet_grid(. ~ `Spanish use`) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  ggtitle("Spanish use") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"),
        plot.title = element_text(hjust = 0.5, size = 8))


ggsave(paste0(figs_path, "/prof_l1_ints.png"), profxl1, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/use__l1_ints.png"), usexl1, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/prof_use_ints.png"), delexuse, width = 200,
       height = 80, units = "mm", dpi = 600)









# Plot GCA 3-way interaction straight away --------------------------------------------------------------------

# Monolingual Spanish speakers

stress_mon <- fits_all_mon %>%
  # mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
  #        condition = fct_relevel(condition, "Present"))%>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin)) + #, fill = condition, color = condition
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) +
  geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  # stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
  #              alpha = 0.5) +
  # geom_point(aes(color = condition), size = 1.3, show.legend = F) +
  # geom_point(aes(color = condition), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_big + legend_adj #+ labs(color = "Condition")


# L2 speakers

base_l2 <- fits_all_l2 %>%
  mutate(Proficiency = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         # Stress = if_else(condition_sum == 1, "Present", "Preterite"),
         # Stress = fct_relevel(Stress, 'Present'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = Proficiency)) + #fill = Stress, color = Stress, 
  facet_grid(`Spanish use` ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  # scale_color_brewer(palette = "Set1", name = "Condition") +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points")) #cm


stress_l2 <- grid.arrange(base_l2,     #en_gca_plot <- 
                      bottom = textGrob('Spanish use', rot = 270,
                                        x = .96, y = 4, gp = gpar(fontsize = 9,
                                                                    fontfamily = 'Times')))


ggsave(paste0(figs_path, "/stress_mon.png"), stress_mon, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_l2.png"), stress_l2, width = 180,
       height = 120, units = "mm", dpi = 600)


# --------------------------------------------------------------------





# All three variables 
base_l2_l1ot2 <- fits_all_l2_l1ot2 %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         # Stress = if_else(condition_sum == 1, "Present", "Preterite"),
         # Stress = fct_relevel(Stress, 'Present'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish proficiency`)) + #fill = Stress, color = Stress, 
  facet_grid(`Spanish use` ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  # scale_color_brewer(palette = "Set1", name = "Condition") +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points")) #cm



stress_l2_l1ot2 <- grid.arrange(base_l2_l1ot2,  
                               bottom = textGrob('Spanish use', rot = 270, x = 0.98, y = 4.3, 
                                                 gp = gpar(fontsize = 9, fontfamily = 'Times')))

ggsave(paste0(figs_path, "/stress_l2_l1ot2.png"), stress_l2_l1ot2, width = 180,
       height = 120, units = "mm", dpi = 600)


# Proficiency and use
prof_plot <- fits_all_l2_l1ot2 %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         # Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
         # Stress = fct_relevel(Stress, 'Paroxytone')
  ) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish proficiency`)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom')

use_plot <- fits_all_l2_ints %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         # Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
         # Stress = fct_relevel(Stress, 'Paroxytone')
  ) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish use`)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom')

ggsave(paste0(figs_path, "/prof_l1ot2.png"), prof_plot, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/use_l1ot2.png"), use_plot, width = 180,
       height = 120, units = "mm", dpi = 600)


# 2-way interactions
profxl1 <- fits_all_l2_l1ot2 %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish proficiency`)) + 
  facet_grid(. ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"))



usexl1 <- fits_all_l2_l1ot2 %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish use`)) + 
  facet_grid(. ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"))


delexuse <- fits_all_l2_l1ot2 %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish proficiency`)) + 
  facet_grid(. ~ `Spanish use`) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  ggtitle("Spanish use") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"),
        plot.title = element_text(hjust = 0.5, size = 8))


ggsave(paste0(figs_path, "/prof_l1_l1ot2.png"), profxl1, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/use__l1_l1ot2.png"), usexl1, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/prof_use_l1ot2.png"), delexuse, width = 200,
       height = 80, units = "mm", dpi = 600)





# Switch prof and use in 3-way plot
fits_all_l2_l1ot2 %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         `Spanish use` = as.factor(use_z),
         # Stress = if_else(condition_sum == 1, "Present", "Preterite"),
         # Stress = fct_relevel(Stress, 'Present'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish use`)) + #fill = Stress, color = Stress, 
  facet_grid(`Spanish proficiency` ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  # scale_color_brewer(palette = "Set1", name = "Condition") +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to final syllable onset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = 'bottom', #c(0.1, 0.9) 
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"))






# overlay model on empirical data
fits_all_l2_overlay <- fits_all_l2_l1ot2

fits_all_l2_overlay$target_prop <- exp(fits_all_l2_overlay$fit)
fits_all_l2_overlay$l1 <- fits_all_l2_overlay$l1_sum

fits_all_l2_overlay %<>% mutate(l1 = if_else(l1 == -1, "English", "Chinese"),
                                l1 = fct_relevel(l1, "English", "Chinese"))


gca <- fits_all_l2_overlay %>%
  # mutate(`Spanish proficiency` = as.factor(DELE_z),
  #        `Spanish use` = as.factor(use_z)) %>%
  ggplot(., aes(x = time_zero, y = target_prop, #ymax = ymax, ymin = ymin,
                lty = l1)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = .5) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                     labels = c("-200", "-100", "0", "100", "200")) +
  labs(x = "Time (ms) relative to verb's final syllable onset",
       y = "Exponential value of\nfixation probability",
       lty = "L1") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom')





tc <- stress50 %>%
  filter(time_zero > -5 & time_zero < 5, l1 != 'es') %>%
  ggplot(., aes(x = time_zero, y = target_count, lty = l1)) +
  geom_vline(xintercept = 4, lty = 3, color = 'grey40') +
  # geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  scale_color_discrete(name="L1",
                       breaks = c('en', 'ma'),
                       labels = c("English", 'Chinese')) +
  scale_x_continuous(breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8),
                     labels = c("-400", "-300", "-200", "-100", "0", "100", "200", "300", "400")) +
  # ggtitle("Time course per verbal tense") +
  # xlab("Time in 50 ms bins (0 = marker time before accounting for 200 ms processing)") +
  labs(y = "Count of fixations\non target", xiintercept = 1.5) +
  # annotate("text", x = 3.65, y = 0.53, label = '200ms',
  #                       angle = 90, size = 3, hjust = 0) +
  theme_grey(base_size = 10, base_family = "Times") +
  theme(legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank())


grid.newpage()
grid.draw(rbind(ggplotGrob(tc), ggplotGrob(gca), size = "last"))






stress50 %>%
  filter(time_zero > -5 & time_zero < 5, l1 != 'es') %>%
  ggplot(., aes(x = time_zero, y = target_count, color = l1)) +
  geom_vline(xintercept = 4, lty = 3) +
  #geom_hline(yintercept = 0.5, lty = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  geom_smooth(method = "loess", se = FALSE) 
  # geom_line(data = fits_all_l2_overlay, 
  #           aes(stat_summary(fun.y = "mean", geom = "line", size = .5)))
  scale_color_discrete(name="L1",
                       breaks = c('en', 'ma'),
                       labels = c("English", 'Chinese')) +
  scale_x_continuous(breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8),
                     labels = c("-400", "-300", "-200", "-100", "0", "100", "200", "300", "400")) +
  # ggtitle("Time course per verbal tense") +
  # xlab("Time in 50 ms bins (0 = marker time before accounting for 200 ms processing)") +
  labs(y = "Proportion of fixations on target") +
  # annotate("text", x = 3.65, y = 0.53, label = '200ms',
  #                       angle = 90, size = 3, hjust = 0) +
  theme_grey(base_size = 10, base_family = "Times") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())


