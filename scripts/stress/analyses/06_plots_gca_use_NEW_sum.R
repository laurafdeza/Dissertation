# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
source(here::here("scripts", "01_helpers.R"))
source(here::here("scripts", "02_load_data.R"))

# Get path to saved models 
gca_mods_path <- here("mods", "stress", "gca") 

# Load models as list and store full mod to global env
# load(paste0(gca_mods_path, "/gca_mon_mods_sum.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_mods_sum.Rdata"))
# load(paste0(gca_mods_path, "/model_preds.Rdata"))
load(paste0(gca_mods_path, "/model_preds_l2_sum.Rdata"))
# list2env(gca_mon_mods, globalenv())
# list2env(gca_l2_mods, globalenv())
list2env(model_preds, globalenv())
# list2env(model_preds_l2, globalenv())

# Set path for saving figs
figs_path <- here("figs", "stress", "gca")

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
# 
# 
# 
# 
# # Plot GCA --------------------------------------------------------------------
# 
# # Monolingual Spanish speakers
# 
# stress_mon <- model_preds$fits_all_mon %>%
#   mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
#          condition = fct_relevel(condition, "Present"))%>%
#   ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
#                 fill = condition)) + #, color = condition
#   geom_hline(yintercept = 0, lty = 3, size = 0.4) +
#   geom_vline(xintercept = 4, lty = 3, size = 0.4) +
#   stat_summary(fun.y = "mean", geom = "line", size = 1) + 
#   # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
#   stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
#                alpha = 0.5) +
#   geom_point(aes(color = condition), size = 1.3, show.legend = F) +
#   geom_point(aes(color = condition), size = 0.85, show.legend = F) +
#   scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
#                      labels = c("-200", "0", "200", "400", "600")) +
#   labs(x = "Time (ms) relative to target syllable offset",
#        y = "Empirical logit of looks to target") +
#   theme_big + legend_adj #+ labs(color = "Condition")


# L2 speakers

# Proficiency
prof_sum <- fits_all_l2_dele %>%
  mutate(Proficiency = as.factor(DELE_z),
         Stress = if_else(condition_sum == 1, "Paroxytone/Present", "Oxytone/Preterite"),
         Stress = fct_relevel(Stress, 'Paroxytone/Present'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = Stress, color = Stress, lty = Proficiency)) +
  facet_grid(. ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "Condition") +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = c(0.1, 0.70))

prof_sum_cond <- fits_all_l2_dele %>%
  mutate(Proficiency = as.factor(DELE_z),
         Stress = if_else(condition_sum == 1, "Present", "Preterite"),
         Stress = fct_relevel(Stress, 'Present'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = L1, color = L1, lty = Proficiency)) +
  facet_grid(. ~ Stress) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "L1") +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = c(0.12, 0.70))


# L2 use
use_sum <- fits_all_l2_use %>%
  mutate(`L2 use` = as.factor(use_z),
         Stress = if_else(condition_sum == 1, "Paroxytone/Present", "Oxytone/Preterite"),
         Stress = fct_relevel(Stress, 'Paroxytone/Present'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = Stress, color = Stress, lty = `L2 use`)) +
  facet_grid(. ~ L1) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "Condition") +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = c(0.1, 0.70))

use_sum_cond <- fits_all_l2_use %>%
  mutate(`L2 use` = as.factor(use_z),
         Stress = if_else(condition_sum == 1, "Present", "Preterite"),
         Stress = fct_relevel(Stress, 'Present'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = L1, color = L1, lty = `L2 use`)) +
  facet_grid(. ~ Stress) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "L1") +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_2 +
  theme(legend.position = c(0.12, 0.70))

ggsave(paste0(figs_path, "/use_sum.png"), use_sum, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/use_sum_cond.png"), use_sum_cond, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/prof_sum.png"), prof_sum, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/prof_sum_cond.png"), prof_sum_cond, width = 180,
       height = 120, units = "mm", dpi = 600)
