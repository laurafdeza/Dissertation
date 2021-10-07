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

# Load models 
load(paste0(gca_mods_path, "/gca_mon_mods_4equal.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_dele_4equal.Rdata")) 
load(paste0(gca_mods_path, "/gca_l2_use_4equal.Rdata"))
load(paste0(gca_mods_path, "/model_preds_4equal.Rdata"))


# Set path for saving figs
figs_path <- here("figs", "stress", "gca", "LL_changes")

# -----------------------------------------------------------------------------






# Plot GCA nested-interactions --------------------------------------------------------------------


# L2 speakers

# All three variables 
dele_plot <- model_preds_4equal$fits_all_dele %>%
  mutate(`Spanish proficiency` = as.factor(DELE_z),
         #`Spanish use` = as.factor(use_z),
         Stress = if_else(condition_sum == 1, "Oxytone", 'Paroxytone'),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish proficiency`, fill = Stress, color = Stress)) +
  facet_grid(. ~ L1) +
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

plogis(-1)
#0.2689414
plogis(1)
#0.7310586
plogis(2)
#0.8807971
plogis(3)
#0.9525741


ggsave(paste0(figs_path, "/dele_plot_4equal.png"), dele_plot, width = 180,
       height = 120, units = "mm", dpi = 600)


use_plot <- model_preds_4equal$fits_all_use %>%
  mutate(`Spanish use` = as.factor(use_z),
         #`Spanish use` = as.factor(use_z),
         Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish use`, fill = Stress, color = Stress)) +
  facet_grid(. ~ L1) +
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

ggsave(paste0(figs_path, "/use_plot_4equal.png"), use_plot, width = 180,
       height = 120, units = "mm", dpi = 600)


model_preds_4equal$fits_all_use %>%
  mutate(`Spanish use` = as.factor(use_z),
         #`Spanish use` = as.factor(use_z),
         Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         L1 = if_else(l1_sum == 1, "Mandarin Chinese", "English"),
         L1 = fct_relevel(L1, "English")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Spanish use`, fill = L1, color = L1)) +
  facet_grid(. ~ Stress) +
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



# Monolingual Spanish speakers

stress_mon <- model_preds_4equal$fits_all_mon %>%
  mutate(Stress = if_else(condition_sum == 1, "Oxytone", "Paroxytone"),
         Stress = fct_relevel(Stress, "Paroxytone"))%>%
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

ggsave(paste0(figs_path, "/mon_plot_4equal.png"), stress_mon, width = 180,
       height = 120, units = "mm", dpi = 600)


