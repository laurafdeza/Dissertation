# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
library(grid)
library(gridExtra)

source(here::here("scripts", "01_helpers.R"))
# source(here::here("scripts", "02_load_data.R"))

# Get path to saved models
gca_mods_path  <- here("mods", "vision", "gca", "cont_speed_verb")

# Load models as lists
load(paste0(gca_mods_path, "/mon_mods_psequal.Rdata"))
load(paste0(gca_mods_path, "/model_preds_psequal.Rdata"))


# Set path for saving figures
figs_path <- here("figs", "vision", "gca", 'cont_speed_verb')

# -----------------------------------------------------------------------------






# Plot raw data ---------------------------------------------------------------

# Merge visuospatial data into main dataframe
wm <- read_csv("./data/clean/wm_processing_speed.csv")
wm_score <- read_csv("./data/clean/ospan_set_z_scores.csv")

vision50 <- left_join(x = stress50, y = wm, by = "participant", all.x=TRUE)
vision50 <- left_join(x = vision50, y = wm_score, by = "participant", all.x=TRUE)






vision50 <- vision50 %>%
  filter(., time_zero >= -5 & time_zero <= 5) %>%
  mutate(., #l1 = fct_relevel(l1, "es", "en", "ma"),
         stress_sum = if_else(cond == "1", -1, 1)) %>%           # 1 = present, 2 = preterit        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

vision50 <- filter(vision50) %>% select(-DELE, -percent_l2_week, 
                                                      -prof, -group)

vision50 <- na.omit(vision50)

plot_data <- vision50

plot_data$cond <- as.factor(as.character(plot_data$cond))
plot_data$l1 <- as.factor(as.character(plot_data$l1))

condition_names <- c(
  `1` = 'Paroxytone\n(CANta)',
  `2` = 'Oxytone\n(canTÓ)'
)


timecourse <- plot_data %>%
    #na.omit(.) %>%
    filter(., time_zero >= -4, time_zero <= 4) %>%
    mutate(., l1 = fct_relevel(l1, "es", "en", "ma")) %>%
    ggplot(., aes(x = time_zero, y = target_prop, fill = l1)) +
    facet_grid(cond ~ ., labeller = as_labeller(condition_names)) +
    geom_hline(yintercept = 0.5, color = 'grey40', lty = 3) +
    geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
    geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
    stat_summary(fun.y = "mean", geom = "line", size = .5) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
                 stroke = 0.5, pch = 21) +
    scale_fill_discrete(name = "L1",
                      labels = c("Spanish", "English", "Mandarin Chinese")) +
    scale_x_continuous(breaks = c(-4, -2, 0, 2, 4),
                       labels = c("-200", "-100", "0", "100", "200")) +
    labs(y = 'Proportion of target fixations',
         x = 'Time relative to target syllable offset (ms)',
         caption = "Mean +/- 95% CI") +
    annotate("text", x = 3.8, y = 0.2, label = '200ms',
             angle = 90, size = 3, hjust = 0, family = 'Times') +
    theme_grey(base_size = 12, base_family = "Times") +
    theme(legend.position = 'bottom')


ggsave(paste0(figs_path, '/timecourse_michele.png'),
       plot = timecourse,
       height = 3.5, width = 4.5, units = 'in')



# -----------------------------------------------------------------------------




# Plot GCA --------------------------------------------------------------------

# Monolinguals
mon_michele <- model_preds_michele$fits_mon_wm %>%
  mutate(`Verbal WM` = as.factor(ospan),
         Stress = if_else(stress_sum == -1, "Paroxytone/present", "Oxytone/preterite"),
         Stress = fct_relevel(Stress, 'Paroxytone/present')) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin, 
                lty = `Verbal WM`)) +
  facet_grid(Stress ~ .) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.1, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       lty = 'Verbal WM') +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))


ggsave(paste0(figs_path, "/mon_michele.png"), mon_michele, width = 180,
       height = 120, units = "mm", dpi = 600)


# English speakers
en_michele <- model_preds_michele$fits_en_wm %>%
  mutate(`Verbal WM` = as.factor(ospan),
         Proficiency = as.factor(DELE_z),
         Stress = if_else(stress_sum == -1, "Paroxytone/present", "Oxytone/preterite"),
         Stress = fct_relevel(Stress, 'Paroxytone/present')) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin, 
                lty = `Verbal WM`)) +
  facet_grid(Stress ~ Proficiency) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.1, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       lty = 'Verbal WM') +
  ggtitle("Spanish proficiency") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))


ggsave(paste0(figs_path, "/en_michele.png"), en_michele, width = 180,
       height = 120, units = "mm", dpi = 600)
