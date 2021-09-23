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
source(here::here("scripts", "02_load_data.R"))

# Get path to saved models
gca_mods_path  <- here("mods", "vision", "gca", "continuous")

# Get path to saved models 
load(paste0(gca_mods_path, "/mon_mods.Rdata"))
load(paste0(gca_mods_path, "/en_mods.Rdata"))
load(paste0(gca_mods_path, "/ma_mods.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))


# Set path for saving figures
figs_path <- here("figs", "vision", "gca")

# -----------------------------------------------------------------------------






# Plot raw data ---------------------------------------------------------------

# Merge visuospatial data into main dataframe
vision <- read_csv("./data/clean/vision_scores.csv")
corsi <- read_csv("./data/clean/corsi_z_scores.csv")

visuospatial_df <- left_join(x = vision, y = corsi, by = "participant", all.x=TRUE)

vision50 <- left_join(x = stress50, y = visuospatial_df, by = "participant", all.x=TRUE)

vision50 <- na.omit(vision50)




vision50 <- vision50 %>%
  filter(., time_zero >= -6 & time_zero <= 8) %>%
  filter(., participant != 'ies04' & participant != 'ies17' & 
           participant != 'ies28' & participant != 'aes32') %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma"),
         stress_sum = if_else(cond == "1", 1, -1)
         # ,
         # car_dev_z = (car_dev - mean(car_dev)) / sd(car_dev)
         ) %>%           # 1 = present/paroxytone, -1/2 = preterit/oxytone        
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")


vision50$cond <- as.factor(as.character(vision50$cond))


condition_names <- c(
  `1` = 'Paroxytone\n(CANta)',
  `2` = 'Oxytone\n(canTÓ)'
)


stress_p1 <- vision50 %>%
    #na.omit(.) %>%
    filter(., time_zero >= -10, time_zero <= 8) %>%
    mutate(., l1 = fct_relevel(l1, "es", "en", "ma")) %>%
    ggplot(., aes(x = time_zero, y = target_prop, fill = l1)) +
    facet_grid(cond ~ ., labeller = as_labeller(condition_names)) +
    geom_hline(yintercept = 0.5, color = 'grey40', lty = 3) +
    geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
    geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
    stat_summary(fun.y = "mean", geom = "line", size = 1) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
                 stroke = 0.5, pch = 21) +
    # scale_fill_brewer(palette = 'Greys', name = "L1",
    #                    labels = c("Spanish", "English", "Mandarin Chinese")) +
    scale_fill_manual(values = c('#000000', '#736F6E', '#FFFFFF'), name = "Tense",
                      labels = c("Spanish", "English", "Mandarin Chinese")) +
    scale_x_continuous(breaks = c(-10, 0, 10),
                       labels = c("-500", "0", "500")) +
    labs(y = 'Proportion of target fixations',
         x = 'Time relative to target syllable offset (ms)',
         caption = "Mean +/- 95% CI") +
    annotate("text", x = 3.3, y = 0.02, label = '200ms',
             angle = 90, size = 3, hjust = 0, family = 'Times') +
    theme_grey(base_size = 12, base_family = "Times") +
    theme(legend.position = 'bottom')

l1_names <- c(
  `es` = 'Spanish\nspeakers',
  `en` = 'English\nspeakers',
  `ma` = 'Mandarin\nChinese\nspeakers'
)

stress_p2 <- vision50 %>%
  #na.omit(.) %>%
  filter(., time_zero >= -6, time_zero <= 8) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma")) %>%
  ggplot(., aes(x = time_zero, y = target_prop, fill = cond)) +
  facet_grid(l1 ~ ., labeller = as_labeller(l1_names)) +
  geom_hline(yintercept = 0.5, color = 'grey40', lty = 3) +
  geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
  geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
  stat_summary(fun.y = "mean", geom = "line", size = 1, alpha = .7) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
               stroke = 0.5, pch = 21, alpha = .8) +
  # scale_fill_brewer(palette = 'Greys', name = "Tense",
  #                   labels = c("Present", "Preterite")) +
  scale_fill_manual(values = c('#000000', '#736F6E'), name = "Tense",
                    labels = c("Paroxytone", "Oxytone")) +
  scale_x_continuous(breaks = c(-10, 0, 10),
                     labels = c("-500", "0", "500")) +
  labs(y = 'Proportion of target fixations',
       x = 'Time relative to target syllable offset (ms)',
       caption = "Mean +/- 95% CI") +
  annotate("text", x = 3.3, y = 0.02, label = '200ms',
           angle = 90, size = 3, hjust = 0, family = 'Times') +
  theme_grey(base_size = 12, base_family = "Times") +
  theme(legend.position = 'bottom')


ggsave('timecourse_l1.png',
       plot = stress_p1, dpi = 600, device = "png",
       path = figs_path,
       height = 3.5, width = 4.5, units = 'in')
ggsave('timecourse_tense.png',
       plot = stress_p2, dpi = 600, device = "png",
       path = figs_path,
       height = 3.5, width = 4.5, units = 'in')



# -----------------------------------------------------------------------------




# Plot GCA --------------------------------------------------------------------

# Spanish speakers

# Car model
mon_car <- model_preds$fits_mon_car %>%
  mutate(`Visuospatial prediction timing` = as.factor(car_dev),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone')) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin, 
                lty = `Visuospatial prediction timing`)) +
  facet_grid(. ~ Stress) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       lty = 'Stress condition') +
  ggtitle('Stress pattern') +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))


# Corsi model
mon_corsi <- model_preds$fits_mon_corsi %>%
  mutate(`Visuospatial WM score` = as.factor(corsi),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone')) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin, 
                lty = `Visuospatial WM score`)) +
  facet_grid(. ~ Stress) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       lty = 'Stress condition') +
  ggtitle('Stress pattern') +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))


# English speakers

# Car
plot1 <- model_preds$fits_en_car %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial prediction timing` = as.factor(car_dev)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Visuospatial prediction timing`)) +
  facet_grid(Proficiency ~ Stress) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       lty = 'Visuospatial prediction timing') +
  ggtitle('Stress condition') +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
        legend.position = 'bottom',
        plot.title = element_text(size = 9, hjust = 0.5),
        plot.margin = margin(t = 5, l = 5, r = 24))

plot2 <- grid.arrange(plot1,     
             bottom = textGrob('Proficiency', rot = 270,
                               x = .97, y = 4.3, gp = gpar(fontsize = 9,
                                                           fontfamily = 'Times')))

## Proficiency
en_prof <- model_preds$fits_en_car %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial prediction timing` = as.factor(car_dev)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = Proficiency)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = 0.8) +
  # stat_summary(fun.data = mean_cl_boot, geom = 'ribbon', size = 0.5,
  #              pch = 21, alpha = 0.5) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))

## Stress
en_stress <- model_preds$fits_en_car %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial prediction timing` = as.factor(car_dev)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = Stress)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = 0.8) +
  # stat_summary(fun.data = mean_cl_boot, geom = 'ribbon', size = 0.5,
  #              pch = 21, alpha = 0.5) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))

## Visuospatial prediction
en_car <- model_preds$fits_en_car %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial prediction timing` = as.factor(car_dev)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Visuospatial prediction timing`)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = 0.8) +
  # stat_summary(fun.data = mean_cl_boot, geom = 'ribbon', size = 0.5,
  #              pch = 21, alpha = 0.5) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", 'dotted', 'dashed')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))


# Corsi
plot3 <- model_preds$fits_en_corsi %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial WM score` = as.factor(corsi)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Visuospatial WM score`)) +
  facet_grid(Proficiency ~ Stress) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  ggtitle('Stress condition') +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))

plot4 <- grid.arrange(plot3,     
                      bottom = textGrob('Proficiency', rot = 270,
                                        x = .97, y = 4.3, gp = gpar(fontsize = 9,
                                                                     fontfamily = 'Times')))


# Mandarin Chinese speakers

# Car
plot5 <- model_preds$fits_ma_car %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial prediction timing` = as.factor(car_dev)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Visuospatial prediction timing`)) +
  facet_grid(Proficiency ~ Stress) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       lty = 'Visuospatial prediction timing') +
  ggtitle('Stress condition') +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))

plot6 <- grid.arrange(plot5,     
                      bottom = textGrob('Proficiency', rot = 270,
                                        x = .97, y = 4.3, gp = gpar(fontsize = 9,
                                                                     fontfamily = 'Times')))



## Proficiency
ma_prof <- model_preds$fits_ma_car %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial prediction timing` = as.factor(car_dev)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = Proficiency)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = 0.8) +
  # stat_summary(fun.data = mean_cl_boot, geom = 'ribbon', size = 0.5,
  #              pch = 21, alpha = 0.5) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))

## Stress
ma_stress <- model_preds$fits_ma_car %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial prediction timing` = as.factor(car_dev)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = Stress)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = 0.8) +
  # stat_summary(fun.data = mean_cl_boot, geom = 'ribbon', size = 0.5,
  #              pch = 21, alpha = 0.5) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))

## Visuospatial prediction
ma_car <- model_preds$fits_ma_car %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial prediction timing` = as.factor(car_dev)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Visuospatial prediction timing`)) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  # geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", size = 0.8) +
  # stat_summary(fun.data = mean_cl_boot, geom = 'ribbon', size = 0.5,
  #              pch = 21, alpha = 0.5) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", 'dotted', 'dashed')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))


# Corsi
plot7 <- model_preds$fits_ma_corsi %>%
  mutate(Proficiency = as.factor(prof_std),
         Proficiency = fct_relevel(Proficiency, "1", "0", "-1"),
         Stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         Stress = fct_relevel(Stress, 'Paroxytone'),
         `Visuospatial WM score` = as.factor(corsi)) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                lty = `Visuospatial WM score`)) +
  facet_grid(Proficiency ~ Stress) +
  geom_hline(yintercept = 0, size = .5, color = "grey40", linetype = 'dotted') +
  geom_vline(xintercept = 4, size = .5, color = "grey40", linetype = 'dotted') +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_linetype_manual(values=c("solid", "dashed", 'dotted')) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  ggtitle('Stress condition') +
  theme_grey(base_size = 10, base_family = "Times") + 
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(t = 5, l = 5, r = 24))

plot8 <- grid.arrange(plot7,     
                      bottom = textGrob('Proficiency', rot = 270,
                                        x = .97, y = 4.3, gp = gpar(fontsize = 9,
                                                                     fontfamily = 'Times')))

# Save plots
ggsave(paste0(figs_path, "/mon_gca_car.png"), mon_car, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/mon_gca_corsi.png"), mon_corsi, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/en_gca_car.png"), plot2, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/en_gca_corsi.png"), plot4, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/en_car.png"), en_car, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/en_prof.png"), en_prof, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/en_stress.png"), en_stress, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/ma_gca_car.png"), plot6, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/ma_gca_corsi.png"), plot8, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/ma_car.png"), ma_car, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/ma_prof.png"), ma_prof, width = 180,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/ma_stress.png"), ma_stress, width = 180,
       height = 120, units = "mm", dpi = 600)
