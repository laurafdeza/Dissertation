# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "02_load_data.R"))

# Get path to saved models
gca_mods_path  <- here("mods", "vision", "gca")

# Load models as list and store full mod to global env
load(paste0(gca_mods_path, "/full_mods.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))

gca_final_model <- full_mods$gca_full_mod_int_1

# Set path for saving figs
figs_path <- here("figs", "vision", "gca")

# -----------------------------------------------------------------------------






# Prep data ---------------------------------------------------------------


# # For raw data plot, need to add non linx data to 1st df (car_dev + corsi)
# vision <- read_csv("./data/clean/vision_scores.csv")
# corsi <- read_csv("./data/clean/corsi_z_scores.csv")
# 
# visuospatial_df <- left_join(x = vision, y = corsi, by = "participant", all.x=TRUE)
# 
# vision50 <- left_join(x = stress50, y = visuospatial_df, by = "participant", all.x=TRUE)
# 
# vision50 <- na.omit(vision50)

# -----------------------------------------------------------------------------




# Plot GCA --------------------------------------------------------------------

# Present group differences in movement anticipation
car_gca <- model_preds$fits_all %>%
  mutate(stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         stress = fct_relevel(stress, "Paroxytone"),
         group = fct_relevel(group, "SS", "AE", "AM", "IE", "IM")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin, color = car_dev)) +
  facet_grid(group ~ stress) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.75) +
  geom_point(aes(color = car_dev), color = "black", size = 1.3, show.legend = F) +
  geom_point(aes(color = car_dev), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big + #legend_adj
  ggtitle('GCA considering visuospatial anticipation performance')

# visuospatial WM
corsi_gca <- model_preds$fits_all %>%
  mutate(stress = if_else(stress_sum == 1, "Paroxytone", "Oxytone"),
         stress = fct_relevel(stress, "Paroxytone"),
         group = fct_relevel(group, "SS", "AE", "AM", "IE", "IM")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin, color = corsi)) + 
  facet_grid(group ~ stress) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.75) +
  geom_point(aes(color = corsi), color = "black", size = 1.3, show.legend = F) +
  geom_point(aes(color = corsi), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big + #legend_adj 
  ggtitle("GCA considering visuospatial WM")




ggsave(paste0(figs_path, "/car_gca.png"), car_gca, width = 150,
       height = 120, units = "mm", dpi = 600)
# ggsave(paste0(figs_path, "/stress_p2.eps"), stress_p2, width = 150,
#        height = 120, units = "mm", dpi = 600, device = cairo_ps)
ggsave(paste0(figs_path, "/corsi_gca.png"), corsi_gca, width = 150,
       height = 120, units = "mm", dpi = 600)
# ggsave(paste0(figs_path, "/stress_p3.eps"), stress_p3, width = 150,
#        height = 120, units = "mm", dpi = 600, device = cairo_ps)

# -----------------------------------------------------------------------------








# Effect size plot ------------------------------------------------------------

bind_rows(
  data.comp %>%
    group_by(participant, group, time_zero, coda) %>%
    summarize(competition = mean(eLog)) %>%
    spread(coda, competition) %>%
    mutate(Effect = `1` - `0`,
           condition_type = "Syllable structure",
           fit_type = "raw") %>%
    select(-`0`, -`1`),
  data.comp %>%
    group_by(participant, group, time_zero, condition) %>%
    summarize(competition = mean(eLog)) %>%
    spread(condition, competition) %>%
    mutate(Effect = unstressed - stressed,
           condition_type = "Lexical stress",
           fit_type = "raw") %>%
    select(-unstressed, -stressed),
  data.comp %>%
    group_by(participant, group, time_zero, coda) %>%
    summarize(competition = mean(GCA_Full)) %>%
    spread(coda, competition) %>%
    mutate(Effect = `1` - `0`,
           condition_type = "Syllable structure",
           fit_type = "model") %>%
    select(-`0`, -`1`),
  data.comp %>%
    group_by(participant, group, time_zero, condition) %>%
    summarize(competition = mean(GCA_Full)) %>%
    spread(condition, competition) %>%
    mutate(Effect = unstressed - stressed,
           condition_type = "Lexical stress",
           fit_type = "model") %>%
    select(-unstressed, -stressed)) %>%
  spread(fit_type, Effect) %>%
  #filter(time_zero < 10) %>%
  ggplot(., aes(x = time_zero, y = raw, color = group)) +
    facet_grid(. ~ condition_type) +
    geom_hline(yintercept = 0, color = "white", size = 3) +
    geom_vline(xintercept = 4, color = "white", size = 3) +
    stat_summary(fun.data = mean_se, geom = "pointrange") +
    #stat_summary(aes(y = model, color = group), fun.y = mean,
    #             geom = 'line', size = 1.4) +
    geom_smooth(method = 'gam', formula = y ~ poly(x, 3), se = F,
                show.legend = FALSE) +
    labs(x = "Time (ms) relative to target syllable offset",
         y = "Effect", caption = "Mean +/- 95% CI") +
    scale_color_brewer(palette = "Set1", name = "", guide = 'legend',
                       labels = c("SS", "LA", "INT")) +
    scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                       labels = c("-200", "0", "200", "400", "600")) +
    theme_grey(base_size = 15, base_family = "Times New Roman")

# -----------------------------------------------------------------------------











