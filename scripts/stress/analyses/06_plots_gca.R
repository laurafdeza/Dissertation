# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "02_load_data.R"))

# Get path to saved models 
gca_mods_path  <- here("reports", "mods", "gca")      # a level might be missing

# Load models as list and store full mod to global env
load(paste0(gca_mods_path, "/full_mods.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))

gca_full_mod_int_3 <- full_mods$gca_full_mod_int_3

# Set path for saving figs
figs_path <- here("figs", "stress", "eye_track")

# -----------------------------------------------------------------------------






# Plot raw data ---------------------------------------------------------------

df_stress_50 <- stress50 %>%
  filter(.,  group %in% c('mon', 'aes', 'ies', 'ams', 'ims'),
            !participant %in% c("", "", "", "", ""))

condition_names <- c(
  `stressed` = 'Paroxytone',
  `unstressed` = 'Oxytone',
#  `0` = 'CV',
#  `1` = 'CVC'
)

stress_p1 <- df_stress_50 %>%        # settings were coda so some more things might need changing
    na.omit(.) %>%
    filter(., time_zero >= -10, time_zero <= 20) %>%
    mutate(., group = fct_relevel(group, "mon", "aes", "ies", "ams", "ims")) %>%
    ggplot(., aes(x = time_zero, y = targetProp, fill = group, shape = group)) +
#    facet_grid(condition ~ coda, labeller = as_labeller(condition_names)) +
    geom_hline(yintercept = 0.5, color = 'white', size = 3) +
    geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
    geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
                 stroke = 0.5, pch = 21) +
    scale_fill_brewer(palette = 'Set1', name = "",
                       labels = c("mon", "aes", "ies", "ams", "ims")) +
    scale_x_continuous(breaks = c(-10, 0, 10, 20),
                       labels = c("-500", "0", "500", "1000")) +
    labs(y = 'Proportion of target fixations',
         x = 'Time relative to target syllable offset (ms)',
         caption = "Mean +/- 95% CI") +
    annotate("text", x = 3.3, y = 0.02, label = '200ms',
             angle = 90, size = 3, hjust = 0) +
    theme_grey(base_size = 12, base_family = "Times")

ggsave_cols(paste0(figs_path, "/stress_p1.png"), stress_p1, cols = 1)
ggsave_cols(paste0(figs_path, "/stress_p1.eps"), stress_p1, cols = 1, device = cairo_ps)

# -----------------------------------------------------------------------------




# Plot GCA --------------------------------------------------------------------

# Within group differences
stress_p2 <- model_preds$fits_all %>%
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = cond, color = cond)) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.75) +
  geom_point(aes(shape = cond), color = "black", size = 1.3, show.legend = F) +
  geom_point(aes(shape = cond), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "Syllable structure") +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big + legend_adj

# Comparisons with natives
stress_p3 <- model_preds$fits_all %>%
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = group, color = group)) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.75) +
  geom_point(aes(shape = group), color = "black", size = 1.3, show.legend = F) +
  geom_point(aes(shape = group), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "Group") +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big + legend_adj_2

ggsave(paste0(figs_path, "/stress_p2.png"), stress_p2, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_p2.eps"), stress_p2, width = 150,
       height = 120, units = "mm", dpi = 600, device = cairo_ps)
ggsave(paste0(figs_path, "/stress_p3.png"), stress_p3, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_p3.eps"), stress_p3, width = 150,
       height = 120, units = "mm", dpi = 600, device = cairo_ps)

# -----------------------------------------------------------------------------








# Effect size plot ------------------------------------------------------------

bind_rows(
  data.comp %>%
    group_by(participant, group, time_zero, condition) %>%
    summarize(competition = mean(eLog)) %>%
    spread(condition, competition) %>%
    mutate(Effect = unstressed - stressed,
           condition_type = "Lexical stress",
           fit_type = "raw") %>%
    select(-unstressed, -stressed),
  data.comp %>%
    group_by(participant, group, time_zero, condition) %>%
    summarize(competition = mean(GCA_Full)) %>%
    spread(condition, competition) %>%
    mutate(Effect = unstressed - stressed,
           condition_type = "Lexical stress",
           fit_type = "model") %>%
    select(-unstressed, -stressed)) %>%
  spread(fit_type, Effect) %>%
  #filter(time_zero < 10) %>%                                          # was already commented by C.
  ggplot(., aes(x = time_zero, y = raw, color = group)) +
#    facet_grid(. ~ condition_type) +
    geom_hline(yintercept = 0, color = "white", size = 3) +
    geom_vline(xintercept = 4, color = "white", size = 3) +
    stat_summary(fun.data = mean_se, geom = "pointrange") +
    #stat_summary(aes(y = model, color = group), fun.y = mean,         # was already commented by C.
    #             geom = 'line', size = 1.4) +                         # was already commented by C.
    geom_smooth(method = 'gam', formula = y ~ poly(x, 3), se = F,
                show.legend = FALSE) +
    labs(x = "Time (ms) relative to target syllable offset",
         y = "Effect", caption = "Mean +/- 95% CI") +
    scale_color_brewer(palette = "Set1", name = "", guide = 'legend',
                       labels = c("mon", 'aes', 'ies', 'ams', 'ims')) +
    scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                       labels = c("-200", "0", "200", "400", "600")) +
    theme_grey(base_size = 15, base_family = "Times New Roman")

# -----------------------------------------------------------------------------











