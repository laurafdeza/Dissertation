# Original script by Joseph for DurStress project
# Modified by Cris
#
# Last update: 06/12/2019
# T-tests ---------------------------------------------------------------------
#
# - Question: can participants predict morphology before onset of target
#   suffix?
# - We test this for each group of participants for each type of word
#   (paroxytone, oxytone)
# - This analysis does not compare groups or conditions
#
# -----------------------------------------------------------------------------


# Load data -------------------------------------------------------------------

source(here::here("scripts", "02_load_data.R"))

# -----------------------------------------------------------------------------



# Data prep -------------------------------------------------------------------

# Create stress subset that includes interpreters, late advanced and natives
# Must check on the participants that are removed
# We want to analyze proportion of target gaze at target onset
# so we need to make a subset of the data that only uses the
# target onset bin (adjusted 200ms for VWP)
# We are using 10ms bins so we want time_zero bin 20 (200/10 = 20)
# We remove participants to make them comparable in terms of WM

df_short <- stress10 %>%
  filter(participant != "MON04" & participant != "MON05" &
           participant != "MON06" &
           participant != "MON11" & participant != "MON016" & 
           participant != "MON18" & participant != "MON20" & 
           participant != "MON21" & participant != "MON23" &
           participant != "MON24" & participant != "MON27" &
           participant != "IES05" & participant != "IES09" & 
           participant != 'IES16' & participant != 'IES21' &
           participant != "IES23" & participant != "IES30" &
           participant != "IES33" & participant != "AES01" &
           participant != 'AES03' & participant != 'AES04' &
           participant != "AES05" & participant != "AES07" &
           participant != "AES11" & participant != "AES14" &
           participant != 'AES16' & participant != 'AES23',
            time_zero == 20)



# -----------------------------------------------------------------------------






# T-tests unrelated ---------------------------------------------------------------------
## Stress


# Quick and dirty mean of target fixations as a function of
# group and condition (stressed, unstressed *1st syllable*)

df_short %>%
  na.omit(.) %>%
  group_by(., group, cond) %>%
  summarise(., meanFix = mean(target_prop))

# We will test this for each group in each condition (stressed, untressed)
# using a one-sided t-test. Specifically, we are testing the
# hypothesis that the proportion of looks is greater than
# chance (50%).
# - H0: u = 0.50
# - Ha: u > 0.50
# The generic code is:
# t.test(myVector, alternative = "greater", my = 0.33, conf.level = 0.95)

stress_unrel_ttest <- df_short %>%
#  filter(., exp == "stress_unrelated") %>%
  na.omit(.) %>%
  group_by(., group, cond, participant) %>%
  summarise(., meanFix = mean(target_prop)) %>%
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.5, conf.level = 0.95)))

# Convert pvalues from scientific notation
stress_unrel_ttest$p.value <- format(stress_unrel_ttest$p.value, scientific = F)
stress_unrel_ttest$sig <- "N.S."
stress_unrel_ttest[stress_unrel_ttest$p.value <= 0.05/6, 'sig'] <- "*"

# Print results
print(as.data.frame(stress_unrel_ttest[, c(1:7, 11)]))

saveRDS(stress_unrel_ttest, "./reports/mods/stress_ttest.rds", compress = "xz")
 
# -----------------------------------------------------------------------------


# Graph for stress_unrelated

# stress_unrel_ttest$group <- factor(stress_unrel_ttest$group, levels = c("aes", "ies", "ams", "ims", "mon"))
# stress_unrel_ttest$cond <- factor(stress_unrel_ttest$cond, levels = c("pres", "pret"))   # Not sure if these are the column names
stress_unrel_ttest$cond <- as.factor(stress_unrel_ttest$cond)

stress_unel_TargetFix <- stress_unrel_ttest %>% 
  ggplot(., aes(x = group, y = estimate, color = cond, 
                group = interaction(group, cond), dodge = cond)) +
  geom_hline(yintercept = 0.5, lty = 3) + 
  geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                 position = position_dodge(width = 0.75), size = 1) +
  geom_point(position = position_dodge(width = 0.75), size = 4) +
  ylim(0, 1.0) + ylab('Target fixations') + xlab('Group') + 
  ggtitle('Mean fixations on target and\nlower-bound 95% confidence interval') + 
  scale_color_brewer(palette = "Set2", name = '', labels = c('Present', 'Preterit')) + # reorder?
#  scale_x_discrete(labels = c("aes", "ams", "ies", "ims", "mon")) +  
  theme_bw(base_size = 16, base_family = 'Times') # -> stress_unel_TargetFix


# Looks good, save as .png file. 
ggsave('stress_unrel_TargetFix.png', plot = stress_unel_TargetFix, dpi = 600, device = "png", path = "./reports/figs")

