

# Load data and models --------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
library("report") 
library(nortest)

# Load data
stress50 <- read_csv(here("data", "clean", "stress_50ms_allparticipants.csv"))


# -----------------------------------------------------------------------------







# Data prep -------------------------------------------------------------------


stress_gc_subset <- stress50 %>%
  filter(., time_zero >= -4 & time_zero <= 12 & l1 != 'ma') %>%
  select(., -prof) %>%
  mutate(., condition_sum = if_else(cond == "1", -1, 1)) %>%       # 1 = present, 2 = past
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")


l2_data <- stress_gc_subset%>%
  filter(., l1 == 'en') %>% 
  mutate(., use_z = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week),
         DELE_z = (DELE - mean(DELE))/sd(DELE))


# -----------------------------------------------------------------------------


# ASSUMPTION 1. Normality
ggplot(l2_data, aes(x=DELE_z)) + 
  geom_density()
# if error for graphics, run dev.off() and then run code again

ggplot(l2_data, aes(x=use_z)) + 
  geom_density()


ad.test(l2_data$DELE_z) # Anderson-Darling normality test for large samples
# A = 147.71, p-value < 2.2e-16

ad.test(l2_data$use_z) # neither normal > Spearman
# A = 584.18, p-value < 2.2e-16


# Spearman correlation between 2 variables
cor(l2_data$DELE_z, l2_data$use_z, method = "spearman")
# 0.3844002


# Regarding the direction of the relationship: On the one hand, 
# a negative correlation implies that the two variables under consideration 
# vary in opposite directions, that is, if a variable increases the other decreases 
# and vice versa. On the other hand, a positive correlation implies 
# that the two variables under consideration vary in the same direction, i.e., 
# if a variable increases the other one increases and if one decreases 
# the other one decreases as well.
# 
# Regarding the strength of the relationship: The more extreme 
# the correlation coefficient (the closer to -1 or 1), the stronger the relationship. 
# This also means that a correlation close to 0 indicates that the two variables 
# are independent, that is, as one variable increases, there is no tendency 
# in the other variable to either decrease or increase.
#
# Source: https://statsandr.com/blog/correlation-coefficient-and-correlation-test-in-r/#data

# Plot correlation
corr_all <- ggplot(l2_data) +
  aes(x = DELE_z, y = use_z) +
  geom_point(size = .8) + #colour = "#0c4c8a"
  xlab("L2 proficiency") + ylab("L2 use") +
  geom_smooth(method=lm) + # for linear
  theme_gray(base_size = 12,
             base_family = "Times") +
  theme(legend.position = "bottom",
        legend.box = "vertical")

figs_path <- here("figs", "use_prof")
ggsave(paste0(figs_path, "/correlation_all.png"), corr_all, width = 180,
       height = 120, units = "mm", dpi = 600)


# a correlation test is used to test whether the correlation (denoted ρ) 
# between 2 variables is significantly different from 0 or not in the population.
# H0:  ρ = 0 (meaning that there is no linear relationship between the two variables)
# H1:  ρ ≠ 0 (meaning that there is a linear relationship between the two variables)
# Example: The p-value of the correlation test between these 2 variables is 0.62. 
# At the 5% significance level, we do not reject the null hypothesis of no correlation. 
# We therefore conclude that we do not reject the hypothesis that there is no linear
# relationship between the 2 variables.
test <- cor.test(l2_data$DELE_z, l2_data$use_z, method = "spearman", exact=FALSE)
test
# At 5% sig level, p < .05 we reject the null hypothesis of no correlation.
# We therefore conclude that we reject the hypothesis that there is no linear relationship

report(test)
# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Spearman's rank correlation rho between l2_data$DELE_z and l2_data$use_z is positive, statistically significant, and large (rho = 0.38, S = 5.51e+11, p < .001)
# ---------------------------------------------------------------------------------


