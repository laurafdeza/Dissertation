

# Load data and models --------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
library("report") 
library(nortest)

# Load data
source(here::here("scripts", "02_load_data.R"))


# -----------------------------------------------------------------------------







# Data prep -------------------------------------------------------------------


stress_gc_subset <- stress50 %>%
  # select(., -WM_set) %>%
  filter(., time_zero >= -4 & time_zero <= 12) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma"),
            condition_sum = if_else(cond == "1", 1, -1)) %>%       # 1 = present, 2 = past
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")


l2_data <- stress_gc_subset%>%
  filter(., l1 != 'es') %>% 
  filter(., participant != 'ies04' & participant != 'ies17' & participant != 'ies28' & participant != 'aes32') %>%
  mutate(., l1_sum = if_else(l1 == 'en', -1, 1))


# -----------------------------------------------------------------------------


# ASSUMPTION 1. Normality
ggplot(l2_data, aes(x=DELE_z)) + 
  geom_density()

ggplot(l2_data, aes(x=use_z)) + 
  geom_density()


ad.test(l2_data$DELE_z) # Anderson-Darling normality test for large samples
# A = 210.38, p-value < 2.2e-16

ad.test(l2_data$use_z) # neither normal > Spearman
# A = 1028.2, p-value < 2.2e-16


# Pearson correlation between 2 variables
cor(l2_data$DELE_z, l2_data$use_z) 
# r = 0.2558775 --> both increase together. Weak relationship bc closer to 0 than to 1

# Spearman correlation between 2 variables
cor(l2_data$DELE_z, l2_data$use_z, method = "spearman")
# 0.2747663

# Kendall correlation between 2 variables
cor(l2_data$DELE_z, l2_data$use_z, method = "kendall") # , digits = 2
# 0.2041101

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
ggplot(l2_data) +
  aes(x = DELE_z, y = use_z, color = l1, shape = l1) +
  geom_point() + #colour = "#0c4c8a"
  xlab("Proficiency") + ylab("L2 use") +
  geom_smooth(size = .5) + # remove for loess method / method=lm for linear
  geom_smooth(method=lm) +
  scale_color_discrete(name = "L1", labels = c("English", "Mandarin Chinese")) +
  guides(shape = FALSE) +
  theme_gray(base_size = 12,
             base_family = "Times") +
  theme(legend.position = "bottom",
        legend.box = "vertical")
  # theme_minimal()

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
# The Pearson's product-moment correlation between l2_data$DELE_z and l2_data$use_z is positive, statistically significant, and medium (r = 0.26, 95% CI [0.25, 0.27], t(33737) = 48.62, p < .001)

# ---------------------------------------------------------------------------------


# ENGLISH SPEAKERS

en_data <- l2_data %>%
  filter(., l1 == 'en') 

# Pearson correlation between 2 variables
cor(en_data$DELE_z, en_data$use_z) 
# 0.2572734 --> both increase together. Weak relationship bc closer to 0 than to 1

# Spearman correlation between 2 variables
cor(en_data$DELE_z, en_data$use_z, method = "spearman")
# 0.3372631

# Kendall correlation between 2 variables
cor(en_data$DELE_z, en_data$use_z, method = "kendall") # , digits = 2
# 0.2519289

# Plot correlation
ggplot(en_data) +
  aes(x = DELE_z, y = use_z) +
  geom_point(colour = "#0c4c8a") +
  xlab("Proficiency") + ylab("L2 use") +
  theme_minimal()

report(cor.test(en_data$DELE_z, en_data$use_z, method = "spearman", exact=FALSE))

# ---------------------------------------------------------------------------------


# MANDARIN CHINESE SPEAKERS

ma_data <- l2_data %>%
  filter(., l1 == 'ma') 

# Pearson correlation between 2 variables
cor(ma_data$DELE_z, ma_data$use_z) 
# 0.2632409 --> both increase together. Weak relationship bc closer to 0 than to 1

# Spearman correlation between 2 variables
cor(ma_data$DELE_z, ma_data$use_z, method = "spearman")
# 0.2554505

# Kendall correlation between 2 variables
cor(ma_data$DELE_z, ma_data$use_z, method = "kendall") # , digits = 2
# 0.1868851

# Plot correlation
ggplot(ma_data) +
  aes(x = DELE_z, y = use_z) +
  geom_point(colour = "#0c4c8a") +
  xlab("Proficiency") + ylab("L2 use") +
  theme_minimal()

report(cor.test(ma_data$DELE_z, ma_data$use_z, method = "spearman", exact=FALSE))


