

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
corr_all <- ggplot(l2_data) +
  aes(x = DELE_z, y = use_z) +
  geom_point(size = .8) + #colour = "#0c4c8a"
  xlab("L2 proficiency") + ylab("L2 use") +
  geom_smooth(method=lm) + # for linear
  theme_gray(base_size = 12,
             base_family = "Times") +
  theme(legend.position = "bottom",
        legend.box = "vertical")

figs_path <- here("figs", "stress", "gca")
ggsave(paste0(figs_path, "/correlation_all.png"), corr_all, width = 180,
       height = 120, units = "mm", dpi = 600)

# Plot correlation by group
corr_plot <- ggplot(l2_data) +
  aes(x = DELE_z, y = use_z, shape = l1, color = l1) +
  geom_point(size = 2) + #colour = "#0c4c8a"
  xlab("L2 proficiency") + ylab("L2 use") +
  #geom_smooth(size = .5) + # remove for loess method 
  geom_smooth(method=lm, size = .7) + # for linear
  scale_color_discrete(name = "L1", labels = c("English", "Mandarin Chinese")) +
  scale_shape_discrete(name = "L1", labels = c("English", "Mandarin Chinese")) +
  # guides(shape = FALSE) +
  # guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_gray(base_size = 12,
             base_family = "Times") +
  theme(legend.position = "bottom",
        legend.box = "vertical")
  # theme_minimal()

figs_path <- here("figs", "stress", "gca")
ggsave(paste0(figs_path, "/correlation_byl1.png"), corr_plot, width = 180,
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







# ---------------------------------------------------------------------------------

# Build up random effects to test time terms
  
  mod_ot0 <- lmer(eLog ~ 1 + (1 | participant),
                  control = lmerControl(optimizer = 'bobyqa'),   # , optCtrl=list(maxfun=2e5)
                  data = l2_data, weights = 1/wts, REML = F)
  
  mod_ot1 <- update(mod_ot0, . ~ . -(1 | participant) +
                      ot1 + (1 + ot1 | participant))
  
  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + ot1 | participant) +
             ot2 + (1 + ot1 + ot2 | participant))
  
  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + ot1 + ot2 | participant) +
             ot3 + (1 + ot1 + ot2 + ot3 | participant))
  
  anova(mod_ot0, mod_ot1, mod_ot2, mod_ot3)
  #         npar    AIC    BIC logLik deviance    Chisq Df Pr(>Chisq)    
  # mod_ot0    3 189913 189938 -94953   189907                           
  # mod_ot1    6 186090 186141 -93039   186078 3828.730  3  < 2.2e-16 ***
  # mod_ot2   10 186018 186102 -92999   185998   80.268  4  < 2.2e-16 ***
  # mod_ot3   15 185604 185730 -92787   185574  424.141  5  < 2.2e-16 ***
    
  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))
  
  mod_ot5 <- update(mod_ot4, . ~ . -(1 | target) + (1 + ot1 | target))
  
  mod_ot6 <- update(mod_ot5, . ~ . -(1 + ot1 | target) +
                       + (1 + ot1 + ot2 | target))
  
  mod_ot7 <- update(mod_ot6, . ~ . -(1 + ot1 + ot2 | target) +
                       + (1 + ot1 + ot2 + ot3 | target))
  
  anova(mod_ot3, mod_ot4, mod_ot5, mod_ot6, mod_ot7)
  #         npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)    
  # mod_ot3   15 185604 185730 -92787   185574                          
  # mod_ot4   16 185123 185258 -92546   185091 482.360  1  < 2.2e-16 ***
  # mod_ot5   18 184755 184907 -92360   184719 372.212  2  < 2.2e-16 ***
  # mod_ot6   21 184623 184800 -92290   184581 138.163  3  < 2.2e-16 ***
  # mod_ot7   25 184569 184780 -92260   184519  61.562  4  1.362e-12 ***
  
  
  
  # Build fixed effects
  
  gca_l2_mod_base <- mod_ot7
  # lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +         
  #        (1 + ot1 + ot2 + ot3 | participant) +
  #        (1 + ot1 + ot2 + ot3 | target), 
  #      control = lmerControl(optimizer = 'bobyqa',
  #                            optCtrl = list(maxfun = 2e4)),
  #      data = l2_data, REML = F)
  
  # add proficiency effect to intercept, linear slope, quadratic, and cubic time terms
  
  gca_l2_mod_dele_0 <- update(gca_l2_mod_base,   . ~ . + DELE_z) 
  gca_l2_mod_dele_1 <- update(gca_l2_mod_dele_0, . ~ . + ot1:DELE_z) 
  gca_l2_mod_dele_2 <- update(gca_l2_mod_dele_1, . ~ . + ot2:DELE_z)
  gca_l2_mod_dele_3 <- update(gca_l2_mod_dele_2, . ~ . + ot3:DELE_z)
  
  dele_anova <-
    anova(gca_l2_mod_base, gca_l2_mod_dele_0, gca_l2_mod_dele_1,
          gca_l2_mod_dele_2, gca_l2_mod_dele_3)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
  # gca_l2_mod_base     25 184569 184780 -92260   184519                        
  # gca_l2_mod_dele_0   26 184566 184785 -92257   184514 4.9477  1   0.026125 * 
  # gca_l2_mod_dele_1   27 184561 184789 -92254   184507 7.1857  1   0.007349 **
  # gca_l2_mod_dele_2   28 184562 184798 -92253   184506 1.3448  1   0.246184   
  # gca_l2_mod_dele_3   29 184563 184807 -92252   184505 1.2716  1   0.259462 
  
  
  
  # add L2 use effect to intercept, linear slope, quadratic, and cubic time terms
  
  gca_l2_mod_use_0 <- update(gca_l2_mod_dele_1,    . ~ . + use_z) 
  gca_l2_mod_use_1 <- update(gca_l2_mod_use_0, . ~ . + ot1:use_z) 
  gca_l2_mod_use_2 <- update(gca_l2_mod_use_1, . ~ . + ot2:use_z)
  gca_l2_mod_use_3 <- update(gca_l2_mod_use_2, . ~ . + ot3:use_z)
  
  use_anova <-
    anova(gca_l2_mod_dele_1, gca_l2_mod_use_0, gca_l2_mod_use_1,
          gca_l2_mod_use_2, gca_l2_mod_use_3)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_l2_mod_dele_1   27 184561 184789 -92254   184507                     
  # gca_l2_mod_use_0    28 184561 184797 -92253   184505 1.8411  1     0.1748
  # gca_l2_mod_use_1    29 184562 184807 -92252   184504 1.0397  1     0.3079
  # gca_l2_mod_use_2    30 184562 184815 -92251   184502 1.9168  1     0.1662
  # gca_l2_mod_use_3    31 184562 184823 -92250   184500 2.2572  1     0.1330
  
  
  
  # add interaction
  
  gca_l2_mod_int_0 <- update(gca_l2_mod_dele_1,    . ~ . + DELE_z:use_z) 
  gca_l2_mod_int_1 <- update(gca_l2_mod_int_0, . ~ . + ot1:DELE_z:use_z) 
  gca_l2_mod_int_2 <- update(gca_l2_mod_int_1, . ~ . + ot2:DELE_z:use_z)
  gca_l2_mod_int_3 <- update(gca_l2_mod_int_2, . ~ . + ot3:DELE_z:use_z)
  
  int_anova <-
    anova(gca_l2_mod_dele_1, gca_l2_mod_int_0, gca_l2_mod_int_1,
          gca_l2_mod_int_2, gca_l2_mod_int_3)
  #                   npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
  # gca_l2_mod_dele_1   27 184561 184789 -92254   184507                     
  # gca_l2_mod_int_0    28 184563 184799 -92254   184507 0.0198  1     0.8881
  # gca_l2_mod_int_1    29 184565 184809 -92254   184507 0.0716  1     0.7890
  # gca_l2_mod_int_2    30 184566 184819 -92253   184506 0.7740  1     0.3790
  # gca_l2_mod_int_3    31 184566 184827 -92252   184504 2.1277  1     0.1447
  
  
  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  gca_l2_mod_l1_0 <- update(gca_l2_mod_dele_1, . ~ . + l1_sum) 
  gca_l2_mod_l1_1 <- update(gca_l2_mod_l1_0, . ~ . + ot1:l1_sum) 
  gca_l2_mod_l1_2 <- update(gca_l2_mod_l1_1, . ~ . + ot2:l1_sum) 
  gca_l2_mod_l1_3 <- update(gca_l2_mod_l1_2, . ~ . + ot3:l1_sum)
  
  l1_anova <-
    anova(gca_l2_mod_dele_1, gca_l2_mod_l1_0, gca_l2_mod_l1_1,
          gca_l2_mod_l1_2, gca_l2_mod_l1_3)
  #                   npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)  
  # gca_l2_mod_dele_1   27 184561 184789 -92254   184507                         
  # gca_l2_mod_l1_0     28 184563 184799 -92253   184507  0.3421  1   0.558645   
  # gca_l2_mod_l1_1     29 184555 184799 -92248   184497 10.0551  1   0.001519 **
  # gca_l2_mod_l1_2     30 184556 184809 -92248   184496  0.3615  1   0.547678   
  # gca_l2_mod_l1_3     31 184556 184817 -92247   184494  2.3920  1   0.121956 
  
  
  # add interaction
  
  gca_l2_mod_l1_int_0 <- update(gca_l2_mod_l1_1,     . ~ . + l1_sum:DELE_z:use_z) 
  gca_l2_mod_l1_int_1 <- update(gca_l2_mod_l1_int_0, . ~ . + ot1:l1_sum:DELE_z:use_z) 
  gca_l2_mod_l1_int_2 <- update(gca_l2_mod_l1_int_1, . ~ . + ot2:l1_sum:DELE_z:use_z)
  gca_l2_mod_l1_int_3 <- update(gca_l2_mod_l1_int_2, . ~ . + ot3:l1_sum:DELE_z:use_z)
  
  l1_int_anova <-
    anova(gca_l2_mod_l1_1, gca_l2_mod_int_0, gca_l2_mod_int_1,
          gca_l2_mod_int_2, gca_l2_mod_int_3)
  # npar    AIC    BIC logLik deviance   Chisq Df Pr(>Chisq)   
  # gca_l2_mod_int_0   28 184563 184799 -92254   184507                         
  # gca_l2_mod_l1_1    29 184555 184799 -92248   184497 10.3773  1   0.001276 **
  # gca_l2_mod_int_1   29 184565 184809 -92254   184507  0.0000  0   1.000000   
  # gca_l2_mod_int_2   30 184566 184819 -92253   184506  0.7740  1   0.378990   
  # gca_l2_mod_int_3   31 184566 184827 -92252   184504  2.1277  1   0.144655 
  
  mod_type <- "gca_l2_mod"
  mod_spec <- c('_base', 
                "_l1_0", "_l1_1", "_l1_2", "_l1_3", 
                "_dele_0", "_dele_1", "_dele_2", "_dele_3", 
                "_int_0", "_int_1", "_int_2", "_int_3",
                "_use_0", "_use_1", "_use_2", "_use_3" 
                )
  
  # Store ind models in list
  gca_corr <- mget(c(paste0(mod_type, mod_spec)))
  
  save(gca_corr,
       file = here("mods", "stress", "gca",
                   "gca_corr.Rdata"))
  