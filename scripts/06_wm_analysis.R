# Working memory analysis -----------------------------------------------------




# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------


# Last reviewed on 2019/06/03


## 1. We check homogeneity for all participants in PSTM and WM
## 2. We removed participants from LAL group because of the DELE 
## (INT had higher L2 proficency, hence we removed the 5 lower proficiency LAL speakers)



wm_all <- read_csv(here("data", "pupurri_analysis.csv"))


# Create a L1 variable from the participant id info
wm_all <- wm_all %>%
  separate(., col = group, into = c('proficiency', 'l1'), sep = 1, remove = FALSE) %>%  # es = EN speaker, ms = MA speaker, on = monolingual ES speaker
  select(., -proficiency)


# Set accuracy as a continues variable
wm_all$WM_set <- as.numeric(wm_all$WM_set)


# Check that all participants per group are there

wm_all %>%
  group_by(group) %>%
  summarise(n_distinct(participant))

# Check means per L1
wm_all %>%
  group_by(l1) %>%
  summarise(mean = round(mean(WM_set),2), sd = round(sd(WM_set), 2), 
            n = length(unique(participant))) %>%
  knitr::kable()
# |l1 | mean|   sd|  n|
# |:--|----:|----:|--:|
# |es | 8.89| 2.11| 65|
# |ms | 7.78| 2.14| 64|
# |on | 6.20| 2.72| 30|

# Check means per group
wm_all %>%
  group_by(group) %>%
  summarise(mean = round(mean(WM_set),2), sd = round(sd(WM_set), 2), 
            n = length(unique(participant))) %>%
  knitr::kable()
# |group | mean|   sd|  n|
#   |:-----|----:|----:|--:|
#   |aes   | 9.38| 2.01| 32|
#   |ams   | 8.00| 2.18| 32|
#   |ies   | 8.42| 2.14| 33|
#   |ims   | 7.56| 2.11| 32|
#   |mon   | 6.20| 2.72| 30|


bartlett.test(WM_set ~ l1, data = wm_all) # variance ok
# Bartlett's K-squared = 3.1316, df = 2, p-value = 0.2089

bartlett.test(WM_set ~ group, data = filter(wm_all, l1 != "es")) # variance ok
# Bartlett's K-squared = 2.3365, df = 2, p-value = 0.3109

bartlett.test(WM_set ~ group, data = filter(wm_all, l1 != "on")) # variance ok
# Bartlett's K-squared = 0.22006, df = 3, p-value = 0.9743

bartlett.test(WM_set ~ group, data = filter(wm_all, l1 != "ms")) # variance ok
# Bartlett's K-squared = 3.1462, df = 2, p-value = 0.2074

library(TOSTER)
# wm toast mon vs int (all good)
TOSTtwo(m1 = 6.20, sd1 = 2.72, n1 = 30, # mon
        m2 = 8.89, sd2 = 2.11, n2 = 65, # EN speakers
        low_eqbound_d = -0.3, 
        high_eqbound_d = 0.3, 
        alpha = 0.05)

# wm toast mon vs lal (all good)
TOSTtwo(m1 = 6.20, sd1 = 2.72, n1 = 30, # mon
        m2 = 7.78, sd2 = 2.14, n2 = 64, # Mandarin speakers
        low_eqbound_d = -0.3, 
        high_eqbound_d = 0.3, 
        alpha = 0.05)

# wm toast int vs lal (all good)
TOSTtwo(m1 = 8.89, sd1 = 2.11, n1 = 65, # EN speakers
        m2 = 7.78, sd2 = 2.14, n2 = 64, # Mandarin speakers
        low_eqbound_d = -0.3, 
        high_eqbound_d = 0.3, 
        alpha = 0.05)

# Check that WM is ok after removing the 5 participants

wm_all %>%
  group_by(l1) %>%
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
           participant != 'AES16' & participant != 'AES23') %>%
  summarise(mean = round(mean(WM_set),2), sd = round(sd(WM_set), 2), 
            n = length(unique(participant))) %>%
  knitr::kable()
# |l1 | mean|   sd|  n|
#   |:--|----:|----:|--:|
#   |es | 8.08| 1.78| 49|
#   |ms | 7.78| 2.14| 64|
#   |on | 7.65| 1.95| 20|


bartlett.test(WM_set ~ l1, data = filter(wm_all,participant != "MON04" & participant != "MON05" &
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
                                              participant != 'AES16' & participant != 'AES23'))

# Bartlett's K-squared = 1.8361, df = 2, p-value = 0.3993

# wm toast mon vs int (all good)
TOSTtwo(m1 = 7.65, sd1 = 1.95, n1 = 20, # mon
        m2 = 8.08, sd2 = 1.78, n2 = 49, # EN
        low_eqbound_d = -0.3, 
        high_eqbound_d = 0.3, 
        alpha = 0.05)


# wm toast mon vs lal (all good)
TOSTtwo(m1 = 7.65, sd1 = 1.95, n1 = 20, # mon
        m2 = 7.78, sd2 = 2.14, n2 = 64, # MA
        low_eqbound_d = -0.3, 
        high_eqbound_d = 0.3, 
        alpha = 0.05)

# wm toast int vs lal (all good)
TOSTtwo(m1 = 8.08, sd1 = 1.78, n1 = 49, # EN
        m2 = 7.78, sd2 = 2.14, n2 = 64, # MA
        low_eqbound_d = -0.3, 
        high_eqbound_d = 0.3, 
        alpha = 0.05)

# Graph to check outliers in wm

wm_all %>% 
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
           participant != 'AES16' & participant != 'AES23') %>%
  ggplot(., aes(x = l1, y = WM_set, label = participant, color = l1)) + 
  geom_text() + 
  stat_summary(fun.data = mean_sdl, geom = "pointrange") + 
  coord_flip()


# -----------------------------------------------------------------------------


# 
# 
# wm_dataset <- read.spss(here("data", "raw", "gating.sav"), to.data.frame = TRUE)
# 
# str(wm_dataset)
# 
# wm_df <- wm_dataset %>%
#       select(Exp, Group, participant = ExperimentName, WM, Condition) %>%
#       filter(., Exp == 'S', Group != 'HS' & Group != 'IN') %>%
#       as.tbl(.)
# 
# # remove unwanted characters from column
# wm_df$participant <- gsub(" ", "", paste(wm_df$participant))
# 
# str(wm_df)
# 
# 
# # Remove participants according to homogeneity of variance
# # test re: working memory
# remove <- c("L20", "L21", "L22", "L23", "L30", "L31", "L02", "L05",
#             "L06", "L08", "L10", "L15")
# wm_df2 <- filter(wm_df, !participant %in% remove)
# 
# 
# str(wm_df)
# str(wm_df2)
# 
# 
# wm_df %>%
#   ggplot(., aes(x = Group, y = WM)) +
#     geom_boxplot()
# 
# wm_df2 %>%
#   ggplot(., aes(x = Group, y = WM)) +
#     geom_boxplot()
# 
# unique(wm_df[wm_df$Group == 'L', 'WM'])
# unique(wm_df2[wm_df2$Group == 'L', 'WM'])
