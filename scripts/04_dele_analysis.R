# DELE analysis ---------------------------------------------------------------



# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------


dem_all <- read_csv("./data/pupurri_analysis.csv")

### DELE

dem_all$DELE <- as.numeric(dem_all$DELE)
dem_all %>% 
  filter(., group %in% c("ies", "aes", "ims", "ams")) %>%
  group_by(., group) %>% 
  summarize(., DELE_mean = mean(DELE), 
            DELE_sd = sd(DELE), 
            n = n())

# results including all participants collected

#   group   DELE_mean  DELE_sd n
# 1 aes        45.4    4.26    32
# 2 ams        45.5    3.97    32
# 3 ies        31.7    4.58    33
# 4 ims        32.8    4.23    32


bartlett.test(DELE ~ group, data = dem_all) # all good
# Bartlett's K-squared = 0.64069, df = 3, p-value = 0.8871

# DELE tost aes vs ams (OK, not significant and between dotted lines so they are not different)
TOSTER::TOSTtwo(m1 = 45.4, sd1 = 4.26, n1 = 32, # aes
        m2 = 45.5, sd2 = 3.97, n2 = 32, # ams
        low_eqbound_d = -0.3, 
        high_eqbound_d = 0.3, 
        alpha = 0.05)
# t(61.69) = -0.0971, p = 0.923

# DELE tost ies vs ims (idem)
TOSTER::TOSTtwo(m1 = 31.7, sd1 = 4.58, n1 = 33, # ies
                m2 = 32.8, sd2 = 4.23, n2 = 32, # ims
                low_eqbound_d = -0.3, 
                high_eqbound_d = 0.3, 
                alpha = 0.05)
# t(62.85) = -1.006, p = 0.318