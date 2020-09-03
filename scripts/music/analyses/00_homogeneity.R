library(tidyverse); library("TOSTER")

# Load data
dem_all <- read_csv("./data/pupurri_analysis.csv")


# Remove participants to make the groups homogenous in L2
# proficiency (DELE)

dem_all <- dem_all %>%
  filter(., participant != "MON04" & participant != "MON05" & participant != "MON06" &
           participant != "MON11" & participant != "MON16" & participant != "MON18" &
           participant != "MON20" & participant != "MON21" & participant != "MON23" &
           participant != "MON24" & participant != "MON27") %>%
  filter(., participant != "IES05" & participant != "IES09" & participant != "IES16" &
           participant != "IES21" & participant != "IES23" & participant != "IES30" &
           participant != "IES33") %>%
  filter(., participant != "AES01" & participant != "AES03" & participant != "AES04" &
           participant != "AES05" & participant != "AES07" & participant != "AES11" &
           participant != "AES14" & participant != "AES16" & participant != "AES23") %>%
  filter(., participant != "IMS02" & participant != "IMS14" & participant != "IMS17" &
           participant != "AMS02" &
           participant != "AMS05" & participant != "AMS07" & participant != "AMS09" & 
           participant != "AMS12" & participant != "AMS13" & participant != "AMS14" &
           participant != "AMS16" & participant != "AMS21" & participant != "AMS24" &
           participant != "AMS28" & participant != "AMS30") %>%
  separate(., col = group,
           into = c("proficiency", "l1"), # es = EN speaker, ms = MA speaker, on = ES speaker
           sep = 1,
           remove = FALSE)

unique(dem_all$participant)

# dem_all %>%
#   group_by(l1) %>%
#   summarize(n = n_distinct(id))



# Check info: id, group, dele, pstm, wm, age, aoa_l2, months_abroad,
# l1_use_week, l2_use_week, years_word_int

glimpse(dem_all)

# We're missing wm for LA07 (after computer crashed and we lost data)

# Set wm as numeric

dem_all$DELE <- as.numeric(dem_all$DELE)

# Create a table with mean + sd for: wm, pstm, dele, aoa_l2, months_abroad,
# l1_use_week, l2_use_week, years_word_int

dem_all %>%
  group_by(., group) %>%
  filter(., group != "mon") %>%
  summarise(., dele = round(mean(DELE),2),
               dele_sd = round(sd(DELE),2),
            n = length(unique(participant))) %>%
  knitr::kable()
#                wm = round(mean(wm, na.rm=TRUE),2),
#                wm_sd = round(sd(wm, na.rm=TRUE),2),
#                pstm = round(mean(pstm),2),
#                pstm_sd = round(sd(pstm),2),
#                ,
#                ) 


# |group |  dele| dele_sd|  n|
# |:-----|-----:|-------:|--:|
# |aes   | 45.48|    4.09| 23|
# |ams   | 45.50|    3.97| 32|
# |ies   | 31.19|    4.78| 26|
# |ims   | 32.47|    4.09| 30|




dem_all %>%
  group_by(., group) %>%
  summarise(aoa_l2 = round(mean(AoA_L2),2),
            aoa_l2_sd = round(sd(AoA_L2),2),
            abroad = round(mean(mo_ES_country),2),
            abroad_sd = round(sd(mo_ES_country),2),
            l1_use = round(mean(percent_l1_week),2),
            l1_use_sd = round(sd(percent_l1_week),2),
            l2_use = round(mean(percent_l2_week),2),
            l2_use_sd = round(sd(percent_l2_week),2),
            n = length(unique(participant))
           # ,
           # n = n_distinct(id)
            ) %>% knitr::kable()
# |group | aoa_l2| aoa_l2_sd| abroad| abroad_sd| l1_use| l1_use_sd| l2_use| l2_use_sd|  n|
# |:-----|------:|---------:|------:|---------:|------:|---------:|------:|---------:|--:|
# |aes   |  15.57|      3.49|  55.00|     37.43|  59.57|     16.09|  40.43|     16.09| 23|
# |ams   |  17.88|      2.83|  45.12|     55.68|  52.81|     22.57|  46.56|     21.76| 32|
# |ies   |  18.04|      6.65|  26.69|     21.03|  71.73|     17.83|  28.27|     17.83| 26|
# |ims   |  19.93|      4.14|  37.30|     33.42|  65.17|     19.63|  34.83|     19.63| 30|
# |mon   |   5.32|      2.26|   0.47|      2.06|  11.95|      9.48|  87.53|      9.64| 19|



# Some SD values don't work (not sure why), but they work with this function
# aggregate(wm ~ group, data = dem_all, FUN = sd, na.rm=TRUE)
# aggregate(pstm ~ group, data = dem_all, FUN = sd)
# aggregate(dele ~ group, data = dem_all, FUN = sd)
# aggregate(aoa_l2 ~ group, data = dem_all, FUN = sd)



## Homogeneity of variances tests
dem_all <- dem_all %>%
  filter(., group != "mon")
## All seem ok
bartlett.test(DELE ~ l1, data = dem_all)
bartlett.test(AoA_L2 ~ l1, data = dem_all) # p-value = 0.002705
bartlett.test(mo_ES_country ~ l1, data = dem_all) # p-value = 0.01553
bartlett.test(percent_l1_week ~ l1, data = dem_all)
bartlett.test(percent_l2_week ~ l1, data = dem_all)



# dele tost intermediate groups
# all good
TOSTtwo(m1 = 31.19, sd1 = 4.78, n1 = 26, # EN
        m2 = 32.47, sd2 = 4.09, n2 = 30, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# dele tost advanced groups
# all good
TOSTtwo(m1 = 45.48, sd1 = 4.09, n1 = 23, # EN
        m2 = 45.50, sd2 = 3.97, n2 = 32, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)



# age of acquistion tost intermediate
# all good
TOSTtwo(m1 = 18.04, sd1 = 6.65, n1 = 26, # EN
        m2 = 19.59, sd2 = 3.75, n2 = 29, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# age of acquistion tost advanced
# all good
TOSTtwo(m1 = 15.57, sd1 = 3.49, n1 = 23, # EN
        m2 = 16.45, sd2 = 2.37, n2 = 20, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)





# time abroad toast la vs in
# all good, but Barlett is not ok, there's a lot of variance in the int group, what do we do here?
TOSTtwo(m1 = 34.14, sd1 = 86.99, n1 = 22, # in
        m2 = 12.68, sd2 = 15.13, n2 = 25, # la
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# # l1 use in a normal week toast la vs in
# # They are different, Interpreters use less their L1
# TOSTtwo(m1 = 64.09, sd1 = 11.51, n1 = 22, # in
#         m2 = 72.76, sd2 = 12.91, n2 = 25, # la
#         low_eqbound_d = -0.3,
#         high_eqbound_d = 0.3,
#         alpha = 0.05)

# l2 use in a normal week toast 
# all good, this might be because some interpreters also use other L2s
# for example, the United Nations interpreters speak French and use it more often than Spanish
TOSTtwo(m1 = 31.59, sd1 = 14.59, n1 = 22, # in
        m2 = 27.24, sd2 = 12.91, n2 = 25, # la
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
