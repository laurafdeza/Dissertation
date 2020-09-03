# Speech rate analysis --------------------------------------------------------
# Last update. 06/03/2019
# Needs to be reviewed

# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------


## Get speech rate

# Load file after running praat script on the audio stimuli

speech_rate <- read_csv("./data/speech_rate.txt")

glimpse(speech_rate)

mean(speech_rate$`speechrate (nsyll/dur)`) # 4.366774
sd(speech_rate$`speechrate (nsyll/dur)`)  # 0.6752352


# -----------------------------------------------------------------------------













# allinfo <- read_csv("./data/potpurri_master.csv")
# 
# glimpse(allinfo)
# 
# test <- "this sentence is fairly short"
# str_count(test, '\\w+')
# sentence_info <- allinfo %>%
#   # filter(., exp == "stress") %>%
#   select(sentence)
# 
# # # Remove duplicate rows
# # sentence_info <- sentence_info[!duplicated(sentence_info$sentence), ]
# 
# as.vector(sentence_info)
# 
# ## Get # of syllables with the function str_count() from the stringr package
# 
# stress_sent <- filter(allinfo, exp == "stress", version == 1 | version == 2)
# 
# max(str_count(stress_sent$sentence, '\\s+')+1)
# min(str_count(stress_sent$sentence, '\\s+')+1)
# mean(str_count(stress_sent$sentence, '\\s+')+1)
# sd(str_count(stress_sent$sentence, '\\s+')+1)
# 
# # > max(str_count(allinfo$sentence, '\\s+')+1)
# # [1] 5
# # > min(str_count(allinfo$sentence, '\\s+')+1)
# # [1] 5
# # > mean(str_count(allinfo$sentence, '\\s+')+1)
# # [1] 5
# # > sd(str_count(allinfo$sentence, '\\s+')+1)
# # [1] 0



# -----------------------------------------------------------------------------

# From Joseph's script:

# sr_df <- read_csv(here("stimuli", "stress", "output", "stress_speech_rate.csv")) %>%
#   rename(., dur = `dur (s)`,
#             phonTime = `phonationtime (s)`,
#             speechRate = `speechrate (nsyll/dur)`,
#             artRate = `articulation rate (nsyll / phonationtime)`,
#             asd = `ASD (speakingtime/nsyll)`)
# 
# glimpse(sr_df)
# 
# mean(sr_df$speechRate)
# sd(sr_df$speechRate)
# 
# mean(sr_df$dur)
# sd(sr_df$dur)
# 
# ggplot(sr_df, aes(x = speechRate, y = asd)) +
#   geom_point() +
#   geom_smooth()
