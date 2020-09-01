# load packages
library(plyr)
library(tidyr)
library(tidyverse)

# load data
stress_50 <- read.delim("./data/stress_50bin.txt")

# Check gaze fixation columns have different values
unique(stress_50$AVERAGE_IA_1_SAMPLE_COUNT)  # looking at target according to IA_#_ID
unique(stress_50$AVERAGE_IA_2_SAMPLE_COUNT)  # looking at distractor
unique(stress_50$AVERAGE_IA_0_SAMPLE_COUNT)  # elsewhere

# create variable group
ss_50 <- stress_50 %>%
  separate(., col = RECORDING_SESSION_LABEL,
           into = c("group", "group_member"),
           sep = 3,
           remove = FALSE) %>%
  
  #select and rename variables of interest
  select(., RECORDING_SESSION_LABEL, TRIAL_INDEX, BIN_INDEX,
         AVERAGE_IA_0_SAMPLE_COUNT, AVERAGE_IA_0_SAMPLE_COUNT_.,
         AVERAGE_IA_1_SAMPLE_COUNT, AVERAGE_IA_1_SAMPLE_COUNT_.,
         AVERAGE_IA_2_SAMPLE_COUNT, AVERAGE_IA_2_SAMPLE_COUNT_.,
         ACCURACY, RT, block, cond, 
         id, lex_freq, phonot_freq, 
         t01, t02, t03, t04, t05, t06, t07, target, version, group) %>%
  dplyr::rename(., participant = RECORDING_SESSION_LABEL,
         trial = TRIAL_INDEX, 
         bin = BIN_INDEX,
         target_count = AVERAGE_IA_1_SAMPLE_COUNT, 
         target_prop = AVERAGE_IA_1_SAMPLE_COUNT_.,
         offset_prev_word = t01,
         onset_v1 = t02,
         onset_c2 = t03,
         onset_c3 = t04,
         onset_v2 = t05,
         offset_target = t06,
         endSentence = t07,
         sentence_id = id) %>%
  
  # remove incorrect
  filter(., ACCURACY == 1) %>%
  
  # drop unused levels of factors
  droplevels(.)

# Adjust processing time in bin columns
ss_50$bin_adj <- ss_50$bin - 4

write.csv(ss_50, "data/clean/stress50clean.csv")

















# center time course
# Center timecourse 
cortaDF <- filter(ss_50, target == "corta", sentence_id == "S01_C1")
cortaDF$binN <- 999
sufONcorta <- round(unique(cortaDF$onset_c3) / 50)
cortaDF[cortaDF$bin == sufONcorta, 'binN'] <- 0
cortaDF$binN <- cortaDF$bin - sufONcorta
cortaDF$condToken <- 1

cortoDF <- filter(ss_50, target == "cortÃ³", sentence_id == "S01_C2")
cortoDF$binN <- 999
sufONcorto <- round(unique(cortoDF$onset_c3) / 50)
cortoDF[cortoDF$bin == sufONcorto, 'binN'] <- 0
cortoDF$binN <- cortoDF$bin - sufONcorto
cortoDF$condToken <- 1

pintaDF <- filter(ss_50, target == "pinta", sentence_id == "S02_C1")
pintaDF$binN <- 999
sufONpinta <- round(unique(pintaDF$onset_c3) / 50)
pintaDF[pintaDF$bin == sufONpinta, 'binN'] <- 0
pintaDF$binN <- pintaDF$bin - sufONpinta
pintaDF$condToken <- 1

pintoDF <- filter(ss_50, target == "pintÃ³", sentence_id == "S02_C2")
pintoDF$binN <- 999
sufONpinto <- round(unique(pintoDF$onset_c3) / 50)
pintoDF[pintoDF$bin == sufONpinto, 'binN'] <- 0
pintoDF$binN <- pintoDF$bin - sufONpinto
pintoDF$condToken <- 1

firmaDF <- filter(ss_50, target == "firma", sentence_id == "S03_C1")
firmaDF$binN <- 999
sufONfirma <- round(unique(firmaDF$onset_c3) / 50)
firmaDF[firmaDF$bin == sufONfirma, 'binN'] <- 0
firmaDF$binN <- firmaDF$bin - sufONfirma
firmaDF$condToken <- 1

firmoDF <- filter(ss_50, target == "firmÃ³", sentence_id == "S03_C2")
firmoDF$binN <- 999
sufONfirmo <- round(unique(firmoDF$onset_c3) / 50)
firmoDF[firmoDF$bin == sufONfirmo, 'binN'] <- 0
firmoDF$binN <- firmoDF$bin - sufONfirmo
firmoDF$condToken <- 1

planchaDF <- filter(ss_50, target == "plancha", sentence_id == "S04_C1")
planchaDF$binN <- 999
sufONplancha <- round(unique(planchaDF$onset_c3) / 50)
planchaDF[planchaDF$bin == sufONplancha, 'binN'] <- 0
planchaDF$binN <- planchaDF$bin - sufONplancha
planchaDF$condToken <- 1

planchoDF <- filter(ss_50, target == "planchÃ³", sentence_id == "S04_C2")
planchoDF$binN <- 999
sufONplancho <- round(unique(planchoDF$onset_c3) / 50)
planchoDF[planchoDF$bin == sufONplancho, 'binN'] <- 0
planchoDF$binN <- planchoDF$bin - sufONplancho
planchoDF$condToken <- 1

cantaDF <- filter(ss_50, target == "canta", sentence_id == "S05_C1")
cantaDF$binN <- 999
sufONcanta <- round(unique(cantaDF$onset_c3) / 50)
cantaDF[cantaDF$bin == sufONcanta, 'binN'] <- 0
cantaDF$binN <- cantaDF$bin - sufONcanta
cantaDF$condToken <- 1

cantoDF <- filter(ss_50, target == "cantÃ³", sentence_id == "S05_C2")
cantoDF$binN <- 999
sufONcanto <- round(unique(cantoDF$onset_c3) / 50)
cantoDF[cantoDF$bin == sufONcanto, 'binN'] <- 0
cantoDF$binN <- cantoDF$bin - sufONcanto
cantoDF$condToken <- 1

buscaDF <- filter(ss_50, target == "busca", sentence_id == "S06_C1")
buscaDF$binN <- 999
sufONbusca <- round(unique(buscaDF$onset_c3) / 50)
buscaDF[buscaDF$bin == sufONbusca, 'binN'] <- 0
buscaDF$binN <- buscaDF$bin - sufONbusca
buscaDF$condToken <- 1

buscoDF <- filter(ss_50, target == "buscÃ³", sentence_id == "S06_C2")
buscoDF$binN <- 999
sufONbusco <- round(unique(buscoDF$onset_c3) / 50)
buscoDF[buscoDF$bin == sufONbusco, 'binN'] <- 0
buscoDF$binN <- buscoDF$bin - sufONbusco
buscoDF$condToken <- 1

compraDF <- filter(ss_50, target == "compra", sentence_id == "S07_C1")
compraDF$binN <- 999
sufONcompra <- round(unique(compraDF$onset_c3) / 50)
compraDF[compraDF$bin == sufONcompra, 'binN'] <- 0
compraDF$binN <- compraDF$bin - sufONcompra
compraDF$condToken <- 1

comproDF <- filter(ss_50, target == "comprÃ³", sentence_id == "S07_C2")
comproDF$binN <- 999
sufONcompro <- round(unique(comproDF$onset_c3) / 50)
comproDF[comproDF$bin == sufONcompro, 'binN'] <- 0
comproDF$binN <- comproDF$bin - sufONcompro
comproDF$condToken <- 1

gastaDF <- filter(ss_50, target == "gasta", sentence_id == "S08_C1")
gastaDF$binN <- 999
sufONgasta <- round(unique(gastaDF$onset_c3) / 50)
gastaDF[gastaDF$bin == sufONgasta, 'binN'] <- 0
gastaDF$binN <- gastaDF$bin - sufONgasta
gastaDF$condToken <- 1

gastoDF <- filter(ss_50, target == "gastÃ³", sentence_id == "S08_C2")
gastoDF$binN <- 999
sufONgasto <- round(unique(gastoDF$onset_c3) / 50)
gastoDF[gastoDF$bin == sufONgasto, 'binN'] <- 0
gastoDF$binN <- gastoDF$bin - sufONgasto
gastoDF$condToken <- 1

juntaDF <- filter(ss_50, target == "junta", sentence_id == "S09_C1")
juntaDF$binN <- 999
sufONjunta <- round(unique(juntaDF$onset_c3) / 50)
juntaDF[juntaDF$bin == sufONjunta, 'binN'] <- 0
juntaDF$binN <- juntaDF$bin - sufONjunta
juntaDF$condToken <- 1

juntoDF <- filter(ss_50, target == "juntÃ³", sentence_id == "S09_C2")
juntoDF$binN <- 999
sufONjunto <- round(unique(juntoDF$onset_c3) / 50)
juntoDF[juntoDF$bin == sufONjunto, 'binN'] <- 0
juntoDF$binN <- juntoDF$bin - sufONjunto
juntoDF$condToken <- 1

salvaDF <- filter(ss_50, target == "salva", sentence_id == "S10_C1")
salvaDF$binN <- 999
sufONsalva <- round(unique(salvaDF$onset_c3) / 50)
salvaDF[salvaDF$bin == sufONsalva, 'binN'] <- 0
salvaDF$binN <- salvaDF$bin - sufONsalva
salvaDF$condToken <- 1

salvoDF <- filter(ss_50, target == "salvÃ³", sentence_id == "S10_C2")
salvoDF$binN <- 999
sufONsalvo <- round(unique(salvoDF$onset_c3) / 50)
salvoDF[salvoDF$bin == sufONsalvo, 'binN'] <- 0
salvoDF$binN <- salvoDF$bin - sufONsalvo
salvoDF$condToken <- 1

lanzaDF <- filter(ss_50, target == "lanza", sentence_id == "S11_C1")
lanzaDF$binN <- 999
sufONlanza <- round(unique(lanzaDF$onset_c3) / 50)
lanzaDF[lanzaDF$bin == sufONlanza, 'binN'] <- 0
lanzaDF$binN <- lanzaDF$bin - sufONlanza
lanzaDF$condToken <- 1

lanzoDF <- filter(ss_50, target == "lanzÃ³", sentence_id == "S11_C2")
lanzoDF$binN <- 999
sufONlanzo <- round(unique(lanzoDF$onset_c3) / 50)
lanzoDF[lanzoDF$bin == sufONlanzo, 'binN'] <- 0
lanzoDF$binN <- lanzoDF$bin - sufONlanzo
lanzoDF$condToken <- 1

pescaDF <- filter(ss_50, target == "pesca", sentence_id == "S12_C1")
pescaDF$binN <- 999
sufONpesca <- round(unique(pescaDF$onset_c3) / 50)
pescaDF[pescaDF$bin == sufONpesca, 'binN'] <- 0
pescaDF$binN <- pescaDF$bin - sufONpesca
pescaDF$condToken <- 1

pescoDF <- filter(ss_50, target == "pescÃ³", sentence_id == "S12_C2")
pescoDF$binN <- 999
sufONpesco <- round(unique(pescoDF$onset_c3) / 50)
pescoDF[pescoDF$bin == sufONpesco, 'binN'] <- 0
pescoDF$binN <- pescoDF$bin - sufONpesco
pescoDF$condToken <- 1

montaDF <- filter(ss_50, target == "monta", sentence_id == "S13_C1")
montaDF$binN <- 999
sufONmonta <- round(unique(montaDF$onset_c3) / 50)
montaDF[montaDF$bin == sufONmonta, 'binN'] <- 0
montaDF$binN <- montaDF$bin - sufONmonta
montaDF$condToken <- 1

montoDF <- filter(ss_50, target == "montÃ³", sentence_id == "S13_C2")
montoDF$binN <- 999
sufONmonto <- round(unique(montoDF$onset_c3) / 50)
montoDF[montoDF$bin == sufONmonto, 'binN'] <- 0
montoDF$binN <- montoDF$bin - sufONmonto
montoDF$condToken <- 1

mandaDF <- filter(ss_50, target == "manda", sentence_id == "S14_C1")
mandaDF$binN <- 999
sufONmanda <- round(unique(mandaDF$onset_c3) / 50)
mandaDF[mandaDF$bin == sufONmanda, 'binN'] <- 0
mandaDF$binN <- mandaDF$bin - sufONmanda
mandaDF$condToken <- 1

mandoDF <- filter(ss_50, target == "mandÃ³", sentence_id == "S14_C2")
mandoDF$binN <- 999
sufONmando <- round(unique(mandoDF$onset_c3) / 50)
mandoDF[mandoDF$bin == sufONmando, 'binN'] <- 0
mandoDF$binN <- mandoDF$bin - sufONmando
mandoDF$condToken <- 1

saltaDF <- filter(ss_50, target == "salta", sentence_id == "S15_C1")
saltaDF$binN <- 999
sufONsalta <- round(unique(saltaDF$onset_c3) / 50)
saltaDF[saltaDF$bin == sufONsalta, 'binN'] <- 0
saltaDF$binN <- saltaDF$bin - sufONsalta
saltaDF$condToken <- 1

saltoDF <- filter(ss_50, target == "saltÃ³", sentence_id == "S15_C2")
saltoDF$binN <- 999
sufONsalto <- round(unique(saltoDF$onset_c3) / 50)
saltoDF[saltoDF$bin == sufONsalto, 'binN'] <- 0
saltoDF$binN <- saltoDF$bin - sufONsalto
saltoDF$condToken <- 1

cargaDF <- filter(ss_50, target == "carga", sentence_id == "S16_C1")
cargaDF$binN <- 999
sufONcarga <- round(unique(cargaDF$onset_c3) / 50)
cargaDF[cargaDF$bin == sufONcarga, 'binN'] <- 0
cargaDF$binN <- cargaDF$bin - sufONcarga
cargaDF$condToken <- 1

cargoDF <- filter(ss_50, target == "cargÃ³", sentence_id == "S16_C2")
cargoDF$binN <- 999
sufONcargo <- round(unique(cargoDF$onset_c3) / 50)
cargoDF[cargoDF$bin == sufONcargo, 'binN'] <- 0
cargoDF$binN <- cargoDF$bin - sufONcargo
cargoDF$condToken <- 1

df_related <- do.call("bind_rows", list(cortaDF, cortoDF, buscaDF, buscoDF,
                                        firmaDF, firmoDF, mandaDF, mandoDF,
                                        montaDF, montoDF, cantaDF, cantoDF,
                                        compraDF, comproDF, cargaDF, cargoDF,
                                        juntaDF, juntoDF, salvaDF, salvoDF,
                                        gastaDF, gastoDF, lanzaDF, lanzoDF,
                                        saltaDF, saltoDF, pescaDF, pescoDF,
                                        pintaDF, pintoDF, planchaDF, planchoDF))


write.csv(df_related,"data/clean/timecourse_50_clean.csv", row.names = FALSE)



# find data row for bin at C3 onset (verb structure = C1V1C2.C3V2)
bin_S01_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S01_C1", "onset_c3"])) / 50) 
bin_S01_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S01_C2", "onset_c3"])) / 50)

st_S01c3 <- filter(ss_50, sentence_id == "S01_C1", bin_adj == as.numeric(bin_S01_stc3))
un_S01c3 <- filter(ss_50, sentence_id == "S01_C2", bin_adj == as.numeric(bin_S01_unc3)) 



bin_S02_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S02_C1", "onset_c3"])) / 50) 
bin_S02_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S02_C2", "onset_c3"])) / 50)

st_S02c3 <- filter(ss_50, sentence_id == "S02_C1", bin_adj == as.numeric(bin_S02_stc3))
un_S02c3 <- filter(ss_50, sentence_id == "S02_C2", bin_adj == as.numeric(bin_S02_unc3))



bin_S03_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S03_C1", "onset_c3"])) / 50) 
bin_S03_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S03_C2", "onset_c3"])) / 50)

st_S03c3 <- filter(ss_50, sentence_id == "S03_C1", bin_adj == as.numeric(bin_S03_stc3))
un_S03c3 <- filter(ss_50, sentence_id == "S03_C2", bin_adj == as.numeric(bin_S03_unc3))



bin_S04_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S04_C1", "onset_c3"])) / 50) 
bin_S04_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S04_C2", "onset_c3"])) / 50)

st_S04c3 <- filter(ss_50, sentence_id == "S04_C1", bin_adj == as.numeric(bin_S04_stc3))
un_S04c3 <- filter(ss_50, sentence_id == "S04_C2", bin_adj == as.numeric(bin_S04_unc3))



bin_S05_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S05_C1", "onset_c3"])) / 50) 
bin_S05_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S05_C2", "onset_c3"])) / 50)

st_S05c3 <- filter(ss_50, sentence_id == "S05_C1", bin_adj == as.numeric(bin_S05_stc3))
un_S05c3 <- filter(ss_50, sentence_id == "S05_C2", bin_adj == as.numeric(bin_S05_unc3))



bin_S06_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S06_C1", "onset_c3"])) / 50) 
bin_S06_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S06_C2", "onset_c3"])) / 50)

st_S06c3 <- filter(ss_50, sentence_id == "S06_C1", bin_adj == as.numeric(bin_S06_stc3))
un_S06c3 <- filter(ss_50, sentence_id == "S06_C2", bin_adj == as.numeric(bin_S06_unc3))



bin_S07_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S07_C1", "onset_c3"])) / 50) 
bin_S07_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S07_C2", "onset_c3"])) / 50)

st_S07c3 <- filter(ss_50, sentence_id == "S07_C1", bin_adj == as.numeric(bin_S07_stc3))
un_S07c3 <- filter(ss_50, sentence_id == "S07_C2", bin_adj == as.numeric(bin_S07_unc3))



bin_S08_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S08_C1", "onset_c3"])) / 50) 
bin_S08_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S08_C2", "onset_c3"])) / 50)

st_S08c3 <- filter(ss_50, sentence_id == "S08_C1", bin_adj == as.numeric(bin_S08_stc3))
un_S08c3 <- filter(ss_50, sentence_id == "S08_C2", bin_adj == as.numeric(bin_S08_unc3))



bin_S09_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S09_C1", "onset_c3"])) / 50) 
bin_S09_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S09_C2", "onset_c3"])) / 50)

st_S09c3 <- filter(ss_50, sentence_id == "S09_C1", bin_adj == as.numeric(bin_S09_stc3))
un_S09c3 <- filter(ss_50, sentence_id == "S09_C2", bin_adj == as.numeric(bin_S09_unc3))



bin_S10_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S10_C1", "onset_c3"])) / 50) 
bin_S10_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S10_C2", "onset_c3"])) / 50)

st_S10c3 <- filter(ss_50, sentence_id == "S10_C1", bin_adj == as.numeric(bin_S10_stc3))
un_S10c3 <- filter(ss_50, sentence_id == "S10_C2", bin_adj == as.numeric(bin_S10_unc3))



bin_S11_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S11_C1", "onset_c3"])) / 50) 
bin_S11_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S11_C2", "onset_c3"])) / 50)

st_S11c3 <- filter(ss_50, sentence_id == "S11_C1", bin_adj == as.numeric(bin_S11_stc3))
un_S11c3 <- filter(ss_50, sentence_id == "S11_C2", bin_adj == as.numeric(bin_S11_unc3))



bin_S12_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S12_C1", "onset_c3"])) / 50) 
bin_S12_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S12_C2", "onset_c3"])) / 50)

st_S12c3 <- filter(ss_50, sentence_id == "S12_C1", bin_adj == as.numeric(bin_S12_stc3))
un_S12c3 <- filter(ss_50, sentence_id == "S12_C2", bin_adj == as.numeric(bin_S12_unc3))



bin_S13_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S13_C1", "onset_c3"])) / 50) 
bin_S13_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S13_C2", "onset_c3"])) / 50)

st_S13c3 <- filter(ss_50, sentence_id == "S13_C1", bin_adj == as.numeric(bin_S13_stc3))
un_S13c3 <- filter(ss_50, sentence_id == "S13_C2", bin_adj == as.numeric(bin_S13_unc3))



bin_S14_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S14_C1", "onset_c3"])) / 50) 
bin_S14_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S14_C2", "onset_c3"])) / 50)

st_S14c3 <- filter(ss_50, sentence_id == "S14_C1", bin_adj == as.numeric(bin_S14_stc3))
un_S14c3 <- filter(ss_50, sentence_id == "S14_C2", bin_adj == as.numeric(bin_S14_unc3))



bin_S15_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S15_C1", "onset_c3"])) / 50) 
bin_S15_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S15_C2", "onset_c3"])) / 50)

st_S15c3 <- filter(ss_50, sentence_id == "S15_C1", bin_adj == as.numeric(bin_S15_stc3))
un_S15c3 <- filter(ss_50, sentence_id == "S15_C2", bin_adj == as.numeric(bin_S15_unc3))



bin_S16_stc3 <- round((unique(ss_50[ss_50$sentence_id == "S16_C1", "onset_c3"])) / 50) 
bin_S16_unc3 <- round((unique(ss_50[ss_50$sentence_id == "S16_C2", "onset_c3"])) / 50)

st_S16c3 <- filter(ss_50, sentence_id == "S16_C1", bin_adj == as.numeric(bin_S16_stc3))
un_S16c3 <- filter(ss_50, sentence_id == "S16_C2", bin_adj == as.numeric(bin_S16_unc3))

# bind individual rows into one single dataframe
st50_onC3 <- rbind(st_S01c3, un_S01c3, st_S02c3, un_S02c3, st_S03c3, un_S03c3, st_S04c3, un_S04c3,
                   st_S05c3, un_S05c3, st_S06c3, un_S06c3, st_S07c3, un_S07c3, st_S08c3, un_S08c3,
                   st_S09c3, un_S09c3, st_S10c3, un_S10c3, st_S12c3, un_S12c3, st_S12c3, un_S12c3,
                   st_S13c3, un_S13c3, st_S14c3, un_S14c3, st_S15c3, un_S15c3, st_S16c3, un_S16c3)


write.csv(st50_onC3, "data/clean/fixations_onC3_50.csv")

# Bonferroni correction for alpha level
# Divided by 5 because using same data for multiple tests (5 different groups)
print(alphaAdj <- 0.05/10)
# 0.005

# One Sample t-test
# proportion of fixations on target by monolinguals when tense = present
mon_pres_c3 <- st50_onC3 %>% 
  filter(., group == "mon", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 4.3768, df = 29, p-value = 7.133e-05

# proportion of fixations on target by monolinguals when tense = preterit
mon_pret_c3 <- st50_onC3 %>% 
  filter(., group == "mon", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 4.9592, df = 29, p-value = 1.42e-05


# proportion of fixations on target by intermediate EN when tense = present
ies_pres_c3 <- st50_onC3 %>% 
  filter(., group == "ies", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.8065, df = 13, p-value = 0.9926


# proportion of fixations on target by intermediate EN when tense = preterit
ies_pret_c3 <- st50_onC3 %>% 
  filter(., group == "ies", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 1.6859, df = 13, p-value = 0.05783


# proportion of fixations on target by advanced EN when tense = present
aes_pres_c3 <- st50_onC3 %>% 
  filter(., group == "aes", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.8624, df = 20, p-value = 0.9952


# proportion of fixations on target by advanced EN when tense = preterit
aes_pret_c3 <- st50_onC3 %>% 
  filter(., group == "aes", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.27751, df = 20, p-value = 0.3921



# proportion of fixations on target by intermediate CH when tense = present
ims_pres_c3 <- st50_onC3 %>% 
  filter(., group == "ims", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.4941, df = 13, p-value = 0.9866


# proportion of fixations on target by intermediate CH when tense = preterit
ims_pret_c3 <- st50_onC3 %>% 
  filter(., group == "ims", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.42373, df = 13, p-value = 0.3393


# proportion of fixations on target by advanced CH when tense = present
ams_pres_c3 <- st50_onC3 %>% 
  filter(., group == "ams", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pres_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = -1.5285, df = 18, p-value = 0.9281


# proportion of fixations on target by advanced CH when tense = preterit
ams_pret_c3 <- st50_onC3 %>% 
  filter(., group == "ams", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pret_c3$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.32532, df = 18, p-value = 0.3743
















################################################################################
s50 <- ss_50

## create function to find time marker
# to debug options(error=recover)
on_c2 <- ldply(s10, function(s10) {
  for (element in s50$sentence_id) {
    bin_mark_i <- round(unique(s10[s10[sentence_id == i], "onset_c3"]) / 10)
    
    bin_sum_i <- filter(s10, s10[sentence_id == i], bin_adj == as.numeric(bin_mark_i))
  }
}
)



function(x, y) {
  
  bin_mark <- round(unique(df[df$x, y]) / 10)
  
  bin_sum_x <- filter(df, x, bin_adj == bin_num)
  
  rbind(bin_sum_x)
}



#####
# Gaze fixations at segmental disambiguator onset
bin_S01_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S01_C1", "onset_v2"])) / 10) 
bin_S01_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S01_C2", "onset_v2"])) / 10)

st_S01v2 <- filter(ss_10, sentence_id == "S01_C1", bin_adj == as.numeric(bin_S01_stv2))
un_S01v2 <- filter(ss_10, sentence_id == "S01_C2", bin_adj == as.numeric(bin_S01_unv2)) 



bin_S02_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S02_C1", "onset_v2"])) / 10) 
bin_S02_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S02_C2", "onset_v2"])) / 10)

st_S02v2 <- filter(ss_10, sentence_id == "S02_C1", bin_adj == as.numeric(bin_S02_stv2))
un_S02v2 <- filter(ss_10, sentence_id == "S02_C2", bin_adj == as.numeric(bin_S02_unv2))



bin_S03_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S03_C1", "onset_v2"])) / 10) 
bin_S03_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S03_C2", "onset_v2"])) / 10)

st_S03v2 <- filter(ss_10, sentence_id == "S03_C1", bin_adj == as.numeric(bin_S03_stv2))
un_S03v2 <- filter(ss_10, sentence_id == "S03_C2", bin_adj == as.numeric(bin_S03_unv2))



bin_S04_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S04_C1", "onset_v2"])) / 10) 
bin_S04_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S04_C2", "onset_v2"])) / 10)

st_S04v2 <- filter(ss_10, sentence_id == "S04_C1", bin_adj == as.numeric(bin_S04_stv2))
un_S04v2 <- filter(ss_10, sentence_id == "S04_C2", bin_adj == as.numeric(bin_S04_unv2))



bin_S05_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S05_C1", "onset_v2"])) / 10) 
bin_S05_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S05_C2", "onset_v2"])) / 10)

st_S05v2 <- filter(ss_10, sentence_id == "S05_C1", bin_adj == as.numeric(bin_S05_stv2))
un_S05v2 <- filter(ss_10, sentence_id == "S05_C2", bin_adj == as.numeric(bin_S05_unv2))



bin_S06_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S06_C1", "onset_v2"])) / 10) 
bin_S06_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S06_C2", "onset_v2"])) / 10)

st_S06v2 <- filter(ss_10, sentence_id == "S06_C1", bin_adj == as.numeric(bin_S06_stv2))
un_S06v2 <- filter(ss_10, sentence_id == "S06_C2", bin_adj == as.numeric(bin_S06_unv2))



bin_S07_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S07_C1", "onset_v2"])) / 10) 
bin_S07_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S07_C2", "onset_v2"])) / 10)

st_S07v2 <- filter(ss_10, sentence_id == "S07_C1", bin_adj == as.numeric(bin_S07_stv2))
un_S07v2 <- filter(ss_10, sentence_id == "S07_C2", bin_adj == as.numeric(bin_S07_unv2))



bin_S08_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S08_C1", "onset_v2"])) / 10) 
bin_S08_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S08_C2", "onset_v2"])) / 10)

st_S08v2 <- filter(ss_10, sentence_id == "S08_C1", bin_adj == as.numeric(bin_S08_stv2))
un_S08v2 <- filter(ss_10, sentence_id == "S08_C2", bin_adj == as.numeric(bin_S08_unv2))



bin_S09_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S09_C1", "onset_v2"])) / 10) 
bin_S09_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S09_C2", "onset_v2"])) / 10)

st_S09v2 <- filter(ss_10, sentence_id == "S09_C1", bin_adj == as.numeric(bin_S09_stv2))
un_S09v2 <- filter(ss_10, sentence_id == "S09_C2", bin_adj == as.numeric(bin_S09_unv2))



bin_S10_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S10_C1", "onset_v2"])) / 10) 
bin_S10_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S10_C2", "onset_v2"])) / 10)

st_S10v2 <- filter(ss_10, sentence_id == "S10_C1", bin_adj == as.numeric(bin_S10_stv2))
un_S10v2 <- filter(ss_10, sentence_id == "S10_C2", bin_adj == as.numeric(bin_S10_unv2))



bin_S11_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S11_C1", "onset_v2"])) / 10) 
bin_S11_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S11_C2", "onset_v2"])) / 10)

st_S11v2 <- filter(ss_10, sentence_id == "S11_C1", bin_adj == as.numeric(bin_S11_stv2))
un_S11v2 <- filter(ss_10, sentence_id == "S11_C2", bin_adj == as.numeric(bin_S11_unv2))



bin_S12_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S12_C1", "onset_v2"])) / 10) 
bin_S12_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S12_C2", "onset_v2"])) / 10)

st_S12v2 <- filter(ss_10, sentence_id == "S12_C1", bin_adj == as.numeric(bin_S12_stv2))
un_S12v2 <- filter(ss_10, sentence_id == "S12_C2", bin_adj == as.numeric(bin_S12_unv2))



bin_S13_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S13_C1", "onset_v2"])) / 10) 
bin_S13_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S13_C2", "onset_v2"])) / 10)

st_S13v2 <- filter(ss_10, sentence_id == "S13_C1", bin_adj == as.numeric(bin_S13_stv2))
un_S13v2 <- filter(ss_10, sentence_id == "S13_C2", bin_adj == as.numeric(bin_S13_unv2))



bin_S14_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S14_C1", "onset_v2"])) / 10) 
bin_S14_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S14_C2", "onset_v2"])) / 10)

st_S14v2 <- filter(ss_10, sentence_id == "S14_C1", bin_adj == as.numeric(bin_S14_stv2))
un_S14v2 <- filter(ss_10, sentence_id == "S14_C2", bin_adj == as.numeric(bin_S14_unv2))



bin_S15_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S15_C1", "onset_v2"])) / 10) 
bin_S15_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S15_C2", "onset_v2"])) / 10)

st_S15v2 <- filter(ss_10, sentence_id == "S15_C1", bin_adj == as.numeric(bin_S15_stv2))
un_S15v2 <- filter(ss_10, sentence_id == "S15_C2", bin_adj == as.numeric(bin_S15_unv2))



bin_S16_stv2 <- round((unique(ss_10[ss_10$sentence_id == "S16_C1", "onset_v2"])) / 10) 
bin_S16_unv2 <- round((unique(ss_10[ss_10$sentence_id == "S16_C2", "onset_v2"])) / 10)

st_S16v2 <- filter(ss_10, sentence_id == "S16_C1", bin_adj == as.numeric(bin_S16_stv2))
un_S16v2 <- filter(ss_10, sentence_id == "S16_C2", bin_adj == as.numeric(bin_S16_unv2))

# bind individual rows into one single dataframe
st10_onv2 <- rbind(st_S01v2, un_S01v2, st_S02v2, un_S02v2, st_S03v2, un_S03v2, st_S04v2, un_S04v2,
                   st_S05v2, un_S05v2, st_S06v2, un_S06v2, st_S07v2, un_S07v2, st_S08v2, un_S08v2,
                   st_S09v2, un_S09v2, st_S10v2, un_S10v2, st_S12v2, un_S12v2, st_S12v2, un_S12v2,
                   st_S13v2, un_S13v2, st_S14v2, un_S14v2, st_S15v2, un_S15v2, st_S16v2, un_S16v2)


# Bonferroni correction for alpha level
# Divided by 5 because using same data for multiple tests (5 different groups)
print(alphaAdj <- 0.05/10)
# 0.005

# One Sample t-test
# proportion of fixations on target by monolinguals when tense = present
mon_pres_v2 <- st10_onv2 %>% 
  filter(., group == "mon", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 7.8217, df = 29, p-value = 6.303e-09

# proportion of fixations on target by monolinguals when tense = preterit
mon_pret_v2 <- st10_onv2 %>% 
  filter(., group == "mon", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(mon_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 5.5193, df = 29, p-value = 2.998e-06


# proportion of fixations on target by intermediate EN when tense = present
ies_pres_v2 <- st10_onv2 %>% 
  filter(., group == "ies", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.6194, df = 13, p-value = 0.9894


# proportion of fixations on target by intermediate EN when tense = preterit
ies_pret_v2 <- st10_onv2 %>% 
  filter(., group == "ies", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ies_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.85843, df = 13, p-value = 0.2031


# proportion of fixations on target by advanced EN when tense = present
aes_pres_v2 <- st10_onv2 %>% 
  filter(., group == "aes", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = -2.4616, df = 20, p-value = 0.9885


# proportion of fixations on target by advanced EN when tense = preterit
aes_pret_v2 <- st10_onv2 %>% 
  filter(., group == "aes", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(aes_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 0.93483, df = 20, p-value = 0.1805



# proportion of fixations on target by intermediate CH when tense = present
ims_pres_v2 <- st10_onv2 %>% 
  filter(., group == "ims", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = -1.569, df = 13, p-value = 0.9297


# proportion of fixations on target by intermediate CH when tense = preterit
ims_pret_v2 <- st10_onv2 %>% 
  filter(., group == "ims", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ims_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 1.3879, df = 13, p-value = 0.09424


# proportion of fixations on target by advanced CH when tense = present
ams_pres_v2 <- st10_onv2 %>% 
  filter(., group == "ams", cond == "1") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pres_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = -0.90233, df = 18, p-value = 0.8106


# proportion of fixations on target by advanced CH when tense = preterit
ams_pret_v2 <- st10_onv2 %>% 
  filter(., group == "ams", cond == "2") %>% 
  group_by(., participant) %>% 
  summarize(., mean_prop = mean(target_prop))

t.test(ams_pret_v2$mean_prop, alternative = "greater", mu = 0.5)
# t = 1.1074, df = 18, p-value = 0.1414







