# load packages
library(Hmisc)
library(plyr)
library(tidyr)
library(tidyverse)

# load data
stress_10 <- read_tsv("data/stress_10.txt")

# create variable group
ss_10 <- stress_10 %>%
  separate(., col = RECORDING_SESSION_LABEL,
           into = c("group", "group_member"),
           sep = 3,
           remove = FALSE) %>%
  
  #select and rename variables of interest
  select(., RECORDING_SESSION_LABEL, TRIAL_INDEX, BIN_INDEX,
         AVERAGE_IA_0_SAMPLE_COUNT, `AVERAGE_IA_0_SAMPLE_COUNT_%`,
         AVERAGE_IA_1_SAMPLE_COUNT, `AVERAGE_IA_1_SAMPLE_COUNT_%`,
         AVERAGE_IA_2_SAMPLE_COUNT, `AVERAGE_IA_2_SAMPLE_COUNT_%`,
         IA_1_ID, IA_2_ID, IA_0_ID, ACCURACY, RT, cond, 
         distractor, id, lex_freq, phonot_freq, sentence, 
         t01, t02, t03, t04, t05, t06, t07, target, version, group) %>%
  dplyr::rename(., participant = RECORDING_SESSION_LABEL,
                trial = TRIAL_INDEX, 
                bin = BIN_INDEX,
                target_count = AVERAGE_IA_1_SAMPLE_COUNT, 
                target_prop = `AVERAGE_IA_1_SAMPLE_COUNT_%`,
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
ss_10$bin_adj <- ss_10$bin - 20

# Center timecourse 
cortaDF <- filter(ss_10, target == "corta", sentence_id == "S01_C1")
cortaDF$binN <- 999
sufONcorta <- round(unique(cortaDF$onset_c3) / 10)
cortaDF[cortaDF$bin == sufONcorta, 'binN'] <- 0
cortaDF$binN <- cortaDF$bin - sufONcorta
cortaDF$condToken <- 1

cortoDF <- filter(ss_10, target == "cortÃ³", sentence_id == "S01_C2")
cortoDF$binN <- 999
sufONcorto <- round(unique(cortoDF$onset_c3) / 10)
cortoDF[cortoDF$bin == sufONcorto, 'binN'] <- 0
cortoDF$binN <- cortoDF$bin - sufONcorto
cortoDF$condToken <- 1

pintaDF <- filter(ss_10, target == "pinta", sentence_id == "S02_C1")
pintaDF$binN <- 999
sufONpinta <- round(unique(pintaDF$onset_c3) / 10)
pintaDF[pintaDF$bin == sufONpinta, 'binN'] <- 0
pintaDF$binN <- pintaDF$bin - sufONpinta
pintaDF$condToken <- 1

pintoDF <- filter(ss_10, target == "pintÃ³", sentence_id == "S02_C2")
pintoDF$binN <- 999
sufONpinto <- round(unique(pintoDF$onset_c3) / 10)
pintoDF[pintoDF$bin == sufONpinto, 'binN'] <- 0
pintoDF$binN <- pintoDF$bin - sufONpinto
pintoDF$condToken <- 1

firmaDF <- filter(ss_10, target == "firma", sentence_id == "S03_C1")
firmaDF$binN <- 999
sufONfirma <- round(unique(firmaDF$onset_c3) / 10)
firmaDF[firmaDF$bin == sufONfirma, 'binN'] <- 0
firmaDF$binN <- firmaDF$bin - sufONfirma
firmaDF$condToken <- 1

firmoDF <- filter(ss_10, target == "firmÃ³", sentence_id == "S03_C2")
firmoDF$binN <- 999
sufONfirmo <- round(unique(firmoDF$onset_c3) / 10)
firmoDF[firmoDF$bin == sufONfirmo, 'binN'] <- 0
firmoDF$binN <- firmoDF$bin - sufONfirmo
firmoDF$condToken <- 1

planchaDF <- filter(ss_10, target == "plancha", sentence_id == "S04_C1")
planchaDF$binN <- 999
sufONplancha <- round(unique(planchaDF$onset_c3) / 10)
planchaDF[planchaDF$bin == sufONplancha, 'binN'] <- 0
planchaDF$binN <- planchaDF$bin - sufONplancha
planchaDF$condToken <- 1

planchoDF <- filter(ss_10, target == "planchÃ³", sentence_id == "S04_C2")
planchoDF$binN <- 999
sufONplancho <- round(unique(planchoDF$onset_c3) / 10)
planchoDF[planchoDF$bin == sufONplancho, 'binN'] <- 0
planchoDF$binN <- planchoDF$bin - sufONplancho
planchoDF$condToken <- 1

cantaDF <- filter(ss_10, target == "canta", sentence_id == "S05_C1")
cantaDF$binN <- 999
sufONcanta <- round(unique(cantaDF$onset_c3) / 10)
cantaDF[cantaDF$bin == sufONcanta, 'binN'] <- 0
cantaDF$binN <- cantaDF$bin - sufONcanta
cantaDF$condToken <- 1

cantoDF <- filter(ss_10, target == "cantÃ³", sentence_id == "S05_C2")
cantoDF$binN <- 999
sufONcanto <- round(unique(cantoDF$onset_c3) / 10)
cantoDF[cantoDF$bin == sufONcanto, 'binN'] <- 0
cantoDF$binN <- cantoDF$bin - sufONcanto
cantoDF$condToken <- 1

buscaDF <- filter(ss_10, target == "busca", sentence_id == "S06_C1")
buscaDF$binN <- 999
sufONbusca <- round(unique(buscaDF$onset_c3) / 10)
buscaDF[buscaDF$bin == sufONbusca, 'binN'] <- 0
buscaDF$binN <- buscaDF$bin - sufONbusca
buscaDF$condToken <- 1

buscoDF <- filter(ss_10, target == "buscÃ³", sentence_id == "S06_C2")
buscoDF$binN <- 999
sufONbusco <- round(unique(buscoDF$onset_c3) / 10)
buscoDF[buscoDF$bin == sufONbusco, 'binN'] <- 0
buscoDF$binN <- buscoDF$bin - sufONbusco
buscoDF$condToken <- 1

compraDF <- filter(ss_10, target == "compra", sentence_id == "S07_C1")
compraDF$binN <- 999
sufONcompra <- round(unique(compraDF$onset_c3) / 10)
compraDF[compraDF$bin == sufONcompra, 'binN'] <- 0
compraDF$binN <- compraDF$bin - sufONcompra
compraDF$condToken <- 1

comproDF <- filter(ss_10, target == "comprÃ³", sentence_id == "S07_C2")
comproDF$binN <- 999
sufONcompro <- round(unique(comproDF$onset_c3) / 10)
comproDF[comproDF$bin == sufONcompro, 'binN'] <- 0
comproDF$binN <- comproDF$bin - sufONcompro
comproDF$condToken <- 1

gastaDF <- filter(ss_10, target == "gasta", sentence_id == "S08_C1")
gastaDF$binN <- 999
sufONgasta <- round(unique(gastaDF$onset_c3) / 10)
gastaDF[gastaDF$bin == sufONgasta, 'binN'] <- 0
gastaDF$binN <- gastaDF$bin - sufONgasta
gastaDF$condToken <- 1

gastoDF <- filter(ss_10, target == "gastÃ³", sentence_id == "S08_C2")
gastoDF$binN <- 999
sufONgasto <- round(unique(gastoDF$onset_c3) / 10)
gastoDF[gastoDF$bin == sufONgasto, 'binN'] <- 0
gastoDF$binN <- gastoDF$bin - sufONgasto
gastoDF$condToken <- 1

juntaDF <- filter(ss_10, target == "junta", sentence_id == "S09_C1")
juntaDF$binN <- 999
sufONjunta <- round(unique(juntaDF$onset_c3) / 10)
juntaDF[juntaDF$bin == sufONjunta, 'binN'] <- 0
juntaDF$binN <- juntaDF$bin - sufONjunta
juntaDF$condToken <- 1

juntoDF <- filter(ss_10, target == "juntÃ³", sentence_id == "S09_C2")
juntoDF$binN <- 999
sufONjunto <- round(unique(juntoDF$onset_c3) / 10)
juntoDF[juntoDF$bin == sufONjunto, 'binN'] <- 0
juntoDF$binN <- juntoDF$bin - sufONjunto
juntoDF$condToken <- 1

salvaDF <- filter(ss_10, target == "salva", sentence_id == "S10_C1")
salvaDF$binN <- 999
sufONsalva <- round(unique(salvaDF$onset_c3) / 10)
salvaDF[salvaDF$bin == sufONsalva, 'binN'] <- 0
salvaDF$binN <- salvaDF$bin - sufONsalva
salvaDF$condToken <- 1

salvoDF <- filter(ss_10, target == "salvÃ³", sentence_id == "S10_C2")
salvoDF$binN <- 999
sufONsalvo <- round(unique(salvoDF$onset_c3) / 10)
salvoDF[salvoDF$bin == sufONsalvo, 'binN'] <- 0
salvoDF$binN <- salvoDF$bin - sufONsalvo
salvoDF$condToken <- 1

lanzaDF <- filter(ss_10, target == "lanza", sentence_id == "S11_C1")
lanzaDF$binN <- 999
sufONlanza <- round(unique(lanzaDF$onset_c3) / 10)
lanzaDF[lanzaDF$bin == sufONlanza, 'binN'] <- 0
lanzaDF$binN <- lanzaDF$bin - sufONlanza
lanzaDF$condToken <- 1

lanzoDF <- filter(ss_10, target == "lanzÃ³", sentence_id == "S11_C2")
lanzoDF$binN <- 999
sufONlanzo <- round(unique(lanzoDF$onset_c3) / 10)
lanzoDF[lanzoDF$bin == sufONlanzo, 'binN'] <- 0
lanzoDF$binN <- lanzoDF$bin - sufONlanzo
lanzoDF$condToken <- 1

pescaDF <- filter(ss_10, target == "pesca", sentence_id == "S12_C1")
pescaDF$binN <- 999
sufONpesca <- round(unique(pescaDF$onset_c3) / 10)
pescaDF[pescaDF$bin == sufONpesca, 'binN'] <- 0
pescaDF$binN <- pescaDF$bin - sufONpesca
pescaDF$condToken <- 1

pescoDF <- filter(ss_10, target == "pescÃ³", sentence_id == "S12_C2")
pescoDF$binN <- 999
sufONpesco <- round(unique(pescoDF$onset_c3) / 10)
pescoDF[pescoDF$bin == sufONpesco, 'binN'] <- 0
pescoDF$binN <- pescoDF$bin - sufONpesco
pescoDF$condToken <- 1

montaDF <- filter(ss_10, target == "monta", sentence_id == "S13_C1")
montaDF$binN <- 999
sufONmonta <- round(unique(montaDF$onset_c3) / 10)
montaDF[montaDF$bin == sufONmonta, 'binN'] <- 0
montaDF$binN <- montaDF$bin - sufONmonta
montaDF$condToken <- 1

montoDF <- filter(ss_10, target == "montÃ³", sentence_id == "S13_C2")
montoDF$binN <- 999
sufONmonto <- round(unique(montoDF$onset_c3) / 10)
montoDF[montoDF$bin == sufONmonto, 'binN'] <- 0
montoDF$binN <- montoDF$bin - sufONmonto
montoDF$condToken <- 1

mandaDF <- filter(ss_10, target == "manda", sentence_id == "S14_C1")
mandaDF$binN <- 999
sufONmanda <- round(unique(mandaDF$onset_c3) / 10)
mandaDF[mandaDF$bin == sufONmanda, 'binN'] <- 0
mandaDF$binN <- mandaDF$bin - sufONmanda
mandaDF$condToken <- 1

mandoDF <- filter(ss_10, target == "mandÃ³", sentence_id == "S14_C2")
mandoDF$binN <- 999
sufONmando <- round(unique(mandoDF$onset_c3) / 10)
mandoDF[mandoDF$bin == sufONmando, 'binN'] <- 0
mandoDF$binN <- mandoDF$bin - sufONmando
mandoDF$condToken <- 1

saltaDF <- filter(ss_10, target == "salta", sentence_id == "S15_C1")
saltaDF$binN <- 999
sufONsalta <- round(unique(saltaDF$onset_c3) / 10)
saltaDF[saltaDF$bin == sufONsalta, 'binN'] <- 0
saltaDF$binN <- saltaDF$bin - sufONsalta
saltaDF$condToken <- 1

saltoDF <- filter(ss_10, target == "saltÃ³", sentence_id == "S15_C2")
saltoDF$binN <- 999
sufONsalto <- round(unique(saltoDF$onset_c3) / 10)
saltoDF[saltoDF$bin == sufONsalto, 'binN'] <- 0
saltoDF$binN <- saltoDF$bin - sufONsalto
saltoDF$condToken <- 1

cargaDF <- filter(ss_10, target == "carga", sentence_id == "S16_C1")
cargaDF$binN <- 999
sufONcarga <- round(unique(cargaDF$onset_c3) / 10)
cargaDF[cargaDF$bin == sufONcarga, 'binN'] <- 0
cargaDF$binN <- cargaDF$bin - sufONcarga
cargaDF$condToken <- 1

cargoDF <- filter(ss_10, target == "cargÃ³", sentence_id == "S16_C2")
cargoDF$binN <- 999
sufONcargo <- round(unique(cargoDF$onset_c3) / 10)
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
                                        

write.csv(df_related,"data/clean/timecourse_10_clean.csv", row.names = FALSE)


# Time course
df_related %>%
  group_by(group) %>%
  summarise(n_distinct(participant))

fig_names <- c(`1` = "Present tense",
               `2` = "Preterit tense")

timecourse_stress <- df_related %>%
  filter(., binN > -50, binN < 80) %>%
  ggplot(., aes(x = binN, y = target_prop, color = group)) +
  facet_grid(. ~ cond, labeller = as_labeller(fig_names)) + 
  geom_hline(yintercept = 0.5, color = "white", size = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', pch = 21, 
               fill = "white") +
  theme_grey()

timecourse_stress + xlab("Time in 10 ms bins") + ylab("Proportion of fixations on target") + 
  guides(color=guide_legend(title="Group")) + #scale_color_discrete(labels=c("Interpreters", "Proficient L2", "Monolinguals")) +
  theme(legend.position="bottom")



























# calculate lowest max
maxes <- ss_10 %>%
  group_by(., target) %>%
  summarize(max = max(bin_adj)) %>%
  as.data.frame(.)
binAdjMaxMin <- min(maxes$max)

# calculate highest low
mins <- ss_10 %>%
  group_by(., target) %>%
  summarize(min = min(bin_adj)) %>%
  as.data.frame(.)
binAdjMinMax <- max(mins$min)

# subset data based on new ranges
df_short <- ss_10 %>% filter(., bin_adj <= binAdjMaxMin & bin_adj >= binAdjMinMax)

# create new adjusted variable that ranges from 1 to max
df_short$binREadj <- (df_short$bin_adj - binAdjMinMax) + 1



## @knitr binAdjustments

# Bin adjustments

# Where does target suffix begin in time course
suffixOnsets <- ss_10 %>%
  group_by(., target) %>%
  summarize(., sufOnset = unique(onset_c3) / 10)

# Center time course so that suffix onset = 0
suffixOnsetAdj <- ss_10 %>%
  group_by(., target) %>%
  summarize(., sufOnsetAdj = ((unique(onset_c3) / 10) -
                                (unique(onset_c3) / 10)))

# Where does the target word begin in the time course?
twOnsets <- ss_10 %>%
  group_by(., target) %>%
  summarize(., twOnset = unique(offset_prev_word) / 10)

# Adjust to centered time course
twOnsetAdj <- cbind(suffixOnsets, twOnsets[, 2])
twOnsetAdj <- mutate(twOnsetAdj, diff = sufOnset - twOnset, twOnsetAdj = 0 - diff)
