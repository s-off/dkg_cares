###### DKG CARES #####

## Packages und Themeeinstellungen

# install.packages("haven")
# install.packages("tidyverse")
# install.packages("naniar")
# install.packages("gt")
# install.packages("gtsummary")
# install.packages("labelled")
# install.packages("survival")
# install.packages("survminer")


library(haven)
library(tidyverse)
library(naniar)
library(gt)
library(gtsummary)
library(kableExtra)
library(labelled)
library(survival)
library(survminer)

## gtsummary Theme
theme_gtsummary_journal(journal = "jama")

theme_gtsummary_language(
  language = c("de"),
  decimal.mark = ",",
  big.mark = ".",
  iqr.sep = NULL,
  ci.sep = NULL,
  set_theme = TRUE
)

## Reset theme
#reset_gtsummary_theme()


# Import Data

MCB <- read_spss("input/SUFRSDLV17MCB.sav")

KOB <- read_spss("input/SUFRSDLV17KOB.sav")

BYB <- read_spss("input/SUFRSDLV17BYB.sav")

BFB <- read_spss("input/SUFRSDLV17BFB.sav")


### Subset mit Cases, die im Berichtzeitraum eine medizinische Rehabilitation auf Grund einer Krebsdiagnose begonnen haben
  ## Es wird nur die erste Episode (erste Rehabilitation auf Grund einer Krebsdiagnose betrachtet)

cancer_CASES_MCB_KOB <- MCB %>% 
  
  # Zusammenfassen der Primärvariablen
  mutate(
    mcdg1_icd = recode(mcdg1_icd,
                       `21` = 20.5,
                       `22` = 20.5,
                       `23` = 20.5) # Value = 22,5 um position in ICD beizubehalten
  ) %>% 
  
  # Lablen der Primärvariablen
  
  add_value_labels(
    mcdg1_icd = c("Bösartige Neubildungen der weiblichen Genitalorgane (C51-58)" = 20.5)
  ) %>% 
  
  # add_value_labels(
  #   mcdg1_icd = 20.5 = "Bösartige Neubildungen der weiblichen Genitalorgane (C51-58)"
  # ) %>%
  
  filter(
    # Bewilligungsdiagnosengruppe: Neubildungen (ICD10-Nr. C00 - D48)
    mcdggr == 26,
    
    # 1. med. Entlassungsdiagnose ICD-10: C00-96
    mcdg1_icd %in% c(7,8,10,12,16,20,20.5,24,25,26,28,29,30,31),
    
    # Bewilligte Maßnahmeart: Ca(Krebs)-Reha-Leistung nach ? 15 SGB VI oder nach ? 31 Abs. 1 Nr. 2 SGB VI
    mcmsat %in% c(31,32),
    
    # Art des Versicherungs- bzw. Rentnerstatus: Pflichtversicherter, freiwillig Versicherter
    mcpsgral %in% c(2,3),
    
    # # Umfang der Datenmeldung medizinische Rehabilitation: Reha-Leistung beendet, Datensatz vollständig // aktuell nicht nötig, weil alle beendet und vollständig
    # mcumdt == 2,
    
    # Entlassungsform: regulär
    mceafo == 1,
    
    # Ausschluss von Wiederholungsleistungen
    mcwhb == 0
    ) %>%
  
  # Erste Episode wird durch RehaStart (Beginn der medizinischen Reha-Leistung) identifiziert
  mutate(
    RehaStart = mcbemsj*12 + mcbemsm) %>%
  arrange(case, RehaStart) %>% 
  distinct(case, .keep_all = TRUE) %>%
  
  # Zusammenführen von MCB- und KOB-Daten
  inner_join(KOB, by = "case") %>%
  
  mutate(
    # Ende der medizinischen Reha-Leistung in Monaten
    RehaEnde = mcenmsj*12 + mcenmsm,
    
    # Alter in Monaten
    GBZeitpkt = gbja*12 + gbmo,
    
    # Rentenbeginn in Monaten
    RTZeitpkt = rtbej*12 + rtbem,
    
    # Todeszeitpunkt in Monaten für Personen, die im Beobachtungszeitraum versterben - Nachmeldungen von Todesfällen werden als Überlebend analysiert
    TDZeitpkt = case_when(
      tddtj %in% c(1:9998) ~ tddtj*12 + tddtm),
    
    TDZeitpkt = ifelse(
      TDZeitpkt > 24216, NA, TDZeitpkt),
    
    # Alter in Monaten bei Beginn der ersten med. Reha
    AlterRehaStart = RehaStart - GBZeitpkt,
    
    # Zeitraum zwischen Beginn der ersten med. Reha und dem Rentenbeginn in Monaten
    ZeitRehaStartRente = case_when(
      rtbej %in% c(1:9998) ~ RTZeitpkt - RehaStart),
    
    # Zeitraum zwischen Beginn der ersten med. Reha und dem Todeszeitpunkt in Monaten
    ZeitRehaStartTod = case_when(
      TDZeitpkt > 1 ~ TDZeitpkt - RehaStart),
    
    # Bezug einer Rente vor der Beginn der ersten med. Reha
    RTvorReha = if_else(
      ZeitRehaStartRente < 0, 1, 0)
    ) %>%
  
  filter(
    # Ausschluss von Randfällen am Ende des Beobachtungszeitraums (nach dem 30. Juni 2016)
    RehaStart < (2016*12+6),
    
    # Keine Rente vor Beginn der ersten med. Reha
    !RTvorReha %in% c(1),
    
    # Alter zwischen 18 und 63
    AlterRehaStart %in% c(216:756),
    
    # Ausschluss von Personen mit einer einer Dauer der medizinischen Reha-Leistung von über 6 Wochen
    mcdams < 43
    )


### Zusammenführen von MCB-/KOB- und BFB-Daten

## Alle Cases im endgültigen Datensatz mit Beginn und Ende der ersten medizinischen Reha-Leistung in Monaten
CasesForJoin <- cancer_CASES_MCB_KOB %>% 
  select(case, RehaStart, RehaEnde)

## Episodenauswahl der BFB-Daten
BFB_Episode <- BFB %>%
  left_join(CasesForJoin, by = "case") %>% 
  
  # Beginn der LTA in Monaten
  mutate(LTAStart = bfbemsj*12 + bfbemsm) %>% 
  
  # Auswahl des Subsets
  filter(RehaStart > 1,
         LTAStart > RehaStart) %>% 
  
  # Erste Episode wird durch LTASTart (Beginn der Leistungen zur Teilhabe am Arbeitsleben) identifiziert
  arrange(case, LTAStart) %>% 
  distinct(case, .keep_all = TRUE) %>% 
  select(-c(RehaStart,RehaEnde))

## Join

cancer_CASES_MCB_KOB_BFB <- cancer_CASES_MCB_KOB %>% 
  left_join(BFB_Episode, by = "case")


### Erstellen/Bearbeiten der Variablen

CARES_rtw <- cancer_CASES_MCB_KOB_BFB %>% 
  mutate(
    Geschlecht = case_when(
      sex == 1 ~ "männlich",
      sex == 2 ~ "weiblich"),
    
    AlterRehaStart_Jahr = case_when(
      gbja > 0 ~ mcbemsj-gbja),
    
    Tod = if_else(
      TDZeitpkt > 1, 1, 0, missing = 0),
    
    #Variable, die für eine Überlebenszeitanalyse genutz werden kann. TimeToEventMortality beschreibt die Zeit vom Ende der ersten med. Reha bis zum Versterben bzw. dem Ende des Beobachtungszeitraums (Dezember 2017)
    TimeToEventMortality = case_when(
      Tod == 1 ~ ZeitRehaStartTod,
      Tod == 0 ~ (2017*12+13)-RehaEnde), #2017*12+13 = Januar 2018 - wer bis 31.12.2017 nicht verstorben ist wird als überlebend/zensiert analysiert
    
    # bd (alte Berichtform) umcodieren in sb (neue Berichtform)
    sb = replace(sb, bd %in% c(3:6), 4),
    
    # Zeitraum zwischen Ende der ersten med. Reha und dem Beginn der Leistungen zur Teilhabe am Arbeitsleben
    ZeitRehaEndeLTAStart = case_when(
      LTAStart > 1 ~ LTAStart - RehaEnde),
    
    
    ## Komorbiditäten
    
    Komorb = case_when(
      mcdg2_icd %in% c(999) & mcdg3_icd %in% c(999) & mcdg4_icd %in% c(999) & mcdg5_icd %in% c(999) ~ 0,
      mcdg2_icd %in% c(2:6) | mcdg3_icd %in% c(2:6) | mcdg4_icd %in% c(2:6) | mcdg5_icd %in% c(2:6) ~ 1,
      
    )
    
#     
#     
#     # Variable: Bestimmte infektiöse und parasitäre Krankheiten (ICD-10 A00.-B99.) bei Beginn der ersten med. Reha
#     KomorbInfekt = if_else(
#       mcdg2_icd %in% c(2:6) | mcdg3_icd %in% c(2:6) | mcdg4_icd %in% c(2:6) | mcdg5_icd %in% c(2:6), 1, 0),
#     
#     # Variable: Weitere bösartige Tumorerkrankung (ICD-10 C00.-C97.) bei Beginn der ersten med. Reha
#     KomorbBN = if_else(
#       mcdg2_icd %in% c(7:31) | mcdg3_icd %in% c(7:31) | mcdg4_icd %in% c(7:31) | mcdg5_icd %in% c(7:31), 1, 0),
#     
#     # Variable: Zusätzliche Neubildung (In-Situ, Gutartig, unsicherem/unbekanntem Verhaltens) (ICD-10 D00.-D89.) bei Beginn der ersten med. Reha
#     KomorbGN = if_else(
#       mcdg2_icd %in% c(32:41) | mcdg3_icd %in% c(32:41) | mcdg4_icd %in% c(32:41) | mcdg5_icd %in% c(32:41), 1, 0),
#     
#     # Variable: Diabetes mellitus (ICD-10 E00.-E14.) bei Beginn der ersten med. Reha
#     KomorbDiabetes = if_else(
#       mcdg2_icd %in% c(42:45) | mcdg3_icd %in% c(42:45) | mcdg4_icd %in% c(42:45) | mcdg5_icd %in% c(42:45), 1, 0),
#     
#     # Variable: Ernährungserkrankung und Mangelernährung (ICD-10 E40.-64.) bei Beginn der ersten med. Reha
#     KomorbErnaerung = if_else(
#       mcdg2_icd %in% c(47) | mcdg3_icd %in% c(47) | mcdg4_icd %in% c(47) | mcdg5_icd %in% c(47), 1, 0),
#     
#     # Variable: Überernährung (ICD-10 E65.-68.) bei Beginn der ersten med. Reha
#     KomorbUeberernaerung = if_else(
#       mcdg2_icd %in% c(48) | mcdg3_icd %in% c(48) | mcdg4_icd %in% c(48) | mcdg5_icd %in% c(48), 1, 0),
#     
#     # Variable: Psychische oder Verhaltensstörung (ICD-10 F00.-F99.) bei Beginn der ersten med. Reha
#     KomorbPsych = if_else(
#       mcdg2_icd %in% c(50:74) | mcdg3_icd %in% c(50:74) | mcdg4_icd %in% c(50:74) | mcdg5_icd %in% c(50:74), 1, 0),
#     
#     
#         # # Variable: Demenz oder Organisches Psychosyndrom (ICD-10 F00.-F09.) bei Beginn der ersten med. Reha
#         # KomorbDemenz = if_else(
#         #   mcdg2_icd %in% c(50:52) | mcdg3_icd %in% c(50:52) | mcdg4_icd %in% c(50:52) | mcdg5_icd %in% c(50:52), 1, 0),
#         # 
#         # # Variable: Psychische oder Verhaltensstörung durch psychotrope Substanzen (ICD-10 F10.-F19.) bei Beginn der ersten med. Reha
#         # KomorbSucht = if_else(
#         #   mcdg2_icd %in% c(53:54) | mcdg3_icd %in% c(53:54) | mcdg4_icd %in% c(53:54) | mcdg5_icd %in% c(53:54), 1, 0),
#         # 
#         # # Variable: Schizophrenie, schizotype oder wahnhafte Störung (ICD-10 F20.-F29.) bei Beginn der ersten med. Reha
#         # KomorbSchizo = if_else(
#         #   mcdg2_icd %in% c(55:56) | mcdg3_icd %in% c(55:56) | mcdg4_icd %in% c(55:56) | mcdg5_icd %in% c(55:56), 1, 0),
#         # 
#         # # Variable: Affektive Störung (ICD-10 F30.-F39.) bei Beginn der ersten med. Reha
#         # KomorbDepri = if_else(
#         #   mcdg2_icd %in% c(57:60) | mcdg3_icd %in% c(57:60) | mcdg4_icd %in% c(57:60) | mcdg5_icd %in% c(57:60), 1, 0),
#         # 
#         # # Variable: Neurotische, Belastungs- oder somatoforme Störung (ICD-10 F40.-F48.) bei Beginn der ersten med. Reha
#         # KomorbAngst = if_else(
#         #   mcdg2_icd %in% c(61:67) | mcdg3_icd %in% c(61:67) | mcdg4_icd %in% c(61:67) | mcdg5_icd %in% c(61:67), 1, 0),
#         # 
#         # # Variable: Verhaltensauffälligkeit mit körperlicher Störung  (ICD-10 F50.-F59.) bei Beginn der ersten med. Reha
#         # KomorbEssstoer = if_else(
#         #   mcdg2_icd %in% c(68:69) | mcdg3_icd %in% c(68:69) | mcdg4_icd %in% c(68:69) | mcdg5_icd %in% c(68:69), 1, 0),
#         # 
#         # # Variable: Persönlichkeits- oder Verhaltensstörung (ICD-10 F60.-F69.) bei Beginn der ersten med. Reha
#         # KomorbPersoest = if_else(
#         #   mcdg2_icd %in% c(70) | mcdg3_icd %in% c(70) | mcdg4_icd %in% c(70) | mcdg5_icd %in% c(70), 1, 0),
#         # 
#         # # Variable: Intelligenzminderung (ICD-10 F70.-F79.) bei Beginn der ersten med. Reha
#         # KomorbIntell = if_else(
#         #   mcdg2_icd %in% c(71) | mcdg3_icd %in% c(71) | mcdg4_icd %in% c(71) | mcdg5_icd %in% c(71), 1, 0),
#         # 
#         # # Variable: Entwicklungsstörung (ICD-10 F80.-F89.) bei Beginn der ersten med. Reha
#         # KomorbEntwick = if_else(
#         #   mcdg2_icd %in% c(72) | mcdg3_icd %in% c(72) | mcdg4_icd %in% c(72) | mcdg5_icd %in% c(72), 1, 0),
#         # 
#         # # Variable: Verhaltens- oder emotionale Störung (ICD-10 F90.-F98.) bei Beginn der ersten med. Reha
#         # KomorbEmo = if_else(
#         #   mcdg2_icd %in% c(73) | mcdg3_icd %in% c(73) | mcdg4_icd %in% c(73) | mcdg5_icd %in% c(73), 1, 0),
#         # 
#         # # Variable: Nicht näher bezeichnete psychische Störungen (ICD-10 F99.) bei Beginn der ersten med. Reha
#         # KomorbsonstPsych = if_else(
#         #   mcdg2_icd %in% c(74) | mcdg3_icd %in% c(74) | mcdg4_icd %in% c(74) | mcdg5_icd %in% c(74), 1, 0),
# 
# 
# 
# 
# 
# # 
# #     # Variable: xx (ICD-10 xx-xx) bei Beginn der ersten med. Reha
# #     KomorbXx = if_else(
# #       mcdg2_icd %in% c(xx) | mcdg3_icd %in% c(xx) | mcdg4_icd %in% c(xx) | mcdg5_icd %in% c(xx), 1, 0),
# # 
# #     # Variable: xx (ICD-10 xx-xx) bei Beginn der ersten med. Reha
# #     KomorbXx = if_else(
# #       mcdg2_icd %in% c(xx) | mcdg3_icd %in% c(xx) | mcdg4_icd %in% c(xx) | mcdg5_icd %in% c(xx), 1, 0),
# # 
# #     # Variable: xx (ICD-10 xx-xx) bei Beginn der ersten med. Reha
# #     KomorbXx = if_else(
# #       mcdg2_icd %in% c(xx) | mcdg3_icd %in% c(xx) | mcdg4_icd %in% c(xx) | mcdg5_icd %in% c(xx), 1, 0),
    
    ) %>%
  
  # Hinzufügen von lables zur gtsummary Ausgabe
  set_variable_labels(
  #   # Komorbiditäten
  #   KomorbInfekt = "Bestimmte infektiöse und parasitäre Krankheiten (ICD-10 A00-B99)",
  #   KomorbBN = "Weitere bösartige Tumorerkrankung (ICD-10 C00.-C97.)",
  #   KomorbGN = "Zusätzliche Neubildung (In-Situ, Gutartig, unsicheren/unbekannten Verhaltens) (ICD-10 D00.-D89.)",
  #   KomorbDiabetes = "Diabetes mellitus (ICD-10 E00.-E14.)",
  #   KomorbErnaerung = "Ernährungserkrankung und Mangelernährung (ICD-10 E40.-64.)",
  #   KomorbUeberernaerung = "Überernährung (ICD-10 E65.-68.)",
  #   KomorbPsych = "Psychische oder Verhaltensstörung (ICD-10 F00.-F99.)",
  #   KomorbDemenz = "Demenz oder Organisches Psychosyndrom (ICD-10 F00.-F09.)",
  #   KomorbSucht = "Psychische oder Verhaltensstörung durch psychotrope Substanzen (ICD-10 F10.-F19.)",
  #   KomorbSchizo = "Schizophrenie, schizotype oder wahnhafte Störung (ICD-10 F20.-F29.)",
  #   KomorbDepri = "Affektive Störung (ICD-10 F30.-F39.)",
  #   KomorbAngst = "Neurotische, Belastungs- oder somatoforme Störung (ICD-10 F40.-F48.)",
  #   KomorbEssstoer = "Verhaltensauffälligkeit mit körperlicher Störung  (ICD-10 F50.-F59.)",
  #   KomorbPersoest = "Persönlichkeits- oder Verhaltensstörung (ICD-10 F60.-F69.)",
  #   KomorbIntell = "Intelligenzminderung (ICD-10 F70.-F79.)",
  #   KomorbEntwick = "Entwicklungsstörung (ICD-10 F80.-F89.)",
  #   KomorbEmo = "Verhaltens- oder emotionale Störung (ICD-10 F90.-F98.)",
  #   KomorbsonstPsych = "Nicht näher bezeichnete psychische Störungen (ICD-10 F99.)"#,
  #   # xx = "xx",
  #   # xx = "xx",
  #   # xx = "xx",
  
    Komorb = "Komorbiditäten",
  ) %>%
  
  # Anpassung der Variablenlabel für die gtsummary Ausgabe
  
  add_value_labels(
    mcdg1_icd = c(
      "BN der weiblichen Genitalorgane (C51-58)"  = 20.5,
      "BN der männlichen Genitalorgane (C60-63)" = 24,
      "BN der Harnorgane (C64-68)" = 25,
      "BN Leber, Gallengänge /-blase, Pankreas, sonstige Verdauungsorgane, mehrere Lokalisationen" = 28
    ),
    tlrt = c(
      "kein Rentenbezug" = 0
    )
  ) %>% 
  
  set_value_labels(
    Komorb = c(
      "Keine Komorbiditäten" = 0,
      "Bestimmte infektiöse und parasitäre Krankheiten (ICD-10 A00.-B99.)" = 1
    )
  ) %>% 
  
  # Entfernen der leeren Diagnoselabels für gtsummary Ausgabe
  drop_unused_value_labels()





### Datensatzbeschreibung

CARES_rtw %>%
  select(Geschlecht,
         AlterRehaStart_Jahr,
         mcfmsd,
         mcdg1_icd,
         Komorb,
         # # Komorbiditäten - Mediator
         # KomorbInfekt,
         # KomorbBN,
         # KomorbGN,
         # KomorbDiabetes,
         # KomorbErnaerung,
         # KomorbUeberernaerung,
         # KomorbDemenz,
         # KomorbSucht,
         # KomorbSchizo,
         # KomorbDepri,
         # KomorbAngst,
         # KomorbEssstoer,
         # KomorbPersoest,
         mcdg1_erg,               # Mediator
         tlrt,
         ZeitRehaStartRente,
         Tod,
         ZeitRehaStartTod,
         sb,                      # Mediator
         mcsvbh,                  # Mediator
         mcbfkl_gr,               # Mediator
         mcaivoaq,
         mcaiufzt,                # Mediator
         mcaift,                  # Mediator
         mcleft_lb,               # Mediator
         mcleft_at,               # Mediator
         mcstbf,                  # Mediator
         mcahb,
         mcakk,                   # Mediator
         mcaqda,
         mc3voncms,
         mc5voncms,
         mc6voncms,
         mc7voncms,
         mc8voncms,
         mc9voncms,
         mc10voncms,
         mc11voncms,
         mc12voncms,
         mc1pqle,
         mc2pqle,
         mc3pqle,
         mc4pqle,
         mc5pqle,
         mc6pqle,
         mc7pqle,
         mcwhot_ostwest,
         mcwhot_skt,
         mcmsot_bland) %>% 
  mutate_at(vars("sb", "mcsvbh", "mcbfkl_gr", "mcaivoaq", "mcaift", "mcahb", "mcakk", "mcfmsd", "mcwhot_ostwest", "mcwhot_skt", "mcmsot_bland", "mcstbf", "mcaiufzt", "mcleft_lb", "mcleft_at", "mcdg1_icd", "Komorb", "mcdg1_erg", "tlrt", "mc3voncms", "mc5voncms", "mc6voncms", "mc7voncms", "mc8voncms", "mc9voncms", "mc10voncms", "mc11voncms", "mc12voncms", "mc1pqle", "mc2pqle", "mc3pqle", "mc4pqle", "mc5pqle", "mc6pqle", "mc7pqle"), haven::as_factor) %>%
  tbl_summary(
    by = Geschlecht,
    statistic = all_continuous() ~ "{mean} ({sd})",
    label = list(
      AlterRehaStart_Jahr ~ "Alter bei Beginn der ersten Rehabilitationsleistung",
      mcdg1_icd ~ "Medizinische Entlassungsdiagnose der Rehabilitation",
      ZeitRehaStartRente ~ "Zeitraum vom Beginn der ersten medizinischen Rehabilitation bis zum Rentenbeginn in Monaten",
      Tod ~ "Versterben im Beobachtungszeitraum",
      ZeitRehaStartTod ~ "Zeitraum vom Beginn der ersten medizinischen Rehabilitation bis zum Todeszeitpunkt in Monaten",
      mcbfkl_gr ~ "Berufsgruppenklassifikation nach dem Statistikband zur Rehabilitation der Rentenversicherung (Grundlage ist die KldB 88)",
      mcaqda ~ "Zeitraum von der Antragstellung eines Antrags auf medizinische Rehabilitation bis zu seiner Bewilligung in Tagen"),
    digits = list(
      AlterRehaStart_Jahr ~ 1,
      ZeitRehaStartTod ~ 0),
    missing = "no") %>% 
  add_overall() %>% 
  modify_header(label = "**Variable**") %>% 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Geschlecht**") %>% 
  bold_labels() #%>% 
  as_kable_extra()



### Survival

# Surv Objekt für Tod und Zeit bis zum Tod oder Ende des Beobachtungszeitraums
Surv_icd <- Surv(CARES_rtw$TimeToEventMortality, CARES_rtw$Tod, type = "right")

# Kaplan-Meier-Schätzer für die verschiedenen Diagnosegruppen
fit_icd <- survfit(Surv_icd~mcdg1_icd, data = CARES_rtw)

# Kaplan-Meier-Kurve für die verschiedenen Diagnosegruppen
plot(fit_icd, xlab = "Überleben in Monaten", ylab = "Prozentuales Überleben", yscale = 100,
     col = c(7:31))
legend("bottomleft", title="ICD-Gruppen", c("C00-14","C15, C16","C17, C18","C19, C20, C21","C30-39","C40-41","C43-44","C45-49","C50","C51-53","C54","C56-58","C60-63","C64-68","C69-72","C73-75","C22-26, C76-80, C97","C81","C82-90","C91-96"),
       fill= c(7:31))

ggsurvplot(fit_icd, data = CARES_rtw, pval = TRUE)



####


# BYB-Auswahl aller Episoden der Personen, die im endgültigen Datensatz zur Analyse eingeschlossen sind

BYB_Episoden <- BYB %>% 
  inner_join(CasesForJoin, by = "case")

# Erstellen des Vektors für die Auswertung der Monate bis zur RTW

ed_data <- readRDS("BYB_Episoden.rds")
fillValueNoWork = 85 
beginYear = 2007
endYear = 2017
prefCaseNum = ed_data[1,1]
prefYear = ed_data[1,3]
currentVec <- c(1:(12*(endYear - beginYear +1)))
currentVec[1:(12*(endYear - beginYear +1))] <- -1
AllCases <- data.frame(test = 1:132)
currentVec[1:12]<- unlist(ed_data[1,5:16],use.names = FALSE)
for (indexDataRows in 2:nrow(ed_data)) {
  print(indexDataRows)
  currCaseNum = ed_data[indexDataRows,1]
  currYear = ed_data[indexDataRows,3]
  if (currCaseNum != prefCaseNum) {
    if(prefYear < endYear){
      for (yearDiff in ((prefYear +1):endYear)) {
        currentVec[((12*(yearDiff - beginYear) +1)):(12*(yearDiff - (beginYear -1)))]<- fillValueNoWork
      }
    }
    nameCol <- sprintf("Case_%d", prefCaseNum)
    AllCases[nameCol] <- currentVec
    currentVec[1:(12*(endYear - beginYear +1))] <- -1
    
    prefYear = beginYear -1
    
  }
  if (currYear > prefYear +1) {
    for (yearDiff in (prefYear +1):(currYear-1)) {
      currentVec[(12*(yearDiff - beginYear) +1):(12*(yearDiff - (beginYear -1)))]<- fillValueNoWork
    }
  }
  #print(unlist(currentVec))
  currentVec[((12*(currYear - beginYear) +1)):(12*(currYear - (beginYear -1)))]<- unlist(ed_data[indexDataRows,(5:16)],use.names = FALSE)
  #print(unlist(currentVec))
  prefCaseNum = currCaseNum
  prefYear = currYear
}














##### Experimental #####


x <- CARES_rtw %>% 
  group_by(mcdg5_icd) %>% 
  summarise(count = n())

view(x)


x <- CARES_rtw %>%
  group_by(Tod, TimeToEventMortality) %>%
  tally() %>%
  spread(Tod, n)

view(x)
