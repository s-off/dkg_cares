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
# install.packages("dagitty")
# install.packages("ggdag")


library(haven)
library(tidyverse)
library(naniar)
library(gt)
library(gtsummary)
library(kableExtra)
library(flextable)
library(labelled)
library(survival)
library(survminer)
library(dagitty)
library(ggdag)


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

theme_gtsummary_compact(set_theme = TRUE, font_size = NULL)

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
  select(case, RehaStart, mcbemsj, mcbemsm, RehaEnde, mcenmsj, mcenmsm, mcdg1_icd)

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
  select(-c(RehaStart, mcbemsj, mcbemsm, RehaEnde, mcenmsj, mcenmsm, mcdg1_icd))

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
    
    # Leistung zur Teilhabe am Arbeitsleben
    LTA = if_else(
      LTAStart > 1, 1, 0, missing = 0),
    
    # Zeitraum zwischen Ende der ersten med. Reha und dem Beginn der Leistungen zur Teilhabe am Arbeitsleben
    ZeitRehaEndeLTAStart = case_when(
      LTAStart > 1 ~ LTAStart - RehaEnde),
    
    # Rentenart
    
    Rente = case_when(
      tlrt == 0 ~ 0,
      tlrt %in% c(10,11,12,13,18) ~ 1,
      tlrt %in% c(20:28) ~ 2),
    
    ZeitRehaStartAltersRente = case_when(
      Rente == 1 ~ RTZeitpkt - RehaStart),
    
    ZeitRehaStartEMRente = case_when(
      Rente == 2 ~ RTZeitpkt - RehaStart),
    
    
    ## Komorbiditäten
    
      # Keine bekannten Komorbiditäten
    Komorb00 = if_else(
      mcdg2_icd %in% c(999) & mcdg3_icd %in% c(999) & mcdg4_icd %in% c(999) & mcdg5_icd %in% c(999), 1, 0, missing = 0),
      
    # Bestimmte infektiöse und parasitäre Krankheiten (A00.-B99.)
    Komorb01 = if_else(
      mcdg2_icd %in% c(2:6) | mcdg3_icd %in% c(2:6) | mcdg4_icd %in% c(2:6) | mcdg5_icd %in% c(2:6), 1, 0, missing = 0),
    
    # Weitere bösartige Tumorerkrankung (C00.-C97.)
    Komorb02 = if_else(
      mcdg2_icd %in% c(7:31) | mcdg3_icd %in% c(7:31) | mcdg4_icd %in% c(7:31) | mcdg5_icd %in% c(7:31), 1, 0, missing = 0),
    
    # Zusätzliche Neubildung (In-Situ, Gutartig, unsicherem/unbekanntem Verhaltens) (D00.-D89.)
    Komorb03 = if_else(
      mcdg2_icd %in% c(32:41) | mcdg3_icd %in% c(32:41) | mcdg4_icd %in% c(32:41) | mcdg5_icd %in% c(32:41), 1, 0, missing = 0),
    
    # Diabetes mellitus (E00.-E14.)
    Komorb04 = if_else(mcdg2_icd %in% c(42:45) | mcdg3_icd %in% c(42:45) | mcdg4_icd %in% c(42:45) | mcdg5_icd %in% c(42:45), 1, 0, missing = 0),
    
    # Ernährungserkrankung und Mangelernährung (E40.-64.)
    Komorb05 = if_else(
      mcdg2_icd %in% c(47) | mcdg3_icd %in% c(47) | mcdg4_icd %in% c(47) | mcdg5_icd %in% c(47), 1, 0, missing = 0),
    
    # Überernährung (ICD-10 E65.-68.)
    Komorb06 = if_else(
      mcdg2_icd %in% c(48) | mcdg3_icd %in% c(48) | mcdg4_icd %in% c(48) | mcdg5_icd %in% c(48), 1, 0, missing = 0),
    
    # Psychische oder Verhaltensstörung (F00.-F99.)
    Komorb07 = if_else(
      mcdg2_icd %in% c(50:74) | mcdg3_icd %in% c(50:74) | mcdg4_icd %in% c(50:74) | mcdg5_icd %in% c(50:74), 1, 0, missing = 0),
    
    # Krankheit des zentralen Nervensystems (G00.-G37.)
    Komorb08 = if_else(
      mcdg2_icd %in% c(75:79) | mcdg3_icd %in% c(75:79) | mcdg4_icd %in% c(75:79) | mcdg5_icd %in% c(75:79), 1, 0, missing = 0),
    
    # Migräne und Schlafstörungen (G43.-44., G47.)
    Komorb09 = if_else(
      mcdg2_icd %in% c(81) | mcdg3_icd %in% c(81) | mcdg4_icd %in% c(81) | mcdg5_icd %in% c(81), 1, 0, missing = 0),
    
    # Krankheit des peripheren Nervensystems (G50.-82., ohne G80., G83.)
    Komorb10 = if_else(
      mcdg2_icd %in% c(83:88) | mcdg3_icd %in% c(83:88) | mcdg4_icd %in% c(83:88) | mcdg5_icd %in% c(83:88), 1, 0, missing = 0),
    
    # Ischämische oder pulmonale Herzkrankheiten (I20.-28.)
    Komorb11 = if_else(
      mcdg2_icd %in% c(98:100) | mcdg3_icd %in% c(98:100) | mcdg4_icd %in% c(98:100) | mcdg5_icd %in% c(98:100), 1, 0, missing = 0),
    
    # Chronische Krankheiten der unteren Atemwege (J43.-47.)
    Komorb12 = if_else(
      mcdg2_icd %in% c(100:111) | mcdg3_icd %in% c(100:111) | mcdg4_icd %in% c(100:111) | mcdg5_icd %in% c(100:111), 1, 0, missing = 0),
    
    # Nichtinfektiöse Enteritis, Kolitis oder Darmerkrankung (K50.-67.)
    Komorb13 = if_else(
      mcdg2_icd %in% c(117:119) | mcdg3_icd %in% c(117:119) | mcdg4_icd %in% c(117:119) | mcdg5_icd %in% c(117:119), 1, 0, missing = 0),
    
    # Krankheiten des Muskel-Skelett-Systems und des Bindegewebes (M00.-99.)
    Komorb14 = if_else(
      mcdg2_icd %in% c(129:153) | mcdg3_icd %in% c(129:153) | mcdg4_icd %in% c(129:153) | mcdg5_icd %in% c(129:153), 1, 0, missing = 0),
    
    # Niereninsuffizienz (N17.-19.)
    Komorb15 = if_else(
      mcdg2_icd %in% c(154) | mcdg3_icd %in% c(154) | mcdg4_icd %in% c(154) | mcdg5_icd %in% c(154), 1, 0, missing = 0)
    ) %>%
  
  # Hinzufügen von lables zur gtsummary Ausgabe
  set_variable_labels(
    AlterRehaStart_Jahr = "Alter bei Beginn der ersten med. Reha",
    Tod = "Versterben im Beobachtungszeitraum",
    LTA = "Leistung zur Teilhabe am Arbeitsleben",
    ZeitRehaEndeLTAStart = "Zeitraum zwischen Ende der ersten med. Reha und dem Beginn der Leistungen zur Teilhabe am Arbeitsleben in Monaten",
    ZeitRehaStartTod = "Zeitraum vom Beginn der ersten med. Reha bis zum Todeszeitpunkt in Monaten",
    ZeitRehaStartAltersRente = "Zeitraum zwischen Ende der med. Reha und dem Bezug einer Altersrente in Monaten",
    ZeitRehaStartEMRente = "Zeitraum zwischen Ende der med. Reha und dem Bezug einer Erwerbsminderungsrente in Monaten",
    Komorb00 = "Keine bekannten Komorbiditäten",
    Komorb01 = "Bestimmte infektiöse und parasitäre Krankheiten (A00.-B99.)",
    Komorb02 = "Weitere bösartige Tumorerkrankung (C00.-C97.)",
    Komorb03 = "Zusätzliche Neubildung (In-Situ, Gutartig, unsicherem/unbekanntem Verhaltens) (D00.-D89.)",
    Komorb04 = "Diabetes mellitus (E00.-E14.)",
    Komorb05 = "Ernährungserkrankung und Mangelernährung (E40.-64.)",
    Komorb06 = "Überernährung (E65.-68.)",
    Komorb07 = "Psychische oder Verhaltensstörung (F00.-F99.)",
    Komorb08 = "Krankheit des zentralen Nervensystems (G00.-G37.)",
    Komorb09 = "Migräne und Schlafstörungen (G43.-44., G47.)",
    Komorb10 = "Krankheit des peripheren Nervensystems (G50.-82., ohne G80., G83.)",
    Komorb11 = "Ischämische oder pulmonale Herzkrankheiten (I20.-28.)",
    Komorb12 = "Chronische Krankheiten der unteren Atemwege (J43.-47.)",
    Komorb13 = "Nichtinfektiöse Enteritis, Kolitis oder Darmerkrankung (K50.-67.)",
    Komorb14 = "Krankheiten des Muskel-Skelett-Systems und des Bindegewebes (M00.-99.)",
    Komorb15 = "Niereninsuffizienz (N17.-19.)"
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
    Rente = c(
      "kein Rentenbezug" = 0,
      "Rente wegen Alters" = 1,
      "Rente wegen verminderter Erwerbsfähigkeit" = 2
    )
  ) %>% 
  
  # Entfernen der leeren Diagnoselabels für gtsummary Ausgabe
  drop_unused_value_labels()


### Tabelle 1

CARES_rtw %>% 
  select(
    AlterRehaStart_Jahr,       # Alter bei Beginn der ersten Rehabilitationsleistung
    Geschlecht,                # Geschlecht
    mcdg1_icd,                 # Medizinische Entlassungsdiagnose der Rehabilitation
    mcdg1_erg,                 # 1. Diagnose: Behandlungsergebnis
    Komorb00,                  # Komorbiditäten
    Komorb01,
    Komorb02,
    Komorb03,
    Komorb04,
    Komorb05,
    Komorb06,
    Komorb07,
    Komorb08,
    Komorb09,
    Komorb10,
    Komorb11,
    Komorb12,
    Komorb13,
    Komorb14,
    Komorb15,
    sb,                        # Schulbildung	
    mcaivoaq,                  # Erwerbsstatus u. -umfang vor Antragsstellung 
    mcbfkl_gr,                 # Berufsgruppenklassifikation
    mcstbf,                    # Stellung im Beruf
    mcaiufzt,                  # Arbeitsunfähigkeit in den letzten 12 Monaten
    #mcaift,                    # Arbeitsfähigkeit
    mcleft_lb,                 # Leistungsfähigkeit im letzten Beruf
    mc8voncms,                   # Empfehlung LTA
    LTA,                       # Leistung zur Teilhabe am Arbeitsleben
    ZeitRehaEndeLTAStart,      # Zeitraum zwischen Ende der ersten med. Reha und dem Beginn der Leistungen zur Teilhabe am Arbeitsleben
    mcfmsd,                    # Familienstand
    mcwhot_ostwest,            # Wohnort: altes/neues Bundesgebiet
    mcwhot_skt,                # Siedlungsstruktureller Kreistyp
    #sa,                        # Staatsangehörigkeit des Antragstellers
    Tod,                       # Versterben im Beobachtungszeitraum
    ZeitRehaStartTod,          # Zeitraum vom Beginn der ersten medizinischen Rehabilitation bis zum Todeszeitpunkt in Monaten
    Rente,                     # Teilrentenkennzeichen in Altersrente/EM-Rente
    ZeitRehaStartAltersRente,  # Zeitraum zwischen Ende der med. Reha und dem Bezug einer Altersrente
    ZeitRehaStartEMRente       # Zeitraum zwischen Ende der med. Reha und dem Bezug einer Erwerbsminderungsrente
  ) %>% 
  mutate_at(vars("mcdg1_icd", "mcdg1_erg", "sb", "mcaivoaq", "mcbfkl_gr", "mcstbf", "mcaiufzt", "mcleft_lb", "mc8voncms", "mcfmsd", "mcwhot_ostwest", "mcwhot_skt", "Rente"), haven::as_factor) %>%
  tbl_summary(
    #by = mcdg1_icd,
    statistic = all_continuous() ~ "{mean} ({sd})",
    label = list(
      mcdg1_icd ~ "Primäre Tumorerkrankung",
      mcbfkl_gr ~ "Berufsgruppenklassifikation nach dem Statistikband zur Rehabilitation der Rentenversicherung (Grundlage ist die KldB 88)",
      mcaivoaq ~ "Erwerbsstatus und -umfang vor Antragsstellung"
    ),
    missing = "no",
    type = list(LTA ~ "dichotomous")
  ) %>%
  modify_header(label = "**Variable**") %>% 
  bold_labels()

### Outcomes

CARES_rtw %>% 
  select(
    mcdg1_icd,                 # Medizinische Entlassungsdiagnose der Rehabilitation
    Tod,                       # Versterben im Beobachtungszeitraum
    ZeitRehaStartTod,          # Zeitraum vom Beginn der ersten medizinischen Rehabilitation bis zum Todeszeitpunkt in Monaten
    Rente,                     # Teilrentenkennzeichen in Altersrente/EM-Rente
    ZeitRehaStartAltersRente,  # Zeitraum zwischen Ende der med. Reha und dem Bezug einer Altersrente
    ZeitRehaStartEMRente       # Zeitraum zwischen Ende der med. Reha und dem Bezug einer Erwerbsminderungsrente
  ) %>% 
  mutate_at(vars("mcdg1_icd", "Rente"), haven::as_factor) %>%
  tbl_summary(
    by = mcdg1_icd,
    statistic = all_continuous() ~ "{mean} ({sd})",
    label = list(
      mcdg1_icd ~ "Primäre Tumorerkrankung"
    ),
    missing = "no"
  ) %>%
  add_p() %>% 
  modify_header(label = "**Variable**") %>% 
  bold_labels()


### ausführliche Datensatzbeschreibung

CARES_rtw %>%
  select(Geschlecht,
         AlterRehaStart_Jahr,
         mcfmsd,
         mcdg1_icd,
         Komorb,
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
      mcsvbh ~ "Sozialversicherungspflichtige Beschäftigung bei Beginn der ersten Rehabilitationsleistung",
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


### RTW - time to event


# BYB-Auswahl aller Episoden der Personen, die im endgültigen Datensatz zur Analyse eingeschlossen sind

BYB_Episoden <- BYB %>% 
  inner_join(CasesForJoin, by = "case")

BYB_1_Jahr_nach_Reha <- BYB_Episoden %>% 
  filter(
    byja == (mcenmsj+1)
  ) %>% 
  mutate(
    rtw182 = if_else(
      bybhbytg1 > 181 | bybhbytg2 > 181, 1, 0)
  ) %>% 
  # Entfernen der leeren Diagnoselabels für gtsummary Ausgabe
  drop_unused_value_labels()

BYB_2_Jahr_nach_Reha <- BYB_Episoden %>% 
  filter(
    byja == (mcenmsj+2)
  ) %>% 
  mutate(
    rtw182 = if_else(
      bybhbytg1 > 181 | bybhbytg2 > 181, 1, 0)
  ) %>% 
  # Entfernen der leeren Diagnoselabels für gtsummary Ausgabe
  drop_unused_value_labels()

BYB_3_Jahr_nach_Reha <- BYB_Episoden %>% 
  filter(
    byja == (mcenmsj+3)
  ) %>% 
  mutate(
    rtw182 = if_else(
      bybhbytg1 > 181 | bybhbytg2 > 181, 1, 0)
  ) %>% 
  # Entfernen der leeren Diagnoselabels für gtsummary Ausgabe
  drop_unused_value_labels()

### RTW 1 Jahr

BYB_1_Jahr_nach_Reha %>% 
  select(
    rtw182,
    mcdg1_icd
  ) %>% 
  mutate_at(vars("mcdg1_icd"), haven::as_factor) %>%
  tbl_summary(
    by = mcdg1_icd,
    statistic = all_continuous() ~ "{mean} ({sd})",
    label = list(
      mcdg1_icd ~ "Primäre Tumorerkrankung"
    ),
    missing = "no"
  ) %>%
  add_p() %>% 
  modify_header(label = "**Variable**") %>% 
  bold_labels()

BYB_1_Jahr_nach_Reha %>% 
  select(
    rtw182,
  ) %>% 
  tbl_summary(
    statistic = all_continuous() ~ "{mean} ({sd})",
    missing = "no"
  ) %>%
  modify_header(label = "**Variable**") %>% 
  bold_labels()



# Erstellen des Vektors für die Auswertung der Monate bis zur RTW

ed_data <- BYB_Episoden
ed_reha_end <- ed_reha
count_in_end_table = 2 
EndOfTherapyTime <- data.frame(Case_831 = ((ed_reha_end[1,2]-2007*12):(ed_reha_end[1,2]-2007*12)))

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
  if (currCaseNum != prefCaseNum) { # Wenn wir zum nächsten Fall kommen, wird der Vektor mit "end values" ergänzt unter AllCases gespeichert
    if(prefYear < endYear){
      for (yearDiff in ((prefYear +1):endYear)) {
        currentVec[((12*(yearDiff - beginYear) +1)):(12*(yearDiff - (beginYear -1)))]<- fillValueNoWork
      }
    }
    nameCol <- sprintf("Case_%d", prefCaseNum)
    AllCases[nameCol] <- currentVec
    limit =0
    while(ed_data[indexDataRows,1] != ed_reha_end[count_in_end_table,1]){ # einige cases sind nur entweder in MCB/KPB oder in BYB - müssen geskippt werden
      count_in_end_table = count_in_end_table +1
      if(limit> 20 ){
        break
      }
      limit = limit +1
      
    }
    if(ed_data[indexDataRows,1] != ed_reha_end[count_in_end_table,1]){
      count_in_end_table = count_in_end_table +1
      if(ed_data[indexDataRows,1] != ed_reha_end[count_in_end_table,1]){
        print(currCaseNum)
        print(ed_reha_end[count_in_end_table,1] )
        print("dsoifjoidjfoidsjfpijdspifjdspijfoidsjfoiWRONG: case numer ",ed_reha_end[count_in_end_table,1] )
      }
    }
    nameColNew <- sprintf("Case_%d", currCaseNum)
    EndOfTherapyTime[nameColNew] <- c((ed_reha_end[count_in_end_table,2]-2007*12):(ed_reha_end[count_in_end_table,2]-2007*12))
    count_in_end_table = count_in_end_table +1
    currentVec[1:(12*(endYear - beginYear +1))] <- -1
    
    prefYear = beginYear -1
    
  }
  if (currYear > prefYear +1) { # wenn keine Daten vorhanden: fill - default value
    for (yearDiff in (prefYear +1):(currYear-1)) {
      currentVec[(12*(yearDiff - beginYear) +1):(12*(yearDiff - (beginYear -1)))]<- fillValueNoWork
    }
  }
  #print(unlist(currentVec))
  # Falls daten vorhanden - fill mit richtigen Werten
  currentVec[((12*(currYear - beginYear) +1)):(12*(currYear - (beginYear -1)))]<- unlist(ed_data[indexDataRows,(5:16)],use.names = FALSE) 
  #print(unlist(currentVec))
  prefCaseNum = currCaseNum
  prefYear = currYear
}

# Berechnung: Erste Veränderung im Zeitverlauf
CaseTmeTillChange <- data.frame(test= 1:1)
for (indexCase in 1:ncol(EndOfTherapyTime)) {
  endOfTherapyfromOneTo132 <- EndOfTherapyTime[1,indexCase] 
  nameCol <- names(EndOfTherapyTime)[indexCase]
  if (endOfTherapyfromOneTo132 >= (12*(endYear - beginYear +1))) {
    CaseTmeTillChange[nameCol] <- -1
  }
  currentVec = AllCases[,indexCase+1]
  prefOccupation = currentVec[endOfTherapyfromOneTo132]
  for (currentMonth in (endOfTherapyfromOneTo132+1):length(currentVec)) {
    if(prefOccupation != currentVec[currentMonth]){
      CaseTmeTillChange[nameCol] <- currentMonth - endOfTherapyfromOneTo132 
      break
    }
    if(currentMonth == length(currentVec)){
      CaseTmeTillChange[nameCol] <- -2 # wenn keine Veränderung bis zum Ende aufgetreten ist
    }
  }
}

