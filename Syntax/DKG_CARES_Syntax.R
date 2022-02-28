###### DKG CARES #####

## Packages und Themeeinstellungen

# install.packages("haven")
# install.packages("tidyverse")
# install.packages("gt")
# install.packages("gtsummary")
# install.packages("labelled")

library(tidyverse)
library(gt)
library(gtsummary)
library(labelled)

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
  filter(
    # Bewilligungsdiagnosengruppe: Neubildungen (ICD10-Nr. C00 - D48)
    mcdggr == 26,
    
    # 1. med. Entlassungsdiagnose ICD-10: C00-96
    mcdg1_icd %in% c(7:31),
    
    # Bewilligte Maßnahmeart: Ca(Krebs)-Reha-Leistung nach ? 15 SGB VI oder nach ? 31 Abs. 1 Nr. 2 SGB VI
    mcmsat %in% c(31,32),
    
    # Art des Versicherungs- bzw. Rentnerstatus: Pflichtversicherter, freiwillig Versicherter
    mcpsgral %in% c(2,3),
    
    # Umfang der Datenmeldung medizinische Rehabilitation: Reha-Leistung beendet, Datensatz vollständig
    mcumdt == 2,
    
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
    # Alter in Monaten
    GBZeitpkt = gbja*12 + gbmo,
    
    # Rentenbeginn in Monaten
    RTZeitpkt = rtbej*12 + rtbem,
    
    # Todeszeitpunkt in Monaten
    TDZeitpkt = tddtj*12 + tddtm,
    
    # Alter in Monaten bei Beginn der ersten med. Reha
    AlterRehaStart = RehaStart - GBZeitpkt,
    
    # Zeitraum zwischen Beginn der ersten med. Reha und dem Rentenbeginn in Monaten
    ZeitRehaStartRente = case_when(
      rtbej %in% c(1:9998) ~ RTZeitpkt - RehaStart),
    
    # Zeitraum zwischen Beginn der ersten med. Reha und dem Todeszeitpunkt in Monaten
    ZeitRehaStartTod = case_when(
      tddtj %in% c(1:9998) ~ TDZeitpkt - RehaStart),
    
    # Bezug einer Rente vor der Beginn der ersten med. Reha
    RTvorReha = if_else(
      ZeitRehaStartRente < 0, 1, 0)
    
    ) %>%
  
  filter(
    # Keine Rente vor Beginn der ersten med. Reha
    !RTvorReha %in% c(1),
    
    # Alter zwischen 18 und 67
    AlterRehaStart %in% c(216:804),
    
    # Ausschluss Artefakt mit Todesjahr 1992
    tddtj != 1992
    
    )



### Erstellen/Bearbeiten der Variablen

CARES_rtw <- cancer_CASES_MCB_KOB %>% 
  mutate(
    Geschlecht = case_when(
      sex == 1 ~ "männlich",
      sex == 2 ~ "weiblich"),
    
    AlterRehaStart_Jahr = case_when(
      gbja > 0 ~ mcbemsj-gbja),
    
    Tod = if_else(
      tddtj > 0, 1, 0),
    
    # bd (alte Berichtform) umcodieren in sb (neue Berichtform)
    sb = replace(sb, bd %in% c(3:6), 4),
    
    # Variable: psychische Erkrankung (ICD-10 F Diagnose) bei Beginn der ersten med. Reha
    psychischeErkrankung = if_else(
      mcdg2_icd %in% c(50:74) | mcdg3_icd %in% c(50:74) | mcdg4_icd %in% c(50:74) | mcdg5_icd %in% c(50:74), 1, 0)
    
    ) %>%
  
  # Entfernen der leeren Diagnoselabels für gtsummary Ausgabe - Behält nur die lables für mmcdg1_icd %in% c(7:31), siehe Subseterstellung
  remove_value_labels(mcdg1_icd = c(1:6, 32:999))



### Datensatzbeschreibung

CARES_rtw %>%
  select(Geschlecht,
         AlterRehaStart_Jahr,
         mcfmsd,
         mcdg1_icd,
         psychischeErkrankung,    # Mediator
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
  mutate_at(vars("sb", "mcsvbh", "mcbfkl_gr", "mcaivoaq", "mcaift", "mcahb", "mcakk", "mcfmsd", "mcwhot_ostwest", "mcwhot_skt", "mcmsot_bland", "mcstbf", "mcaiufzt", "mcleft_lb", "mcleft_at", "mcdg1_icd", "mcdg1_erg", "tlrt", "mc3voncms", "mc5voncms", "mc6voncms", "mc7voncms", "mc8voncms", "mc9voncms", "mc10voncms", "mc11voncms", "mc12voncms", "mc1pqle", "mc2pqle", "mc3pqle", "mc4pqle", "mc5pqle", "mc6pqle", "mc7pqle"), haven::as_factor) %>%
  tbl_summary(
    by = Geschlecht,
    statistic = all_continuous() ~ "{mean} ({sd})",
    label = list(
      AlterRehaStart_Jahr ~ "Alter bei Beginn der ersten Rehabilitationsleistung",
      psychischeErkrankung ~ "Psychische oder Verhaltensstörungen nach ICD-10 als weitere medizinische Entlassungsdiagnose",
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
  bold_labels()











##### Experimental #####


x <- CARES_rtw %>% 
  group_by(mc2voncms) %>% 
  summarise(count = n())

view(x)


x <- CARES_rtw %>%
  group_by(x, y) %>%
  tally() %>%
  spread(x, n)

view(x)