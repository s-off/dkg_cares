###### DKG CARES #####

# Packages

# install.packages("haven")
# install.packages("tidyverse")
# install.packages("gt")
# install.packages("gtsummary")
# install.packages("labelled")

library(haven)
library(tidyverse)
library(gt)
library(gtsummary)
library(labelled)

# Import Data

MCB <- read_spss("input/SUFRSDLV17MCB.sav")

KOB <- read_spss("input/SUFRSDLV17KOB.sav")

BYB <- read_spss("input/SUFRSDLV17BYB.sav")

### Subset mit Cases, die im Berichtzeitraum eine medizinische Rehabilitation auf Grund einer Krebsdiagnose begonnen haben
  ## Es wird nur die erste Episode (erste Rehabilitation auf Grund einer Krebsdiagnose betrachtet)

cancer_CASES_MCB_KOB <- MCB %>% 
  filter(mcdggr == 26,                                  # Bewilligungsdiagnosengruppe: Neubildungen (ICD10-Nr. C00 - D48)
         mcdg1_icd %in% c(7:31),                        # 1. med. Entlassungsdiagnose ICD-10: C00-96
         mcmsat %in% c(31,32),                          # Bewilligte Maßnahmeart: Ca(Krebs)-Reha-Leistung nach ? 15 SGB VI oder nach ? 31 Abs. 1 Nr. 2 SGB VI
         mcpsgral %in% c(2,3),                          # Art des Versicherungs- bzw. Rentnerstatus: Pflichtversicherter, freiwillig Versicherter
         mcumdt == 2,                                   # Umfang der Datenmeldung medizinische Rehabilitation: Reha-Leistung beendet, Datensatz vollständig
         mceafo == 1,                                   # Entlassungsform: regulär
         mcwhb == 0) %>%                                # Ausschluss von Wiederholungsleistungen
  mutate(RehaStart = mcbemsj*12 + mcbemsm) %>%          # Erste Episode wird durch RehaStart (Beginn der medizinischen Reha-Leistung) identifiziert
  arrange(case, RehaStart) %>% 
  distinct(case, .keep_all = TRUE) %>%
  inner_join(KOB, by = "case") %>%                      # Zusammenführen von MCB- und KOB-Daten
  mutate(GBZeitpkt = gbja*12 + gbmo,                    # Alter in Monaten
         RTZeitpkt = rtbej*12 + rtbem,                  # Rentenbeginn in Monaten
         TDZeitpkt = tddtj*12 + tddtm,                  # Todeszeitpunkt in Monaten
         AlterRehaStart = RehaStart - GBZeitpkt,        # Alter in Monaten bei Beginn der ersten med. Reha
         ZeitRehaStartRente = case_when(                # Zeitraum zwischen Beginn der ersten med. Reha und dem Rentenbeginn in Monaten
           rtbej %in% c(1:9998) ~ RTZeitpkt - RehaStart),
         ZeitRehaStartTod = case_when(                  # Zeitraum zwischen Beginn der ersten med. Reha und dem Todeszeitpunkt in Monaten
           tddtj %in% c(1:9998) ~ TDZeitpkt - RehaStart),
         RTvorReha = if_else(
           ZeitRehaStartRente < 0, 1, 0)) %>%
  filter(!RTvorReha %in% c(1),                          # Keine Rente vor Beginn der ersten med. Reha
         AlterRehaStart %in% c(216:804),                # Alter zwischen 18 und 67
         tddtj != 1992)                                 # Ausschluss Artefakt mit Todesjahr 1992



### Erstellen/Bearbeiten der Variablen

CARES_rtw <- cancer_CASES_MCB_KOB %>% 
  mutate(
    Geschlecht = case_when(
      sex == 1 ~ "männlich",
      sex == 2 ~ "weiblich"),
    AlterRehaStart_Jahr = case_when(
      gbja > 0 ~ mcbemsj-gbja),
    Tod = if_else(
      tddtj > 0, 1, 0)
    ) %>%
  remove_value_labels(mcdg1_icd = c(1:6, 32:999)) # Entfernen der leeren Diagnoselabels für gtsummary Ausgabe - Behält nur die lables für mmcdg1_icd %in% c(7:31), siehe Subseterstellung



### Datensatzbeschreibung

CARES_rtw %>%
  select(Geschlecht,
         AlterRehaStart_Jahr,
         mcfmsd,
         mcdg1_icd,
         tlrt,
         ZeitRehaStartRente,
         Tod,
         ZeitRehaStartTod,
         mcaiufzt,
         mcleft_lb,
         mcleft_at,
         bd,
         sb,
         bb,
         mcahb,
         mcaqda,
         mcstbf,
         mcwhot_ostwest,
         mcwhot_skt,
         mcmsot_bland) %>% 
  mutate_at(vars("bd", "sb", "bb", "mcahb", "mcfmsd", "mcwhot_ostwest", "mcwhot_skt", "mcmsot_bland", "mcstbf", "mcaiufzt", "mcleft_lb", "mcleft_at", "mcdg1_icd", "tlrt"), haven::as_factor) %>%
  tbl_summary(
    by = Geschlecht,
    statistic = all_continuous() ~ "{mean} ({sd})",
    label = list(
      AlterRehaStart_Jahr ~ "Alter bei Beginn der ersten Rehabilitationsleistung",
      mcdg1_icd ~ "Medizinische Entlassungsdiagnose der Rehabilitation",
      ZeitRehaStartRente ~ "Zeitraum vom Beginn der ersten medizinischen Rehabilitation bis zum Rentenbeginn in Monaten",
      Tod ~ "Versterben im Beobachtungszeitraum",
      ZeitRehaStartTod ~ "Zeitraum vom Beginn der ersten medizinischen Rehabilitation bis zum Todeszeitpunkt in Monaten",
      mcaqda ~ "Zeitraum von der Antragstellung eines Antrags auf medizinische Rehabilitation bis zu seiner Bewilligung in Tagen"),
    digits = list(
      AlterRehaStart_Jahr ~ 1,
      ZeitRehaStartTod ~ 0),
    missing = "no") %>% 
  add_overall() %>% 
  modify_header(label = "**Variable**") %>% 
  bold_labels()






##### Experimental #####


x <- CARES_rtw %>% 
  group_by(breha_3) %>% 
  summarise(count = n())

view(x)
