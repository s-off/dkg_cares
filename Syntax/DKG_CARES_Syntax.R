###### DKG CARES #####

# Packages

# install.packages("haven")
# install.packages("tidyverse")
# install.packages("stargazer")
# install.packages("ggplot2")
# install.packages("knitr")
# install.packages("gt")
# install.packages("gtsummary")

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
  inner_join(KOB, by = "case") %>%                      # Zusammenführung von MCB- und KOB-Daten
  mutate(GBZeitpkt = gbja*12 + gbmo,                    # Alter in Monaten
         RTZeitpkt = rtbej*12 + rtbem,                  # Rentenbeginn in Monaten
         TDZeitpkt = tddtj*12 + tddtm,                  # Todeszeitpunkt in Monaten
         AlterRehaStart = RehaStart - GBZeitpkt,        # Alter in Monaten bei Beginn der ersten med. Reha
         ZeitRehaStartTod = case_when(                  # Zeitraum zwischen Beginn der ersten med. Reha und dem Todeszeitpunkt in Monaten
           tddtj %in% c(1:9998) ~ TDZeitpkt - RehaStart)) %>%
  filter(RTZeitpkt < RehaStart,                         # Keine Rente vor Beginn der ersten med. Reha
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
    Rente = if_else(
      rtbej > 0, 1, 0),
    Tod = if_else(
      tddtj > 0, 1, 0)
    ) %>%
  remove_value_labels(mcdg1_icd = c(1:6, 32:999)) # Entfernen der leeren Diagnoselabels für gtsummary Ausgabe - Behält nur die lables für mmcdg1_icd %in% c(7:31), siehe Subseterstellung



### Datensatzbeschreibung

CARES_rtw %>%
  select(Geschlecht,
         mcdg1_icd,
         AlterRehaStart_Jahr,
         mcfmsd,
         Rente,
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
  mutate_at(vars("bd", "sb", "bb", "mcahb", "mcfmsd", "mcwhot_ostwest", "mcwhot_skt", "mcmsot_bland", "mcstbf", "mcaiufzt", "mcleft_lb", "mcleft_at", "mcdg1_icd"), haven::as_factor) %>%
  tbl_summary(
    by = Geschlecht,
    missing_text = "Fehlende Beobachtungen",
    statistic = all_continuous() ~ "{mean} ({sd})",
    label = list(
      AlterRehaStart_Jahr ~ "Alter bei Beginn der ersten Rehabilitationsleistung",
      Rente ~ "Versichertenrente",
      Tod ~ "Versterben im Beobachtungszeitraum",
      ZeitRehaStartTod ~ "Zeitraum vom Beginn der ersten med. Reha und dem Todeszeitpunkt (Monate)",
      mcaqda ~ "Zeitraum von der Antragstellung eines Antrags auf medizinische Rehabilitation bis zu seiner Bewilligung in Tagen"),
    digits = list(
      AlterRehaStart_Jahr ~ 1,
      ZeitRehaStartTod ~ 0)) %>% 
  add_overall() %>% 
  modify_header(label = "**Variable**") %>% 
  bold_labels()





##### Count #####


x <- CARES_rtw %>% 
  group_by(mcdg1_icd) %>% 
  summarise(count = n())

view(x)




x <- CARES_rtw %>% 
  group_by(mcdg1_icd) %>% 
  summarise(count = n())

view(x)








### evtl. noch nützlich

## Geburtsjahr Plot
# Geburtsjahr <- cancer_CASES_MCB_KOB %>%
#   filter(gbja>0)
# 
# hist(Geburtsjahr$gbja, breaks=100)



# dataset with few variables to summarize
x1 <- final_CARES %>% 
  select(AlterRehaStart, Geschlecht)

table1 <- tbl_summary(x1)



# simple count
CARES_rtw %>% 
  count(mcaift)

# simple mean
final_CARES %>%
  group_by(GroupeVariable) %>%
  summarize(Mean = mean(MeanVariable, na.rm=TRUE))



# () tlrt
x <- final_CARES %>%
  group_by(tlrt) %>% 
  summarise(mean = round(mean(AlterRehaStart, na.rm=TRUE), 2),
            min = min(AlterRehaStart, na.rm=TRUE),
            max = max(AlterRehaStart, na.rm=TRUE),
            count = n())


##### Notes #####

## Rehabeginn abhängig vom Monat. Um Neujahr weniger Rehaantritte (Teilweise bis zu 50 %).
## breaks=96 weil 8 Jahre = 96 Monate
# hist(cancer_CASES_MCB$RehaStart, breaks=96)
# 
# RehaBeginnMonat <- cancer_CASES_MCB_KOB %>% 
#   count(RehaStart, mcbemsm, wt=korrektur)
# 
# ggplot(RehaBeginnMonat, aes( RehaStart, n, colour = mcbemsm) ) +
#   geom_point( ) + geom_line( )



