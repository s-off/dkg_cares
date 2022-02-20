#################
### DKG CARES ###
#################

# Packages #

# install.packages("haven")
# install.packages("tidyverse")
# install.packages("stargazer")
# install.packages("ggplot2")
# install.packages("knitr")
# install.packages("gtsummary")

library(haven)
library(tidyverse)
library(stargazer)
library(knitr)
library(gtsummary)

# Import Data

MCB <- read_spss("input/SUFRSDLV17MCB.sav")

KOB <- read_spss("input/SUFRSDLV17KOB.sav")

BYB <- read_spss("input/SUFRSDLV17BYB.sav")

### Subset mit Cases, die im Berichtzeitraum eine medizinische Rehabilitation auf Grund einer Krebsdiagnose begonnen haben
  ## Es wird nur die erste Episode (erste Rehabilitation auf Grund einer Krebsdiagnose betrachtet)

cancer_CASES_MCB <- MCB %>% 
  filter(mcdggr == 26,                           # Bewilligungsdiagnosengruppe: Neubildungen (ICD10-Nr. C00 - D48)
         mcdg1_icd %in% c(7:31),                 # 1. med. Entlassungsdiagnose ICD-10: C00-96
         mcmsat %in% c(31,32),                   # Bewilligte Maßnahmeart: Ca(Krebs)-Reha-Leistung nach § 15 SGB VI oder nach § 31 Abs. 1 Nr. 2 SGB VI
         mceafo == 1) %>%                        # Entlassungsform: regulär
  mutate(RehaStart=mcbemsj*12+mcbemsm) %>%       # Erste Episode wird durch RehaStart (Beginn der medizinischen Reha-Leistung) identifiziert
  arrange(case, RehaStart) %>% 
  distinct(case, .keep_all = TRUE)


### Zusammenführung von MCB- und KOB-Daten für alle Krebspatienten: Inner join cancer_IDs (MCB selection) with KOB

cancer_CASES_MCB_KOB <- inner_join(cancer_CASES_MCB, KOB, by = "case")

### Erstellen der Variablen

final_CARES <- cancer_CASES_MCB_KOB %>% 
  mutate(
    AlterRehaStart = case_when(
      gbja > 0 ~ mcbemsj-gbja),
    Geschlecht = case_when(
      sex == 1 ~ "m",
      sex == 2 ~ "f")
    )


### Anwendung der Ein- und Ausschlusskriterien
  # ??? Alter < 18 Jahre? Kinder ausschließen?
  # ??? mcwhb==1
    # erste Episode = Wiederholungsleistung ausschließen (n=13, meist zu Beginn des Beobachtungszeitraums)
  # ??? mcrar
    # Reha-Leistung aus dem Rentenverfahren: Umdeutung eines Rentenantrags in einen
      # Reha-Antrag handelt (§ 116 SGB VI). D. h. die Reha-Leistung erfolgt entweder nach
      # Rentenantraganstellung oder nachdem ein Rentenantrag abgelehnt wurde oder bei
      # laufendem Bezug einer Rente wegen verminderter Erwerbsfähigkeit.
  # (x) MCEAFO: ggf. nach Entlassungsform Unterscheiden
  # () gbja: n=14 haben "0"
  # () tlrt
                x <- final_CARES %>%
                  group_by(tlrt) %>% 
                  summarise(mean = round(mean(AlterRehaStart, na.rm=TRUE), 2),
                            min = min(AlterRehaStart, na.rm=TRUE),
                            max = max(AlterRehaStart, na.rm=TRUE),
                            count = n())
  # () Was mach ich mit Verstorbenen?
        # siehe Notion
  # ()
                
              


                

x <- final_CARES %>% 
  group_by(mceafo) %>% 
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
final_CARES %>% 
  count(mcaift)

# simple mean
final_CARES %>%
  group_by(GroupeVariable) %>%
  summarize(Mean = mean(MeanVariable, na.rm=TRUE))





### 1. Zusammensetzung der Stichprobe

Diagnose <- cancer_CASES_MCB_KOB %>% 
  count(sex, mcdg1_icd, wt = korrektur)

stargazer(Diagnose, type = "html", summary = FALSE, out = "Krebsdiagonsen nach Geschlecht")

mcaqat <- cancer_IDs_KOB %>% 
  count(mcdg1_icd, mcaqat) %>% # not weighted (, wt = korrektur)
  spread(mcaqat, n, fill = 0)

stargazer(mcaqat, type = "html", summary = FALSE, out = "Art der beantragten medizinischen Reha-Leistung nach Diagnoseschlüssel")


mcmsat <- cancer_IDs_KOB %>% 
  count(mcdg1_icd, mcmsat) %>% # not weighted (, wt = korrektur)
  spread(mcmsat, n, fill = 0)

stargazer(mcmsat, type = "html", summary = FALSE, out = "Bewilligte Maßnahmeart nach Diagnoseschlüssel")



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



