################################################################################
# Chapter 4: Clean Data - GRISELDA Network
################################################################################

# Griselda dataset with doses -------------------------------------------------

library(dplyr)
library(stringr)

# Drug lookup table: class + fluoxetine-equivalent conversion ratio -----------
drug_info <- tibble::tribble(
  ~agent,             ~Class,     ~ratio,
  "amitriptyline",    "TCA",      0.327,
  "clomipramine",     "TCA",      0.344,
  "nefazodone",       "SARI",     0.075,
  "trazodone",        "SARI",     0.100,
  "citalopram",       "SSRI",     1.000,
  "escitalopram",     "SSRI",     2.222,
  "fluoxetine",       "SSRI",     1.000,
  "fluvoxamine",      "SSRI",     0.279,
  "paroxetine",       "SSRI",     1.176,
  "sertraline",       "SSRI",     0.406,
  "desvenlafaxine",   "SNRI",     0.400,
  "duloxetine",       "SNRI",     0.333,
  "levomilnacipran",  "SNRI",     0.200,
  "milnacipran",      "SNRI",     0.200,
  "venlafaxine",      "SNRI",     0.268,
  "agomelatine",      "Atypical", 0.752,
  "bupropion",        "Atypical", 0.115,
  "mirtazapine",      "Atypical", 0.784,
  "reboxetine",       "Atypical", 3.448,
  "vilazodone",       "Atypical", 1.000,
  "vortioxetine",     "Atypical", 2.000,
  "placebo",          "Placebo",  0
)

# Single-arm studies to exclude -----------------------------------------------
single_arm_studies <- str_trim(c(
  "Blacker1988",
  "CL3-20098-036",
  "CL3-20098-062",
  "DeRonchi1998",
  "Dierick1996",
  "Judd1993",
  "Kennedy2014 (EudraCT 2009-011238-84, CL3-20098-069)",
  "Khan2007 (SCT-MD-23) (NCT00108979)",
  "Mao2015 (NCT01098318)",
  "Marchesi1998",
  "Noguera1991",
  "OntiverosSanchez1998",
  "Ventura2007 (SCT-MD-18)",
  "Versiani1999"
))

# Initial cleaning: fixed dosing only, rename, fix types ----------------------
Griselda <- griselda_full %>%
  filter(Dosing == "Fixed") %>%
  rename(
    mean_base = `Mean...17`,
    mean_fu   = `Mean...20`
  ) %>%
  mutate(
    # Strip dose range suffix (e.g. "20-40" -> "20")
    Dose = sub("-.*", "", Dose),
    # "*" is used as missing
    across(everything(), ~ if_else(.x == "*", NA, .x)),
    across(c(Dose, mean_base, mean_fu, Responders, No_randomised),
           as.numeric),
    response_rate  = Responders / No_randomised,
    remission_rate = as.numeric(Remitters) / No_randomised
  ) %>%
  # Drop rows with NAs in the key fields used downstream
  filter(if_all(c(Drug, Dose, No_randomised, Responders, Weeks), ~ !is.na(.x)))


# Build Griselda_clean for MBNMAdose ------------------------------------------
Griselda_clean <- Griselda %>%
  rename(
    studyID = StudyID,
    agent   = Drug,
    dose    = Dose,
    n       = No_randomised,
    r       = Responders,
    weeks   = Weeks
  ) %>%
  mutate(
    studyID = str_trim(studyID),
    across(c(dose, n, r, weeks), as.numeric)
  ) %>%
  # Drop unpublished studies
  filter(!str_detect(str_to_lower(Year_Published), "unpublish|in press|n/?a")) %>%
  # Drop single-arm studies
  filter(!studyID %in% single_arm_studies) %>%
  # Placebo arms always get dose 0
  group_by(studyID) %>%
  mutate(dose = if_else(agent == "placebo" & any(agent != "placebo"), 0, dose)) %>%
  ungroup() %>%
  # Attach Class and fluoxetine-equivalent dose from the lookup table
  left_join(drug_info, by = "agent") %>%
  mutate(
    dose_FE = if_else(agent == "placebo", 0, dose * ratio),
    r = round(r),
    n = round(n)
  ) %>%
  select(studyID, Class, agent, dose, dose_FE, n, r, weeks) %>%
  filter(if_all(c(agent, dose, n, r, weeks), ~ !is.na(.x)))

# Quick check
str(Griselda_clean)
summary(Griselda_clean)