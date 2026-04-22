# Griselda dataset with doses

library(dplyr)
library(writexl)
library(stringr)

# Subset for Fixed dosing only
Griselda <- na.omit(griselda_full[griselda_full$Dosing == "Fixed",])

# Remove hyphen and everything after it from Dose column
Griselda$Dose <- sub("-.*", "", Griselda$Dose)

# Rename variables
Griselda <- Griselda %>%
  rename(mean_base = Mean...17)

Griselda <- Griselda %>%
  rename(mean_fu = Mean...20)

# Convert both columns to numeric
Griselda$mean_fu <- as.numeric(Griselda$mean_fu)

Griselda$mean_base <- as.numeric(Griselda$mean_base)

# Set missing values * to missing
Griselda[Griselda == "*"] <- NA

# Calculate response rate
Griselda$Responders <- as.numeric(Griselda$Responders)

Griselda$No_randomised <- as.numeric(Griselda$No_randomised)

Griselda$response_rate <- Griselda$Responders / Griselda$No_randomised

# Calculate remission rate
Griselda$Remitters <- as.numeric(Griselda$Remitters)

Griselda$remission_rate <- Griselda$Remitters / Griselda$No_randomised


### Try to detemine whether endpoint is follow-up or change score ###

table(is.na(Griselda$mean_base))

# Create a variable to flag the outcome type
Griselda$outcome_type <- ifelse(Griselda$mean_fu < 0, "change", "unknown")

# Look at the distribution
table(Griselda$outcome_type, useNA = "always")

# For HAMD studies, if calculated_endpoint > 54, mean_fu must be endpoint

Griselda$calculated_endpoint <- Griselda$mean_base + Griselda$mean_fu

Griselda$outcome_type <- ifelse(
  Griselda$outcome_type == "unknown" &
    grepl("HAMD", Griselda$`Definition of response`, ignore.case = TRUE) & 
    !is.na(Griselda$calculated_endpoint) &
    Griselda$calculated_endpoint > 54,
  "end",
  Griselda$outcome_type
)

# For MADRS studies, if calculated_endpoint > 60, mean_fu must be endpoint
Griselda$outcome_type <- ifelse(
  Griselda$outcome_type == "unknown" &
    grepl("MADRS", Griselda$`Definition of response`, ignore.case = TRUE) & 
    !is.na(Griselda$calculated_endpoint) &
    Griselda$calculated_endpoint > 60,
  "end",
  Griselda$outcome_type
)

# If response rate is high (say >40%) but mean_fu is large (>15), 
# then mean_fu is probably an endpoint (low final score = good)

# If response rate is high (>40%) and mean_fu is negative/small,
# then mean_fu is probably a change score (big reduction = good)
Griselda$outcome_type <- ifelse(
  Griselda$outcome_type == "unknown" &
    !is.na(Griselda$response_rate) &
    !is.na(Griselda$mean_fu) &
    Griselda$response_rate > 0.3 &  # Decent response rate
    Griselda$mean_fu > 12,           # But mean_fu is moderate-high
  "end",                           # Must be endpoint (improved TO ~12)
  Griselda$outcome_type
)

Griselda$outcome_type <- ifelse(
  Griselda$outcome_type == "unknown" &
    !is.na(Griselda$response_rate) &
    !is.na(Griselda$mean_fu) &
    Griselda$response_rate > 0.3 &   # Decent response rate
    Griselda$mean_fu > 0 & Griselda$mean_fu < 10,  # But mean_fu is small positive
  "change",                        # Must be change (improved BY ~8)
  Griselda$outcome_type
)

# High remission rate + moderate mean_fu → endpoint
Griselda$outcome_type <- ifelse(
  Griselda$outcome_type == "unknown" &
    !is.na(Griselda$remission_rate) &
    Griselda$remission_rate > 0.25 &  # Decent remission rate
    Griselda$mean_fu > 8 & Griselda$mean_fu < 15,  # Moderate values
  "end",
  Griselda$outcome_type
)

# Low/moderate remission + higher mean_fu → endpoint
Griselda$outcome_type <- ifelse(
  Griselda$outcome_type == "unknown" &
    !is.na(Griselda$remission_rate) &
    Griselda$remission_rate < 0.30 &
    Griselda$mean_fu > 15 & Griselda$mean_fu < 30,
  "end",
  Griselda$outcome_type
)

# For each study, find the most common outcome_type (excluding "unknown")
Griselda <- Griselda %>%
  group_by(StudyID) %>%
  mutate(
    # Get the known outcome type from other arms in this study
    study_outcome_type = {
      known_types <- outcome_type[outcome_type != "unknown"]
      if(length(known_types) > 0) known_types[1] else NA
    }
  ) %>%
  ungroup()

# Apply that type to the unknowns in the same study
Griselda$outcome_type <- ifelse(
  Griselda$outcome_type == "unknown" & !is.na(Griselda$study_outcome_type),
  Griselda$study_outcome_type,
  Griselda$outcome_type
)

# See how many left are unknown
table(Griselda$outcome_type, useNA = "always")

# 2 Left, inspected and almost certainly endpoints
Griselda$outcome_type <- ifelse(
  Griselda$StudyID == "Goodarzi 2015(IRCT2012101811155N1)",
  "end",
  Griselda$outcome_type
)

# Check it worked
table(Griselda$outcome_type, useNA = "always")

# Rename mean_fu to be clearer
Griselda$mean_reported <- Griselda$mean_fu

# Create separate change and endpoint variables
Griselda$mean_change <- ifelse(Griselda$outcome_type == "change", 
                               Griselda$mean_reported,
                               Griselda$mean_base - Griselda$mean_reported)

Griselda$mean_end <- ifelse(Griselda$outcome_type == "end",
                            Griselda$mean_reported,
                            Griselda$mean_base + Griselda$mean_reported)

# Same for SDs if you need them
Griselda$sd_change <- Griselda$  # Assuming sd_fu corresponds to reported outcome
Griselda$sd_end <- Griselda$sd_fu

# View the result
head(Griselda)

# Create a lookup vector for drug classes
drug_classes <- c(
  "amitriptyline"    = "TCA",
  "clomipramine"     = "TCA",
  "nefazodone"       = "SARI",
  "trazodone"        = "SARI",
  "citalopram"       = "SSRI",
  "escitalopram"     = "SSRI",
  "fluoxetine"       = "SSRI",
  "fluvoxamine"      = "SSRI",
  "paroxetine"       = "SSRI",
  "sertraline"       = "SSRI",
  "desvenlafaxine"   = "SNRI",
  "duloxetine"       = "SNRI",
  "levomilnacipran"  = "SNRI",
  "milnacipran"      = "SNRI",
  "venlafaxine"      = "SNRI",
  "agomelatine"      = "Atypical",
  "bupropion"        = "Atypical",
  "mirtazapine"      = "Atypical",
  "reboxetine"       = "Atypical",
  "vilazodone"       = "Atypical",
  "vortioxetine"     = "Atypical",
  "placebo"          = "Placebo"
)

# Add Class variable to Griselda
Griselda$Class <- drug_classes[Griselda$Drug]

# Check the result
head(Griselda)

# Subset using bracket notation
Griselda_SSRI <- Griselda[Griselda$Class == "SSRI", ]

# View the result
head(Griselda_SSRI)   # shows first few rows

# Compare Griselda to SSRI, eventually change names in griselda to match 
print(ssri)
print(Griselda_SSRI, width = Inf)

# Check column names
names(Griselda_SSRI)
names(ssri)

# Check classes/types
sapply(Griselda_SSRI[, c("Drug", "Dose", "No_randomised", "Responders", "Weeks")], class)
sapply(ssri[, c("agent", "dose", "n", "r", "weeks")], class)

summary(Griselda_SSRI[, c("Dose", "No_randomised", "Responders", "Weeks")])
summary(ssri[, c("dose", "n", "r", "weeks")])

sort(unique(Griselda_SSRI$Drug))
sort(unique(ssri$agent))

library(dplyr)

# Sum of participants and responders by drug/agent
Griselda_SSRI %>% group_by(Drug) %>% summarise(total_n = sum(No_randomised, na.rm = TRUE),
                                               total_r = sum(Responders, na.rm = TRUE))

ssri %>% group_by(agent) %>% summarise(total_n = sum(n, na.rm = TRUE),
                                       total_r = sum(r, na.rm = TRUE))

# Looks good! Clean Griselda so it matches more
library(dplyr)
library(stringr)

# 1️⃣ Start with a clean copy and rename variables
Griselda_clean <- Griselda %>%
  rename(
    agent = Drug,
    dose = Dose,
    n = No_randomised,
    r = Responders,
    weeks = Weeks
  ) %>%
  select(StudyID, Class, agent, dose, n, r, weeks, Year_Published) %>%
  # Convert numeric variables
  mutate(
    dose = as.numeric(dose),
    n = as.numeric(n),
    r = as.numeric(r),
    weeks = as.numeric(weeks)
  ) %>%
  # 2️⃣ Filter out unpublished studies
  mutate(Year_Published_clean = str_to_lower(Year_Published)) %>%
  filter(
    !str_detect(Year_Published_clean, "unpublish|in press|n/?a")
  ) %>%
  select(-Year_Published_clean, -Year_Published) %>%
  
  group_by(StudyID) %>%
  mutate(
    dose = if_else(agent == "Placebo" & any(agent != "Placebo"), 0, dose)
  ) %>%
  ungroup()  %>%
  
  
  # 3️⃣ Remove any rows with NA in key variables
  filter(
    !is.na(agent) & !is.na(dose) & !is.na(n) & !is.na(r) & !is.na(weeks)
  )

# 4️⃣ Quick check
str(Griselda_clean)
summary(Griselda_clean)

#Format  for MBNMA dose:
# Remove names from Class
Griselda_clean$Class <- unname(Griselda_clean$Class)

# Rename StudyID -> studyID
Griselda_clean <- Griselda_clean %>%
  rename(studyID = StudyID)

# Check column names
colnames(Griselda_clean)
# Should include: studyID, agent, dose, r, n, weeks, Class

# Vector of studyIDs to remove
single_arm_studies <- c(
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
)

# Trim whitespace from studyID
Griselda_clean <- Griselda_clean %>%
  mutate(studyID = str_trim(studyID))

# Also remove any problematic single-arm studies
single_arm_studies <- str_trim(single_arm_studies)

Griselda_clean <- Griselda_clean %>%
  filter(!studyID %in% single_arm_studies)

# Conversion ratios: multiply raw mg/day by ratio to get fluoxetine-equivalent
drug_ratios <- c(
  "amitriptyline"    = 0.327,
  "clomipramine"     = 0.344,
  "nefazodone"       = 0.075,
  "trazodone"        = 0.100,
  "citalopram"       = 1.000,
  "escitalopram"     = 2.222,
  "fluoxetine"       = 1.000,
  "fluvoxamine"      = 0.279,
  "paroxetine"       = 1.176,
  "sertraline"       = 0.406,
  "desvenlafaxine"   = 0.400,
  "duloxetine"       = 0.333,
  "levomilnacipran"  = 0.200,
  "milnacipran"      = 0.200,
  "venlafaxine"      = 0.268,
  "agomelatine"      = 0.752,
  "bupropion"        = 0.115,
  "mirtazapine"      = 0.784,
  "reboxetine"       = 3.448,
  "vilazodone"       = 1.000,
  "vortioxetine"     = 2.000
)

Griselda_clean$dose_FE <- Griselda_clean$dose * drug_ratios[Griselda_clean$agent]

Griselda_clean <- Griselda_clean %>%
  mutate(r = round(r),
         n = round(n),
         dose_FE = ifelse(agent == "placebo", 0, dose_FE))
############TIAN########################

Tian <- Tian %>%
  mutate(Dose_minutes = Min_day_avg * Frequency # Total min/week (sensitivity)
  )

