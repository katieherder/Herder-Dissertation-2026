# Main project directory
setwd("C:/Users/katie/Desktop/Depression-NMA")

# === Quick check ===
head(Griselda_clean)
str(Griselda_clean)

# Quick summary
nrow(Griselda_clean)                         # total arms
length(unique(Griselda_clean$studyID))       # total studies
table(Griselda_clean$agent)                  # arms per agent
sort(unique(Griselda_clean$dose))            # all dose levels

Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(agent) %>%
  summarise(
    n_arms = n(),
    n_studies = n_distinct(studyID),
    n_doses = n_distinct(dose),
    doses = paste(sort(unique(dose)), collapse = ", ")
  ) %>% print(n = Inf)

# Check for head-to-head trials (no placebo arm in the study)
Griselda_clean %>%
  group_by(studyID) %>%
  filter(all(agent != "placebo")) %>%
  distinct(studyID, agent)

# === Overall summary ===
cat("Total studies:", n_distinct(Griselda_clean$studyID), "\n")
cat("Total arms:", nrow(Griselda_clean), "\n")
cat("Total participants:", sum(Griselda_clean$n), "\n")
cat("Median sample size per arm:", median(Griselda_clean$n), "\n")

# === Arms per study ===
cat("\n=== Arms per study ===\n")
Griselda_clean %>% group_by(studyID) %>% summarise(k = n()) %>% pull(k) %>% table()

# === Agent summary ===
cat("\n=== By agent ===\n")
Griselda_clean %>%
  group_by(agent) %>%
  summarise(
    n_arms = n(),
    n_studies = n_distinct(studyID),
    n_participants = sum(n),
    median_n = median(n),
    n_doses = n_distinct(dose),
    dose_range = paste0(min(dose), "--", max(dose)),
    doses = paste(sort(unique(dose)), collapse = ", "),
    overall_response = round(sum(r) / sum(n), 3)
  ) %>%
  print(width = Inf)

# === Studies without placebo ===
cat("\n=== Studies without placebo ===\n")
Griselda_clean %>%
  group_by(studyID) %>%
  filter(all(agent != "placebo")) %>%
  print(n = Inf)

# === Studies comparing different active agents ===
cat("\n=== Studies comparing different active agents ===\n")
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(studyID) %>%
  filter(n_distinct(agent) > 1) %>%
  arrange(studyID) %>%
  select(studyID, agent, dose) %>%
  print(n = Inf)

# === Follow-up duration (weeks) ===
cat("\n=== Weeks ===\n")
table(Griselda_clean$weeks)

# === Response rates by agent and dose (exclude placebo) ===
cat("\n=== Response rate by agent and dose ===\n")
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(agent, dose) %>%
  summarise(
    k = n(),
    total_n = sum(n),
    total_r = sum(r),
    response_rate = round(total_r / total_n, 3),
    .groups = "drop"
  ) %>%
  arrange(agent, dose) %>%
  print(n = Inf)

# === Placebo response rate ===
cat("\n=== Placebo response rate ===\n")
Griselda_clean %>%
  filter(agent == "placebo") %>%
  summarise(
    k = n(),
    total_n = sum(n),
    total_r = sum(r),
    response_rate = round(total_r / total_n, 3)
  )

# === Placebo response rate === # CLASS LEVEL! 
# Class-level summary table for main text
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(Class) %>%
  summarise(
    n_agents = n_distinct(agent),
    n_studies = n_distinct(studyID),
    n_arms = n(),
    n_participants = sum(n),
    n_dose_levels = n_distinct(paste(agent, dose)),
    median_arm_n = median(n),
    overall_response = round(sum(r) / sum(n), 3)
  )

# Then a separate summary flagging data-sparse agents
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(agent, Class) %>%
  summarise(
    n_arms = n(),
    n_studies = n_distinct(studyID),
    n_doses = n_distinct(dose),
    .groups = "drop"
  ) %>%
  filter(n_arms <= 3 | n_doses <= 1) %>%
  arrange(n_arms)

# Head-to-head structure by class pair
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(studyID) %>%
  filter(n_distinct(Class) > 1) %>%
  summarise(classes = paste(sort(unique(Class)), collapse = " vs ")) %>%
  count(classes, sort = TRUE)

# 1. Cross-class comparison counts (motivates why class pooling matters)
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(studyID) %>%
  filter(n_distinct(Class) > 1) %>%
  ungroup() %>%
  summarise(
    n_cross_class_studies = n_distinct(studyID),
    n_cross_class_arms = n()
  )

# 2. Within-class comparison counts
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(studyID) %>%
  filter(n_distinct(agent) > 1, n_distinct(Class) == 1) %>%
  ungroup() %>%
  count(Class)

# 3. Dose coverage by class (are some classes too sparse for Emax?)
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(Class, agent) %>%
  summarise(
    n_doses = n_distinct(dose),
    dose_range_FE = paste0(
      round(min(dose) * drug_ratios[agent[1]], 1), "--",
      round(max(dose) * drug_ratios[agent[1]], 1)
    ),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# === Network plot ===
griselda_network <- mbnma.network(data.ab = Griselda_clean)
jpeg("figures/Figure_Griselda_network.jpeg", width = 1800, height = 1500, res = 300)
par(mar = c(2, 2, 2, 2))
plot(griselda_network, main = "Griselda Network")
dev.off()

# Tian data
library(tidyverse)

# 1. Verify number of studies and comparisons
Tian %>%
  summarise(
    n_studies = n_distinct(Study),
    n_comparisons = n(),
    n_control_comparisons = sum(Agent == "CON"),
    n_exercise_comparisons = sum(Agent != "CON")
  )



# 2. How many studies have more than 2 arms?
multi_arm_count <- Tian %>%
  group_by(Study) %>%
  summarise(n_arms = n()) %>%
  summarise(
    total_studies = n(),
    two_arm_studies = sum(n_arms == 2),
    multi_arm_studies = sum(n_arms > 2)
  )

print(multi_arm_count)

# 3. Verify duration
Tian %>%
  filter(Agent != "CON") %>%
  summarise(
    min_freq = min(Exact_Length, na.rm = TRUE),
    max_freq = max(Exact_Length, na.rm = TRUE),
    median_freq = median(Exact_Length, na.rm = TRUE)
  )

# 4. Verify frequency
Tian %>%
  filter(Agent != "CON") %>%
  summarise(
    min_freq = min(Frequency, na.rm = TRUE),
    max_freq = max(Frequency, na.rm = TRUE),
    median_freq = median(Frequency, na.rm = TRUE)
  )

# 5. Verify minutes per day
Tian %>%
  filter(Agent != "CON") %>%
  summarise(
    min_freq = min(Min_day_avg, na.rm = TRUE),
    max_freq = max(Min_day_avg, na.rm = TRUE),
    median_freq = median(Min_day_avg, na.rm = TRUE)
  )

# 5. Count comparisons by exercise modality
Tian %>%
  count(Agent) %>%
  arrange(desc(n))

######TABLE 3.1###############
# 6. Verify dose distribution
Tian %>%
  filter(Agent != "CON") %>%
  summarise(
    min_dose = min(Dose, na.rm = TRUE),
    max_dose = max(Dose, na.rm = TRUE),
    median_dose = median(Dose, na.rm = TRUE),
    q25_dose = quantile(Dose, 0.25, na.rm = TRUE),
    q75_dose = quantile(Dose, 0.75, na.rm = TRUE)
  )

# 7. Dose distribution by modality
Tian %>%
  group_by(Agent) %>%
  summarise(
    n = n(),
    total_participants = sum(N, na.rm = TRUE),
    median_N = median(N, na.rm = TRUE),
    median_dose = median(Dose, na.rm = TRUE),
    min_dose = min(Dose, na.rm = TRUE),
    max_dose = max(Dose, na.rm = TRUE),
    q25_dose = quantile(Dose, 0.25, na.rm = TRUE),
    q75_dose = quantile(Dose, 0.75, na.rm = TRUE)
  )

###############################################
# 8. Verify sample sizes
Tian %>%
  summarise(
    total_participants = sum(N, na.rm = TRUE),
    median_n_per_comparison = median(N, na.rm = TRUE)
  )

# 9. Check for head-to-head trials (studies with no control arm)
Tian %>%
  group_by(Study) %>%
  summarise(has_control = any(Agent == "CON")) %>%
  filter(!has_control) %>%
  nrow()

# 10. Check for dose-comparison studies (same study, same agent, different doses)
dose_comparison_studies <- Tian %>%
  group_by(Study, Agent) %>%
  filter(Agent != "CON") %>%
  summarise(
    n_arms = n(),
    n_unique_doses = n_distinct(Dose),
    .groups = "drop"
  ) %>%
  filter(n_unique_doses > 1)  # Same agent, but different doses
nrow(dose_comparison_studies)

# 11.  How many unique scales?
Tian %>%
  summarise(n_unique_scales = n_distinct(Scale, na.rm = TRUE))

## Distribution of scales across studies
scale_distribution <- Tian %>%
  group_by(Scale) %>%
  summarise(
    n_arms = n(),
    n_studies = n_distinct(Study)
  ) %>%
  arrange(desc(n_studies))

print(scale_distribution)

# 3. How many studies used each scale?
studies_per_scale <- Tian %>%
  distinct(Study, Scale) %>%
  count(Scale, name = "n_studies") %>%
  arrange(desc(n_studies))

print(studies_per_scale)

# 4. Check if any studies used multiple scales (should be 0 or very few)
studies_multiple_scales <- Tian %>%
  group_by(Study) %>%
  summarise(n_scales = n_distinct(Scale, na.rm = TRUE)) %>%
  filter(n_scales > 1)

print(studies_multiple_scales)
nrow(studies_multiple_scales)


############################################
##########SSRI DATA#########################
library(MBNMAdose)
library(gemtc)
library(magick)
library(dplyr)

head(ssri)
str(ssri)

# Quick summary
nrow(ssri)                          # total arms
length(unique(ssri$studyID))        # total studies
table(ssri$agent)                   # arms per agent
sort(unique(ssri$dose))             # all dose levels

# Dose ranges by agent
ssri %>%
  filter(agent != "Placebo") %>%
  group_by(agent) %>%
  summarise(
    n_arms = n(),
    n_studies = n_distinct(studyID),
    n_doses = n_distinct(dose),
    doses = paste(sort(unique(dose)), collapse = ", ")
  )

# Check for head-to-head trials (no placebo arm in the study)
ssri %>%
  group_by(studyID) %>%
  filter(all(agent != "placebo")) %>%
  distinct(studyID, agent)

ssri %>% filter(studyID == 161001)


# === Overall summary ===
cat("Total studies:", n_distinct(ssri$studyID), "\n")
cat("Total arms:", nrow(ssri), "\n")
cat("Total participants:", sum(ssri$n), "\n")
cat("Median sample size per arm:", median(ssri$n), "\n")

# === Arms per study ===
cat("\n=== Arms per study ===\n")
ssri %>% group_by(studyID) %>% summarise(k = n()) %>% pull(k) %>% table()

# === Agent summary ===
cat("\n=== By agent ===\n")
ssri %>%
  group_by(agent) %>%
  summarise(
    n_arms = n(),
    n_studies = n_distinct(studyID),
    n_participants = sum(n),
    median_n = median(n),
    n_doses = n_distinct(dose),
    dose_range = paste0(min(dose), "--", max(dose)),
    doses = paste(sort(unique(dose)), collapse = ", "),
    overall_response = round(sum(r) / sum(n), 3)
  ) %>%
  print(width = Inf)

# === Head-to-head check ===
cat("\n=== Studies without placebo ===\n")
ssri %>%
  group_by(studyID) %>%
  filter(all(agent != "placebo")) %>%
  print(n = Inf)

# === Studies with multiple active agents ===
cat("\n=== Studies comparing different SSRIs ===\n")
ssri %>%
  filter(agent != "placebo") %>%
  group_by(studyID) %>%
  filter(n_distinct(agent) > 1) %>%
  arrange(studyID) %>%
  select(studyID, agent, dose) %>%
  print(n = Inf)

# === Follow-up duration ===
cat("\n=== Weeks ===\n")
table(ssri$weeks)

# === Risk of bias ===
cat("\n=== Risk of bias ===\n")
ssri %>% distinct(studyID, bias) %>% pull(bias) %>% table()

# === Response rates by agent (for context) ===
cat("\n=== Response rate by agent and dose ===\n")
ssri %>%
  filter(agent != "placebo") %>%
  group_by(agent, dose) %>%
  summarise(
    k = n(),
    total_n = sum(n),
    total_r = sum(r),
    response_rate = round(total_r / total_n, 3),
    .groups = "drop"
  ) %>%
  arrange(agent, dose) %>%
  print(n = Inf)

# === Placebo response rate ===
cat("\n=== Placebo response rate ===\n")
ssri %>%
  filter(agent == "placebo") %>%
  summarise(
    k = n(),
    total_n = sum(n),
    total_r = sum(r),
    response_rate = round(total_r / total_n, 3)
  )

# === Network Plot ===
ssri_network <- mbnma.network(data.ab = ssri)

jpeg("figures/Figure 4.1.jpeg", width = 1800, height = 1500, res = 300)
par(mar = c(2, 2, 2, 2))
plot(ssri_network, main = "")
dev.off()