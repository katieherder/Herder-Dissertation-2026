################################################################################
# Chapter 3: Network characteristics — Tian Dataset
#
# Produces all descriptive numbers reported in the "Network Characteristics
# and Connectivity" section: overall totals, dose distribution by modality,
# arms per dose category, study design breakdown, intervention parameters,
# outcome scales, and network structure.
#
# First run:
#   - 01_load_data
#   - 02_clean_data
################################################################################

library(dplyr)

# --- Overall totals --------------------------------------------------------
# Reported in text: "45 RCTs ... 97 total arms ... 2,964 participants ...
# 58 exercise intervention arms and 39 control arms"
cat("=== Overall ===\n")
cat("Studies:       ", n_distinct(Tian$Study), "\n")
cat("Arms:          ", nrow(Tian), "\n")
cat("Participants:  ", sum(Tian$N), "\n")
cat("Exercise arms: ", sum(Tian$Agent != "CON"), "\n")
cat("Control arms:  ", sum(Tian$Agent == "CON"), "\n")
cat("Median N/arm:  ", median(Tian$N), "\n")

# --- Arms per study --------------------------------------------------------
# Reported in text: "39 studies were two-arm trials, and six studies
# included three or more arms"
cat("\n=== Arms per study ===\n")
Tian %>%
  group_by(Study) %>%
  summarise(k = n(), .groups = "drop") %>%
  summarise(
    total     = n(),
    two_arm   = sum(k == 2),
    multi_arm = sum(k > 2)
  ) %>%
  print()

# --- Arms per modality ------------------------------------------------------
# Reported in text: "AE: k = 23 arms ... RE: k = 6 ... ME: k = 9 ...
# MBE: k = 20 ... CON: k = 39"
cat("\n=== Arms per modality ===\n")
Tian %>% count(Agent) %>% arrange(desc(n)) %>% print()

# --- Intervention parameters ------------------------------------------------
# Reported in text: "duration ranged from 4 to 32 weeks (median = 12),
# frequency ranged from 1 to 7 (median = 3), session time ranged from
# 20 to 90 minutes (median = 60)"
cat("\n=== Intervention parameters (exercise arms only) ===\n")

Tian %>%
  filter(Agent != "CON") %>%
  summarise(
    duration_min    = min(Exact_Length, na.rm = TRUE),
    duration_max    = max(Exact_Length, na.rm = TRUE),
    duration_median = median(Exact_Length, na.rm = TRUE),
    freq_min        = min(Frequency, na.rm = TRUE),
    freq_max        = max(Frequency, na.rm = TRUE),
    freq_median     = median(Frequency, na.rm = TRUE),
    minutes_min     = min(Min_day_avg, na.rm = TRUE),
    minutes_max     = max(Min_day_avg, na.rm = TRUE),
    minutes_median  = median(Min_day_avg, na.rm = TRUE)
  ) %>%
  print()

# --- Evidence structure -----------------------------------------------------
# Reported in text: "Six trials were head-to-head comparisons between
# different exercise types (no control group), and six trials tested
# different doses of the same exercise intervention"
cat("\n=== Head-to-head trials (no control) ===\n")
h2h <- Tian %>%
  group_by(Study) %>%
  summarise(has_control = any(Agent == "CON"), .groups = "drop") %>%
  filter(!has_control)
cat("Studies without control:", nrow(h2h), "\n")

cat("\n=== Within-modality dose-comparison studies ===\n")
dose_comp <- Tian %>%
  filter(Agent != "CON") %>%
  group_by(Study, Agent) %>%
  summarise(n_doses = n_distinct(Dose), .groups = "drop") %>%
  filter(n_doses > 1)
cat("Studies comparing doses within same modality:", nrow(dose_comp), "\n")

# --- Dose distribution (Table 3.1) ------------------------------------------
# Reported in Table 3.1: dose distribution by modality
cat("\n=== Dose distribution by modality (Table 3.1) ===\n")
Tian %>%
  group_by(Agent) %>%
  summarise(
    n_arms         = n(),
    n_participants = sum(N, na.rm = TRUE),
    median_n       = median(N, na.rm = TRUE),
    median_dose    = median(Exact_dose, na.rm = TRUE),
    min_dose       = min(Exact_dose, na.rm = TRUE),
    max_dose       = max(Exact_dose, na.rm = TRUE),
    q1_dose        = quantile(Exact_dose, 0.25, na.rm = TRUE),
    q3_dose        = quantile(Exact_dose, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# Overall dose summary (exercise arms only)
cat("\n=== Overall dose distribution ===\n")
Tian %>%
  filter(Agent != "CON") %>%
  summarise(
    median = median(Exact_dose, na.rm = TRUE),
    min    = min(Exact_dose, na.rm = TRUE),
    max    = max(Exact_dose, na.rm = TRUE),
    q1     = quantile(Exact_dose, 0.25, na.rm = TRUE),
    q3     = quantile(Exact_dose, 0.75, na.rm = TRUE)
  ) %>%
  print()

# --- Arms by dose category (Table 3.2) --------------------------------------
# Reported in Table 3.2: arms per modality-dose node in splitting approach
cat("\n=== Arms by modality and dose category (Table 3.2) ===\n")
Tian %>%
  filter(Agent != "CON") %>%
  count(Agent, Dose) %>%
  pivot_wider(names_from = Dose, values_from = n, values_fill = 0) %>%
  print()

# --- Outcome scales ----------------------------------------------------------
# Reported in text: "Eight scales were used ... HAM-D (n = 16),
# BDI (n = 13), GDS (n = 6) ..."
cat("\n=== Depression scales ===\n")
cat("Unique scales:", n_distinct(Tian$Scale, na.rm = TRUE), "\n")

Tian %>%
  distinct(Study, Scale) %>%
  count(Scale, name = "n_studies") %>%
  arrange(desc(n_studies)) %>%
  print()

# Check for studies using multiple scales (should be 0)
multi_scale <- Tian %>%
  group_by(Study) %>%
  summarise(n_scales = n_distinct(Scale, na.rm = TRUE), .groups = "drop") %>%
  filter(n_scales > 1)
cat("Studies with multiple scales:", nrow(multi_scale), "\n")