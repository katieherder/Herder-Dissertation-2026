################################################################################
# Chapter 4: Network characteristics — SSRI and GRISELDA datasets
#
# This script produces all the descriptive numbers reported in Sections
# "SSRI Network" and "GRISELDA Network" of the results chapter, plus the
# network diagrams.
#
# First run:
#   - 01_load_data
#   - 02_clean_data
################################################################################

library(MBNMAdose)
library(dplyr)

setwd("C:/Users/katie/Desktop/Depression-NMA/Chapter 4")


################################################################################
# 1. SSRI NETWORK
#
# Feeds the "SSRI Network" subsection: overall totals, Table of network
# characteristics by agent, arm-count breakdown, follow-up distribution,
# placebo response rate, and Figure 4.1 (dose-level network diagram).
################################################################################

# --- Overall totals -----------------------------------------------------------
# Reported in text: "60 RCTs ... 145 arms ... 15,174 participants"
# and "median sample size per arm was 99 (range: 10–332)".
cat("=== SSRI: overall ===\n")
cat("Studies:     ", n_distinct(ssri$studyID), "\n")
cat("Arms:        ", nrow(ssri), "\n")
cat("Participants:", sum(ssri$n), "\n")
cat("Median arm n:", median(ssri$n),
    " (range ", min(ssri$n), "–", max(ssri$n), ")\n", sep = "")

# --- Arms-per-study breakdown -------------------------------------------------
# Reported in text: "43 two-arm trials, 11 three-arm, 4 four-arm, 2 five-arm".
cat("\n=== SSRI: arms per study ===\n")
ssri %>%
  group_by(studyID) %>%
  summarise(k = n(), .groups = "drop") %>%
  pull(k) %>%
  table()

# --- Follow-up duration -------------------------------------------------------
# Reported in text: "most commonly 8 weeks (49%) or 6 weeks (42%)".
cat("\n=== SSRI: follow-up weeks ===\n")
ssri %>%
  distinct(studyID, weeks) %>%
  count(weeks) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# --- Table: network characteristics by agent ----------------------------------
# Produces Table "Network characteristics of the SSRI dataset by agent":
# studies, arms, participants, dose levels, dose range, response rate.
cat("\n=== SSRI: Table by agent ===\n")
ssri %>%
  group_by(agent) %>%
  summarise(
    n_studies      = n_distinct(studyID),
    n_arms         = n(),
    n_participants = sum(n),
    n_doses        = n_distinct(dose),
    dose_range     = if (all(dose == 0)) "--"
    else paste0(min(dose[dose > 0]), "–", max(dose[dose > 0])),
    response_rate  = round(sum(r) / sum(n), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(n_studies)) %>%
  print(width = Inf)

# --- Placebo pooled response rate ---------------------------------------------
# Reported in text: "pooled placebo response rate was 39.3% (2,183/5,556)".
cat("\n=== SSRI: placebo pooled rate ===\n")
ssri %>%
  filter(agent == "placebo") %>%
  summarise(k = n(), total_n = sum(n), total_r = sum(r),
            response_rate = round(total_r / total_n, 3))

# --- Head-to-head / within-agent dose-comparison checks -----------------------
# Supports claim: "no trial directly compared two different SSRIs;
# one trial compared two doses of citalopram without a placebo arm".
cat("\n=== SSRI: studies without placebo ===\n")
ssri %>%
  group_by(studyID) %>%
  filter(all(agent != "placebo")) %>%
  print(n = Inf)

cat("\n=== SSRI: studies comparing different SSRIs ===\n")
ssri %>%
  filter(agent != "placebo") %>%
  group_by(studyID) %>%
  filter(n_distinct(agent) > 1) %>%
  select(studyID, agent, dose) %>%
  print(n = Inf)

# --- Figure 4.1: dose-level network diagram -----------------------------------
ssri_network <- mbnma.network(data.ab = ssri)

jpeg("figures/Figure 4.1.jpeg", width = 1800, height = 1500, res = 300)
par(mar = c(2, 2, 2, 2))
plot(ssri_network, main = "")
dev.off()


################################################################################
# 2. GRISELDA NETWORK
#
# Feeds the "GRISELDA Network" subsection: overall totals, class-level table,
# agent-level table with raw and FE dose ranges, head-to-head structure,
# sparse-agent flagging, and the network diagram.
################################################################################

# --- Overall totals -----------------------------------------------------------
# Reported in text: "133 fixed-dose RCTs ... 330 arms ... 40,956 participants"
# and "median sample size per arm was 126 (range: 12–303)".
cat("\n=== GRISELDA: overall ===\n")
cat("Studies:     ", n_distinct(Griselda_clean$studyID), "\n")
cat("Arms:        ", nrow(Griselda_clean), "\n")
cat("Participants:", sum(Griselda_clean$n), "\n")
cat("Median arm n:", median(Griselda_clean$n),
    " (range ", min(Griselda_clean$n), "–", max(Griselda_clean$n), ")\n", sep = "")

# --- Arms-per-study breakdown -------------------------------------------------
# Reported in text: "89 two-arm, 24 three-arm, 20 four-arm".
cat("\n=== GRISELDA: arms per study ===\n")
Griselda_clean %>%
  group_by(studyID) %>%
  summarise(k = n(), .groups = "drop") %>%
  pull(k) %>%
  table()

# --- Follow-up duration -------------------------------------------------------
# Reported in text: "most commonly 8 weeks (54%) or 6 weeks (32%)".
cat("\n=== GRISELDA: follow-up weeks ===\n")
Griselda_clean %>%
  distinct(studyID, weeks) %>%
  count(weeks) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# --- Placebo pooled response rate ---------------------------------------------
# Reported in text: "pooled placebo response rate was 37.3% (4,159/11,140)".
cat("\n=== GRISELDA: placebo pooled rate ===\n")
Griselda_clean %>%
  filter(agent == "placebo") %>%
  summarise(k = n(), total_n = sum(n), total_r = sum(r),
            response_rate = round(total_r / total_n, 3))

# --- Evidence structure: placebo-free and cross-class studies -----------------
# Reported in text: "45 (34%) included no placebo arm" and
# "64 trials (48%) compared two or more different active agents".
cat("\n=== GRISELDA: evidence structure ===\n")
Griselda_clean %>%
  group_by(studyID) %>%
  summarise(
    has_placebo   = any(agent == "placebo"),
    n_active      = n_distinct(agent[agent != "placebo"]),
    n_classes_act = n_distinct(Class[agent != "placebo"]),
    .groups = "drop"
  ) %>%
  summarise(
    total_studies       = n(),
    no_placebo          = sum(!has_placebo),
    pct_no_placebo      = round(100 * mean(!has_placebo), 1),
    multi_active        = sum(n_active > 1),
    pct_multi_active    = round(100 * mean(n_active > 1), 1),
    cross_class         = sum(n_classes_act > 1),
    pct_cross_class     = round(100 * mean(n_classes_act > 1), 1)
  ) %>%
  print()

# --- Cross-class comparison counts --------------------------------------------
# Supports the claim "SNRI vs SSRI (18 studies) and SSRI vs TCA (10 studies)
# the most common".
cat("\n=== GRISELDA: cross-class comparisons ===\n")
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(studyID) %>%
  filter(n_distinct(Class) > 1) %>%
  summarise(classes = paste(sort(unique(Class)), collapse = " vs "),
            .groups = "drop") %>%
  count(classes, sort = TRUE) %>%
  print()

# --- Table: class-level network summary ---------------------------------------
# Produces Table "Network characteristics of the GRISELDA dataset by
# pharmacologic class": agents, studies, arms, N, dose levels, response rate.
# Note: placebo row is reported separately in the table; kept here for totals.
cat("\n=== GRISELDA: Table by class ===\n")
Griselda_clean %>%
  group_by(Class) %>%
  summarise(
    n_agents       = n_distinct(agent[agent != "placebo"]),
    n_studies      = n_distinct(studyID),
    n_arms         = n(),
    n_participants = sum(n),
    n_dose_levels  = n_distinct(paste(agent, dose)[agent != "placebo"]),
    response_rate  = round(sum(r) / sum(n), 3),
    .groups = "drop"
  ) %>%
  print(width = Inf)

# --- Table: agent-level detail with raw and FE dose ranges --------------------
# Produces Table "Network characteristics of the GRISELDA dataset by agent":
# studies, arms, N, dose levels, raw range, FE range, response rate.
# Uses dose_FE already computed in Griselda_clean — no re-conversion needed.
cat("\n=== GRISELDA: Table by agent ===\n")
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(Class, agent) %>%
  summarise(
    n_studies      = n_distinct(studyID),
    n_arms         = n(),
    n_participants = sum(n),
    n_doses        = n_distinct(dose),
    raw_range      = paste0(min(dose), "–", max(dose)),
    fe_range       = paste0(round(min(dose_FE), 1), "–", round(max(dose_FE), 1)),
    response_rate  = round(sum(r) / sum(n), 3),
    .groups = "drop"
  ) %>%
  arrange(Class, agent) %>%
  print(n = Inf, width = Inf)

# --- Sparse agents ------------------------------------------------------------
# Supports claim: "reboxetine, trazodone, nefazodone, and bupropion each
# had ≤3 arms and 1–2 dose levels".
cat("\n=== GRISELDA: sparse agents (≤3 arms or ≤1 dose) ===\n")
Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(Class, agent) %>%
  summarise(n_arms    = n(),
            n_studies = n_distinct(studyID),
            n_doses   = n_distinct(dose),
            .groups = "drop") %>%
  filter(n_arms <= 3 | n_doses <= 1) %>%
  arrange(n_arms) %>%
  print()

# --- GRISELDA network diagram -------------------------------------------------
griselda_network <- mbnma.network(data.ab = Griselda_clean)

jpeg("figures/Figure_Griselda_network.jpeg",
     width = 1800, height = 1500, res = 300)
par(mar = c(2, 2, 2, 2))
plot(griselda_network, main = "Griselda Network")
dev.off()


################################################################################
# 3. OVERLAP BETWEEN SSRI AND GRISELDA DATASETS
#
# Feeds the closing paragraph of the Extended GRISELDA Network subsection:
# "40 exact matches, 15 partial matches, 5 with no matching arms ...
#  32 studies retained after cleaning".
#
# Matching key: agent + dose + n + r. A study is an EXACT match if all its
# arms match arms in GRISELDA; PARTIAL if some but not all arms match;
# NONE if no arms match.
################################################################################

# Match against the full fixed-dose GRISELDA (pre-cleaning), so that studies
# dropped by the cleaning filters still get a chance to match.
griselda_keys <- paste(Griselda$Drug, Griselda$Dose,
                       Griselda$No_randomised, Griselda$Responders, sep = "_")

ssri_overlap <- ssri %>%
  mutate(key     = paste(agent, dose, n, r, sep = "_"),
         matched = key %in% griselda_keys) %>%
  group_by(studyID) %>%
  summarise(
    n_arms   = n(),
    n_match  = sum(matched),
    status   = case_when(
      n_match == n_arms ~ "exact",
      n_match > 0       ~ "partial",
      TRUE              ~ "none"
    ),
    .groups = "drop"
  )

cat("\n=== SSRI↔GRISELDA study-level overlap ===\n")
ssri_overlap %>% count(status) %>% print()

# Also: how many of the exact-match studies survive the cleaning filters
# that produced Griselda_clean. This is the "32 studies retained" number.
exact_ids <- ssri_overlap %>% filter(status == "exact") %>% pull(studyID)
# (Matching SSRI numeric studyIDs to Griselda_clean string studyIDs requires
# a crosswalk; if none exists, report based on the fixed-dose match above
# and note that the reduction from 40 → 32 reflects downstream exclusions.)
cat("Exact-match studies (pre-cleaning):", length(exact_ids), "\n")

# Build arm keys for Griselda_clean (post-cleaning)
griselda_clean_keys <- paste(Griselda_clean$agent, Griselda_clean$dose,
                             Griselda_clean$n, Griselda_clean$r, sep = "_")

# For each SSRI study, check how many of its arms survive into Griselda_clean
ssri_retained <- ssri %>%
  mutate(key = paste(agent, dose, n, r, sep = "_"),
         matched_pre  = key %in% griselda_keys,
         matched_post = key %in% griselda_clean_keys) %>%
  group_by(studyID) %>%
  summarise(
    n_arms       = n(),
    n_match_pre  = sum(matched_pre),
    n_match_post = sum(matched_post),
    status_pre = case_when(
      n_match_pre == n_arms ~ "exact",
      n_match_pre > 0       ~ "partial",
      TRUE                  ~ "none"
    ),
    status_post = case_when(
      n_match_post == n_arms ~ "exact",
      n_match_post > 0       ~ "partial",
      TRUE                   ~ "none"
    ),
    .groups = "drop"
  )

cat("\n=== SSRI studies retained through GRISELDA cleaning ===\n")
cat("Exact match, pre-cleaning: ",
    sum(ssri_retained$status_pre  == "exact"), "\n")
cat("Exact match, post-cleaning:",
    sum(ssri_retained$status_post == "exact"), "\n")

ssri_retained %>% count(status_pre, status_post) %>% print()