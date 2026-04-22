################################################################################
# Chapter 3: Clean Data - Tian Network
#
# Tian dataset is loaded from the supplemental file provided by
# Tian et al. (2024) with minimal cleaning already applied.
# This script adds the sensitivity dose metric (total min/week).
################################################################################

library(dplyr)

Tian <- Tian %>%
  mutate(Dose_minutes = Min_day_avg * Frequency)