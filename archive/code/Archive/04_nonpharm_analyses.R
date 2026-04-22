# Dose-Response Methods Analysis for Exercise Interventions
# Implementation of methods described in Chapter 3 of the dissertation
# Author: Katie Herder
# Date: February 17, 2026

# Load required packages
library(tidyverse)
library(MBNMAdose)
library(gemtc)
library(coda)
library(ggplot2)
library(corrplot)
library(knitr)
library(kableExtra)
library(grid)
library(gridExtra)
library(magick)

# Set seed for reproducibility
set.seed(12345)

# Main project directory
setwd("C:/Users/katie/Desktop/Depression-NMA")

head(Tian)
#==============================================================================
# 3.3.1 NETWORK CHARACTERISTICS AND CONNECTIVITY
#==============================================================================

# Basic network statistics
network_stats <- Tian %>%
  summarise(
    n_studies = length(unique(Study)),
    n_arms = nrow(.),
    n_participants = sum(N),
    n_exercise_arms = sum(Agent != "CON"),
    n_control_arms = sum(Agent == "CON"),
  )

print("Network Statistics:")
print(network_stats)

dose_distribution_no_con <- Tian %>%
  filter(Agent != "CON") %>%
  summarise(
    n_arms = n(),
    n_participants = sum(N, na.rm = TRUE),
    n_median_participants = median(N, na.rm = TRUE),
    median_dose_met = median(Exact_dose, na.rm = TRUE),
    min_dose_met = min(Exact_dose, na.rm = TRUE),
    max_dose_met = max(Exact_dose, na.rm = TRUE),
    q1_dose_met = quantile(Exact_dose, 0.25, na.rm = TRUE),
    q3_dose_met = quantile(Exact_dose, 0.75, na.rm = TRUE)
  )

print("Overall Dose Distribution:")
print(dose_distribution_no_con)

# Create dose distribution summary by agent (for Table 3.1)
dose_distribution <- Tian %>%
  group_by(Agent) %>%
  summarise(
    n_arms = n(),
    n_participants = sum(N),
    n_median_participants = median(N),
    median_dose_met = median(Exact_dose),
    min_dose_met = min(Exact_dose),
    max_dose_met = max(Exact_dose),
    q1_dose_met = quantile(Exact_dose, 0.25, na.rm = TRUE),
    q3_dose_met = quantile(Exact_dose, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

print(dose_distribution)

# Check dose range and disagreement for each agent-dose combination
Tian %>%
  filter(Agent != "CON") %>%
  group_by(Agent, Dose) %>%
  summarise(
    k = n(),
    min_exact = min(Exact_dose),
    max_exact = max(Exact_dose),
    range_exact = max(Exact_dose) - min(Exact_dose),
    .groups = "drop"
  ) %>%
  arrange(Agent, Dose)

#==============================================================================
# DATA PREPARATION FOR DIFFERENT METHODS
#==============================================================================

# For MBNMAdose (splitting)
split_data <- Tian %>%
  transmute(
    studyID = Study,
    agent   = Agent,
    dose    = Dose,      # categorized MET-min/week
    y       = y,
    se      = SE,
    n       = N
  )

# For MBNMAdose (MBNMA methods)
mbnma_data <- Tian %>%
  transmute(
    studyID = Study,
    agent   = Agent,
    dose    = Exact_dose / 1000,      # primary metric, raw MET-min/week, rescaled to thousands of MET-min/week
    y       = y,
    se      = SE,
    n       = N
  )

# Alternative dose metric for sensitivity analysis
mbnma_data_min <- mbnma_data %>%
  mutate(dose = Tian$Dose_minutes)

# For gemtc (lumping) - aggregate duplicate treatment arms within studies
gemtc_data <- Tian %>%
  group_by(Study, Agent) %>%
  summarise(
    mean = weighted.mean(y, w = N),          # weighted mean of change scores
    std.err = sqrt(1 / sum(1 / SE^2)),       # pooled SE (inverse-variance)
    sampleSize = sum(N),                      # total sample size
    .groups = "drop"
  ) %>%
  rename(study = Study, treatment = Agent)

# Verify no duplicates remain
any(duplicated(gemtc_data[, c("study", "treatment")]))  # should be FALSE

print(paste("MBNMAdose rows:", nrow(mbnma_data)))
print(paste("gemtc rows:", nrow(gemtc_data)))


lumped_data <- Tian %>%
  filter(!Study %in% c("Hanssen et al, 2018", 
                        "Krogh et al,2012", 
                        "Scott et al, 2019", 
                        "Streeter et al,2017", 
                        "Streeter, 2020")) %>%
  transmute(
    studyID = Study,
    agent = Agent,
    dose = ifelse(Agent == "CON", 0, 1),
    y = y,
    se = SE,
    n = N
  ) %>%
  group_by(studyID, agent, dose) %>%
  summarise(
    y = weighted.mean(y, w = n),
    se = sqrt(1 / sum(1 / se^2)),
    n = sum(n),
    .groups = "drop"
  )

network_lumped <- mbnma.network(data.ab = lumped_data)
result_lumped <- nma.run(network_lumped, method = "random", link = "smd")
#==============================================================================
# 3.3.1 NETWORK CHARACTERISTICS AND CONNECTIVITY
#==============================================================================

print("TABLE 3.2 - Network connectivity and dose distribution:")

# Network connectivity matrix
connectivity_matrix <- Tian %>%
  select(Study, Agent, Exact_dose) %>%
  mutate(dose_cat = case_when(
    Agent == "CON" ~ "Control",
    Exact_dose <= 400 ~ "Low (≤400)",
    Exact_dose <= 800 ~ "Medium (401-800)",
    TRUE ~ "High (>800)"
  )) %>%
  group_by(Study) %>%
  summarise(connections = paste(paste(Agent, dose_cat, sep="-"), collapse = " vs "), .groups = 'drop')

print("Network connections by study:")
print(connectivity_matrix, n=45)

dose_counts <- Tian %>%
  group_by(Agent, Dose) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(Agent, Dose)

print(dose_counts)

#------------------------------------------------------------------------------
# NETWORK OBJECT AND PLOTS
#------------------------------------------------------------------------------

# --- Save Panel A alone ---
jpeg("figures/tmp_panelA.jpeg", width = 1800, height = 1500, res = 300)
par(mar = c(2, 2, 2, 2))  # fixed: was missing 4th value
network_gemtc <- mtc.network(
  data.ab = gemtc_data  
)
plot(network_gemtc, main = "")
mtext("A", side = 3, line = 0.5, 
      at = par("usr")[1] + 0.02*(par("usr")[2]-par("usr")[1]),
      cex = 2, font = 2)
dev.off()

# --- Save Panel B alone ---
jpeg("figures/tmp_panelB.jpeg", width = 1600, height = 1500, res = 300)
par(mar = c(2, 2, 2, 2))
mbnma_network <- mbnma.network(data.ab = split_data)
plot(mbnma_network, main = "")
mtext("B", side = 3, line = 0.5, 
      at = par("usr")[1] + 0.02*(par("usr")[2]-par("usr")[1]),
      cex = 2, font = 2)
dev.off()

# --- Combine side by side with magick ---
img_a <- image_read("figures/tmp_panelA.jpeg")
img_b <- image_read("figures/tmp_panelB.jpeg")

combined <- image_append(c(img_a, img_b), stack = FALSE)  # FALSE = side by side
image_write(combined, "figures/Figure 3.1.jpeg")

# Clean up temp files
file.remove("figures/tmp_panelA.jpeg", "figures/tmp_panelB.jpeg")

#==============================================================================
# 3.3.2 PRIMARY ANALYSIS: MET-MINUTES/WEEK
#==============================================================================

# Assumes data preparation already complete:
# - split_data: studyID, agent, dose (categorized), y, se, n
# - mbnma_data: studyID, agent, dose (exact MET-min/week), y, se, n
# - gemtc_data: study, treatment, mean, std.err, sampleSize

set.seed(12345)

# Storage for results
results_all <- list()
convergence_all <- data.frame()
#==============================================================================
# METHOD 1: LUMPING (gemtc)
# Standard NMA ignoring dose - each modality is one node
#==============================================================================

cat("\n========== METHOD 1: LUMPING ==========\n")

# Create gemtc network
network_lumping <- mtc.network(data.ab = gemtc_data)

# Fit random-effects NMA
model_lumping <- mtc.model(network_lumping,
                           linearModel = "random",
                           n.chain = 4)

results_lumping <- mtc.run(model_lumping,
                           n.adapt = 5000,
                           n.iter = 20000,
                           thin = 10)

# --- Convergence diagnostics (on original model) ---
gd_lumping <- gelman.diag(results_lumping, multivariate = FALSE)
cat("\nGelman-Rubin diagnostics:\n")
print(gd_lumping)

ess_lumping <- effectiveSize(as.mcmc.list(results_lumping))
cat("\nEffective sample sizes:\n")
print(ess_lumping)

# --- Effects relative to CON ---
re_vs_con <- relative.effect(results_lumping, t1 = "CON")
cat("\nEffects vs CON:\n")
print(summary(re_vs_con))

# --- DIC ---
cat("\nDIC:\n")
print(summary(results_lumping)$DIC)

# Store convergence info
convergence_all <- rbind(convergence_all, data.frame(
  Method = "Lumping",
  MaxRhat = max(gd_lumping$psrf[, "Point est."], na.rm = TRUE),
  MinESS = min(ess_lumping, na.rm = TRUE),
  Converged = max(gd_lumping$psrf[, "Point est."], na.rm = TRUE) < 1.05
))

# Store results (keep original model + CON-referenced effects)
results_all$lumping <- list(
  model = results_lumping,
  re_vs_con = re_vs_con,
  network = network_lumping
)

# --- Plots ---
jpeg("figures/trace_lumping.jpeg", width = 2400, height = 1800, res = 300)
plot(results_lumping)
dev.off()

jpeg("figures/forest_lumping.jpeg", width = 2400, height = 1200, res = 300)
forest(relative.effect.table(results_lumping), t1 = "CON")
dev.off()

cat("\n--- Lumping complete ---\n")

#==============================================================================
# METHOD 2: SPLITTING (MBNMAdose nma.run)
# Each agent-dose combination is a separate node
# Uses the categorized dose (0, 300, 600, 900, 1200)
#==============================================================================

cat("\n========== METHOD 2: SPLITTING ==========\n")

# Create network using categorized doses
# split_data has: studyID, agent, dose (categorized), y, se, n
network_split <- mbnma.network(data.ab = split_data)

# nma.run fits a standard NMA treating each agent-dose as a separate treatment
# This is the "splitting" approach - no dose-response assumption
result_split <- nma.run(network_split, 
                        method = "random",
                        n.iter = 20000,
                        link = "smd")

cat("\nSplitting NMA summary:\n")
print(result_split)

result_split$jagsresult$BUGSoutput$DIC
result_split$jagsresult$BUGSoutput$summary["sd", "50%"]
result_split$jagsresult$BUGSoutput$summary["totresdev", "50%"]

# Store results
results_all$splitting <- list(
  model = result_split,
  network = network_split
)

split_bugs <- results_all$splitting$model$jagsresult$BUGSoutput
split_ess <- split_bugs$summary[, "n.eff"]
min(split_ess[split_ess > 1])

cat("\n--- Splitting complete ---\n")


#==============================================================================
# METHOD 3: LINEAR MBNMA
# f(x) = beta_1 * x  (degree=1 polynomial)
#==============================================================================

cat("\n========== METHOD 3: LINEAR MBNMA ==========\n")

# Create network using exact (continuous) doses
network_mbnma <- mbnma.network(data.ab = mbnma_data)

# Linear: dpoly(degree=1) with relative effects, random between-study effects
result_linear <- mbnma.run(network_mbnma,
                           fun = dpoly(degree = 1, beta.1 = "rel"),
                           method = "random",
                           link = "smd",
                           n.iter = 20000,
                           n.chains = 3,
                           pD = TRUE)

cat("\nLinear MBNMA summary:\n")
summary(result_linear)

# Store
results_all$linear <- list(model = result_linear, network = network_mbnma)

cat("\n--- Linear MBNMA complete ---\n")


#==============================================================================
# METHOD 4: EMAX MBNMA
# f(x) = Emax * x / (ED50 + x)
#==============================================================================

cat("\n========== METHOD 4: EMAX MBNMA ==========\n")

# Emax with relative effects on both parameters
# Note: ED50 is modeled on log-scale internally
result_emax <- mbnma.run(network_mbnma,
                         fun = demax(emax = "rel", ed50 = "rel"),
                         method = "random",
                         link = "smd",
                         n.iter = 20000,
                         n.chains = 3,
                         pD = TRUE)

cat("\nEmax MBNMA summary:\n")
summary(result_emax)

results_all$emax <- list(model = result_emax, network = network_mbnma)

cat("\n--- Emax MBNMA complete ---\n")


#==============================================================================
# METHOD 5: QUADRATIC MBNMA
#==============================================================================

cat("\n========== METHOD 5: QUADRATIC MBNMA ==========\n")

# Quadratic using duser() with rescaled doses (already in thousands via mbnma_data)
quad_fun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))

result_quad <- mbnma.run(network_mbnma,
                         fun = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
                         method = "random",
                         link = "smd",
                         n.iter = 20000,
                         n.chains = 3,
                         pD = TRUE)

# Check ESS
cat("\nQuadratic parameter ESS:\n")
print(result_quad$BUGSoutput$summary[, c("mean", "Rhat", "n.eff")])

cat("\nQuadratic MBNMA summary:\n")
summary(result_quad)

results_all$quadratic <- list(model = result_quad, network = network_mbnma)

cat("\n--- Quadratic MBNMA complete ---\n")


#==============================================================================
# METHOD 6: NATURAL SPLINE (NS) MBNMA
# f(x) = sum of spline basis functions with knots at 10th, 50th, 90th pctiles
#==============================================================================

cat("\n========== METHOD 6: NS MBNMA ==========\n")

# Calculate knot positions from observed doses (excluding controls)
exercise_doses <- mbnma_data$dose[mbnma_data$agent != "CON" & mbnma_data$dose > 0]
knot_positions <- quantile(exercise_doses, probs = c(0.10, 0.50, 0.90))
cat("RCS knot positions (MET-min/week):", knot_positions, "\n")

# NS spline with 3 knots
result_rcs <- mbnma.run(network_mbnma,
                        fun = dspline(type = "ns", 
                                      knots = c(0.10, 0.50, 0.90),
                                      beta.1 = "rel", 
                                      beta.2 = "rel"),
                        method = "random",
                        link = "smd",
                        n.iter = 20000,
                        n.chains = 3,
                        pD = TRUE)

cat("\nRCS MBNMA summary:\n")
summary(result_rcs)

results_all$rcs <- list(model = result_rcs, network = network_mbnma)

cat("\n--- RCS MBNMA complete ---\n")


#==============================================================================
# CONVERGENCE SUMMARY TABLE
#==============================================================================
cat("\n\n========== CONVERGENCE SUMMARY ==========\n")

# Start fresh
convergence_all <- data.frame()

# --- Lumping (already extracted above, re-add) ---
convergence_all <- rbind(convergence_all, data.frame(
  Method = "Lumping",
  MaxRhat = max(gd_lumping$psrf[, "Point est."], na.rm = TRUE),
  MinESS = min(ess_lumping, na.rm = TRUE),
  Converged = max(gd_lumping$psrf[, "Point est."], na.rm = TRUE) < 1.05
))

# --- Splitting (jagsresult$BUGSoutput) ---
split_bugs <- results_all$splitting$model$jagsresult$BUGSoutput
split_rhat <- split_bugs$summary[, "Rhat"]
split_ess <- split_bugs$summary[, "n.eff"]
convergence_all <- rbind(convergence_all, data.frame(
  Method = "Splitting",
  MaxRhat = max(split_rhat, na.rm = TRUE),
  MinESS = min(split_ess[split_ess > 0], na.rm = TRUE),
  Converged = max(split_rhat, na.rm = TRUE) < 1.05
))

# --- Linear, Emax, Quadratic, NCS (standard BUGSoutput) ---
for (nm in c("linear", "emax", "quadratic", "rcs")) {
  label <- c(linear = "Linear MBNMA", emax = "Emax MBNMA",
             quadratic = "Quadratic MBNMA", rcs = "NCS MBNMA")[nm]
  bugs <- results_all[[nm]]$model$BUGSoutput
  rhat_vals <- bugs$summary[, "Rhat"]
  ess_vals <- bugs$summary[, "n.eff"]
  convergence_all <- rbind(convergence_all, data.frame(
    Method = label,
    MaxRhat = max(rhat_vals, na.rm = TRUE),
    MinESS = min(ess_vals[ess_vals > 0], na.rm = TRUE),
    Converged = max(rhat_vals, na.rm = TRUE) < 1.05
  ))
}

print(convergence_all)

#==============================================================================
# MODEL COMPARISON TABLE (Table 3.3)
#==============================================================================
cat("\n\n========== MODEL COMPARISON (Table 3.3) ==========\n")

# Function to extract DIC components from MBNMAdose models
extract_fit_mbnma <- function(model) {
  list(
    DIC = model$BUGSoutput$DIC,
    pD = model$BUGSoutput$pD,
    resdev = model$BUGSoutput$summary["totresdev", "mean"],
    tau = tryCatch(
      model$BUGSoutput$summary["sd", "50%"],
      error = function(e) NA
    )
  )
}

# Build comparison table
model_names <- c("Lumping", "Splitting", "Linear MBNMA", 
                 "Emax MBNMA", "Quadratic MBNMA", "NCS MBNMA")
model_comparison <- data.frame(Method = model_names)

# --- Lumping (gemtc stores DIC differently) ---
lumping_summary <- summary(results_all$lumping$model)
model_comparison$DIC[1] <- lumping_summary$DIC[["DIC"]]
model_comparison$pD[1] <- lumping_summary$DIC[["pD"]]
model_comparison$ResidDev[1] <- lumping_summary$DIC[["Dbar"]]
re_summary <- summary(results_all$lumping$model)$summaries$quantiles
tau_row <- grep("sd.d", rownames(re_summary))
model_comparison$tau[1] <- ifelse(length(tau_row) > 0, 
                                  re_summary[tau_row, "50%"], NA)

# --- Splitting (different structure: jagsresult$BUGSoutput) ---
split_bugs <- results_all$splitting$model$jagsresult$BUGSoutput
model_comparison$DIC[2] <- split_bugs$DIC
model_comparison$pD[2] <- split_bugs$pV
model_comparison$ResidDev[2] <- split_bugs$summary["totresdev", "mean"]
model_comparison$tau[2] <- split_bugs$summary["sd", "50%"]

# --- Linear ---
fit <- extract_fit_mbnma(results_all$linear$model)
model_comparison[3, c("DIC","pD","ResidDev","tau")] <- c(fit$DIC, fit$pD, fit$resdev, fit$tau)

# --- Emax ---
fit <- extract_fit_mbnma(results_all$emax$model)
model_comparison[4, c("DIC","pD","ResidDev","tau")] <- c(fit$DIC, fit$pD, fit$resdev, fit$tau)

# --- Quadratic ---
quad_bugs <- results_all$quadratic$model$BUGSoutput
model_comparison$DIC[5] <- quad_bugs$DIC
model_comparison$pD[5] <- quad_bugs$pD
model_comparison$ResidDev[5] <- quad_bugs$summary["totresdev", "mean"]
model_comparison$tau[5] <- quad_bugs$summary["sd", "50%"]

# --- NCS ---
fit <- extract_fit_mbnma(results_all$rcs$model)
model_comparison[6, c("DIC","pD","ResidDev","tau")] <- c(fit$DIC, fit$pD, fit$resdev, fit$tau)

# Format and print
model_comparison$DIC <- round(model_comparison$DIC, 1)
model_comparison$pD <- round(model_comparison$pD, 1)
model_comparison$ResidDev <- round(model_comparison$ResidDev, 1)
model_comparison$tau <- round(model_comparison$tau, 3)

# Add delta DIC (relative to best among MBNMAdose models only)
mbnma_rows <- 2:6
model_comparison$DeltaDIC <- NA
model_comparison$DeltaDIC[mbnma_rows] <- round(
  model_comparison$DIC[mbnma_rows] - min(model_comparison$DIC[mbnma_rows], na.rm = TRUE), 1)

cat("\nTable 3.3: Model Performance Summary\n")
print(model_comparison)

# Save results
save(results_all, convergence_all, model_comparison, 
     file = "results/primary_analysis_results.RData")
names(results_all)
print(convergence_all)
print(model_comparison)

#==============================================================================
# EFFECT ESTIMATE COMPARISON - SPLITTING vs EMAX
# Continues after primary analysis (all models in results_all)
#
# This script:
#   1. Extracts treatment effects from splitting vs Placebo_0
#   2. Generates dose-response predictions from Emax model
#   3. Compares SMDs at standard doses (300, 600, 900, 1200 MET-min/week)
#   4. Identifies MED (minimum effective dose, SMD ≤ -0.5) from Emax
#   5. Compares agent rankings: splitting vs Emax
#   6. Generates dose-response curve plot (Figure 3.2)
#==============================================================================

# Standard doses for comparison
standard_doses_raw <- c(300, 600, 900, 1200)
agents <- c("AE", "MBE", "ME", "RE")

#==============================================================================
# 1. SPLITTING: Extract treatment effects vs Placebo_0
#    These are on the SMD scale (link = "smd")
#==============================================================================

cat("\n========== SPLITTING EFFECTS ==========\n")

split_labs <- results_all$splitting$model$trt.labs
split_summary <- results_all$splitting$model$jagsresult$BUGSoutput$summary

cat("Splitting treatment labels:\n")
print(split_labs)

# Build data frame: for each agent-dose node, extract SMD vs Placebo
splitting_effects <- data.frame()

for (j in 2:length(split_labs)) {
  lab <- split_labs[j]
  parts <- strsplit(lab, "_")[[1]]
  agent <- parts[1]
  dose_cat <- as.numeric(parts[2])
  param_name <- paste0("d[", j, "]")
  
  splitting_effects <- rbind(splitting_effects, data.frame(
    Agent = agent,
    Model = "Splitting",
    Dose = dose_cat,
    Median = split_summary[param_name, "50%"],
    Lower = split_summary[param_name, "2.5%"],
    Upper = split_summary[param_name, "97.5%"]
  ))
}

cat("\nSplitting effects (SMD vs Placebo):\n")
print(splitting_effects)

#==============================================================================
# 2. EMAX: Generate dose-response predictions
#    predict(model, E0 = 0) gives relative effects vs placebo (= SMD)
#    Doses must be on the SCALED metric (thousands of MET-min/week)
#==============================================================================

cat("\n========== EMAX PREDICTIONS ==========\n")

# First, test predict() to see its output structure
cat("\nTesting predict() structure...\n")
test_pred <- predict(results_all$emax$model, E0 = 0, n.dose = 5)
cat("Class:", class(test_pred), "\n")
cat("Names:", names(test_pred), "\n")
str(test_pred, max.level = 2)
cat("\nFirst element of predicts:\n")
print(head(test_pred$predicts[[1]]))

#==============================================================================
# NOTE: Run everything above first, then check the predict() output.
# The column names printed above will tell you what to use below.
# You may need to adjust "50%", "2.5%", "97.5%" to match actual names.
# Once confirmed, continue with the rest of the script.
#==============================================================================

# Fine grid for dose-response curves (scaled: 0 to 1.2 = 0 to 1200 MET-min/wk)
dose_grid_scaled <- seq(0, 1.2, by = 0.01)

# exact.doses must be a named list: one entry per agent
dose_list <- list(
  AE = dose_grid_scaled,
  MBE = dose_grid_scaled,
  ME = dose_grid_scaled,
  RE = dose_grid_scaled
)

pred_emax <- predict(results_all$emax$model, E0 = 0,
                     exact.doses = dose_list)

# Extract predictions into a data frame
emax_preds <- data.frame()

for (ag_name in names(pred_emax$predicts)) {
  if (grepl("Placebo", ag_name)) next
  
  ag_list <- pred_emax$predicts[[ag_name]]
  
  for (dose_name in names(ag_list)) {
    samples <- ag_list[[dose_name]]
    dose_val <- as.numeric(dose_name)
    
    emax_preds <- rbind(emax_preds, data.frame(
      Agent = ag_name,
      Model = "Emax",
      Dose_scaled = dose_val,
      Dose_raw = dose_val * 1000,
      Median = median(samples),
      Lower = quantile(samples, 0.025),
      Upper = quantile(samples, 0.975)
    ))
  }
}

row.names(emax_preds) <- NULL

cat("Emax predictions extracted. Rows:", nrow(emax_preds), "\n")
cat("\nSample (AE at standard doses):\n")
print(emax_preds %>% filter(Agent == "AE", Dose_raw %in% c(300, 600, 900, 1200)))

#==============================================================================
# 3. COMPARE SMDs AT STANDARD DOSES: SPLITTING vs EMAX
#==============================================================================

cat("\n========== COMPARISON AT STANDARD DOSES ==========\n")

# Emax predictions at standard doses
emax_at_standard <- emax_preds %>%
  filter(abs(Dose_raw - 300) < 10 |
           abs(Dose_raw - 600) < 10 |
           abs(Dose_raw - 900) < 10 |
           abs(Dose_raw - 1200) < 10) %>%
  mutate(Dose = round(Dose_raw / 300) * 300) %>%
  select(Agent, Model, Dose, Median, Lower, Upper)

# Splitting at standard doses
splitting_at_standard <- splitting_effects %>%
  filter(Dose %in% standard_doses_raw)

# Combine
comparison <- bind_rows(splitting_at_standard, emax_at_standard) %>%
  arrange(Agent, Dose, Model)

# Print side-by-side
cat("\nSplitting vs Emax at standard doses:\n")
comparison_wide <- comparison %>%
  mutate(Effect = sprintf("%.2f (%.2f, %.2f)", Median, Lower, Upper)) %>%
  select(Agent, Dose, Model, Effect) %>%
  pivot_wider(names_from = Model, values_from = Effect)
print(comparison_wide, n = 50)

# Calculate disagreement (difference in medians)
disagreement <- comparison %>%
  select(Agent, Dose, Model, Median) %>%
  pivot_wider(names_from = Model, values_from = Median) %>%
  mutate(
    Diff = Emax - Splitting,
    Disagreement = abs(Diff) > 0.3
  )

cat("\nDisagreement (|Emax - Splitting| > 0.3 SMD flagged):\n")
print(disagreement)

#==============================================================================
# 4. MINIMUM EFFECTIVE DOSE (MED) FROM EMAX
#    Defined as lowest dose where median SMD ≤ -0.5
#==============================================================================

cat("\n========== MINIMUM EFFECTIVE DOSE ==========\n")

med_emax <- emax_preds %>%
  filter(Median <= -0.5) %>%
  group_by(Agent) %>%
  summarise(
    MED = min(Dose_raw),
    SMD_at_MED = Median[which.min(Dose_raw)],
    .groups = "drop"
  )

cat("\nMinimum Effective Dose (Emax model, SMD ≤ -0.5):\n")
print(med_emax)

# Agents that never reach -0.5
agents_no_med <- setdiff(agents, med_emax$Agent)
if (length(agents_no_med) > 0) {
  cat("\nAgents that do not reach SMD ≤ -0.5 within observed range:",
      paste(agents_no_med, collapse = ", "), "\n")
}

#==============================================================================
# 5. AGENT RANKINGS: SPLITTING vs EMAX
#    Rank agents by median SMD (most negative = rank 1 = best)
#==============================================================================

cat("\n========== AGENT RANKINGS ==========\n")

rankings <- comparison %>%
  group_by(Model, Dose) %>%
  mutate(Rank = base::rank(Median)) %>%
  ungroup()

rankings_wide <- rankings %>%
  select(Agent, Dose, Model, Rank) %>%
  pivot_wider(names_from = Model, values_from = Rank) %>%
  mutate(Rank_Shift = abs(Emax - Splitting))

cat("\nAgent rankings (Splitting vs Emax):\n")
print(rankings_wide, n = 50)

cat("\nRanking shifts > 0:\n")
print(rankings_wide %>% filter(Rank_Shift > 0))

#==============================================================================
# 6. DOSE-RESPONSE PLOT (Figure 3.2)
#    Emax curves per agent with splitting point estimates overlaid
#==============================================================================

cat("\n========== GENERATING FIGURE 3.2 ==========\n")

fig_3.2 <- ggplot() +
  # Emax credible bands
  geom_ribbon(data = emax_preds,
              aes(x = Dose_raw, ymin = Lower, ymax = Upper),
              fill = "steelblue", alpha = 0.15) +
  # Emax dose-response curves
  geom_line(data = emax_preds,
            aes(x = Dose_raw, y = Median),
            color = "steelblue", linewidth = 1) +
  # Splitting point estimates with 95% CrI
  geom_point(data = splitting_effects,
             aes(x = Dose, y = Median),
             shape = 16, size = 2.5, color = "black") +
  geom_errorbar(data = splitting_effects,
                aes(x = Dose, ymin = Lower, ymax = Upper),
                width = 30, color = "black", linewidth = 0.5) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -0.5, linetype = "dotted", color = "red", alpha = 0.6) +
  # Facet by agent
  facet_wrap(~ Agent, ncol = 2,
             labeller = labeller(Agent = c(
               "AE" = "Aerobic Exercise",
               "MBE" = "Mind-Body Exercise",
               "ME" = "Mixed Exercise",
               "RE" = "Resistance Exercise"
             ))) +
  labs(
    x = "Exercise Dose (MET-min/week)",
    y = "SMD vs Control (Hedges' g)",
    caption = "Blue line/band: Emax model (median, 95% CrI). Black points: Splitting NMA estimates (95% CrI).\nDotted red line: minimum clinically important difference (SMD = -0.5)."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggsave("figures/Figure_3.2.jpeg", fig_3.2,
       width = 10, height = 8, dpi = 300)

cat("Figure 3.2 saved.\n")

#==============================================================================
# 7. SAVE OUTPUTS
#==============================================================================

save(comparison, splitting_effects, emax_preds,
     med_emax, rankings, disagreement,
     file = "results/effect_comparison_results.RData")

print(comparison_wide, n = 50)
print(disagreement)
print(med_emax) 
print(rankings_wide, n = 50) 

cat("\n========== EFFECT COMPARISON COMPLETE ==========\n")

cat("\n\n========== PRIMARY ANALYSIS COMPLETE ==========\n")


#==============================================================================
# 3.3.3 SENSITIVITY ANALYSIS: TOTAL MINUTES/WEEK
# Repeats primary analysis using alternative dose metric
# Only splitting + 4 MBNMA models (lumping is dose-agnostic, unchanged)
#==============================================================================

set.seed(12345)

# Storage
results_sens <- list()

#==============================================================================
# PREPARE SENSITIVITY NETWORK
#==============================================================================

cat("\n========== SENSITIVITY ANALYSIS: MINUTES/WEEK ==========\n")

# Check dose range
cat("Minutes/week dose range:\n")
print(summary(Tian$Dose_minutes[Tian$Agent != "CON"]))

# Split data for sensitivity (categorized doses stay the same)
# NOTE: splitting uses categorized doses which are based on MET-min/week,
# not minutes/week. For a true sensitivity analysis of the splitting model,
# you would need to re-categorize based on minutes/week.
# Here we only repeat MBNMA models with the continuous minutes/week metric.

# MBNMA data with minutes/week
mbnma_data_min <- Tian %>%
  transmute(
    studyID = Study,
    agent = Agent,
    dose = Dose_minutes,  # total minutes/week (NOT rescaled)
    y = y,
    se = SE,
    n = N
  )

# Check if rescaling is needed (same logic as primary: if max dose^2 is huge)
cat("\nDose range (min/week):", range(mbnma_data_min$dose[mbnma_data_min$agent != "CON"]), "\n")
cat("Max dose^2:", max(mbnma_data_min$dose[mbnma_data_min$agent != "CON"])^2, "\n")

# Rescale to hundreds of minutes/week if needed for quadratic stability
# Adjust this based on the actual range of Dose_minutes
mbnma_data_min_scaled <- mbnma_data_min %>%
  mutate(dose = dose / 100)  # adjust divisor based on actual range

cat("Rescaled dose range:", range(mbnma_data_min_scaled$dose[mbnma_data_min_scaled$agent != "CON"]), "\n")

mbnma_data_min_scaled <- mbnma_data_min_scaled %>%
  filter(!studyID %in% c("Hanssen et al, 2018", "Krogh et al,2012"))

network_sens <- mbnma.network(data.ab = mbnma_data_min_scaled)

#==============================================================================
# METHOD 1: LINEAR MBNMA (SENSITIVITY)
#==============================================================================

cat("\n========== SENSITIVITY: LINEAR ==========\n")

result_linear_sens <- mbnma.run(network_sens,
                                fun = dpoly(degree = 1, beta.1 = "rel"),
                                method = "random",
                                link = "smd",
                                n.iter = 20000,
                                n.chains = 3,
                                pD = TRUE)

cat("\nLinear sensitivity summary:\n")
summary(result_linear_sens)

results_sens$linear <- list(model = result_linear_sens, network = network_sens)

#==============================================================================
# METHOD 2: EMAX MBNMA (SENSITIVITY)
#==============================================================================

cat("\n========== SENSITIVITY: EMAX ==========\n")

result_emax_sens <- mbnma.run(network_sens,
                              fun = demax(emax = "rel", ed50 = "rel"),
                              method = "random",
                              link = "smd",
                              n.iter = 20000,
                              n.chains = 3,
                              pD = TRUE)

cat("\nEmax sensitivity summary:\n")
summary(result_emax_sens)

results_sens$emax <- list(model = result_emax_sens, network = network_sens)

#==============================================================================
# METHOD 3: QUADRATIC MBNMA (SENSITIVITY)
#==============================================================================

cat("\n========== SENSITIVITY: QUADRATIC ==========\n")

quad_fun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))

result_quad_sens <- mbnma.run(network_sens,
                              fun = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
                              method = "random",
                              link = "smd",
                              n.iter = 20000,
                              n.chains = 3,
                              pD = TRUE)

cat("\nQuadratic sensitivity ESS:\n")
print(result_quad_sens$BUGSoutput$summary[, c("mean", "Rhat", "n.eff")])

cat("\nQuadratic sensitivity summary:\n")
summary(result_quad_sens)

results_sens$quadratic <- list(model = result_quad_sens, network = network_sens)

#==============================================================================
# METHOD 4: NCS MBNMA (SENSITIVITY)
#==============================================================================

cat("\n========== SENSITIVITY: NCS ==========\n")

result_rcs_sens <- mbnma.run(network_sens,
                             fun = dspline(type = "ns",
                                           knots = c(0.10, 0.50, 0.90),
                                           beta.1 = "rel",
                                           beta.2 = "rel"),
                             method = "random",
                             link = "smd",
                             n.iter = 20000,
                             n.chains = 3,
                             pD = TRUE)

cat("\nNCS sensitivity summary:\n")
summary(result_rcs_sens)

results_sens$rcs <- list(model = result_rcs_sens, network = network_sens)

#==============================================================================
# SENSITIVITY CONVERGENCE
#==============================================================================

cat("\n========== SENSITIVITY CONVERGENCE ==========\n")

convergence_sens <- data.frame()

for (nm in c("linear", "emax", "quadratic", "rcs")) {
  label <- c(linear = "Linear MBNMA", emax = "Emax MBNMA",
             quadratic = "Quadratic MBNMA", rcs = "NCS MBNMA")[nm]
  bugs <- results_sens[[nm]]$model$BUGSoutput
  rhat_vals <- bugs$summary[, "Rhat"]
  ess_vals <- bugs$summary[, "n.eff"]
  convergence_sens <- rbind(convergence_sens, data.frame(
    Method = label,
    MaxRhat = max(rhat_vals, na.rm = TRUE),
    MinESS = min(ess_vals[ess_vals > 1], na.rm = TRUE),
    Converged = max(rhat_vals, na.rm = TRUE) < 1.05
  ))
}

cat("\nSensitivity convergence:\n")
print(convergence_sens)

#==============================================================================
# SENSITIVITY MODEL COMPARISON
#==============================================================================

cat("\n========== SENSITIVITY MODEL COMPARISON ==========\n")

model_names_sens <- c("Linear MBNMA", "Emax MBNMA", 
                      "Quadratic MBNMA", "NCS MBNMA")
model_comp_sens <- data.frame(Method = model_names_sens)

for (i in 1:4) {
  nm <- c("linear", "emax", "quadratic", "rcs")[i]
  bugs <- results_sens[[nm]]$model$BUGSoutput
  model_comp_sens$DIC[i] <- bugs$DIC
  model_comp_sens$pD[i] <- bugs$pD
  model_comp_sens$ResidDev[i] <- bugs$summary["totresdev", "mean"]
  model_comp_sens$tau[i] <- tryCatch(
    bugs$summary["sd", "50%"],
    error = function(e) NA
  )
}

model_comp_sens$DIC <- round(model_comp_sens$DIC, 1)
model_comp_sens$pD <- round(model_comp_sens$pD, 1)
model_comp_sens$ResidDev <- round(model_comp_sens$ResidDev, 1)
model_comp_sens$tau <- round(model_comp_sens$tau, 3)
model_comp_sens$DeltaDIC <- round(
  model_comp_sens$DIC - min(model_comp_sens$DIC, na.rm = TRUE), 1)

cat("\nSensitivity Model Comparison (minutes/week):\n")
print(model_comp_sens)


#==============================================================================
# SENSITIVITY EMAX PREDICTIONS (if Emax is still best)
#==============================================================================

cat("\n========== SENSITIVITY EMAX PREDICTIONS ==========\n")

# Dose grid — adjust based on rescaled range
# If you divided by 100, and raw range is e.g. 60-450 min/week,
# then scaled range is 0.6-4.5
dose_range_sens <- range(mbnma_data_min_scaled$dose[mbnma_data_min_scaled$agent != "CON"])
dose_grid_sens <- seq(0, max(dose_range_sens), by = 0.01)

dose_list_sens <- list(
  AE = dose_grid_sens,
  MBE = dose_grid_sens,
  ME = dose_grid_sens,
  RE = dose_grid_sens
)

pred_emax_sens <- predict(results_sens$emax$model, E0 = 0,
                          exact.doses = dose_list_sens)

# Extract predictions
emax_preds_sens <- data.frame()

for (ag_name in names(pred_emax_sens$predicts)) {
  if (grepl("Placebo", ag_name)) next
  
  ag_list <- pred_emax_sens$predicts[[ag_name]]
  
  for (dose_name in names(ag_list)) {
    samples <- ag_list[[dose_name]]
    dose_val <- as.numeric(dose_name)
    
    emax_preds_sens <- rbind(emax_preds_sens, data.frame(
      Agent = ag_name,
      Model = "Emax (min/week)",
      Dose_scaled = dose_val,
      Dose_raw = dose_val * 100,  # back-transform to raw minutes/week
      Median = median(samples),
      Lower = quantile(samples, 0.025),
      Upper = quantile(samples, 0.975)
    ))
  }
}

row.names(emax_preds_sens) <- NULL

# MED from sensitivity
med_emax_sens <- emax_preds_sens %>%
  filter(Median <= -0.5) %>%
  group_by(Agent) %>%
  summarise(
    MED_min = min(Dose_raw),
    SMD_at_MED = Median[which.min(Dose_raw)],
    .groups = "drop"
  )

cat("\nSensitivity MED (Emax, minutes/week):\n")
print(med_emax_sens)
#==============================================================================
# SAVE SENSITIVITY RESULTS
#==============================================================================
result_emax$BUGSoutput$summary[grep("emax|ed50", 
                                    rownames(result_emax$BUGSoutput$summary)), c("50%", "2.5%", "97.5%")]

save(results_sens, convergence_sens, model_comp_sens,
     emax_preds_sens, med_emax_sens,
     file = "results/sensitivity_analysis_results.RData")

cat("\n========== SENSITIVITY ANALYSIS COMPLETE ==========\n")
#==============================================================================
# SAVE SENSITIVITY RESULTS
#==============================================================================

save(results_sens, convergence_sens, model_comp_sens,
     emax_preds_sens, med_emax_sens,
     file = "results/sensitivity_analysis_results.RData")

cat("\n========== SENSITIVITY ANALYSIS COMPLETE ==========\n")

#==============================================================================
# CONSISTENCY ASSESSMENT
#==============================================================================

cat("\n========== CONSISTENCY ASSESSMENT ==========\n")

# UME model for splitting
ume_split <- nma.run(network_split, method = "random", link = "smd",
                     n.iter = 20000, UME = TRUE)

cat("Splitting DIC:", results_all$splitting$model$jagsresult$BUGSoutput$DIC, "\n")
cat("UME DIC:", ume_split$jagsresult$BUGSoutput$DIC, "\n")
cat("Difference:", ume_split$jagsresult$BUGSoutput$DIC - 
      results_all$splitting$model$jagsresult$BUGSoutput$DIC, "\n")

# Node-splitting for Emax
emax_nodesplit <- mbnma.nodesplit(network_mbnma,
                                  fun = demax(emax = "rel", ed50 = "rel"),
                                  method = "random", link = "smd")

print(emax_nodesplit)

# Significant inconsistencies (p < 0.05)
ns_summary <- summary(emax_nodesplit)
significant <- ns_summary[ns_summary$Evidence == "Direct" & 
                            ns_summary$p.value < 0.05, ]
cat("\nSignificant inconsistencies (p < 0.05):\n")
print(significant)

save(ume_split, emax_nodesplit,
     file = "results/consistency_results.RData")

cat("\n========== CONSISTENCY ASSESSMENT COMPLETE ==========\n")

#==============================================================================
# DEVIANCE AND LEVERAGE PLOTS (Supplementary)
#==============================================================================

cat("\n========== DEVIANCE AND LEVERAGE PLOTS ==========\n")

# Emax (primary model)
jpeg("figures/deviance_emax.jpeg", width = 2400, height = 1800, res = 300)
devplot(results_all$emax$model)
dev.off()


cat("Deviance and leverage plots saved.\n")

#==============================================================================
# FINAL SAVE: ALL RESULTS
#==============================================================================

save(results_all, convergence_all, model_comparison,
     comparison, splitting_effects, emax_preds,
     med_emax, rankings, disagreement,
     results_sens, convergence_sens, model_comp_sens,
     emax_preds_sens, med_emax_sens,
     ume_split, emax_nodesplit,
     file = "results/all_chapter3_results.RData")

cat("\n========== ALL ANALYSES COMPLETE ==========\n")
cat("\nOutputs saved:\n")
cat("  results/primary_analysis_results.RData\n")
cat("  results/effect_comparison_results.RData\n")
cat("  results/sensitivity_analysis_results.RData\n")
cat("  results/consistency_results.RData\n")
cat("  results/all_chapter3_results.RData\n")
cat("\nFigures:\n")
cat("  figures/Figure 3.1.jpeg (network diagrams)\n")
cat("  figures/Figure_3.2.jpeg (Emax dose-response curves)\n")
cat("  figures/trace_lumping.jpeg\n")
cat("  figures/forest_lumping.jpeg\n")
