################################################################################
# Chapter 3: Evaluating Dose-Response Methods for Exercise Interventions
#
# This script implements:
#   1. Network characteristics and connectivity
#   2. Data preparation (lumping, splitting, MBNMA)
#   3. Network plots
#   4. Primary analysis: 6 models (lumping, splitting, linear, Emax,
#      quadratic, NCS) using MET-min/week
#   5. Effect comparison: splitting vs Emax at standard doses
#   6. Minimum effective dose and agent rankings
#   7. Consistency assessment (UME, node-splitting)
#   8. Sensitivity analysis: minutes/week
#   9. Cross-method comparison: lumping vs Emax at median dose
#  10. Trace plots and diagnostics
#  11. Summary and output
#
# PREREQUISITE: Tian dataset must be loaded with columns:
#   Study, Agent, Dose (categorized), Exact_dose (MET-min/week),
#   Dose_minutes (min/week), y, SE, N
#
# Dependencies: MBNMAdose (>= 0.5.0), gemtc, tidyverse, coda,
#               ggplot2, magick
################################################################################

# ==============================================================================
# 0. Setup
# ==============================================================================

library(tidyverse)
library(MBNMAdose)
library(gemtc)
library(coda)
library(ggplot2)
library(magick)

set.seed(12345)
setwd("C:/Users/katie/Desktop/Depression-NMA/Chapter 3")

theme_diss <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )
theme_set(theme_diss)

# MCMC settings (consistent across all models)
n_iter   <- 20000
n_chains <- 3
n_burnin <- floor(n_iter / 2)
n_thin   <- max(1, floor((n_iter - n_burnin) / 1000))

cat("MCMC settings:\n",
    "  Iterations:", n_iter, "\n",
    "  Burn-in:", n_burnin, "\n",
    "  Thinning:", n_thin, "\n",
    "  Chains:", n_chains, "\n",
    "  Posterior samples:", n_chains * ((n_iter - n_burnin) / n_thin), "\n")


# ==============================================================================
# 1. Network Characteristics
# ==============================================================================

cat("\n=== Network characteristics ===\n")

# Overall
cat("Studies:", n_distinct(Tian$Study), "\n")
cat("Arms:", nrow(Tian), "\n")
cat("Participants:", sum(Tian$N), "\n")

# Dose distribution (exercise arms only)
Tian %>%
  filter(Agent != "CON") %>%
  summarise(
    n_arms   = n(),
    median   = median(Exact_dose, na.rm = TRUE),
    min      = min(Exact_dose, na.rm = TRUE),
    max      = max(Exact_dose, na.rm = TRUE),
    q1       = quantile(Exact_dose, 0.25, na.rm = TRUE),
    q3       = quantile(Exact_dose, 0.75, na.rm = TRUE)
  ) %>%
  print()

# By modality
Tian %>%
  group_by(Agent) %>%
  summarise(
    n_arms     = n(),
    n_participants = sum(N),
    median_n   = median(N),
    median_dose = median(Exact_dose),
    min_dose   = min(Exact_dose),
    max_dose   = max(Exact_dose),
    q1_dose    = quantile(Exact_dose, 0.25, na.rm = TRUE),
    q3_dose    = quantile(Exact_dose, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# Arms per dose category by modality
cat("\n=== Arms by modality and dose category ===\n")
Tian %>%
  filter(Agent != "CON") %>%
  group_by(Agent, Dose) %>%
  summarise(
    k         = n(),
    min_exact = min(Exact_dose),
    max_exact = max(Exact_dose),
    range     = max(Exact_dose) - min(Exact_dose),
    .groups = "drop"
  ) %>%
  arrange(Agent, Dose) %>%
  print(n = Inf)


# ==============================================================================
# 2. Data Preparation
# ==============================================================================

# ---- Splitting: categorized doses ----
split_data <- Tian %>%
  transmute(
    studyID = Study,
    agent   = Agent,
    dose    = Dose,
    y       = y,
    se      = SE,
    n       = N
  )

# ---- MBNMA: continuous doses (rescaled to thousands of MET-min/week) ----
mbnma_data <- Tian %>%
  transmute(
    studyID = Study,
    agent   = Agent,
    dose    = Exact_dose / 1000,
    y       = y,
    se      = SE,
    n       = N
  )

# ---- Lumping: all non-zero doses collapsed to 1 ----
# Excludes 5 studies that compare doses within the same modality
# (no estimable contrast after collapsing)
lumped_data <- Tian %>%
  filter(!Study %in% c("Hanssen et al, 2018",
                       "Krogh et al,2012",
                       "Scott et al, 2019",
                       "Streeter et al,2017",
                       "Streeter, 2020")) %>%
  transmute(
    studyID = Study,
    agent   = Agent,
    dose    = ifelse(Agent == "CON", 0, 1),
    y       = y,
    se      = SE,
    n       = N
  ) %>%
  group_by(studyID, agent, dose) %>%
  summarise(
    y  = weighted.mean(y, w = n),
    se = sqrt(1 / sum(1 / se^2)),
    n  = sum(n),
    .groups = "drop"
  )

# ---- gemtc format (for modality-level network plot only) ----
gemtc_data <- Tian %>%
  group_by(Study, Agent) %>%
  summarise(
    mean       = weighted.mean(y, w = N),
    std.err    = sqrt(1 / sum(1 / SE^2)),
    sampleSize = sum(N),
    .groups = "drop"
  ) %>%
  rename(study = Study, treatment = Agent)

cat("\n=== Data preparation summary ===\n")
cat("Splitting studies:", n_distinct(split_data$studyID), "| arms:", nrow(split_data), "\n")
cat("MBNMA studies:   ", n_distinct(mbnma_data$studyID), "| arms:", nrow(mbnma_data), "\n")
cat("Lumping studies: ", n_distinct(lumped_data$studyID), "| arms:", nrow(lumped_data), "\n")


# ==============================================================================
# 3. Network Plots (Figure 3.5)
# ==============================================================================

# Panel A: modality-level (gemtc)
jpeg("figures/tmp_panelA.jpeg", width = 1800, height = 1500, res = 300)
par(mar = c(2, 2, 2, 2))
network_gemtc <- mtc.network(data.ab = gemtc_data)
plot(network_gemtc, main = "")
mtext("A", side = 3, line = 0.5,
      at = par("usr")[1] + 0.02 * (par("usr")[2] - par("usr")[1]),
      cex = 2, font = 2)
dev.off()

# Panel B: dose-level (MBNMAdose)
jpeg("figures/tmp_panelB.jpeg", width = 1600, height = 1500, res = 300)
par(mar = c(2, 2, 2, 2))
mbnma_network <- mbnma.network(data.ab = split_data)
plot(mbnma_network, main = "")
mtext("B", side = 3, line = 0.5,
      at = par("usr")[1] + 0.02 * (par("usr")[2] - par("usr")[1]),
      cex = 2, font = 2)
dev.off()

# Combine
img_a <- image_read("figures/tmp_panelA.jpeg")
img_b <- image_read("figures/tmp_panelB.jpeg")
combined <- image_append(c(img_a, img_b), stack = FALSE)
image_write(combined, "figures/Figure 3.5.jpeg")
file.remove("figures/tmp_panelA.jpeg", "figures/tmp_panelB.jpeg")

cat("Network plots saved to figures/Figure 3.5.jpeg\n")


# ==============================================================================
# 4. Primary Analysis: Model Fitting
# ==============================================================================

results_all <- list()

# ---- 4a. Lumping (40 studies, dose-agnostic benchmark) ----
cat("\n--- Fitting: Lumping ---\n")
network_lumped <- mbnma.network(data.ab = lumped_data)
result_lumped <- nma.run(network_lumped,
                         method = "random",
                         link = "smd",
                         n.iter = n_iter,
                         n.chains = n_chains,
                         pD = TRUE,
                         jags.seed = 12345)
results_all$lumping <- list(model = result_lumped, network = network_lumped)

# ---- 4b. Splitting (45 studies, categorized doses) ----
cat("\n--- Fitting: Splitting ---\n")
network_split <- mbnma.network(data.ab = split_data)
result_split <- nma.run(network_split,
                        method = "random",
                        link = "smd",
                        n.iter = n_iter,
                        n.chains = n_chains,
                        pD = TRUE,
                        jags.seed = 12345)
results_all$splitting <- list(model = result_split, network = network_split)

# ---- 4c. Linear MBNMA ----
cat("\n--- Fitting: Linear MBNMA ---\n")
network_mbnma <- mbnma.network(data.ab = mbnma_data)
result_linear <- mbnma.run(network_mbnma,
                           fun = dpoly(degree = 1, beta.1 = "rel"),
                           method = "random",
                           link = "smd",
                           n.iter = n_iter,
                           n.chains = n_chains,
                           pD = TRUE,
                           jags.seed = 1234)
results_all$linear <- list(model = result_linear, network = network_mbnma)

# ---- 4d. Emax MBNMA ----
cat("\n--- Fitting: Emax MBNMA ---\n")
result_emax <- mbnma.run(network_mbnma,
                         fun = demax(emax = "rel", ed50 = "rel"),
                         method = "random",
                         link = "smd",
                         n.iter = n_iter,
                         n.chains = n_chains,
                         pD = TRUE,
                         jags.seed = 12345)
results_all$emax <- list(model = result_emax, network = network_mbnma)

# ---- 4e. Quadratic MBNMA ----
cat("\n--- Fitting: Quadratic MBNMA ---\n")
quad_fun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
result_quad <- mbnma.run(network_mbnma,
                         fun = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
                         method = "random",
                         link = "smd",
                         n.iter = n_iter,
                         n.chains = n_chains,
                         pD = TRUE,
                         jags.seed = 12345)
results_all$quadratic <- list(model = result_quad, network = network_mbnma)

# ---- 4f. Natural cubic spline MBNMA ----
cat("\n--- Fitting: NCS MBNMA ---\n")
exercise_doses <- mbnma_data$dose[mbnma_data$agent != "CON" & mbnma_data$dose > 0]
knot_positions <- quantile(exercise_doses, probs = c(0.10, 0.50, 0.90))
cat("NCS knot positions (rescaled):", round(knot_positions, 3), "\n")

result_rcs <- mbnma.run(network_mbnma,
                        fun = dspline(type = "ns",
                                      knots = c(0.10, 0.50, 0.90),
                                      beta.1 = "rel",
                                      beta.2 = "rel"),
                        method = "random",
                        link = "smd",
                        n.iter = n_iter,
                        n.chains = n_chains,
                        pD = TRUE,
                        jags.seed = 12345)
results_all$rcs <- list(model = result_rcs, network = network_mbnma)

cat("\n=== All 6 models fitted ===\n")


# ==============================================================================
# 5. Convergence and Model Comparison
# ==============================================================================

# ---- 5a. Convergence ----
cat("\n=== Convergence summary ===\n")

convergence_all <- data.frame()

# Lumping and splitting (jagsresult$BUGSoutput)
for (nm in c("lumping", "splitting")) {
  bugs <- results_all[[nm]]$model$jagsresult$BUGSoutput
  rhat <- bugs$summary[, "Rhat"]
  ess  <- bugs$summary[, "n.eff"]
  convergence_all <- rbind(convergence_all, data.frame(
    Method    = c(lumping = "Lumping", splitting = "Splitting")[nm],
    MaxRhat   = max(rhat, na.rm = TRUE),
    MinESS    = min(ess[ess > 1], na.rm = TRUE),
    Converged = max(rhat, na.rm = TRUE) < 1.05
  ))
}

# MBNMA models (BUGSoutput directly)
for (nm in c("linear", "emax", "quadratic", "rcs")) {
  label <- c(linear = "Linear MBNMA", emax = "Emax MBNMA",
             quadratic = "Quadratic MBNMA", rcs = "NCS MBNMA")[nm]
  bugs <- results_all[[nm]]$model$BUGSoutput
  rhat <- bugs$summary[, "Rhat"]
  ess  <- bugs$summary[, "n.eff"]
  convergence_all <- rbind(convergence_all, data.frame(
    Method    = label,
    MaxRhat   = max(rhat, na.rm = TRUE),
    MinESS    = min(ess[ess > 1], na.rm = TRUE),
    Converged = max(rhat, na.rm = TRUE) < 1.05
  ))
}

print(convergence_all)


# ---- 5b. Model comparison (5 primary dose-aware models) ----
cat("\n=== Model comparison (45 studies, SMD scale) ===\n")

extract_fit <- function(model) {
  list(
    DIC      = model$BUGSoutput$DIC,
    pD       = model$BUGSoutput$pD,
    resdev   = model$BUGSoutput$summary["totresdev", "mean"],
    tau      = tryCatch(model$BUGSoutput$summary["sd", "50%"],
                        error = function(e) NA)
  )
}

model_comparison <- data.frame(
  Method = c("Splitting", "Linear MBNMA", "Emax MBNMA",
             "Quadratic MBNMA", "NCS MBNMA")
)

# Splitting (jagsresult structure)
split_bugs <- results_all$splitting$model$jagsresult$BUGSoutput
model_comparison[1, c("DIC", "pD", "ResidDev", "tau")] <- c(
  split_bugs$DIC, split_bugs$pV,
  split_bugs$summary["totresdev", "mean"],
  split_bugs$summary["sd", "50%"]
)

# MBNMA models
for (i in 2:5) {
  nm <- c("linear", "emax", "quadratic", "rcs")[i - 1]
  fit <- extract_fit(results_all[[nm]]$model)
  model_comparison[i, c("DIC", "pD", "ResidDev", "tau")] <- c(
    fit$DIC, fit$pD, fit$resdev, fit$tau
  )
}

model_comparison <- model_comparison %>%
  mutate(across(c(DIC, pD, ResidDev), ~ round(.x, 1)),
         tau = round(tau, 3),
         DeltaDIC = round(DIC - min(DIC, na.rm = TRUE), 1))

cat("\nPrimary model comparison:\n")
print(model_comparison)


# ---- 5c. Lumped benchmark ----
cat("\n=== Lumped NMA benchmark (40 studies) ===\n")

lumped_bugs <- results_all$lumping$model$jagsresult$BUGSoutput
lumped_summ <- lumped_bugs$summary

cat("tau:", round(lumped_summ["sd", "50%"], 3), "\n")

d_rows <- lumped_summ[grep("^d\\[", rownames(lumped_summ)), ]
d_rows <- d_rows[-1, ]  # drop d[1] = CON (reference)

results_lumped <- data.frame(
  Modality = c("AE", "MBE", "ME", "RE"),
  SMD      = round(d_rows[, "50%"], 2),
  Lower    = round(d_rows[, "2.5%"], 2),
  Upper    = round(d_rows[, "97.5%"], 2)
) %>%
  arrange(SMD)

cat("\nLumped treatment effects (SMD vs control):\n")
print(results_lumped)


# ==============================================================================
# 6. Effect Comparison: Splitting vs Emax
# ==============================================================================

agents <- c("AE", "MBE", "ME", "RE")
standard_doses <- c(300, 600, 900, 1200)

# ---- 6a. Extract splitting effects ----
split_labs    <- results_all$splitting$model$trt.labs
split_summary <- results_all$splitting$model$jagsresult$BUGSoutput$summary

splitting_effects <- data.frame()
for (j in 2:length(split_labs)) {
  parts <- strsplit(split_labs[j], "_")[[1]]
  param <- paste0("d[", j, "]")
  splitting_effects <- rbind(splitting_effects, data.frame(
    Agent  = parts[1],
    Model  = "Splitting",
    Dose   = as.numeric(parts[2]),
    Median = split_summary[param, "50%"],
    Lower  = split_summary[param, "2.5%"],
    Upper  = split_summary[param, "97.5%"]
  ))
}

# ---- 6b. Emax predictions on fine dose grid ----
dose_grid_scaled <- seq(0, 1.2, by = 0.01)
dose_list <- setNames(lapply(agents, function(a) dose_grid_scaled), agents)

pred_emax <- predict(results_all$emax$model, E0 = 0,
                     exact.doses = dose_list)

emax_preds <- data.frame()
for (ag_name in names(pred_emax$predicts)) {
  if (grepl("Placebo", ag_name)) next
  ag_list <- pred_emax$predicts[[ag_name]]
  for (dose_name in names(ag_list)) {
    samples  <- ag_list[[dose_name]]
    dose_val <- as.numeric(dose_name)
    emax_preds <- rbind(emax_preds, data.frame(
      Agent      = ag_name,
      Model      = "Emax",
      Dose_scaled = dose_val,
      Dose_raw    = dose_val * 1000,
      Median      = median(samples),
      Lower       = quantile(samples, 0.025),
      Upper       = quantile(samples, 0.975)
    ))
  }
}
row.names(emax_preds) <- NULL

# ---- 6c. Compare at standard doses ----
emax_at_standard <- emax_preds %>%
  filter(abs(Dose_raw - 300) < 10 |
           abs(Dose_raw - 600) < 10 |
           abs(Dose_raw - 900) < 10 |
           abs(Dose_raw - 1200) < 10) %>%
  mutate(Dose = round(Dose_raw / 300) * 300) %>%
  dplyr::select(Agent, Model, Dose, Median, Lower, Upper)

splitting_at_standard <- splitting_effects %>%
  filter(Dose %in% standard_doses)

comparison <- bind_rows(splitting_at_standard, emax_at_standard) %>%
  arrange(Agent, Dose, Model)

comparison_wide <- comparison %>%
  mutate(Effect = sprintf("%.2f (%.2f, %.2f)", Median, Lower, Upper)) %>%
  dplyr::select(Agent, Dose, Model, Effect) %>%
  pivot_wider(names_from = Model, values_from = Effect)

cat("\n=== Splitting vs Emax at standard doses ===\n")
print(comparison_wide, n = 50)

# Disagreement
disagreement <- comparison %>%
  dplyr::select(Agent, Dose, Model, Median) %>%
  pivot_wider(names_from = Model, values_from = Median) %>%
  mutate(Diff = Emax - Splitting,
         Flagged = abs(Diff) > 0.3)

cat("\nDisagreement (|diff| > 0.3 SMD flagged):\n")
print(disagreement)


# ---- 6d. Minimum effective dose (SMD ≤ -0.5) ----
med_emax <- emax_preds %>%
  filter(Median <= -0.5) %>%
  group_by(Agent) %>%
  summarise(MED = min(Dose_raw),
            SMD_at_MED = Median[which.min(Dose_raw)],
            .groups = "drop")

cat("\n=== Minimum effective dose (Emax, SMD ≤ -0.5) ===\n")
print(med_emax)

agents_no_med <- setdiff(agents, med_emax$Agent)
if (length(agents_no_med) > 0) {
  cat("Agents not reaching threshold:", paste(agents_no_med, collapse = ", "), "\n")
}


# ---- 6e. Agent rankings at standard doses ----
rankings <- comparison %>%
  group_by(Model, Dose) %>%
  mutate(Rank = base::rank(Median)) %>%
  ungroup()

rankings_wide <- rankings %>%
  dplyr::select(Agent, Dose, Model, Rank) %>%
  pivot_wider(names_from = Model, values_from = Rank) %>%
  mutate(Rank_Shift = abs(Emax - Splitting))

cat("\n=== Agent rankings (Splitting vs Emax) ===\n")
print(rankings_wide, n = 50)


# ---- 6f. Dose-response plot (Figure 3.6) ----
fig_dr <- ggplot() +
  geom_ribbon(data = emax_preds,
              aes(x = Dose_raw, ymin = Lower, ymax = Upper),
              fill = "steelblue", alpha = 0.15) +
  geom_line(data = emax_preds,
            aes(x = Dose_raw, y = Median),
            color = "steelblue", linewidth = 1) +
  geom_point(data = splitting_effects,
             aes(x = Dose, y = Median),
             shape = 16, size = 2.5, color = "black") +
  geom_errorbar(data = splitting_effects,
                aes(x = Dose, ymin = Lower, ymax = Upper),
                width = 30, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -0.5, linetype = "dotted", color = "red", alpha = 0.6) +
  facet_wrap(~ Agent, ncol = 2,
             labeller = labeller(Agent = c(
               "AE" = "Aerobic Exercise", "MBE" = "Mind-Body Exercise",
               "ME" = "Mixed Exercise", "RE" = "Resistance Exercise"
             ))) +
  labs(x = "Exercise Dose (MET-min/week)",
       y = "SMD vs Control (Hedges' g)",
       caption = paste("Blue line/band: Emax model (median, 95% CrI).",
                       "Black points: Splitting NMA (95% CrI).",
                       "Dotted red line: MCID (SMD = -0.5).")) +
  theme_diss +
  theme(strip.text = element_text(face = "bold"),
        plot.caption = element_text(size = 8, hjust = 0))

ggsave("figures/Figure 3.6.jpeg", fig_dr,
       width = 10, height = 8, dpi = 300)
cat("Dose-response plot saved.\n")


# ==============================================================================
# 7. Consistency Assessment
# ==============================================================================

cat("\n=== Consistency assessment ===\n")

# UME model for splitting
ume_split <- nma.run(network_split, method = "random", link = "smd",
                     n.iter = n_iter, UME = TRUE)

cat("Splitting DIC:", results_all$splitting$model$jagsresult$BUGSoutput$DIC, "\n")
cat("UME DIC:     ", ume_split$jagsresult$BUGSoutput$DIC, "\n")
cat("Difference:  ",
    ume_split$jagsresult$BUGSoutput$DIC -
      results_all$splitting$model$jagsresult$BUGSoutput$DIC, "\n")

# Node-splitting for Emax
emax_nodesplit <- mbnma.nodesplit(network_mbnma,
                                  fun = demax(emax = "rel", ed50 = "rel"),
                                  method = "random", link = "smd")

print(emax_nodesplit)

ns_summary <- summary(emax_nodesplit)
significant <- ns_summary[ns_summary$Evidence == "Direct" &
                            ns_summary$p.value < 0.05, ]
cat("\nSignificant inconsistencies (p < 0.05):\n")
print(significant)


# ==============================================================================
# 8. Sensitivity Analysis: Minutes/Week
# ==============================================================================

cat("\n=== Sensitivity analysis: minutes/week ===\n")

results_sens <- list()

# Prepare data (rescaled to hundreds of min/week)
# Exclude 2 studies with identical min/week across intensity-varying arms
mbnma_data_min <- Tian %>%
  filter(!Study %in% c("Hanssen et al, 2018", "Krogh et al,2012")) %>%
  transmute(
    studyID = Study,
    agent   = Agent,
    dose    = Dose_minutes / 100,
    y       = y,
    se      = SE,
    n       = N
  )

cat("Sensitivity studies:", n_distinct(mbnma_data_min$studyID), "\n")
cat("Dose range (rescaled):",
    range(mbnma_data_min$dose[mbnma_data_min$agent != "CON"]), "\n")

network_sens <- mbnma.network(data.ab = mbnma_data_min)

# ---- 8a. Linear ----
cat("\n--- Sensitivity: Linear ---\n")
result_linear_sens <- mbnma.run(network_sens,
                                fun = dpoly(degree = 1, beta.1 = "rel"),
                                method = "random", link = "smd",
                                n.iter = n_iter, n.chains = n_chains,
                                pD = TRUE, jags.seed = 12345)
results_sens$linear <- list(model = result_linear_sens, network = network_sens)

# ---- 8b. Emax ----
cat("\n--- Sensitivity: Emax ---\n")
result_emax_sens <- mbnma.run(network_sens,
                              fun = demax(emax = "rel", ed50 = "rel"),
                              method = "random", link = "smd",
                              n.iter = n_iter, n.chains = n_chains,
                              pD = TRUE, jags.seed = 12345)
results_sens$emax <- list(model = result_emax_sens, network = network_sens)

# ---- 8c. Quadratic ----
cat("\n--- Sensitivity: Quadratic ---\n")
result_quad_sens <- mbnma.run(network_sens,
                              fun = duser(fun = quad_fun, beta.1 = "rel",
                                          beta.2 = "rel"),
                              method = "random", link = "smd",
                              n.iter = n_iter, n.chains = n_chains,
                              pD = TRUE, jags.seed = 12345)
results_sens$quadratic <- list(model = result_quad_sens, network = network_sens)

# ---- 8d. NCS ----
cat("\n--- Sensitivity: NCS ---\n")
result_rcs_sens <- mbnma.run(network_sens,
                             fun = dspline(type = "ns",
                                           knots = c(0.10, 0.50, 0.90),
                                           beta.1 = "rel", beta.2 = "rel"),
                             method = "random", link = "smd",
                             n.iter = n_iter, n.chains = n_chains,
                             pD = TRUE, jags.seed = 12345)
results_sens$rcs <- list(model = result_rcs_sens, network = network_sens)


# ---- 8e. Sensitivity convergence ----
convergence_sens <- data.frame()
for (nm in c("linear", "emax", "quadratic", "rcs")) {
  label <- c(linear = "Linear MBNMA", emax = "Emax MBNMA",
             quadratic = "Quadratic MBNMA", rcs = "NCS MBNMA")[nm]
  bugs <- results_sens[[nm]]$model$BUGSoutput
  rhat <- bugs$summary[, "Rhat"]
  ess  <- bugs$summary[, "n.eff"]
  convergence_sens <- rbind(convergence_sens, data.frame(
    Method    = label,
    MaxRhat   = max(rhat, na.rm = TRUE),
    MinESS    = min(ess[ess > 1], na.rm = TRUE),
    Converged = max(rhat, na.rm = TRUE) < 1.05
  ))
}

cat("\nSensitivity convergence:\n")
print(convergence_sens)


# ---- 8f. Sensitivity model comparison ----
model_comp_sens <- data.frame(
  Method = c("Linear MBNMA", "Emax MBNMA", "Quadratic MBNMA", "NCS MBNMA")
)

for (i in 1:4) {
  nm <- c("linear", "emax", "quadratic", "rcs")[i]
  fit <- extract_fit(results_sens[[nm]]$model)
  model_comp_sens[i, c("DIC", "pD", "ResidDev", "tau")] <- c(
    fit$DIC, fit$pD, fit$resdev, fit$tau
  )
}

model_comp_sens <- model_comp_sens %>%
  mutate(across(c(DIC, pD, ResidDev), ~ round(.x, 1)),
         tau = round(tau, 3),
         DeltaDIC = round(DIC - min(DIC, na.rm = TRUE), 1))

cat("\nSensitivity model comparison (min/week):\n")
print(model_comp_sens)


# ---- 8g. Sensitivity Emax predictions and MED ----
dose_range_sens <- range(mbnma_data_min$dose[mbnma_data_min$agent != "CON"])
dose_grid_sens  <- seq(0, max(dose_range_sens), by = 0.01)
dose_list_sens  <- setNames(lapply(agents, function(a) dose_grid_sens), agents)

pred_emax_sens <- predict(results_sens$emax$model, E0 = 0,
                          exact.doses = dose_list_sens)

emax_preds_sens <- data.frame()
for (ag_name in names(pred_emax_sens$predicts)) {
  if (grepl("Placebo", ag_name)) next
  ag_list <- pred_emax_sens$predicts[[ag_name]]
  for (dose_name in names(ag_list)) {
    samples  <- ag_list[[dose_name]]
    dose_val <- as.numeric(dose_name)
    emax_preds_sens <- rbind(emax_preds_sens, data.frame(
      Agent       = ag_name,
      Model       = "Emax (min/week)",
      Dose_scaled = dose_val,
      Dose_raw    = dose_val * 100,
      Median      = median(samples),
      Lower       = quantile(samples, 0.025),
      Upper       = quantile(samples, 0.975)
    ))
  }
}
row.names(emax_preds_sens) <- NULL

med_emax_sens <- emax_preds_sens %>%
  filter(Median <= -0.5) %>%
  group_by(Agent) %>%
  summarise(MED_min = min(Dose_raw),
            SMD_at_MED = Median[which.min(Dose_raw)],
            .groups = "drop")

cat("\nSensitivity MED (Emax, min/week):\n")
print(med_emax_sens)


# ==============================================================================
# 9. Cross-Method Comparison: Lumping vs Emax at Median Dose
# ==============================================================================

cat("\n=== Lumping vs Emax at median dose ===\n")

# Primary: median MET-min/week
median_dose_raw <- median(Tian$Exact_dose[Tian$Agent != "CON"], na.rm = TRUE)
cat("Median exercise dose (MET-min/week):", median_dose_raw, "\n")

emax_at_median <- emax_preds %>%
  group_by(Agent) %>%
  slice_min(abs(Dose_raw - median_dose_raw), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Agent, Dose_raw, Median, Lower, Upper)

# Sensitivity: median min/week
median_dose_min <- median(Tian$Dose_minutes[Tian$Agent != "CON"], na.rm = TRUE)
cat("Median exercise dose (min/week):", median_dose_min, "\n")

emax_at_median_sens <- emax_preds_sens %>%
  group_by(Agent) %>%
  slice_min(abs(Dose_raw - median_dose_min), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Agent, Dose_raw, Median, Lower, Upper)

# Combined table
combined_comparison <- emax_at_median %>%
  dplyr::select(Agent, SMD_Emax = Median, Lower_Emax = Lower, Upper_Emax = Upper) %>%
  left_join(
    emax_at_median_sens %>%
      dplyr::select(Agent, SMD_Sens = Median, Lower_Sens = Lower, Upper_Sens = Upper),
    by = "Agent"
  ) %>%
  left_join(
    results_lumped %>%
      dplyr::select(Modality, SMD_Lumped = SMD, Lower_Lumped = Lower, Upper_Lumped = Upper),
    by = c("Agent" = "Modality")
  ) %>%
  mutate(
    Emax_MET = sprintf("%.2f (%.2f, %.2f)", SMD_Emax, Lower_Emax, Upper_Emax),
    Emax_Min = sprintf("%.2f (%.2f, %.2f)", SMD_Sens, Lower_Sens, Upper_Sens),
    Lumped   = sprintf("%.2f (%.2f, %.2f)", SMD_Lumped, Lower_Lumped, Upper_Lumped)
  ) %>%
  arrange(SMD_Emax) %>%
  dplyr::select(Agent, Emax_MET, Emax_Min, Lumped)

cat("\n--- Comparison table ---\n")
cat("Emax_MET = Emax at median", round(median_dose_raw), "MET-min/week\n")
cat("Emax_Min = Emax at median", round(median_dose_min), "min/week\n")
cat("Lumped   = Lumped NMA (dose-agnostic)\n\n")
print(combined_comparison)


# ==============================================================================
# 10. Trace Plots (Emax model diagnostics)
# ==============================================================================

sims <- results_all$emax$model$BUGSoutput$sims.array

# Emax parameters
params_emax <- c("emax[2]", "emax[3]", "emax[4]", "emax[5]")
mcmc_emax <- mcmc.list(
  mcmc(sims[, 1, params_emax]),
  mcmc(sims[, 2, params_emax]),
  mcmc(sims[, 3, params_emax])
)

jpeg("figures/Appendix G_1.jpeg", width = 3000, height = 3200, res = 300)
par(mar = c(2, 2, 2, 1))
plot(mcmc_emax, ask = FALSE)
dev.off()

# ED50 parameters
params_ed50 <- c("ed50[2]", "ed50[3]", "ed50[4]", "ed50[5]")
mcmc_ed50 <- mcmc.list(
  mcmc(sims[, 1, params_ed50]),
  mcmc(sims[, 2, params_ed50]),
  mcmc(sims[, 3, params_ed50])
)

jpeg("figures/Appendix G_2.jpeg", width = 3000, height = 3200, res = 300)
par(mar = c(2, 2, 2, 1))
plot(mcmc_ed50, ask = FALSE)
dev.off()

# SD and total residual deviance
params_other <- c("sd", "totresdev")
mcmc_other <- mcmc.list(
  mcmc(sims[, 1, params_other]),
  mcmc(sims[, 2, params_other]),
  mcmc(sims[, 3, params_other])
)

jpeg("figures/Appendix G_3.jpeg", width = 3000, height = 1800, res = 300)
par(mar = c(2, 2, 2, 1))
plot(mcmc_other, ask = FALSE)
dev.off()

# Deviance plot
jpeg("figures/Appendix H.jpeg", width = 2400, height = 1800, res = 300)
devplot(results_all$emax$model)
dev.off()

cat("Trace and deviance plots saved.\n")


# ==============================================================================
# 11. Emax Parameter Summary
# ==============================================================================

cat("\n=== Emax parameter estimates ===\n")
emax_params <- result_emax$BUGSoutput$summary[
  grep("emax|ed50", rownames(result_emax$BUGSoutput$summary)),
  c("50%", "2.5%", "97.5%")
]
print(round(emax_params, 2))

cat("\nSensitivity Emax parameters:\n")
sens_params <- result_emax_sens$BUGSoutput$summary[
  grep("emax|ed50", rownames(result_emax_sens$BUGSoutput$summary)),
  c("50%", "2.5%", "97.5%")
]
print(round(sens_params, 2))


# ==============================================================================
# 12. Summary and Save
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  RESULTS SUMMARY \u2014 Chapter 3\n")
cat("================================================================\n\n")

cat("--- Model comparison (primary) ---\n")
print(model_comparison)

cat("\n--- Model comparison (sensitivity) ---\n")
print(model_comp_sens)

cat("\n--- Lumped benchmark ---\n")
print(results_lumped)

cat("\n--- Minimum effective dose (primary) ---\n")
print(med_emax)

cat("\n--- Minimum effective dose (sensitivity) ---\n")
print(med_emax_sens)

cat("\n--- Lumping vs Emax at median dose ---\n")
print(combined_comparison)

cat("\n================================================================\n")
cat("  Analysis complete.\n")
cat("================================================================\n")

# ---- Save ----
save(results_all, convergence_all, model_comparison,
     comparison, splitting_effects, emax_preds,
     med_emax, rankings, disagreement,
     results_lumped, results_sens, convergence_sens, model_comp_sens,
     emax_preds_sens, med_emax_sens,
     ume_split, emax_nodesplit,
     combined_comparison,
     file = "results/all_chapter3_results.RData")

cat("\n=== All results saved ===\n")
cat("  Figures: figures/\n")
cat("  Results: results/\n")

writeLines(capture.output(sessionInfo()), "results/session_info.txt")