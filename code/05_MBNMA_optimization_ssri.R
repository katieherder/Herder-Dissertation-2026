################################################################################
# Chapter 4: SSRI Empirical Analysis
# Optimizing Model-Based Network Meta-Analysis for Pharmacologic Networks
#
# This script implements:
#   1. Network setup
#   2. Model fitting (Models A, B, Ref)
#   3. Model comparison (DIC, LPML)
#   4. Dose-response estimation and plotting
#   5. Agent rankings (SUCRA) at clinically equivalent doses
#   6. Prior sensitivity analysis on tau
#
# Dependencies: MBNMAdose (>= 0.5.0), R2jags, ggplot2, dplyr, tidyr, 
#               kableExtra, patchwork
################################################################################

# ==============================================================================
# 0. Setup
# ==============================================================================

library(MBNMAdose)
library(R2jags)
library(ggplot2)
library(dplyr)
library(tidyr)
library(kableExtra)
library(patchwork)
library(coda)        # for MCMC diagnostics
library(stringr)  # for str_replace_all and str_to_title

set.seed(12345)   # reproducibility
setwd("C:/Users/katie/Desktop/Depression-NMA")

# ggplot theme for dissertation figures
theme_diss <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )
theme_set(theme_diss)


# ==============================================================================
# 1. Create Network Object
# ==============================================================================
# Dose equivalence table (fluoxetine 20mg = 1 equivalent unit)
# From Hayasaka et al. (2015) / Furukawa et al. (2019)
dose_equiv <- data.frame(
  agent       = c("citalopram", "escitalopram", "fluoxetine", "paroxetine", "sertraline"),
  equiv_mg    = c(20, 9, 20, 17, 49.3),  # mg/day equivalent to fluoxetine 20mg
  stringsAsFactors = FALSE
)

ssri_equiv <- ssri %>%
  left_join(dose_equiv, by = "agent") %>%
  mutate(dose = ifelse(agent == "placebo", 0, dose / equiv_mg * 20)) %>%
  select(-equiv_mg)

ssri_net <- mbnma.network(ssri_equiv)

# ==============================================================================
# 2. Model Fitting
# ==============================================================================

# MCMC settings (consistent with Chapter 3)
n_iter   <- 20000
n_chains <- 3
n_burnin <- floor(n_iter / 2)   # 10,000
n_thin   <- max(1, floor((n_iter - n_burnin) / 1000))  # 10

cat("MCMC settings:\n",
    "  Iterations:", n_iter, "\n",
    "  Burn-in:", n_burnin, "\n",
    "  Thinning:", n_thin, "\n",
    "  Chains:", n_chains, "\n",
    "  Posterior samples:", n_chains * ((n_iter - n_burnin) / n_thin), "\n")

# ---- Model A: Shared Emax across all agents ----
# All agents share common Emax and ED50 parameters (pooled)
# emax = "rel" (relative effects pooled across agents)
# ed50 = "common" (single ED50 estimated for all agents)
# This is the standard MBNMA approach from Chapter 2's review
cat("\n--- Fitting Model A: Shared Emax ---\n")
model_A <- mbnma.run(
  network  = ssri_net,
  fun      = demax(emax = "rel", ed50 = "common"),
  method   = "random",
  likelihood = "binomial",
  link     = "logit",
  parameters.to.save = c("resdev", "emax", "ed50", "sd", "totresdev", "theta"),
  n.iter   = n_iter,
  n.chains = n_chains,
  n.burnin = n_burnin,
  n.thin   = n_thin,
  jags.seed = 123456
)

print(model_A)
cat("Model A priors:\n")

# ---- Model B: Agent-specific Emax (dmulti) ----
# Each agent gets its own Emax curve — no borrowing of strength
# This isolates the effect of parameter pooling vs. functional form
cat("\n--- Fitting Model B: Agent-specific Emax (dmulti) ---\n")
model_B <- mbnma.run(
  network  = ssri_net,
  fun      = dmulti(
    list(
      "citalopram"   = demax(),
      "escitalopram" = demax(),
      "fluoxetine"   = demax(),
      "paroxetine"   = demax(),
      "sertraline"   = demax()
    )
  ),
  method   = "random",
  likelihood = "binomial",
  link     = "logit",
  parameters.to.save = c("resdev", "emax", "ed50", "sd", "totresdev", "theta"),
  n.iter   = n_iter,
  n.chains = n_chains,
  n.burnin = n_burnin,
  n.thin   = n_thin,
  jags.seed = 123456
)

print(model_B)


# ---- Reference Model (Ref): Shared Linear ----
# Deliberate misspecification benchmark
# Linear = dpoly(degree=1) is equivalent to standard NMA
cat("\n--- Fitting Reference Model: Shared Linear ---\n")
model_Ref <- mbnma.run(
  network  = ssri_net,
  fun      = dpoly(degree = 1),
  method   = "random",
  likelihood = "binomial",
  link     = "logit",
  parameters.to.save = c("resdev", "beta.1", "sd", "totresdev", "theta"),
  n.iter   = n_iter,
  n.chains = n_chains,
  n.burnin = n_burnin,
  n.thin   = n_thin,
  jags.seed = 123456
)

print(model_Ref)

# ---- Lumped NMA: collapse all non-zero doses ----
cat("\n--- Fitting Lumped NMA ---\n")

ssri_lumped <- ssri_equiv
ssri_lumped$dose[ssri_lumped$dose > 0] <- 1

# Drop studies with duplicate agent-dose arms after collapsing
ssri_lumped <- ssri_lumped %>%
  group_by(studyID, agent, dose) %>%
  mutate(dup = n() > 1) %>%
  group_by(studyID) %>%
  filter(!any(dup & dose > 0)) %>%
  ungroup() %>%
  dplyr::select(-dup) %>%
  group_by(studyID) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  as.data.frame()

ssri_net_lumped <- mbnma.network(ssri_lumped)

model_Lumped <- nma.run(
  network    = ssri_net_lumped,
  method     = "random",
  likelihood = "binomial",
  link       = "logit",
  n.iter     = n_iter,
  n.chains   = n_chains,
  n.burnin   = n_burnin,
  n.thin     = n_thin,
  jags.seed  = 123456
)

# ==============================================================================
# 3. Convergence Assessment
# ==============================================================================

check_convergence <- function(model, model_name, is_nma = FALSE) {
  cat("\n=== Convergence diagnostics for", model_name, "===\n")
  
  # nma objects store BUGS output differently
  if (is_nma) {
    bugs <- model$jagsresult$BUGSoutput
  } else {
    bugs <- model$BUGSoutput
  }
  
  if (!is.null(bugs)) {
    rhat_vals <- bugs$summary[, "Rhat"]
    rhat_vals <- rhat_vals[!is.na(rhat_vals)]
    
    cat("  Max Rhat:", round(max(rhat_vals), 4), "\n")
    cat("  Parameters with Rhat > 1.05:", 
        sum(rhat_vals > 1.05), "of", length(rhat_vals), "\n")
    cat("  Parameters with Rhat > 1.10:", 
        sum(rhat_vals > 1.10), "of", length(rhat_vals), "\n")
    
    neff_vals <- bugs$summary[, "n.eff"]
    neff_vals <- neff_vals[!is.na(neff_vals) & neff_vals > 0]
    cat("  Min n.eff:", min(neff_vals), "\n")
    cat("  Median n.eff:", median(neff_vals), "\n")
  }
  
  invisible(NULL)
}

check_convergence(model_A, "Model A (Shared Emax)")
check_convergence(model_B, "Model B (Agent-specific Emax)")
check_convergence(model_Ref, "Ref (Shared Linear)")
check_convergence(model_Lumped, "Lumped NMA", is_nma = TRUE)


# ==============================================================================
# 4. Model Comparison: DIC and LPML
# ==============================================================================

# ---- 4a. DIC Comparison ----
lumped_pD <- model_Lumped$jagsresult$BUGSoutput$DIC - model_Lumped$jagsresult$BUGSoutput$mean$deviance

dic_table <- data.frame(
  Model = c("A: Shared Emax", "B: Agent-specific Emax", "Ref: Shared Linear", "Lumped NMA"),
  DIC   = c(model_A$BUGSoutput$DIC,
            model_B$BUGSoutput$DIC,
            model_Ref$BUGSoutput$DIC,
            model_Lumped$jagsresult$BUGSoutput$DIC),
  pD    = c(model_A$BUGSoutput$pD,
            model_B$BUGSoutput$pD,
            model_Ref$BUGSoutput$pD,
            lumped_pD)
)
dic_table$Dbar <- dic_table$DIC - dic_table$pD

cat("\n=== DIC Comparison ===\n")
print(dic_table, digits = 2)

# ---- 4b. LPML via CPO (Conditional Predictive Ordinate) ----
# 
# MBNMAdose does not compute LPML natively. We compute it from the
# arm-level deviance contributions stored in the model output.
#
# For a binomial outcome with logit link:
#   dev_ik = -2 * [r_ik * log(p_ik) + (n_ik - r_ik) * log(1 - p_ik)]
#            + 2 * [r_ik * log(r_ik/n_ik) + (n_ik - r_ik) * log(1 - r_ik/n_ik)]
#
# The CPO for observation i is estimated via the harmonic mean:
#   CPO_i = 1 / E[1/L(y_i | theta)]
# where the expectation is over posterior draws.
#
# From the deviance residual at each iteration:
#   log L_i = -dev_i / 2  (up to the saturated model constant)
#
# LPML = sum(log(CPO_i))
#
# Since MBNMAdose monitors 'dev' (per-arm deviance contributions),
# we can extract these from the MCMC output.

compute_lpml <- function(model, model_name = "") {
  # Extract per-arm deviance contributions across MCMC iterations
  # 'dev' should be monitored by default in MBNMAdose
  
  mcmc_array <- model$BUGSoutput$sims.array
  param_names <- dimnames(mcmc_array)[[3]]
  
  # Find dev parameters
  dev_idx <- grep("^resdev\\[", param_names)
  
  if (length(dev_idx) == 0) {
    cat("Warning: No 'resdev' parameters found for", model_name, "\n")
    cat("Available parameters:", head(param_names, 20), "\n")
    return(list(lpml = NA, cpo = NA))
  }
  
  # Combine chains: dimensions [iterations, chains, parameters] -> [total_iter, n_arms]
  n_iter_saved <- dim(mcmc_array)[1]
  n_chains_out <- dim(mcmc_array)[2]
  n_arms <- length(dev_idx)
  
  # Stack chains
  dev_matrix <- matrix(NA, nrow = n_iter_saved * n_chains_out, ncol = n_arms)
  for (ch in 1:n_chains_out) {
    rows <- ((ch - 1) * n_iter_saved + 1):(ch * n_iter_saved)
    dev_matrix[rows, ] <- mcmc_array[, ch, dev_idx]
  }
  
  # CPO via harmonic mean estimator
  # dev_ik = -2 * log(L_ik / L_saturated)
  # => log(L_ik) = -dev_ik / 2 + constant (the constant cancels in CPO)
  # CPO_i = 1 / mean(1 / L_i) = 1 / mean(exp(dev_i / 2))
  # log(CPO_i) = -log(mean(exp(dev_i / 2)))
  
  log_cpo <- numeric(n_arms)
  for (j in 1:n_arms) {
    # Numerically stable computation using log-sum-exp
    half_dev <- dev_matrix[, j] / 2
    max_hd <- max(half_dev)
    log_cpo[j] <- -(max_hd + log(mean(exp(half_dev - max_hd))))
  }
  
  lpml <- sum(log_cpo)
  
  cat(sprintf("LPML for %s: %.2f (based on %d arms)\n", model_name, lpml, n_arms))
  
  return(list(lpml = lpml, cpo = log_cpo))
}

lpml_A   <- compute_lpml(model_A, "Model A")
lpml_B   <- compute_lpml(model_B, "Model B")
lpml_Ref <- compute_lpml(model_Ref, "Ref")


# ---- 4c. Combined Model Comparison Table ----
comparison_table <- data.frame(
  Model = c("A: Shared Emax", "B: Agent-specific Emax", "Ref: Shared Linear", "Lumped NMA"),
  DIC   = round(c(model_A$BUGSoutput$DIC,
                  model_B$BUGSoutput$DIC,
                  model_Ref$BUGSoutput$DIC,
                  model_Lumped$jagsresult$BUGSoutput$DIC), 1),
  pD    = round(c(model_A$BUGSoutput$pD,
                  model_B$BUGSoutput$pD,
                  model_Ref$BUGSoutput$pD,
                  lumped_pD), 1),
  LPML  = round(c(lpml_A$lpml, lpml_B$lpml, lpml_Ref$lpml, NA), 1)
)

cat("\n=== Combined Model Comparison ===\n")
print(comparison_table)

jpeg("figures/Chapter 4 DevDev.jpeg", width = 8, height = 6, units = "in", res = 300)
devdev(model_A, model_Ref, dev.type = "resdev")
dev.off()

# ==============================================================================
# 5. Dose-Response Estimation and Visualization
# ==============================================================================

# ---- 5a. Predicted dose-response curves ----

# Define placebo response for predictions
# Use the observed pooled placebo response rate
placebo_rate <- sum(ssri$r[ssri$agent == "placebo"]) / sum(ssri$n[ssri$agent == "placebo"])
cat("\nPooled placebo response rate:", round(placebo_rate, 3), "\n")

# Predict from each model
# E0 = placebo response on log-odds scale for logit link
E0_logit <- log(placebo_rate / (1 - placebo_rate))

# Model A predictions
pred_A <- predict(model_A, E0 = placebo_rate, n.doses = 30)
plot(pred_A) + ggtitle("Model A: Shared Emax — Predicted Dose-Response")

# Model B predictions
pred_B <- predict(model_B, E0 = placebo_rate, n.doses = 30)
plot(pred_B) + ggtitle("Model B: Agent-specific Emax — Predicted Dose-Response")

# Reference model predictions
pred_Ref <- predict(model_Ref, E0 = placebo_rate, n.doses = 30)
plot(pred_Ref) + ggtitle("Ref: Shared Linear — Predicted Dose-Response")


# ---- 5b. Extract and plot dose-response at equivalent doses ----
# Predict at clinically meaningful equivalent doses 
# (0.5x, 1x, 1.5x, 2x the standard equivalent dose)

equiv_doses <- list(
  "citalopram"   = c(0, 10, 20, 30, 40, 60),
  "escitalopram" = c(0, 4, 9, 14, 18, 28),
  "fluoxetine"   = c(0, 5, 10, 20, 40, 60),
  "paroxetine"   = c(0, 10, 17, 20, 30, 40),
  "sertraline"   = c(0, 50, 100, 150, 200)
)

pred_A_equiv <- predict(model_A, E0 = placebo_rate, exact.doses = equiv_doses)
pred_B_equiv <- predict(model_B, E0 = placebo_rate, exact.doses = equiv_doses)


# ---- 5c. Dose-response curve plots ----
plot(pred_A)
plot(pred_B)
plot(pred_Ref)

jpeg("figures/ssri/Figure 4.5.jpeg", width = 10, height = 6, units = "in", res = 300)
plot(pred_A) + ggtitle("Model A: Shared Emax")
dev.off()

jpeg("figures/ssri/Figure 4.6.jpeg", width = 10, height = 6, units = "in", res = 300)
plot(pred_B) + ggtitle("Model B: Agent-specific Emax")
dev.off()

jpeg("figures/ssri/Appendix N.jpeg", width = 10, height = 6, units = "in", res = 300)
plot(pred_Ref) + ggtitle("Ref: Shared Linear")
dev.off()

# ---- 5d. Fitted values overlay (goodness-of-fit) ----
fitplot(model_A)
fitplot(model_B)
fitplot(model_Ref)

jpeg("figures/ssri/Figure 4.5.jpeg", width = 10, height = 8, units = "in", res = 300)
fitplot(model_A)
dev.off()

jpeg("figures/ssri/Chapter 4 Fitplot_B.jpeg", width = 10, height = 8, units = "in", res = 300)
fitplot(model_B)
dev.off()

jpeg("figures/ssri/Chapter 4 Fitplot_Ref.jpeg", width = 10, height = 8, units = "in", res = 300)
fitplot(model_Ref)
dev.off()

# ==============================================================================
# 6. Agent Rankings
# ==============================================================================

# ---- Max observed dose per agent (primary ranking comparison) ----
agent_names <- c("citalopram", "escitalopram", "fluoxetine", "paroxetine", "sertraline")

max_rank_doses <- list(
  "citalopram"   = max(ssri_equiv$dose[ssri_equiv$agent == "citalopram"]),
  "escitalopram" = max(ssri_equiv$dose[ssri_equiv$agent == "escitalopram"]),
  "fluoxetine"   = max(ssri_equiv$dose[ssri_equiv$agent == "fluoxetine"]),
  "paroxetine"   = max(ssri_equiv$dose[ssri_equiv$agent == "paroxetine"]),
  "sertraline"   = max(ssri_equiv$dose[ssri_equiv$agent == "sertraline"])
)

cat("\n=== Max observed equivalent doses for ranking ===\n")
print(data.frame(agent = names(max_rank_doses), max_dose = unlist(max_rank_doses)))

# Also keep 1x equivalent doses for secondary comparison
rank_doses <- list(
  "citalopram" = 20, "escitalopram" = 20, "fluoxetine" = 20,
  "paroxetine" = 20, "sertraline" = 20
)

rank_doses_half <- list(
  "citalopram" = 10, "escitalopram" = 10, "fluoxetine" = 10,
  "paroxetine" = 10, "sertraline" = 10
)

rank_doses_double <- list(
  "citalopram" = 40, "escitalopram" = 40, "fluoxetine" = 40,
  "paroxetine" = 40, "sertraline" = 40
)

# ---- 6a. Rankings at max observed dose ----
cat("\n=== Rankings from Model A at max observed dose ===\n")
pred_A_max <- predict(model_A, E0 = placebo_rate, exact.doses = max_rank_doses)
rank_A_max <- rank(pred_A_max, lower_better = FALSE)
print(rank_A_max)
summary(rank_A_max)

cat("\n=== Rankings from Model B at max observed dose ===\n")
pred_B_max <- predict(model_B, E0 = placebo_rate, exact.doses = max_rank_doses)
rank_B_max <- rank(pred_B_max, lower_better = FALSE)
print(rank_B_max)
summary(rank_B_max)

cat("\n=== Rankings from Ref at max observed dose ===\n")
pred_Ref_max <- predict(model_Ref, E0 = placebo_rate, exact.doses = max_rank_doses)
rank_Ref_max <- rank(pred_Ref_max, lower_better = FALSE)
print(rank_Ref_max)
summary(rank_Ref_max)

# ---- 6b. Lumped NMA rankings (directly from posterior) ----
cat("\n=== Rankings from Lumped NMA ===\n")
lumped_sims <- model_Lumped$jagsresult$BUGSoutput$sims.list$d
# Column 1 is placebo, columns 2-6 are agents
d_agents <- lumped_sims[, 2:6]
colnames(d_agents) <- agent_names

# Posterior mean ranks
rank_mat_lumped <- t(apply(d_agents, 1, function(row) base::rank(-row)))
lumped_mean_ranks <- colMeans(rank_mat_lumped)
lumped_median_d <- apply(d_agents, 2, median)

cat("Lumped NMA posterior median treatment effects:\n")
print(round(lumped_median_d, 3))
cat("\nLumped NMA mean ranks:\n")
print(round(lumped_mean_ranks, 2))

# ---- 6c. Rankings at 1x equivalent dose (secondary comparison) ----
cat("\n=== Rankings from Model A at 1x equivalent dose ===\n")
pred_A_rank <- predict(model_A, E0 = placebo_rate, exact.doses = rank_doses)
rank_A <- rank(pred_A_rank, lower_better = FALSE)
print(rank_A)
summary(rank_A)

cat("\n=== Rankings from Model B at 1x equivalent dose ===\n")
pred_B_rank <- predict(model_B, E0 = placebo_rate, exact.doses = rank_doses)
rank_B <- rank(pred_B_rank, lower_better = FALSE)
print(rank_B)
summary(rank_B)

cat("\n=== Rankings from Ref at 1x equivalent dose ===\n")
pred_Ref_rank <- predict(model_Ref, E0 = placebo_rate, exact.doses = rank_doses)
rank_Ref <- rank(pred_Ref_rank, lower_better = FALSE)
print(rank_Ref)
summary(rank_Ref)

# Rankings at 0.5x and 2x
cat("\n=== Rankings at 0.5x equivalent dose (Model A) ===\n")
pred_A_half <- predict(model_A, E0 = placebo_rate, exact.doses = rank_doses_half)
rank_A_half <- rank(pred_A_half, lower_better = FALSE)
print(rank_A_half)

cat("\n=== Rankings at 2x equivalent dose (Model A) ===\n")
pred_A_double <- predict(model_A, E0 = placebo_rate, exact.doses = rank_doses_double)
rank_A_double <- rank(pred_A_double, lower_better = FALSE)
print(rank_A_double)

cat("\n=== Rankings at 0.5x equivalent dose (Model B) ===\n")
pred_B_half <- predict(model_B, E0 = placebo_rate, exact.doses = rank_doses_half)
rank_B_half <- rank(pred_B_half, lower_better = FALSE)
print(rank_B_half)

cat("\n=== Rankings at 2x equivalent dose (Model B) ===\n")
pred_B_double <- predict(model_B, E0 = placebo_rate, exact.doses = rank_doses_double)
rank_B_double <- rank(pred_B_double, lower_better = FALSE)
print(rank_B_double)

# ---- 6d. Compile SUCRA comparison across models ----
compute_sucra_from_summary <- function(rank_obj, model_name) {
  summ <- summary(rank_obj)$Predictions
  
  treat_names  <- as.character(summ$rank.param)
  mean_ranks   <- summ$mean
  median_ranks <- summ$`50%`
  
  K <- nrow(summ)
  sucra <- sapply(1:K, function(i) (K - mean_ranks[i]) / (K - 1))
  
  df <- data.frame(
    agent       = treat_names,
    median_rank = median_ranks,
    mean_rank   = mean_ranks,
    sucra       = sucra,
    model       = model_name,
    stringsAsFactors = FALSE
  )
  
  df <- df[!grepl("placebo", df$agent, ignore.case = TRUE), ]
  return(df)
}

# SUCRA from max-dose rankings
ranking_results <- bind_rows(
  compute_sucra_from_summary(rank_A_max, "Model A"),
  compute_sucra_from_summary(rank_B_max, "Model B"),
  compute_sucra_from_summary(rank_Ref_max, "Ref: Linear")
)

# Add Lumped SUCRA manually
n_agents <- length(agent_names)
lumped_sucra <- (n_agents - lumped_mean_ranks) / (n_agents - 1)
lumped_sucra_df <- data.frame(
  agent       = paste0(agent_names, "_", unlist(max_rank_doses)),
  median_rank = apply(rank_mat_lumped, 2, median),
  mean_rank   = lumped_mean_ranks,
  sucra       = lumped_sucra,
  model       = "Lumped NMA",
  stringsAsFactors = FALSE
)
ranking_results <- bind_rows(ranking_results, lumped_sucra_df)

sucra_wide <- ranking_results %>%
  dplyr::select(agent, model, sucra) %>%
  pivot_wider(names_from = model, values_from = sucra) %>%
  arrange(desc(`Model A`))

cat("\n=== SUCRA Comparison Across Models (max observed dose) ===\n")
print(sucra_wide, digits = 3)

# ---- 6e. Ranking comparison plot ----
ranking_results$agent_clean <- ranking_results$agent %>%
  str_replace_all("_[0-9.]+$", "") %>%
  str_replace_all("_", " ") %>%
  str_to_title()

p_sucra <- ggplot(ranking_results, aes(x = reorder(agent_clean, sucra), y = sucra, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), name = "SUCRA") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
  labs(title = "Agent Rankings at Max Observed Equivalent Dose",
       subtitle = "SUCRA values across model specifications",
       x = NULL, fill = "Model") +
  theme_minimal(base_size = 14)

print(p_sucra)
ggsave("figures/ssri/Figure 4.7.jpeg", p_sucra, width = 8, height = 5, dpi = 300)

# ---- 6f. Save rank plots ----
jpeg("figures/ssri/Chapter 4 Rank A.jpeg", width = 8, height = 6, units = "in", res = 300)
plot(rank_A_max)
dev.off()

jpeg("figures/ssri/Chapter 4 Rank B.jpeg", width = 8, height = 6, units = "in", res = 300)
plot(rank_B_max)
dev.off()

jpeg("figures/ssri/Chapter 4 Rank C.jpeg", width = 8, height = 6, units = "in", res = 300)
plot(rank_Ref_max)
dev.off()

# ==============================================================================
# 7. Prior Sensitivity Analysis
# ==============================================================================
# 
# Systematically vary the prior on tau (between-study heterogeneity SD)
# across the specifications in Table 4 of the methods.
#
# In MBNMAdose, the heterogeneity SD is named "sd" in the priors list.
# JAGS parameterization notes:
#   - dnorm(mu, tau) where tau = precision = 1/variance
#   - Half-Normal(0, sigma): dnorm(0, 1/sigma^2) T(0,)
#   - Half-Cauchy(0, s):     dt(0, 1/s^2, 1) T(0,)
#   - Log-Normal on log(tau^2): requires custom JAGS code
#     Turner et al: log(tau^2) ~ N(-2.34, 1.62^2)
#     => tau^2 ~ dlnorm(-2.34, 1/1.62^2) which gives prior on variance
#     => For SD: sd = sqrt(tau^2), so we use dlnorm on sd^2
#     Note: JAGS dlnorm(meanlog, preclog) where preclog = 1/sdlog^2
# ==============================================================================

# Define prior specifications
# Each entry: list(label, jags_string_for_sd, family, rationale)
#
# IMPORTANT: MBNMAdose's "sd" parameter IS the between-study heterogeneity SD (tau)

prior_specs <- list(
  
  # Uniform priors
  list(
    label     = "U(0, 2)",
    family    = "Uniform",
    jags_sd   = "dunif(0, 2)",
    rationale = "Moderately constrained"
  ),
  list(
    label     = "U(0, 5)",
    family    = "Uniform",
    jags_sd   = "dunif(0, 5)",
    rationale = "MBNMAdose default"
  ),
  list(
    label     = "U(0, 10)",
    family    = "Uniform",
    jags_sd   = "dunif(0, 10)",
    rationale = "Diffuse"
  ),
  
  # Half-Normal priors: dnorm(0, precision) T(0,)
  # HN(0, sigma) => precision = 1/sigma^2
  list(
    label     = "HN(0, 0.5)",
    family    = "Half-Normal",
    jags_sd   = "dnorm(0, 4) T(0,)",       # 1/0.5^2 = 4
    rationale = "Tightly constrained"
  ),
  list(
    label     = "HN(0, 1)",
    family    = "Half-Normal",
    jags_sd   = "dnorm(0, 1) T(0,)",       # 1/1^2 = 1
    rationale = "Rosenberger et al. (2021)"
  ),
  list(
    label     = "HN(0, 2)",
    family    = "Half-Normal",
    jags_sd   = "dnorm(0, 0.25) T(0,)",    # 1/2^2 = 0.25
    rationale = "Moderately diffuse"
  ),
  
  # Half-Cauchy: dt(0, precision, df=1) T(0,)
  # HC(0, s) => precision = 1/s^2
  list(
    label     = "HC(0, 0.5)",
    family    = "Half-Cauchy",
    jags_sd   = "dt(0, 4, 1) T(0,)",       # 1/0.5^2 = 4
    rationale = "Gelman (2006)"
  ),
  
  # Turner et al. (2012) Log-Normal
  # Prior on log(tau^2) ~ N(-2.34, 1.62^2)
  # In JAGS: tau^2 ~ dlnorm(-2.34, 1/1.62^2)
  # => sd = sqrt(tau^2), so sd^2 ~ dlnorm(-2.34, 1/1.62^2)
  # But MBNMAdose takes priors on sd directly.
  # We can use: sd ~ dlnorm(mu_sd, prec_sd) where
  #   if log(sd^2) = log(var) ~ N(-2.34, 1.62^2)
  #   then log(sd) = log(var)/2 ~ N(-2.34/2, 1.62^2/4)
  #   => log(sd) ~ N(-1.17, 0.81^2)
  #   => sd ~ dlnorm(-1.17, 1/0.81^2)
  #   => precision for log(sd) = 1/0.6561 ≈ 1.524
  list(
    label     = "LN(-2.34, 1.62²)",
    family    = "Log-Normal",
    jags_sd   = "dlnorm(-1.17, 1.524)",    # Prior on sd (not variance)
    rationale = "Turner et al. (2012)"
  )
)

# ---- 7a. Fit Model A under each prior specification ----
cat("\n=== Prior Sensitivity Analysis ===\n")
cat("Fitting Model A (Shared Emax) under", length(prior_specs), "prior specifications...\n\n")

prior_results <- list()

for (i in seq_along(prior_specs)) {
  spec <- prior_specs[[i]]
  cat(sprintf("  [%d/%d] Fitting with tau prior: %s (%s)...\n", 
              i, length(prior_specs), spec$label, spec$family))
  
  # Set the prior for tau
  new_priors <- list(sd = spec$jags_sd)
  
  tryCatch({
    # Fit the model
    model_i <- mbnma.run(
      network  = ssri_net,
      fun      = demax(emax = "rel", ed50 = "common"),
      method   = "random",
      likelihood = "binomial",
      link     = "logit",
      parameters.to.save = c("resdev", "emax", "ed50", "sd", "totresdev"),
      priors   = new_priors,
      n.iter   = n_iter,
      n.chains = n_chains,
      n.burnin = n_burnin,
      n.thin   = n_thin,
      jags.seed  = 123456
    )
    
    # Extract tau posterior
    tau_samples <- model_i$BUGSoutput$sims.list$sd
    if (is.matrix(tau_samples)) tau_samples <- tau_samples[, 1]
    
    # Compute LPML using your pre-defined function
    lpml_i <- compute_lpml(model_i, spec$label)$lpml
    
    # Predictions and rankings
    pred_i <- predict(model_i, E0 = placebo_rate, exact.doses = rank_doses)
    rank_i <- rank(pred_i, lower_better = FALSE)
    
    # SUCRA computation using working function
    sucra_df_i <- compute_sucra_from_summary(rank_i, spec$label)
    sucra_i <- setNames(sucra_df_i$sucra, sucra_df_i$agent)
    
    # Convergence check
    rhat_vals <- model_i$BUGSoutput$summary[, "Rhat"]
    rhat_vals <- rhat_vals[!is.na(rhat_vals)]
    
    # Store results
    prior_results[[i]] <- list(
      label      = spec$label,
      family     = spec$family,
      rationale  = spec$rationale,
      model      = model_i,
      dic        = model_i$BUGSoutput$DIC,
      pD         = model_i$BUGSoutput$pD,
      lpml       = lpml_i,
      tau_median = median(tau_samples),
      tau_mean   = mean(tau_samples),
      tau_lower  = quantile(tau_samples, 0.025),
      tau_upper  = quantile(tau_samples, 0.975),
      tau_sd     = sd(tau_samples),
      sucra      = sucra_i,
      max_rhat   = max(rhat_vals),
      converged  = max(rhat_vals) < 1.05
    )
    
    # Print progress
    cat(sprintf("    tau = %.3f (%.3f, %.3f) | DIC = %.1f | LPML = %.1f | Rhat_max = %.3f\n",
                prior_results[[i]]$tau_median,
                prior_results[[i]]$tau_lower,
                prior_results[[i]]$tau_upper,
                prior_results[[i]]$dic,
                prior_results[[i]]$lpml,
                prior_results[[i]]$max_rhat))
    
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
    prior_results[[i]] <<- list(
      label     = spec$label,
      family    = spec$family,
      error     = e$message,
      converged = FALSE
    )
  })
}


# ---- 7b. Prior Sensitivity Summary Table ----
prior_summary <- do.call(rbind, lapply(prior_results, function(res) {
  if (!is.null(res$error)) {
    return(data.frame(
      Prior     = res$label,
      Family    = res$family,
      tau_med   = NA,
      tau_95CI  = NA,
      DIC       = NA,
      LPML      = NA,
      Converged = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    Prior     = res$label,
    Family    = res$family,
    tau_med   = round(res$tau_median, 3),
    tau_95CI  = sprintf("(%.3f, %.3f)", res$tau_lower, res$tau_upper),
    DIC       = round(res$dic, 1),
    LPML      = round(res$lpml, 1),
    Converged = res$converged,
    stringsAsFactors = FALSE
  )
}))

cat("\n=== Prior Sensitivity: Summary Table ===\n")
print(prior_summary)


# ---- 7c. Prior Sensitivity: SUCRA Stability ----
# Check whether agent rankings change across priors
sucra_by_prior <- do.call(rbind, lapply(prior_results, function(res) {
  if (is.null(res$sucra)) return(NULL)
  data.frame(
    prior = res$label,
    agent = names(res$sucra),
    sucra = as.numeric(res$sucra),
    stringsAsFactors = FALSE
  )
})) %>%
  filter(!grepl("placebo|Placebo", agent, ignore.case = TRUE))

cat("\n=== SUCRA Stability Across Priors ===\n")
sucra_prior_wide <- sucra_by_prior %>%
  pivot_wider(names_from = prior, values_from = sucra)
print(sucra_prior_wide, digits = 3)


# ---- 7d. Prior sensitivity plots ----

# Plot 1: Tau posterior across priors
tau_posterior_df <- do.call(rbind, lapply(prior_results, function(res) {
  if (is.null(res$model)) return(NULL)
  tau_samps <- res$model$BUGSoutput$sims.list$sd
  if (is.matrix(tau_samps)) tau_samps <- tau_samps[, 1]
  data.frame(
    prior = res$label,
    tau   = tau_samps,
    stringsAsFactors = FALSE
  )
}))

p_tau_prior <- ggplot(tau_posterior_df, aes(x = tau, fill = prior, color = prior)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 2), name = expression(tau)) +
  labs(title = "Posterior Distribution of Between-Study Heterogeneity",
       subtitle = "Across prior specifications for Model A (Shared Emax)",
       y = "Density", fill = "Prior on τ", color = "Prior on τ") +
  theme(legend.position = "right")

print(p_tau_prior)

ggsave("figures/ssri/Figure 4.10.png", p_tau_prior, 
       width = 9, height = 5, dpi = 300)


# Plot 2: SUCRA stability
p_sucra_prior <- ggplot(sucra_by_prior, aes(x = prior, y = sucra, color = agent, group = agent)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 1), name = "SUCRA") +
  labs(title = "Ranking Stability Across Prior Specifications",
       subtitle = "SUCRA at 1x equivalent dose, Model A",
       x = "Prior on τ", color = "Agent") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_sucra_prior)

ggsave("figures/ssri/Appendix P.png", p_sucra_prior, 
       width = 10, height = 5, dpi = 300)


# Plot 3: DIC and LPML across priors
fit_by_prior <- do.call(rbind, lapply(prior_results, function(res) {
  if (is.null(res$dic)) return(NULL)
  data.frame(
    prior = res$label,
    DIC   = res$dic,
    LPML  = res$lpml,
    stringsAsFactors = FALSE
  )
}))

p_dic_prior <- ggplot(fit_by_prior, aes(x = prior)) +
  geom_point(aes(y = DIC), color = "#E41A1C", size = 3) +
  geom_line(aes(y = DIC, group = 1), color = "#E41A1C", linewidth = 0.6) +
  labs(title = "Model Fit Across Prior Specifications",
       x = "Prior on τ", y = "DIC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_lpml_prior <- ggplot(fit_by_prior, aes(x = prior)) +
  geom_point(aes(y = LPML), color = "#377EB8", size = 3) +
  geom_line(aes(y = LPML, group = 1), color = "#377EB8", linewidth = 0.6) +
  labs(x = "Prior on τ", y = "LPML") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_fit_combined <- p_dic_prior / p_lpml_prior
print(p_fit_combined)

ggsave("figures/Appendix Q.png", p_fit_combined, 
       width = 9, height = 8, dpi = 300)


# ==============================================================================
# 8. Summary Output for Results Section
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  RESULTS SUMMARY — Chapter 4 SSRI Empirical Analysis\n")
cat("================================================================\n\n")

cat("--- 4.3.1.2 Model Comparison ---\n")
print(comparison_table)

cat("\n--- 4.3.1.3 Dose-Response Estimation ---\n")
cat("(See dose-response curve plots saved to working directory)\n")

cat("\n--- 4.3.1.4 Agent Rankings ---\n")
cat("SUCRA at 1x equivalent dose:\n")
print(sucra_wide, digits = 3)

cat("\n--- 4.3.1.5 Prior Sensitivity ---\n")
cat("Tau posteriors across priors:\n")
print(prior_summary)

cat("\nSUCRA stability across priors:\n")
print(sucra_prior_wide, digits = 3)

# Quantify sensitivity: range of SUCRA values per agent across priors
if (nrow(sucra_by_prior) > 0) {
  sucra_range <- sucra_by_prior %>%
    group_by(agent) %>%
    summarise(
      min_sucra  = min(sucra),
      max_sucra  = max(sucra),
      range      = max(sucra) - min(sucra),
      .groups    = "drop"
    ) %>%
    arrange(desc(range))
  
  cat("\nSUCRA range across priors (max - min):\n")
  print(sucra_range, digits = 3)
  
  cat(sprintf("\nMaximum SUCRA variation across priors: %.3f (%s)\n",
              max(sucra_range$range), sucra_range$agent[which.max(sucra_range$range)]))
}

cat("\n================================================================\n")
cat("  Analysis complete. Figures saved to working directory.\n")
cat("================================================================\n")

# ==============================================================================
# 9. Save Results
# ==============================================================================

# Create results directory
dir.create("results/ssri", recursive = TRUE, showWarnings = FALSE)

# ---- Model comparison table (DIC, pD, LPML for all 4 models) ----
write.csv(comparison_table, "results/ssri/model_comparison.csv", row.names = FALSE)

# ---- Convergence diagnostics for all models ----
convergence_df <- data.frame(
  Model = c("Model A", "Model B", "Ref: Linear", "Lumped NMA"),
  Max_Rhat = c(
    max(model_A$BUGSoutput$summary[, "Rhat"], na.rm = TRUE),
    max(model_B$BUGSoutput$summary[, "Rhat"], na.rm = TRUE),
    max(model_Ref$BUGSoutput$summary[, "Rhat"], na.rm = TRUE),
    max(model_Lumped$jagsresult$BUGSoutput$summary[, "Rhat"], na.rm = TRUE)
  ),
  Min_neff = c(
    min(model_A$BUGSoutput$summary[, "n.eff"][model_A$BUGSoutput$summary[, "n.eff"] > 0], na.rm = TRUE),
    min(model_B$BUGSoutput$summary[, "n.eff"][model_B$BUGSoutput$summary[, "n.eff"] > 0], na.rm = TRUE),
    min(model_Ref$BUGSoutput$summary[, "n.eff"][model_Ref$BUGSoutput$summary[, "n.eff"] > 0], na.rm = TRUE),
    min(model_Lumped$jagsresult$BUGSoutput$summary[, "n.eff"][model_Lumped$jagsresult$BUGSoutput$summary[, "n.eff"] > 0], na.rm = TRUE)
  )
)
write.csv(convergence_df, "results/ssri/convergence_diagnostics.csv", row.names = FALSE)

# ---- SUCRA comparison across models at max observed dose ----
write.csv(sucra_wide, "results/ssri/sucra_comparison_max_dose.csv", row.names = FALSE)

# ---- Full ranking results (mean rank, median rank, SUCRA by model) ----
write.csv(ranking_results, "results/ssri/ranking_results_full.csv", row.names = FALSE)

# ---- Max observed equivalent doses used for ranking ----
max_dose_df <- data.frame(
  agent = names(max_rank_doses),
  max_equiv_dose = unlist(max_rank_doses)
)
write.csv(max_dose_df, "results/ssri/max_equiv_doses.csv", row.names = FALSE)

# ---- Lumped NMA treatment effects and ranks ----
lumped_results <- data.frame(
  agent = agent_names,
  d_median = lumped_median_d,
  d_lower = apply(d_agents, 2, quantile, 0.025),
  d_upper = apply(d_agents, 2, quantile, 0.975),
  mean_rank = lumped_mean_ranks,
  sucra = lumped_sucra
)
write.csv(lumped_results, "results/ssri/lumped_nma_results.csv", row.names = FALSE)

# ---- Dose-response parameter estimates (Model B posteriors) ----
# These are the agent-specific Emax/ED50 that calibrate the simulation DGP
dr_params <- data.frame(
  agent = agent_names,
  emax_median = apply(model_B$BUGSoutput$sims.list$emax, 2, median),
  emax_lower  = apply(model_B$BUGSoutput$sims.list$emax, 2, quantile, 0.025),
  emax_upper  = apply(model_B$BUGSoutput$sims.list$emax, 2, quantile, 0.975),
  ed50_median = apply(model_B$BUGSoutput$sims.list$ed50, 2, median),
  ed50_lower  = apply(model_B$BUGSoutput$sims.list$ed50, 2, quantile, 0.025),
  ed50_upper  = apply(model_B$BUGSoutput$sims.list$ed50, 2, quantile, 0.975)
)
write.csv(dr_params, "results/ssri/model_B_dr_parameters.csv", row.names = FALSE)

# ---- Heterogeneity (tau) posterior from each model ----
tau_summary <- data.frame(
  Model = c("Model A", "Model B", "Ref: Linear", "Lumped NMA"),
  tau_median = c(
    median(model_A$BUGSoutput$sims.list$sd),
    median(model_B$BUGSoutput$sims.list$sd),
    median(model_Ref$BUGSoutput$sims.list$sd),
    median(model_Lumped$jagsresult$BUGSoutput$sims.list$sd)
  ),
  tau_lower = c(
    quantile(model_A$BUGSoutput$sims.list$sd, 0.025),
    quantile(model_B$BUGSoutput$sims.list$sd, 0.025),
    quantile(model_Ref$BUGSoutput$sims.list$sd, 0.025),
    quantile(model_Lumped$jagsresult$BUGSoutput$sims.list$sd, 0.025)
  ),
  tau_upper = c(
    quantile(model_A$BUGSoutput$sims.list$sd, 0.975),
    quantile(model_B$BUGSoutput$sims.list$sd, 0.975),
    quantile(model_Ref$BUGSoutput$sims.list$sd, 0.975),
    quantile(model_Lumped$jagsresult$BUGSoutput$sims.list$sd, 0.975)
  )
)
write.csv(tau_summary, "results/ssri/heterogeneity_tau.csv", row.names = FALSE)

# ---- Prior sensitivity results ----
write.csv(prior_summary, "results/ssri/prior_sensitivity_summary.csv", row.names = FALSE)
write.csv(sucra_by_prior, "results/ssri/prior_sensitivity_sucra.csv", row.names = FALSE)

# ---- Full model objects (for reproducibility / further analysis) ----
save(model_A, model_B, model_Ref, model_Lumped,
     comparison_table, ranking_results, sucra_wide,
     dr_params, tau_summary, prior_results,
     file = "results/ssri/ssri_all_models.RData")

cat("\n=== All results saved to results/ssri/ ===\n")
cat("  model_comparison.csv        — DIC, pD, LPML for all 4 models\n")
cat("  convergence_diagnostics.csv — Max Rhat and min n.eff per model\n")
cat("  sucra_comparison_max_dose.csv — SUCRA at max observed dose, wide format\n")
cat("  ranking_results_full.csv    — Full ranking details (all models)\n")
cat("  max_equiv_doses.csv         — Max equivalent dose per agent used for ranking\n")
cat("  lumped_nma_results.csv      — Lumped NMA treatment effects and ranks\n")
cat("  model_B_dr_parameters.csv   — Agent-specific Emax/ED50 posteriors (DGP calibration)\n")
cat("  heterogeneity_tau.csv       — Tau posterior summaries per model\n")
cat("  prior_sensitivity_summary.csv — Tau and fit across prior specs\n")
cat("  prior_sensitivity_sucra.csv — SUCRA stability across priors\n")
cat("  ssri_all_models.RData       — Full model objects for reproducibility\n")