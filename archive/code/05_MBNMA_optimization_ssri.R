################################################################################
# Chapter 4: SSRI Empirical Analysis
# Optimizing Model-Based Network Meta-Analysis for Pharmacologic Networks
#
# CORRECTED VERSION — March 2026
# FIX: Unified ranking approach across all models. Rankings are now computed
# from d_hat (posterior samples of relative treatment effects on log-odds
# scale) for ALL models (A, B, Ref, Lumped), ensuring cross-model SUCRA
# comparisons are on the same scale.
#
# This script implements:
#   1. Network setup
#   2. Model fitting (Models A, B, Ref, Lumped)
#   3. Model comparison (DIC, LPML)
#   4. Dose-response estimation and plotting
#   5. Agent rankings (SUCRA) — unified d_hat approach
#   6. Prior sensitivity analysis on tau
#
# Dependencies: MBNMAdose (>= 0.5.0), R2jags, ggplot2, dplyr, tidyr, 
#               kableExtra, patchwork, coda, stringr
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
library(coda)
library(stringr)

set.seed(12345)
setwd("C:/Users/katie/Desktop/Depression-NMA")

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

dose_equiv <- data.frame(
  agent       = c("citalopram", "escitalopram", "fluoxetine", "paroxetine", "sertraline"),
  equiv_mg    = c(20, 9, 20, 17, 49.3),
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

# ---- Model A: Shared Emax ----
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

# ---- Model B: Agent-specific Emax ----
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

# ---- Reference Model: Shared Linear ----
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

# ---- Lumped NMA ----
cat("\n--- Fitting Lumped NMA ---\n")
ssri_lumped <- ssri_equiv
ssri_lumped$dose[ssri_lumped$dose > 0] <- 1

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

# ---- LPML ----
compute_lpml <- function(model, model_name = "") {
  mcmc_array <- model$BUGSoutput$sims.array
  param_names <- dimnames(mcmc_array)[[3]]
  dev_idx <- grep("^resdev\\[", param_names)
  if (length(dev_idx) == 0) {
    cat("Warning: No 'resdev' parameters found for", model_name, "\n")
    return(list(lpml = NA, cpo = NA))
  }
  n_iter_saved <- dim(mcmc_array)[1]
  n_chains_out <- dim(mcmc_array)[2]
  n_arms <- length(dev_idx)
  dev_matrix <- matrix(NA, nrow = n_iter_saved * n_chains_out, ncol = n_arms)
  for (ch in 1:n_chains_out) {
    rows <- ((ch - 1) * n_iter_saved + 1):(ch * n_iter_saved)
    dev_matrix[rows, ] <- mcmc_array[, ch, dev_idx]
  }
  log_cpo <- numeric(n_arms)
  for (j in 1:n_arms) {
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

comparison_table <- data.frame(
  Model = c("A: Shared Emax", "B: Agent-specific Emax", "Ref: Shared Linear", "Lumped NMA"),
  DIC   = round(c(model_A$BUGSoutput$DIC, model_B$BUGSoutput$DIC,
                  model_Ref$BUGSoutput$DIC, model_Lumped$jagsresult$BUGSoutput$DIC), 1),
  pD    = round(c(model_A$BUGSoutput$pD, model_B$BUGSoutput$pD,
                  model_Ref$BUGSoutput$pD, lumped_pD), 1),
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

placebo_rate <- sum(ssri$r[ssri$agent == "placebo"]) / sum(ssri$n[ssri$agent == "placebo"])
cat("\nPooled placebo response rate:", round(placebo_rate, 3), "\n")

E0_logit <- log(placebo_rate / (1 - placebo_rate))

pred_A <- predict(model_A, E0 = placebo_rate, n.doses = 30)
pred_B <- predict(model_B, E0 = placebo_rate, n.doses = 30)
pred_Ref <- predict(model_Ref, E0 = placebo_rate, n.doses = 30)

# Dose-response curve plots
jpeg("figures/ssri/Figure 4.5.jpeg", width = 10, height = 6, units = "in", res = 300)
plot(pred_A) + ggtitle("Model A: Shared Emax")
dev.off()

jpeg("figures/ssri/Figure 4.6.jpeg", width = 10, height = 6, units = "in", res = 300)
plot(pred_B) + ggtitle("Model B: Agent-specific Emax")
dev.off()

jpeg("figures/ssri/Appendix N.jpeg", width = 10, height = 6, units = "in", res = 300)
plot(pred_Ref) + ggtitle("Ref: Shared Linear")
dev.off()

# Fitted values
jpeg("figures/ssri/Figure 4.4.jpeg", width = 10, height = 8, units = "in", res = 300)
fitplot(model_A)
dev.off()

jpeg("figures/ssri/Chapter 4 Fitplot_B.jpeg", width = 10, height = 8, units = "in", res = 300)
fitplot(model_B)
dev.off()

jpeg("figures/ssri/Chapter 4 Fitplot_Ref.jpeg", width = 10, height = 8, units = "in", res = 300)
fitplot(model_Ref)
dev.off()


# ==============================================================================
# 6. Agent Rankings — Unified d_hat Approach
# ==============================================================================
#
# RANKING METHOD (applied uniformly to all models):
#   1. Extract posterior samples of relative treatment effects (d_hat)
#      on the log-odds scale at specified doses
#   2. For each posterior draw, rank agents (rank 1 = largest effect = best)
#   3. Compute mean ranks and SUCRA from the rank matrix
#
# This replaces the previous approach where dose-response models used
# predict() → rank() (absolute probabilities) while Lumped ranked from
# relative effects. The unified approach ensures cross-model comparisons
# are on the same scale.
# ==============================================================================

agent_names <- c("citalopram", "escitalopram", "fluoxetine", "paroxetine", "sertraline")

# ---- Helper: Extract d_hat matrix from any model at given doses ----
extract_d_hat <- function(model, model_type, agent_names, doses,
                          is_nma = FALSE) {
  #' Extract posterior samples of relative treatment effects at specified doses.
  #'
  #' @param model       Fitted model object
  #' @param model_type  One of "emax_shared", "emax_agent", "linear", "lumped"
  #' @param agent_names Character vector of agent names (in column order)
  #' @param doses       Named list or vector: dose per agent to evaluate at
  #' @param is_nma      TRUE if model is from nma.run() (different object structure)
  #' @return Matrix [n_sims x n_agents] of relative treatment effects
  
  if (is_nma) {
    sims <- model$jagsresult$BUGSoutput$sims.list
    n_sims <- model$jagsresult$BUGSoutput$n.sims
  } else {
    sims <- model$BUGSoutput$sims.list
    n_sims <- model$BUGSoutput$n.sims
  }
  
  n_agents <- length(agent_names)
  d_hat <- matrix(NA, nrow = n_sims, ncol = n_agents)
  colnames(d_hat) <- agent_names
  
  dose_vec <- if (is.list(doses)) unlist(doses[agent_names]) else doses
  
  if (model_type == "emax_shared") {
    # Model A: agent-specific emax, shared ed50
    ed50_samps <- if (is.matrix(sims$ed50)) sims$ed50[, 1] else as.numeric(sims$ed50)
    for (a in seq_along(agent_names)) {
      emax_samps <- sims$emax[, a]
      d_hat[, a] <- emax_samps * dose_vec[a] / (ed50_samps + dose_vec[a])
    }
    
  } else if (model_type == "emax_agent") {
    # Model B: agent-specific emax and ed50
    for (a in seq_along(agent_names)) {
      emax_samps <- sims$emax[, a]
      ed50_samps <- sims$ed50[, a]
      d_hat[, a] <- emax_samps * dose_vec[a] / (ed50_samps + dose_vec[a])
    }
    
  } else if (model_type == "linear") {
    # Ref: agent-specific slope
    for (a in seq_along(agent_names)) {
      beta_samps <- sims$beta.1[, a]
      d_hat[, a] <- beta_samps * dose_vec[a]
    }
    
  } else if (model_type == "lumped") {
    # Lumped NMA: direct treatment effects (column 1 = placebo)
    for (a in seq_along(agent_names)) {
      d_hat[, a] <- sims$d[, a + 1]
    }
    
  } else {
    stop("Unknown model_type: ", model_type)
  }
  
  return(d_hat)
}

# ---- Helper: Compute ranks and SUCRA from d_hat ----
compute_ranks_from_d_hat <- function(d_hat, model_name) {
  #' Compute mean ranks and SUCRA from posterior d_hat matrix.
  #' Rank 1 = largest effect = best.
  #'
  #' @param d_hat      Matrix [n_sims x n_agents]
  #' @param model_name Label for the model
  #' @return Data frame with agent, mean_rank, median_rank, sucra, model
  
  agent_names <- colnames(d_hat)
  n_agents <- ncol(d_hat)
  
  # Rank each posterior draw
  rank_mat <- t(apply(d_hat, 1, function(row) base::rank(-row)))
  
  mean_ranks   <- colMeans(rank_mat)
  median_ranks <- apply(rank_mat, 2, median)
  sucra        <- (n_agents - mean_ranks) / (n_agents - 1)
  
  data.frame(
    agent       = agent_names,
    median_rank = median_ranks,
    mean_rank   = mean_ranks,
    sucra       = sucra,
    d_median    = apply(d_hat, 2, median),
    d_lower     = apply(d_hat, 2, quantile, 0.025),
    d_upper     = apply(d_hat, 2, quantile, 0.975),
    model       = model_name,
    stringsAsFactors = FALSE
  )
}


# ---- Max observed dose per agent ----
max_rank_doses <- list(
  "citalopram"   = max(ssri_equiv$dose[ssri_equiv$agent == "citalopram"]),
  "escitalopram" = max(ssri_equiv$dose[ssri_equiv$agent == "escitalopram"]),
  "fluoxetine"   = max(ssri_equiv$dose[ssri_equiv$agent == "fluoxetine"]),
  "paroxetine"   = max(ssri_equiv$dose[ssri_equiv$agent == "paroxetine"]),
  "sertraline"   = max(ssri_equiv$dose[ssri_equiv$agent == "sertraline"])
)

cat("\n=== Max observed equivalent doses for ranking ===\n")
print(data.frame(agent = names(max_rank_doses), max_dose = unlist(max_rank_doses)))

# Fixed dose levels for comparisons
rank_doses        <- setNames(rep(20, 5), agent_names)
rank_doses_half   <- setNames(rep(10, 5), agent_names)
rank_doses_double <- setNames(rep(40, 5), agent_names)


# ---- 6a. Rankings at max observed dose ----
cat("\n=== Rankings at max observed dose (all models, unified d_hat) ===\n")

d_hat_A_max   <- extract_d_hat(model_A,      "emax_shared", agent_names, max_rank_doses)
d_hat_B_max   <- extract_d_hat(model_B,      "emax_agent",  agent_names, max_rank_doses)
d_hat_Ref_max <- extract_d_hat(model_Ref,    "linear",      agent_names, max_rank_doses)
d_hat_Lmp_max <- extract_d_hat(model_Lumped, "lumped",       agent_names, max_rank_doses,
                               is_nma = TRUE)

ranks_A_max   <- compute_ranks_from_d_hat(d_hat_A_max,   "Model A")
ranks_B_max   <- compute_ranks_from_d_hat(d_hat_B_max,   "Model B")
ranks_Ref_max <- compute_ranks_from_d_hat(d_hat_Ref_max, "Ref: Linear")
ranks_Lmp_max <- compute_ranks_from_d_hat(d_hat_Lmp_max, "Lumped NMA")

ranking_results_max <- bind_rows(ranks_A_max, ranks_B_max, ranks_Ref_max, ranks_Lmp_max)

cat("\nModel A:\n")
print(ranks_A_max[order(ranks_A_max$mean_rank), c("agent", "mean_rank", "sucra", "d_median")])
cat("\nModel B:\n")
print(ranks_B_max[order(ranks_B_max$mean_rank), c("agent", "mean_rank", "sucra", "d_median")])
cat("\nRef: Linear:\n")
print(ranks_Ref_max[order(ranks_Ref_max$mean_rank), c("agent", "mean_rank", "sucra", "d_median")])
cat("\nLumped NMA:\n")
print(ranks_Lmp_max[order(ranks_Lmp_max$mean_rank), c("agent", "mean_rank", "sucra", "d_median")])


# ---- 6b. Rankings at 1x equivalent dose ----
cat("\n=== Rankings at 1x equivalent dose (20mg FE) ===\n")

d_hat_A_1x   <- extract_d_hat(model_A,   "emax_shared", agent_names, rank_doses)
d_hat_B_1x   <- extract_d_hat(model_B,   "emax_agent",  agent_names, rank_doses)
d_hat_Ref_1x <- extract_d_hat(model_Ref, "linear",      agent_names, rank_doses)

ranks_A_1x   <- compute_ranks_from_d_hat(d_hat_A_1x,   "Model A")
ranks_B_1x   <- compute_ranks_from_d_hat(d_hat_B_1x,   "Model B")
ranks_Ref_1x <- compute_ranks_from_d_hat(d_hat_Ref_1x, "Ref: Linear")


# ---- 6c. Rankings at 0.5x and 2x (dose-level stability) ----
cat("\n=== Ranking stability across dose levels ===\n")

d_hat_A_half   <- extract_d_hat(model_A, "emax_shared", agent_names, rank_doses_half)
d_hat_A_double <- extract_d_hat(model_A, "emax_shared", agent_names, rank_doses_double)
d_hat_B_half   <- extract_d_hat(model_B, "emax_agent",  agent_names, rank_doses_half)
d_hat_B_double <- extract_d_hat(model_B, "emax_agent",  agent_names, rank_doses_double)

ranks_A_half   <- compute_ranks_from_d_hat(d_hat_A_half,   "Model A 0.5x")
ranks_A_double <- compute_ranks_from_d_hat(d_hat_A_double, "Model A 2x")
ranks_B_half   <- compute_ranks_from_d_hat(d_hat_B_half,   "Model B 0.5x")
ranks_B_double <- compute_ranks_from_d_hat(d_hat_B_double, "Model B 2x")

# Dose-level comparison table for Model B
rank_stability_B <- data.frame(
  agent = agent_names,
  `0.5x` = ranks_B_half$mean_rank,
  `1x`   = ranks_B_1x$mean_rank,
  `2x`   = ranks_B_double$mean_rank,
  Max    = ranks_B_max$mean_rank,
  check.names = FALSE
)
rank_stability_B <- rank_stability_B[order(rank_stability_B$`1x`), ]
cat("\nModel B mean ranks across dose levels:\n")
print(rank_stability_B, row.names = FALSE, digits = 2)

# Same for Model A
rank_stability_A <- data.frame(
  agent = agent_names,
  `0.5x` = ranks_A_half$mean_rank,
  `1x`   = ranks_A_1x$mean_rank,
  `2x`   = ranks_A_double$mean_rank,
  Max    = ranks_A_max$mean_rank,
  check.names = FALSE
)
rank_stability_A <- rank_stability_A[order(rank_stability_A$`1x`), ]
cat("\nModel A mean ranks across dose levels:\n")
print(rank_stability_A, row.names = FALSE, digits = 2)


# ---- 6d. SUCRA comparison table (wide format) ----
sucra_wide <- ranking_results_max %>%
  dplyr::select(agent, model, sucra) %>%
  pivot_wider(names_from = model, values_from = sucra) %>%
  arrange(desc(`Model A`))

cat("\n=== SUCRA Comparison Across Models (max observed dose) ===\n")
print(sucra_wide, digits = 3)


# ---- 6e. Ranking comparison plot ----
# Create Model B ordering explicitly
model_b_ranks <- ranking_results_max %>%
  filter(model == "Model B") %>%
  arrange(desc(sucra))

# Set factor levels: highest SUCRA first (will appear at top with coord_flip)
agent_order <- model_b_ranks$agent_clean

ranking_results_max$agent_clean <- factor(
  ranking_results_max$agent_clean,
  levels = rev(agent_order)  # rev because coord_flip puts last level on top
)

p_sucra <- ggplot(ranking_results_max,
                  aes(x = agent_clean, y = sucra, fill = model)) +
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
# ==============================================================================
# 7. Prior Sensitivity Analysis
# ==============================================================================

prior_specs <- list(
  list(label = "U(0, 2)",           family = "Uniform",     jags_sd = "dunif(0, 2)",
       rationale = "Moderately constrained"),
  list(label = "U(0, 5)",           family = "Uniform",     jags_sd = "dunif(0, 5)",
       rationale = "MBNMAdose default"),
  list(label = "U(0, 10)",          family = "Uniform",     jags_sd = "dunif(0, 10)",
       rationale = "Diffuse"),
  list(label = "HN(0, 0.5)",        family = "Half-Normal", jags_sd = "dnorm(0, 4) T(0,)",
       rationale = "Tightly constrained"),
  list(label = "HN(0, 1)",          family = "Half-Normal", jags_sd = "dnorm(0, 1) T(0,)",
       rationale = "Rosenberger et al. (2021)"),
  list(label = "HN(0, 2)",          family = "Half-Normal", jags_sd = "dnorm(0, 0.25) T(0,)",
       rationale = "Moderately diffuse"),
  list(label = "HC(0, 0.5)",        family = "Half-Cauchy", jags_sd = "dt(0, 4, 1) T(0,)",
       rationale = "Gelman (2006)"),
  list(label = "LN(-2.34, 1.62\u00B2)", family = "Log-Normal",
       jags_sd = "dlnorm(-1.17, 1.524)",
       rationale = "Turner et al. (2012)")
)

cat("\n=== Prior Sensitivity Analysis ===\n")
cat("Fitting Model A under", length(prior_specs), "prior specifications...\n\n")

prior_results <- list()

for (i in seq_along(prior_specs)) {
  spec <- prior_specs[[i]]
  cat(sprintf("  [%d/%d] Fitting with tau prior: %s (%s)...\n", 
              i, length(prior_specs), spec$label, spec$family))
  
  new_priors <- list(sd = spec$jags_sd)
  
  tryCatch({
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
    
    tau_samples <- model_i$BUGSoutput$sims.list$sd
    if (is.matrix(tau_samples)) tau_samples <- tau_samples[, 1]
    
    lpml_i <- compute_lpml(model_i, spec$label)$lpml
    
    # Rankings using unified d_hat approach
    d_hat_i <- extract_d_hat(model_i, "emax_shared", agent_names, rank_doses)
    ranks_i <- compute_ranks_from_d_hat(d_hat_i, spec$label)
    sucra_i <- setNames(ranks_i$sucra, ranks_i$agent)
    
    rhat_vals <- model_i$BUGSoutput$summary[, "Rhat"]
    rhat_vals <- rhat_vals[!is.na(rhat_vals)]
    
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


# ---- Prior sensitivity summary table ----
prior_summary <- do.call(rbind, lapply(prior_results, function(res) {
  if (!is.null(res$error)) {
    return(data.frame(Prior = res$label, Family = res$family,
                      tau_med = NA, tau_95CI = NA, DIC = NA, LPML = NA,
                      Converged = FALSE, stringsAsFactors = FALSE))
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


# ---- SUCRA stability across priors ----
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


# ---- Prior sensitivity plots ----

# Tau posterior across priors
tau_posterior_df <- do.call(rbind, lapply(prior_results, function(res) {
  if (is.null(res$model)) return(NULL)
  tau_samps <- res$model$BUGSoutput$sims.list$sd
  if (is.matrix(tau_samps)) tau_samps <- tau_samps[, 1]
  data.frame(prior = res$label, tau = tau_samps, stringsAsFactors = FALSE)
}))

p_tau_prior <- ggplot(tau_posterior_df, aes(x = tau, fill = prior, color = prior)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 2), name = expression(tau)) +
  labs(title = "Posterior Distribution of Between-Study Heterogeneity",
       subtitle = "Across prior specifications for Model A (Shared Emax)",
       y = "Density", fill = "Prior on \u03C4", color = "Prior on \u03C4") +
  theme(legend.position = "right")

print(p_tau_prior)
ggsave("figures/ssri/ssri_tau_prior.png", p_tau_prior, width = 9, height = 5, dpi = 300)

# SUCRA stability
p_sucra_prior <- ggplot(sucra_by_prior, aes(x = prior, y = sucra, color = agent, group = agent)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 1), name = "SUCRA") +
  labs(title = "Ranking Stability Across Prior Specifications",
       subtitle = "SUCRA at 1x equivalent dose, Model A (unified d_hat ranking)",
       x = "Prior on \u03C4", color = "Agent") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_sucra_prior)
ggsave("figures/ssri/Figure 4.11.png", p_sucra_prior, width = 10, height = 5, dpi = 300)

# DIC and LPML across priors
fit_by_prior <- do.call(rbind, lapply(prior_results, function(res) {
  if (is.null(res$dic)) return(NULL)
  data.frame(prior = res$label, DIC = res$dic, LPML = res$lpml,
             stringsAsFactors = FALSE)
}))

p_dic_prior <- ggplot(fit_by_prior, aes(x = prior)) +
  geom_point(aes(y = DIC), color = "#E41A1C", size = 3) +
  geom_line(aes(y = DIC, group = 1), color = "#E41A1C", linewidth = 0.6) +
  labs(title = "Model Fit Across Prior Specifications",
       x = "Prior on \u03C4", y = "DIC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_lpml_prior <- ggplot(fit_by_prior, aes(x = prior)) +
  geom_point(aes(y = LPML), color = "#377EB8", size = 3) +
  geom_line(aes(y = LPML, group = 1), color = "#377EB8", linewidth = 0.6) +
  labs(x = "Prior on \u03C4", y = "LPML") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_fit_combined <- p_dic_prior / p_lpml_prior
print(p_fit_combined)
ggsave("figures/Appendix Q.png", p_fit_combined, width = 9, height = 8, dpi = 300)


# ==============================================================================
# 8. Summary Output
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  RESULTS SUMMARY — Chapter 4 SSRI Empirical Analysis\n")
cat("================================================================\n\n")

cat("--- Model Comparison ---\n")
print(comparison_table)

cat("\n--- Agent Rankings (SUCRA at max observed dose) ---\n")
print(sucra_wide, digits = 3)

cat("\n--- Prior Sensitivity ---\n")
print(prior_summary)

cat("\nSUCRA stability across priors:\n")
print(sucra_prior_wide, digits = 3)

if (nrow(sucra_by_prior) > 0) {
  sucra_range <- sucra_by_prior %>%
    group_by(agent) %>%
    summarise(min_sucra = min(sucra), max_sucra = max(sucra),
              range = max(sucra) - min(sucra), .groups = "drop") %>%
    arrange(desc(range))
  cat("\nSUCRA range across priors:\n")
  print(sucra_range, digits = 3)
}

cat("\n================================================================\n")
cat("  Analysis complete.\n")
cat("================================================================\n")


# ==============================================================================
# 9. Save Results
# ==============================================================================

dir.create("results/ssri", recursive = TRUE, showWarnings = FALSE)

write.csv(comparison_table, "results/ssri/model_comparison.csv", row.names = FALSE)

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

write.csv(sucra_wide, "results/ssri/sucra_comparison_max_dose.csv", row.names = FALSE)
write.csv(ranking_results_max, "results/ssri/ranking_results_full.csv", row.names = FALSE)

max_dose_df <- data.frame(agent = names(max_rank_doses), max_equiv_dose = unlist(max_rank_doses))
write.csv(max_dose_df, "results/ssri/max_equiv_doses.csv", row.names = FALSE)

# Lumped NMA results
lumped_results <- ranks_Lmp_max[, c("agent", "d_median", "d_lower", "d_upper",
                                    "mean_rank", "sucra")]
write.csv(lumped_results, "results/ssri/lumped_nma_results.csv", row.names = FALSE)

# Model B dose-response parameters
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

# Tau posteriors
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

write.csv(prior_summary, "results/ssri/prior_sensitivity_summary.csv", row.names = FALSE)
write.csv(sucra_by_prior, "results/ssri/prior_sensitivity_sucra.csv", row.names = FALSE)

# Dose-level stability tables
write.csv(rank_stability_A, "results/ssri/rank_stability_model_A.csv", row.names = FALSE)
write.csv(rank_stability_B, "results/ssri/rank_stability_model_B.csv", row.names = FALSE)

save(model_A, model_B, model_Ref, model_Lumped,
     comparison_table, ranking_results_max, sucra_wide,
     dr_params, tau_summary, prior_results,
     extract_d_hat, compute_ranks_from_d_hat,
     file = "results/ssri/ssri_all_models.RData")

cat("\n=== All results saved to results/ssri/ ===\n")

cat("\n=== Session Info ===\n")
writeLines(capture.output(sessionInfo()), "results/ssri/session_info.txt")