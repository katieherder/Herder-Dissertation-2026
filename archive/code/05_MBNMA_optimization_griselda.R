################################################################################
# Chapter 4: GRISELDA Multi-Class Empirical Analysis
# Optimizing Model-Based Network Meta-Analysis for Pharmacologic Networks
#
# CORRECTED VERSION — March 2026
# FIX: Unified ranking approach across all models using d_hat (posterior
# samples of relative treatment effects on log-odds scale).
#
# This script implements:
#   1. Network setup (using pre-cleaned Griselda_clean with dose_FE)
#   2. Model fitting (Models A, B, Ref, Lumped)
#   3. Model comparison (DIC, LPML)
#   4. Dose-response estimation and plotting
#   5. Agent rankings (SUCRA) — unified d_hat approach
#   6. Prior sensitivity analysis on tau
#   7. Leave-one-out cross-validation (LOO)
#
# PREREQUISITE: Griselda_clean must be loaded with columns:
#   studyID, Class, agent, dose, n, r, weeks, dose_FE
#
# Dependencies: MBNMAdose (>= 0.5.0), R2jags, ggplot2, dplyr, tidyr,
#               patchwork, coda, stringr, parallel
################################################################################

# ==============================================================================
# 0. Setup
# ==============================================================================

library(MBNMAdose)
library(R2jags)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(coda)
library(stringr)
library(parallel)

set.seed(12345)
setwd("C:/Users/katie/Desktop/Depression-NMA")

theme_diss <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )
theme_set(theme_diss)

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
# 1. Create Network Object
# ==============================================================================

griselda_mbnma <- Griselda_clean %>%
  mutate(dose = ifelse(agent == "placebo", 0, dose_FE)) %>%
  select(studyID, agent, dose, n, r) %>%
  as.data.frame()

stopifnot(sum(is.na(griselda_mbnma$dose)) == 0)

cat("\n=== Creating GRISELDA network ===\n")
cat("Studies:", n_distinct(griselda_mbnma$studyID), "\n")
cat("Arms:", nrow(griselda_mbnma), "\n")
cat("Agents:", n_distinct(griselda_mbnma$agent[griselda_mbnma$agent != "placebo"]), "\n")

griselda_net <- mbnma.network(griselda_mbnma)

placebo_rate_g <- sum(Griselda_clean$r[Griselda_clean$agent == "placebo"]) /
  sum(Griselda_clean$n[Griselda_clean$agent == "placebo"])
cat("Pooled placebo response rate:", round(placebo_rate_g, 3), "\n")


# ==============================================================================
# 2. Reusable Functions
# ==============================================================================

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
    neff_vals <- bugs$summary[, "n.eff"]
    neff_vals <- neff_vals[!is.na(neff_vals) & neff_vals > 0]
    cat("  Min n.eff:", min(neff_vals), "\n")
    cat("  Median n.eff:", median(neff_vals), "\n")
  }
  invisible(NULL)
}

extract_d_hat <- function(model, model_type, agent_names, doses,
                          is_nma = FALSE) {
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
    ed50_samps <- if (is.matrix(sims$ed50)) sims$ed50[, 1] else as.numeric(sims$ed50)
    for (a in seq_along(agent_names)) {
      d_hat[, a] <- sims$emax[, a] * dose_vec[a] / (ed50_samps + dose_vec[a])
    }
  } else if (model_type == "emax_agent") {
    for (a in seq_along(agent_names)) {
      d_hat[, a] <- sims$emax[, a] * dose_vec[a] / (sims$ed50[, a] + dose_vec[a])
    }
  } else if (model_type == "linear") {
    for (a in seq_along(agent_names)) {
      d_hat[, a] <- sims$beta.1[, a] * dose_vec[a]
    }
  } else if (model_type == "lumped") {
    for (a in seq_along(agent_names)) {
      d_hat[, a] <- sims$d[, a + 1]
    }
  } else {
    stop("Unknown model_type: ", model_type)
  }
  return(d_hat)
}

compute_ranks_from_d_hat <- function(d_hat, model_name) {
  agent_names <- colnames(d_hat)
  n_agents <- ncol(d_hat)
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


# ==============================================================================
# 3. Model Fitting
# ==============================================================================

active_agents <- sort(unique(griselda_mbnma$agent[griselda_mbnma$agent != "placebo"]))
cat("\nActive agents (", length(active_agents), "):", paste(active_agents, collapse = ", "), "\n")

# ---- Model A: Shared Emax ----
cat("\n--- Fitting Model A: Shared Emax ---\n")
model_A_g <- mbnma.run(
  network  = griselda_net,
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
print(model_A_g)
check_convergence(model_A_g, "Model A (Shared Emax)")

# ---- Model B: Agent-specific Emax ----
cat("\n--- Fitting Model B: Agent-specific Emax (dmulti) ---\n")
dmulti_list <- setNames(lapply(active_agents, function(a) demax()), active_agents)

model_B_g <- mbnma.run(
  network  = griselda_net,
  fun      = dmulti(dmulti_list),
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
print(model_B_g)
check_convergence(model_B_g, "Model B (Agent-specific Emax)")

# ---- Reference Model: Shared Linear ----
cat("\n--- Fitting Reference Model: Shared Linear ---\n")
model_Ref_g <- mbnma.run(
  network  = griselda_net,
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
print(model_Ref_g)
check_convergence(model_Ref_g, "Ref (Shared Linear)")

# ---- Lumped NMA ----
cat("\n--- Fitting Lumped NMA ---\n")
griselda_lumped <- griselda_mbnma
griselda_lumped$dose[griselda_lumped$dose > 0] <- 1

griselda_lumped <- griselda_lumped %>%
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

griselda_net_lumped <- mbnma.network(griselda_lumped)

model_Lumped_g <- nma.run(
  network    = griselda_net_lumped,
  method     = "random",
  likelihood = "binomial",
  link       = "logit",
  n.iter     = n_iter,
  n.chains   = n_chains,
  n.burnin   = n_burnin,
  n.thin     = n_thin,
  jags.seed  = 123456
)
check_convergence(model_Lumped_g, "Lumped NMA", is_nma = TRUE)


# ==============================================================================
# 4. Model Comparison: DIC and LPML
# ==============================================================================

lpml_A_g   <- compute_lpml(model_A_g, "Model A")
lpml_B_g   <- compute_lpml(model_B_g, "Model B")
lpml_Ref_g <- compute_lpml(model_Ref_g, "Ref")

lumped_pD_g <- as.numeric(model_Lumped_g$jagsresult$BUGSoutput$DIC)[1] - 
  model_Lumped_g$jagsresult$BUGSoutput$mean$deviance

comparison_table_g <- data.frame(
  Model = c("A: Shared Emax", "B: Agent-specific Emax", "Ref: Shared Linear", "Lumped NMA"),
  DIC   = round(c(model_A_g$BUGSoutput$DIC, model_B_g$BUGSoutput$DIC,
                  model_Ref_g$BUGSoutput$DIC,
                  as.numeric(model_Lumped_g$jagsresult$BUGSoutput$DIC)[1]), 1),
  pD    = round(c(model_A_g$BUGSoutput$pD, model_B_g$BUGSoutput$pD,
                  model_Ref_g$BUGSoutput$pD, lumped_pD_g), 1),
  LPML  = round(c(lpml_A_g$lpml, lpml_B_g$lpml, lpml_Ref_g$lpml, NA), 1)
)
comparison_table_g$Dbar <- comparison_table_g$DIC - comparison_table_g$pD

cat("\n=== GRISELDA Model Comparison ===\n")
print(comparison_table_g)

jpeg("figures/griselda/GRISELDA_DevDev_A_vs_Ref.jpeg", width = 8, height = 6, units = "in", res = 300)
devdev(model_A_g, model_Ref_g, dev.type = "resdev")
dev.off()


# ==============================================================================
# 5. Dose-Response Estimation and Visualization
# ==============================================================================

pred_A_g   <- predict(model_A_g, E0 = placebo_rate_g, n.doses = 30)
pred_B_g   <- predict(model_B_g, E0 = placebo_rate_g, n.doses = 30)
pred_Ref_g <- predict(model_Ref_g, E0 = placebo_rate_g, n.doses = 30)

jpeg("figures/griselda/Figure 4.9.jpeg", width = 14, height = 10, units = "in", res = 300)
plot(pred_A_g) + ggtitle("GRISELDA Model A: Shared Emax")
dev.off()

jpeg("figures/griselda/GRISELDA_DR_ModelB.jpeg", width = 14, height = 10, units = "in", res = 300)
plot(pred_B_g) + ggtitle("GRISELDA Model B: Agent-specific Emax")
dev.off()

jpeg("figures/griselda/GRISELDA_DR_Ref.jpeg", width = 14, height = 10, units = "in", res = 300)
plot(pred_Ref_g) + ggtitle("GRISELDA Ref: Shared Linear")
dev.off()

jpeg("figures/griselda/Figure 4.8.jpeg", width = 14, height = 12, units = "in", res = 300)
fitplot(model_A_g)
dev.off()

jpeg("figures/griselda/GRISELDA_Fitplot_B.jpeg", width = 14, height = 12, units = "in", res = 300)
fitplot(model_B_g)
dev.off()

jpeg("figures/griselda/GRISELDA_Fitplot_Ref.jpeg", width = 14, height = 12, units = "in", res = 300)
fitplot(model_Ref_g)
dev.off()


# ==============================================================================
# 6. Agent Rankings — Unified d_hat Approach
# ==============================================================================

# ---- Max observed dose per agent ----
max_rank_doses_g <- griselda_mbnma %>%
  filter(agent != "placebo") %>%
  group_by(agent) %>%
  summarise(max_dose = max(dose), .groups = "drop")

cat("\n=== Max observed FE doses for ranking ===\n")
print(max_rank_doses_g, n = Inf)

max_rank_doses_g_list <- as.list(setNames(max_rank_doses_g$max_dose, max_rank_doses_g$agent))

rank_doses_g <- setNames(as.list(rep(20, length(active_agents))), active_agents)

# ---- 6a. Rankings at max observed dose (unified d_hat) ----
cat("\n=== Rankings at max observed dose (all models, unified d_hat) ===\n")

d_hat_A_max_g   <- extract_d_hat(model_A_g,   "emax_shared", active_agents, max_rank_doses_g_list)
d_hat_B_max_g   <- extract_d_hat(model_B_g,   "emax_agent",  active_agents, max_rank_doses_g_list)
d_hat_Ref_max_g <- extract_d_hat(model_Ref_g, "linear",      active_agents, max_rank_doses_g_list)

lumped_trt_labs <- model_Lumped_g$trt.labs
lumped_active_labs <- lumped_trt_labs[!grepl("Placebo|placebo", lumped_trt_labs)]
lumped_agent_names <- gsub("_1$", "", lumped_active_labs)
d_hat_Lmp_max_g <- extract_d_hat(model_Lumped_g, "lumped", lumped_agent_names,
                                 max_rank_doses_g_list, is_nma = TRUE)

ranks_A_max_g   <- compute_ranks_from_d_hat(d_hat_A_max_g,   "Model A")
ranks_B_max_g   <- compute_ranks_from_d_hat(d_hat_B_max_g,   "Model B")
ranks_Ref_max_g <- compute_ranks_from_d_hat(d_hat_Ref_max_g, "Ref: Linear")
ranks_Lmp_max_g <- compute_ranks_from_d_hat(d_hat_Lmp_max_g, "Lumped NMA")

ranking_results_g <- bind_rows(ranks_A_max_g, ranks_B_max_g,
                               ranks_Ref_max_g, ranks_Lmp_max_g)

# ---- 6b. Rankings at 1x equivalent dose (20mg FE) ----
cat("\n=== Rankings at 20mg FE (all models, unified d_hat) ===\n")

d_hat_A_20_g   <- extract_d_hat(model_A_g,   "emax_shared", active_agents, rank_doses_g)
d_hat_B_20_g   <- extract_d_hat(model_B_g,   "emax_agent",  active_agents, rank_doses_g)
d_hat_Ref_20_g <- extract_d_hat(model_Ref_g, "linear",      active_agents, rank_doses_g)

ranks_A_20_g   <- compute_ranks_from_d_hat(d_hat_A_20_g,   "Model A")
ranks_B_20_g   <- compute_ranks_from_d_hat(d_hat_B_20_g,   "Model B")
ranks_Ref_20_g <- compute_ranks_from_d_hat(d_hat_Ref_20_g, "Ref: Linear")

# ---- 6c. SUCRA comparison (wide format) ----
class_lookup <- Griselda_clean %>%
  filter(agent != "placebo") %>%
  distinct(agent, Class)

ranking_results_g <- ranking_results_g %>%
  left_join(class_lookup, by = "agent")

sucra_wide_g <- ranking_results_g %>%
  dplyr::select(agent, Class, model, sucra) %>%
  pivot_wider(names_from = model, values_from = sucra) %>%
  arrange(Class, desc(`Model A`))

cat("\n=== SUCRA Comparison Across Models (max observed dose) ===\n")
print(sucra_wide_g, n = Inf, width = Inf)

# ---- 6d. SUCRA plot ----
ranking_results_g$agent_clean <- ranking_results_g$agent %>%
  str_replace_all("_", " ") %>%
  str_to_title()

model_b_order <- ranking_results_g %>%
  filter(model == "Model B") %>%
  arrange(sucra) %>%
  pull(agent_clean)

ranking_results_g$agent_clean <- factor(ranking_results_g$agent_clean,
                                        levels = unique(model_b_order))

p_sucra_g <- ggplot(ranking_results_g,
                    aes(x = agent_clean, y = sucra, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), name = "SUCRA") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
  labs(title = "GRISELDA: Agent Rankings at Max Observed FE Dose",
       subtitle = "SUCRA values across model specifications",
       x = NULL, fill = "Model") +
  theme_minimal(base_size = 12)

print(p_sucra_g)
ggsave("figures/griselda/Figure 4.10.jpeg", p_sucra_g, width = 10, height = 12, dpi = 300)

# ---- 6e. Model B rank stability across dose levels ----
rank_doses_10_g <- setNames(as.list(rep(10, length(active_agents))), active_agents)
rank_doses_40_g <- setNames(as.list(rep(40, length(active_agents))), active_agents)

d_hat_B_10_g <- extract_d_hat(model_B_g, "emax_agent", active_agents, rank_doses_10_g)
d_hat_B_40_g <- extract_d_hat(model_B_g, "emax_agent", active_agents, rank_doses_40_g)

ranks_B_10_g <- compute_ranks_from_d_hat(d_hat_B_10_g, "Model B 10mg")
ranks_B_40_g <- compute_ranks_from_d_hat(d_hat_B_40_g, "Model B 40mg")

rank_dose_comparison <- data.frame(
  Agent = active_agents,
  `10mg` = round(ranks_B_10_g$mean_rank, 2),
  `20mg` = round(ranks_B_20_g$mean_rank, 2),
  `40mg` = round(ranks_B_40_g$mean_rank, 2),
  Max    = round(ranks_B_max_g$mean_rank, 2),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rank_dose_comparison <- rank_dose_comparison[order(rank_dose_comparison$`20mg`), ]

cat("\n=== Model B Mean Ranks Across Dose Levels ===\n")
print(rank_dose_comparison, row.names = FALSE)


# ==============================================================================
# 7. Prior Sensitivity Analysis
# ==============================================================================

prior_specs <- list(
  list(label = "U(0, 2)",           family = "Uniform",     jags_sd = "dunif(0, 2)"),
  list(label = "U(0, 5)",           family = "Uniform",     jags_sd = "dunif(0, 5)"),
  list(label = "U(0, 10)",          family = "Uniform",     jags_sd = "dunif(0, 10)"),
  list(label = "HN(0, 0.5)",        family = "Half-Normal", jags_sd = "dnorm(0, 4) T(0,)"),
  list(label = "HN(0, 1)",          family = "Half-Normal", jags_sd = "dnorm(0, 1) T(0,)"),
  list(label = "HN(0, 2)",          family = "Half-Normal", jags_sd = "dnorm(0, 0.25) T(0,)"),
  list(label = "HC(0, 0.5)",        family = "Half-Cauchy", jags_sd = "dt(0, 4, 1) T(0,)"),
  list(label = "LN(-2.34, 1.62^2)", family = "Log-Normal",  jags_sd = "dlnorm(-1.17, 1.524)")
)

cat("\n=== GRISELDA Prior Sensitivity Analysis ===\n")
cat("Fitting Model A under", length(prior_specs), "prior specifications...\n\n")

prior_results_g <- list()

for (i in seq_along(prior_specs)) {
  spec <- prior_specs[[i]]
  cat(sprintf("  [%d/%d] Fitting with tau prior: %s (%s)...\n",
              i, length(prior_specs), spec$label, spec$family))
  
  new_priors <- list(sd = spec$jags_sd)
  
  tryCatch({
    model_i <- mbnma.run(
      network  = griselda_net,
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
      jags.seed = 123456
    )
    
    tau_samples <- model_i$BUGSoutput$sims.list$sd
    if (is.matrix(tau_samples)) tau_samples <- tau_samples[, 1]
    
    lpml_i <- compute_lpml(model_i, spec$label)$lpml
    
    # Rankings using unified d_hat approach
    d_hat_i <- extract_d_hat(model_i, "emax_shared", active_agents, rank_doses_g)
    ranks_i <- compute_ranks_from_d_hat(d_hat_i, spec$label)
    sucra_i <- setNames(ranks_i$sucra, ranks_i$agent)
    
    rhat_vals <- model_i$BUGSoutput$summary[, "Rhat"]
    rhat_vals <- rhat_vals[!is.na(rhat_vals)]
    
    prior_results_g[[i]] <- list(
      label      = spec$label,
      family     = spec$family,
      model      = model_i,
      dic        = model_i$BUGSoutput$DIC,
      lpml       = lpml_i,
      tau_median = median(tau_samples),
      tau_lower  = quantile(tau_samples, 0.025),
      tau_upper  = quantile(tau_samples, 0.975),
      sucra      = sucra_i,
      max_rhat   = max(rhat_vals),
      converged  = max(rhat_vals) < 1.05
    )
    
    cat(sprintf("    tau = %.3f (%.3f, %.3f) | DIC = %.1f | LPML = %.1f | Rhat_max = %.3f\n",
                prior_results_g[[i]]$tau_median,
                prior_results_g[[i]]$tau_lower,
                prior_results_g[[i]]$tau_upper,
                prior_results_g[[i]]$dic,
                prior_results_g[[i]]$lpml,
                prior_results_g[[i]]$max_rhat))
    
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
    prior_results_g[[i]] <<- list(label = spec$label, family = spec$family,
                                  error = e$message, converged = FALSE)
  })
}

# ---- Prior sensitivity summary table ----
prior_summary_g <- do.call(rbind, lapply(prior_results_g, function(res) {
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

cat("\n=== GRISELDA Prior Sensitivity: Summary ===\n")
print(prior_summary_g)

# ---- Tau posterior density plot ----
tau_posterior_g <- do.call(rbind, lapply(prior_results_g, function(res) {
  if (is.null(res$model)) return(NULL)
  tau_samps <- res$model$BUGSoutput$sims.list$sd
  if (is.matrix(tau_samps)) tau_samps <- tau_samps[, 1]
  data.frame(prior = res$label, tau = tau_samps, stringsAsFactors = FALSE)
}))

p_tau_g <- ggplot(tau_posterior_g, aes(x = tau, fill = prior, color = prior)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 2), name = expression(tau)) +
  labs(title = "GRISELDA: Posterior of Between-Study Heterogeneity",
       subtitle = "Across prior specifications for Model A",
       y = "Density", fill = "Prior on \u03C4", color = "Prior on \u03C4") +
  theme(legend.position = "right")

ggsave("figures/griselda/GRISELDA_tau_prior_sensitivity.png", p_tau_g,
       width = 9, height = 5, dpi = 300)

# ---- SUCRA stability across priors ----
sucra_by_prior_g <- do.call(rbind, lapply(prior_results_g, function(res) {
  if (is.null(res$sucra)) return(NULL)
  data.frame(prior = res$label, agent = names(res$sucra),
             sucra = as.numeric(res$sucra), stringsAsFactors = FALSE)
})) %>%
  filter(!grepl("placebo", agent, ignore.case = TRUE))

sucra_range_g <- sucra_by_prior_g %>%
  group_by(agent) %>%
  summarise(min_sucra = min(sucra), max_sucra = max(sucra),
            range = max(sucra) - min(sucra), .groups = "drop") %>%
  arrange(desc(range))

cat("\n=== GRISELDA SUCRA range across priors ===\n")
print(sucra_range_g, n = Inf)


# ==============================================================================
# 8. Summary Output
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  RESULTS SUMMARY — Chapter 4 GRISELDA Empirical Analysis\n")
cat("================================================================\n\n")

cat("--- Model Comparison ---\n")
print(comparison_table_g)

cat("\n--- Agent Rankings (SUCRA at max observed dose) ---\n")
print(sucra_wide_g, n = Inf, width = Inf)

cat("\n--- Prior Sensitivity ---\n")
print(prior_summary_g)

cat("\n--- SUCRA Range Across Priors ---\n")
print(sucra_range_g, n = Inf)

cat("\n================================================================\n")
cat("  Analysis complete.\n")
cat("================================================================\n")


# ==============================================================================
# 9. Save Results
# ==============================================================================

dir.create("results/griselda", recursive = TRUE, showWarnings = FALSE)

write.csv(comparison_table_g, "results/griselda/model_comparison.csv", row.names = FALSE)

convergence_df_g <- data.frame(
  Model = c("Model A", "Model B", "Ref: Linear", "Lumped NMA"),
  Max_Rhat = c(
    max(model_A_g$BUGSoutput$summary[, "Rhat"], na.rm = TRUE),
    max(model_B_g$BUGSoutput$summary[, "Rhat"], na.rm = TRUE),
    max(model_Ref_g$BUGSoutput$summary[, "Rhat"], na.rm = TRUE),
    max(model_Lumped_g$jagsresult$BUGSoutput$summary[, "Rhat"], na.rm = TRUE)
  ),
  Min_neff = c(
    min(model_A_g$BUGSoutput$summary[, "n.eff"][model_A_g$BUGSoutput$summary[, "n.eff"] > 0], na.rm = TRUE),
    min(model_B_g$BUGSoutput$summary[, "n.eff"][model_B_g$BUGSoutput$summary[, "n.eff"] > 0], na.rm = TRUE),
    min(model_Ref_g$BUGSoutput$summary[, "n.eff"][model_Ref_g$BUGSoutput$summary[, "n.eff"] > 0], na.rm = TRUE),
    min(model_Lumped_g$jagsresult$BUGSoutput$summary[, "n.eff"][model_Lumped_g$jagsresult$BUGSoutput$summary[, "n.eff"] > 0], na.rm = TRUE)
  )
)
write.csv(convergence_df_g, "results/griselda/convergence_diagnostics.csv", row.names = FALSE)

write.csv(sucra_wide_g, "results/griselda/sucra_comparison_max_dose.csv", row.names = FALSE)
write.csv(ranking_results_g, "results/griselda/ranking_results_full.csv", row.names = FALSE)
write.csv(max_rank_doses_g, "results/griselda/max_equiv_doses.csv", row.names = FALSE)

lumped_results_g <- ranks_Lmp_max_g[, c("agent", "d_median", "d_lower", "d_upper",
                                        "mean_rank", "sucra")]
write.csv(lumped_results_g, "results/griselda/lumped_nma_results.csv", row.names = FALSE)

tau_summary_g <- data.frame(
  Model = c("Model A", "Model B", "Ref: Linear", "Lumped NMA"),
  tau_median = c(
    median(model_A_g$BUGSoutput$sims.list$sd),
    median(model_B_g$BUGSoutput$sims.list$sd),
    median(model_Ref_g$BUGSoutput$sims.list$sd),
    median(model_Lumped_g$jagsresult$BUGSoutput$sims.list$sd)
  ),
  tau_lower = c(
    quantile(model_A_g$BUGSoutput$sims.list$sd, 0.025),
    quantile(model_B_g$BUGSoutput$sims.list$sd, 0.025),
    quantile(model_Ref_g$BUGSoutput$sims.list$sd, 0.025),
    quantile(model_Lumped_g$jagsresult$BUGSoutput$sims.list$sd, 0.025)
  ),
  tau_upper = c(
    quantile(model_A_g$BUGSoutput$sims.list$sd, 0.975),
    quantile(model_B_g$BUGSoutput$sims.list$sd, 0.975),
    quantile(model_Ref_g$BUGSoutput$sims.list$sd, 0.975),
    quantile(model_Lumped_g$jagsresult$BUGSoutput$sims.list$sd, 0.975)
  )
)
write.csv(tau_summary_g, "results/griselda/heterogeneity_tau.csv", row.names = FALSE)

write.csv(prior_summary_g, "results/griselda/prior_sensitivity_summary.csv", row.names = FALSE)
write.csv(sucra_by_prior_g, "results/griselda/prior_sensitivity_sucra.csv", row.names = FALSE)

write.csv(rank_dose_comparison, "results/griselda/rank_dose_comparison_modelB.csv", row.names = FALSE)

save(model_A_g, model_B_g, model_Ref_g, model_Lumped_g,
     comparison_table_g, ranking_results_g, sucra_wide_g,
     tau_summary_g, prior_results_g,
     extract_d_hat, compute_ranks_from_d_hat,
     file = "results/griselda/griselda_all_models.RData")

cat("\n=== All results saved to results/griselda/ ===\n")


# ==============================================================================
# 11. LOO Cross-Validation
# ==============================================================================

loo_candidates <- Griselda_clean %>%
  filter(agent != "placebo") %>%
  group_by(studyID) %>%
  filter(n_distinct(agent) > 1) %>%
  summarise(
    n_agents = n_distinct(agent),
    n_classes = n_distinct(Class),
    agents = paste(sort(unique(agent)), collapse = " vs "),
    .groups = "drop"
  )

cat("Total studies with cross-agent comparisons:", nrow(loo_candidates), "\n")
cat("  Of which cross-class:", sum(loo_candidates$n_classes > 1), "\n")

loo_studies <- loo_candidates %>%
  filter(n_classes > 1) %>%
  pull(studyID)

cat("LOO studies selected:", length(loo_studies), "\n\n")

run_loo_single <- function(held_out_study, model_fun, model_name,
                           full_data, n_iter, n_chains, n_burnin, n_thin) {
  reduced_data <- full_data %>% filter(studyID != held_out_study)
  held_out <- full_data %>% filter(studyID == held_out_study)
  
  tryCatch({
    reduced_mbnma <- reduced_data %>%
      mutate(dose = ifelse(agent == "placebo", 0, dose_FE)) %>%
      select(studyID, agent, dose, n, r) %>%
      as.data.frame()
    
    reduced_net <- mbnma.network(reduced_mbnma)
    
    reduced_agents <- unique(reduced_mbnma$agent[reduced_mbnma$agent != "placebo"])
    held_agents <- held_out %>%
      filter(agent != "placebo") %>% pull(agent) %>% unique()
    removed_agents <- setdiff(held_agents, reduced_agents)
    
    model_loo <- mbnma.run(
      network    = reduced_net,
      fun        = model_fun,
      method     = "random",
      likelihood = "binomial",
      link       = "logit",
      parameters.to.save = c("resdev", "sd", "totresdev"),
      n.iter     = n_iter,
      n.chains   = n_chains,
      n.burnin   = n_burnin,
      n.thin     = n_thin,
      jags.seed  = 123456
    )
    
    rhat_vals <- model_loo$BUGSoutput$summary[, "Rhat"]
    rhat_vals <- rhat_vals[!is.na(rhat_vals)]
    
    return(list(
      study = held_out_study, model = model_name,
      dic = model_loo$BUGSoutput$DIC, pD = model_loo$BUGSoutput$pD,
      totresdev = model_loo$BUGSoutput$summary["totresdev", "50%"],
      tau = model_loo$BUGSoutput$summary["sd", "50%"],
      max_rhat = max(rhat_vals), converged = max(rhat_vals) < 1.05,
      removed_agents = removed_agents, n_arms_removed = nrow(held_out),
      error = NULL
    ))
  }, error = function(e) {
    return(list(
      study = held_out_study, model = model_name,
      dic = NA, pD = NA, totresdev = NA, tau = NA,
      max_rhat = NA, converged = FALSE,
      removed_agents = NULL, n_arms_removed = NA, error = e$message
    ))
  })
}

cat("Starting LOO cross-validation (Models A and Ref)...\n")

loo_results <- list()
loo_counter <- 0

for (study_id in loo_studies) {
  loo_counter <- loo_counter + 1
  cat(sprintf("\n[%d/%d] Holding out: %s\n", loo_counter, length(loo_studies), study_id))
  
  cat("  Fitting Model A...")
  t0 <- Sys.time()
  res_A <- run_loo_single(study_id, demax(emax = "rel", ed50 = "common"), "Model A",
                          Griselda_clean, n_iter, n_chains, n_burnin, n_thin)
  t1 <- Sys.time()
  if (is.null(res_A$error)) {
    cat(sprintf(" done (%.1f min, DIC=%.1f, converged=%s)\n",
                difftime(t1, t0, units = "mins"), res_A$dic, res_A$converged))
  } else {
    cat(sprintf(" ERROR: %s\n", res_A$error))
  }
  
  cat("  Fitting Ref...")
  t0 <- Sys.time()
  res_Ref <- run_loo_single(study_id, dpoly(degree = 1), "Ref",
                            Griselda_clean, n_iter, n_chains, n_burnin, n_thin)
  t1 <- Sys.time()
  if (is.null(res_Ref$error)) {
    cat(sprintf(" done (%.1f min, DIC=%.1f, converged=%s)\n",
                difftime(t1, t0, units = "mins"), res_Ref$dic, res_Ref$converged))
  } else {
    cat(sprintf(" ERROR: %s\n", res_Ref$error))
  }
  
  loo_results[[study_id]] <- list(A = res_A, Ref = res_Ref)
  saveRDS(loo_results, "results/griselda/GRISELDA_LOO_intermediate.rds")
}

cat("\n=== LOO Cross-Validation Complete ===\n")

# ---- LOO Summary ----
loo_summary <- do.call(rbind, lapply(names(loo_results), function(sid) {
  res <- loo_results[[sid]]
  rbind(
    data.frame(study = sid, model = "Model A",
               dic = res$A$dic, totresdev = res$A$totresdev,
               tau = res$A$tau, converged = res$A$converged,
               error = ifelse(is.null(res$A$error), NA, res$A$error),
               removed = paste(res$A$removed_agents, collapse = ", "),
               stringsAsFactors = FALSE),
    data.frame(study = sid, model = "Ref",
               dic = res$Ref$dic, totresdev = res$Ref$totresdev,
               tau = res$Ref$tau, converged = res$Ref$converged,
               error = ifelse(is.null(res$Ref$error), NA, res$Ref$error),
               removed = paste(res$Ref$removed_agents, collapse = ", "),
               stringsAsFactors = FALSE)
  )
}))

loo_compare <- loo_summary %>%
  filter(!is.na(dic)) %>%
  select(study, model, dic) %>%
  pivot_wider(names_from = model, values_from = dic) %>%
  mutate(diff = `Model A` - Ref, A_better = diff < 0)

cat("\nModel A better in", sum(loo_compare$A_better, na.rm = TRUE),
    "of", nrow(loo_compare), "studies\n")
cat("Mean DIC difference (A - Ref):", round(mean(loo_compare$diff, na.rm = TRUE), 2), "\n")

# ---- Save LOO results ----
write.csv(loo_summary, 
          "C:/Users/katie/Desktop/Depression-NMA/results/griselda/loo_summary.csv", 
          row.names = FALSE)
write.csv(loo_compare, 
          "C:/Users/katie/Desktop/Depression-NMA/results/griselda/loo_dic_comparison.csv", 
          row.names = FALSE)
saveRDS(loo_results, 
        "C:/Users/katie/Desktop/Depression-NMA/results/griselda/GRISELDA_LOO_final.rds")

cat("\n=== Session Info ===\n")
writeLines(capture.output(sessionInfo()), "results/griselda/session_info.txt")