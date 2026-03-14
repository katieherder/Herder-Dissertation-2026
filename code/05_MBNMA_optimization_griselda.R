################################################################################
# Chapter 4: GRISELDA Multi-Class Empirical Analysis
# Optimizing Model-Based Network Meta-Analysis for Pharmacologic Networks
#
# This script implements:
#   1. Network setup (using pre-cleaned Griselda_clean with dose_FE)
#   2. Model fitting (Models A, B, Ref) — parallel to SSRI analysis
#   3. Model comparison (DIC, LPML)
#   4. Dose-response estimation and plotting (class-level panels)
#   5. Agent rankings (SUCRA) at 1x fluoxetine-equivalent dose
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

# Set working directory
setwd("C:/Users/katie/Desktop/Depression-NMA")

# MCMC settings (consistent with SSRI analysis and Chapter 3)
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


# ==============================================================================
# 1. Create Network Object
# ==============================================================================

# Prepare data for MBNMAdose: use dose_FE (fluoxetine-equivalent doses)
# MBNMAdose expects: studyID, agent, dose, n, r
griselda_mbnma <- Griselda_clean %>%
  mutate(
    dose = ifelse(agent == "placebo", 0, dose_FE)
  ) %>%
  select(studyID, agent, dose, n, r) %>%
  as.data.frame()

# Quick check: no NAs in dose for active agents
stopifnot(sum(is.na(griselda_mbnma$dose)) == 0)

cat("\n=== Creating GRISELDA network ===\n")
cat("Studies:", n_distinct(griselda_mbnma$studyID), "\n")
cat("Arms:", nrow(griselda_mbnma), "\n")
cat("Agents:", n_distinct(griselda_mbnma$agent[griselda_mbnma$agent != "placebo"]), "\n")

griselda_net <- mbnma.network(griselda_mbnma)

# Placebo response rate
placebo_rate_g <- sum(Griselda_clean$r[Griselda_clean$agent == "placebo"]) /
  sum(Griselda_clean$n[Griselda_clean$agent == "placebo"])
cat("Pooled placebo response rate:", round(placebo_rate_g, 3), "\n")


# ==============================================================================
# 2. LPML Function (reuse from SSRI script)
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


# ==============================================================================
# 3. Convergence Check Function
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
    neff_vals <- bugs$summary[, "n.eff"]
    neff_vals <- neff_vals[!is.na(neff_vals) & neff_vals > 0]
    cat("  Min n.eff:", min(neff_vals), "\n")
    cat("  Median n.eff:", median(neff_vals), "\n")
  }
  invisible(NULL)
}

# ==============================================================================
# 4. Model Fitting
# ==============================================================================

# Get list of active agents for dmulti specification
active_agents <- sort(unique(griselda_mbnma$agent[griselda_mbnma$agent != "placebo"]))
cat("\nActive agents (", length(active_agents), "):", paste(active_agents, collapse = ", "), "\n")

# ---- Model A: Shared Emax (common ED50 across ALL agents) ----
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


# ---- Model B: Agent-specific Emax (dmulti) ----
cat("\n--- Fitting Model B: Agent-specific Emax (dmulti) ---\n")
# Build dmulti specification for all 21 agents
dmulti_list <- setNames(
  lapply(active_agents, function(a) demax()),
  active_agents
)

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

# ---- Lumped NMA: collapse all non-zero doses ----
cat("\n--- Fitting Lumped NMA ---\n")

griselda_lumped <- griselda_mbnma
griselda_lumped$dose[griselda_lumped$dose > 0] <- 1

# Drop studies with duplicate agent-dose arms after collapsing
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
# 5. Model Comparison: DIC and LPML
# ==============================================================================

lpml_A_g   <- compute_lpml(model_A_g, "Model A")
lpml_B_g   <- compute_lpml(model_B_g, "Model B")
lpml_Ref_g <- compute_lpml(model_Ref_g, "Ref")

lumped_pD_g <- as.numeric(model_Lumped_g$jagsresult$BUGSoutput$DIC)[1] - 
  model_Lumped_g$jagsresult$BUGSoutput$mean$deviance

comparison_table_g <- data.frame(
  Model = c("A: Shared Emax", "B: Agent-specific Emax", "Ref: Shared Linear", "Lumped NMA"),
  DIC   = round(c(model_A_g$BUGSoutput$DIC,
                  model_B_g$BUGSoutput$DIC,
                  model_Ref_g$BUGSoutput$DIC,
                  as.numeric(model_Lumped_g$jagsresult$BUGSoutput$DIC)[1]), 1),
  pD    = round(c(model_A_g$BUGSoutput$pD,
                  model_B_g$BUGSoutput$pD,
                  model_Ref_g$BUGSoutput$pD,
                  lumped_pD_g), 1),
  LPML  = round(c(lpml_A_g$lpml, lpml_B_g$lpml, lpml_Ref_g$lpml, NA), 1)
)
comparison_table_g$Dbar <- comparison_table_g$DIC - comparison_table_g$pD

cat("\n=== GRISELDA Model Comparison ===\n")
print(comparison_table_g)

jpeg("figures/griselda/GRISELDA_DevDev_A_vs_Ref.jpeg", width = 8, height = 6, units = "in", res = 300)
devdev(model_A_g, model_Ref_g, dev.type = "resdev")
dev.off()

# ==============================================================================
# 6. Dose-Response Estimation and Visualization
# ==============================================================================

pred_A_g   <- predict(model_A_g, E0 = placebo_rate_g, n.doses = 30)
pred_B_g   <- predict(model_B_g, E0 = placebo_rate_g, n.doses = 30)
pred_Ref_g <- predict(model_Ref_g, E0 = placebo_rate_g, n.doses = 30)

jpeg("figures/griselda/GRISELDA_DR_ModelA.jpeg", width = 14, height = 10, units = "in", res = 300)
plot(pred_A_g) + ggtitle("GRISELDA Model A: Shared Emax")
dev.off()

jpeg("figures/griselda/GRISELDA_DR_ModelB.jpeg", width = 14, height = 10, units = "in", res = 300)
plot(pred_B_g) + ggtitle("GRISELDA Model B: Agent-specific Emax")
dev.off()

jpeg("figures/griselda/GRISELDA_DR_Ref.jpeg", width = 14, height = 10, units = "in", res = 300)
plot(pred_Ref_g) + ggtitle("GRISELDA Ref: Shared Linear")
dev.off()

jpeg("figures/griselda/GRISELDA_Fitplot_A.jpeg", width = 14, height = 12, units = "in", res = 300)
fitplot(model_A_g)
dev.off()

jpeg("figures/griselda/GRISELDA_Fitplot_B.jpeg", width = 14, height = 12, units = "in", res = 300)
fitplot(model_B_g)
dev.off()

jpeg("figures/griselda/GRISELDA_Fitplot_Ref.jpeg", width = 14, height = 12, units = "in", res = 300)
fitplot(model_Ref_g)
dev.off()


# ==============================================================================
# 7. Agent Rankings
# ==============================================================================

# ---- Max observed dose per agent (primary ranking comparison) ----
max_rank_doses_g <- griselda_mbnma %>%
  filter(agent != "placebo") %>%
  group_by(agent) %>%
  summarise(max_dose = max(dose), .groups = "drop")

cat("\n=== Max observed FE doses for ranking ===\n")
print(max_rank_doses_g, n = Inf)

max_rank_doses_g_list <- as.list(setNames(max_rank_doses_g$max_dose, max_rank_doses_g$agent))

# Also keep 1x equivalent doses for secondary comparison
rank_doses_g <- setNames(
  as.list(rep(20, length(active_agents))),
  active_agents
)

# SUCRA computation function
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

# ---- 7a. Rankings at max observed dose ----
cat("\n=== Rankings from Model A at max observed dose ===\n")
pred_A_max_g <- predict(model_A_g, E0 = placebo_rate_g, exact.doses = max_rank_doses_g_list)
rank_A_max_g <- rank(pred_A_max_g, lower_better = FALSE)
print(rank_A_max_g)

cat("\n=== Rankings from Model B at max observed dose ===\n")
pred_B_max_g <- predict(model_B_g, E0 = placebo_rate_g, exact.doses = max_rank_doses_g_list)
rank_B_max_g <- rank(pred_B_max_g, lower_better = FALSE)
print(rank_B_max_g)

cat("\n=== Rankings from Ref at max observed dose ===\n")
pred_Ref_max_g <- predict(model_Ref_g, E0 = placebo_rate_g, exact.doses = max_rank_doses_g_list)
rank_Ref_max_g <- rank(pred_Ref_max_g, lower_better = FALSE)
print(rank_Ref_max_g)

# ---- 7b. Lumped NMA rankings (directly from posterior) ----
cat("\n=== Rankings from Lumped NMA ===\n")
lumped_sims_g <- model_Lumped_g$jagsresult$BUGSoutput$sims.list$d
# Column 1 is placebo, remaining columns are active agents
# Need to figure out agent ordering from trt.labs
lumped_trt_labs <- model_Lumped_g$trt.labs
lumped_active_labs <- lumped_trt_labs[!grepl("Placebo|placebo", lumped_trt_labs)]
lumped_agent_names <- gsub("_1$", "", lumped_active_labs)

d_agents_g <- lumped_sims_g[, 2:ncol(lumped_sims_g)]
colnames(d_agents_g) <- lumped_agent_names

rank_mat_lumped_g <- t(apply(d_agents_g, 1, function(row) base::rank(-row)))
lumped_mean_ranks_g <- colMeans(rank_mat_lumped_g)
lumped_median_d_g <- apply(d_agents_g, 2, median)

cat("Lumped NMA posterior median treatment effects:\n")
print(round(lumped_median_d_g, 3))
cat("\nLumped NMA mean ranks:\n")
print(round(lumped_mean_ranks_g, 2))

# ---- 7c. Rankings at 1x equivalent dose (secondary comparison) ----
cat("\n=== Rankings from Model A at 20mg FE ===\n")
pred_A_rank_g <- predict(model_A_g, E0 = placebo_rate_g, exact.doses = rank_doses_g)
rank_A_g <- rank(pred_A_rank_g, lower_better = FALSE)
print(rank_A_g)

cat("\n=== Rankings from Model B at 20mg FE ===\n")
pred_B_rank_g <- predict(model_B_g, E0 = placebo_rate_g, exact.doses = rank_doses_g)
rank_B_g <- rank(pred_B_rank_g, lower_better = FALSE)
print(rank_B_g)

cat("\n=== Rankings from Ref at 20mg FE ===\n")
pred_Ref_rank_g <- predict(model_Ref_g, E0 = placebo_rate_g, exact.doses = rank_doses_g)
rank_Ref_g <- rank(pred_Ref_rank_g, lower_better = FALSE)
print(rank_Ref_g)

# ---- 7d. Compile SUCRA across models (max dose) ----
ranking_results_g <- bind_rows(
  compute_sucra_from_summary(rank_A_max_g, "Model A"),
  compute_sucra_from_summary(rank_B_max_g, "Model B"),
  compute_sucra_from_summary(rank_Ref_max_g, "Ref: Linear")
)

# Add Lumped SUCRA manually
n_agents_g <- length(lumped_agent_names)
lumped_sucra_g <- (n_agents_g - lumped_mean_ranks_g) / (n_agents_g - 1)
lumped_sucra_df_g <- data.frame(
  agent       = paste0(lumped_agent_names, "_", unlist(max_rank_doses_g_list[lumped_agent_names])),
  median_rank = apply(rank_mat_lumped_g, 2, median),
  mean_rank   = lumped_mean_ranks_g,
  sucra       = lumped_sucra_g,
  model       = "Lumped NMA",
  stringsAsFactors = FALSE
)
ranking_results_g <- bind_rows(ranking_results_g, lumped_sucra_df_g)

# Add class information
class_lookup <- Griselda_clean %>%
  filter(agent != "placebo") %>%
  distinct(agent, Class)

ranking_results_g <- ranking_results_g %>%
  mutate(agent_base = gsub("_[0-9.]+$", "", agent)) %>%
  left_join(class_lookup, by = c("agent_base" = "agent")) %>%
  dplyr::select(-agent_base)

sucra_wide_g <- ranking_results_g %>%
  dplyr::select(agent, Class, model, sucra) %>%
  pivot_wider(names_from = model, values_from = sucra) %>%
  arrange(Class, desc(`Model A`))

cat("\n=== SUCRA Comparison Across Models (max observed dose) ===\n")
print(sucra_wide_g, n = Inf, width = Inf)

# ---- 7e. SUCRA plot by class ----
ranking_results_g$agent_clean <- ranking_results_g$agent %>%
  str_replace_all("_[0-9.]+$", "") %>%
  str_replace_all("_", " ") %>%
  str_to_title()

# Get Model B ordering (highest SUCRA first at top after coord_flip)
model_b_order <- ranking_results_g %>%
  filter(model == "Model B") %>%
  arrange(sucra) %>%
  pull(agent_clean)

ranking_results_g$agent_clean <- factor(ranking_results_g$agent_clean, levels = unique(model_b_order))

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

# Save rank plots
jpeg("figures/griselda/GRISELDA_Rank_A.jpeg", width = 10, height = 8, units = "in", res = 300)
plot(rank_A_max_g)
dev.off()

jpeg("figures/griselda/GRISELDA_Rank_B.jpeg", width = 10, height = 8, units = "in", res = 300)
plot(rank_B_max_g)
dev.off()

jpeg("figures/griselda/GRISELDA_Rank_Ref.jpeg", width = 10, height = 8, units = "in", res = 300)
plot(rank_Ref_max_g)
dev.off()

# ==============================================================================
# 8. Prior Sensitivity Analysis
# ==============================================================================

# Same prior specifications as SSRI analysis
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
    
    # Rankings
    pred_i <- predict(model_i, E0 = placebo_rate_g, exact.doses = rank_doses_g)
    rank_i <- rank(pred_i, lower_better = FALSE)
    sucra_df_i <- compute_sucra_from_summary(rank_i, spec$label)
    sucra_i <- setNames(sucra_df_i$sucra, sucra_df_i$agent)
    
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
       y = "Density", fill = "Prior on τ", color = "Prior on τ") +
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


################################################################################
# Chapter 4: GRISELDA LOO Cross-Validation
# 
# Run this AFTER Models A, B, Ref have been fit and saved.
# Only Models A and Ref are included in LOO (Model B has convergence issues
# on full data and will be worse on reduced data).
################################################################################

# ==============================================================================
# Identify LOO candidate studies
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
cat("  Of which within-class:", sum(loo_candidates$n_classes == 1), "\n")

loo_studies <- loo_candidates %>%
  filter(n_classes > 1) %>%
  pull(studyID)

cat("LOO studies selected:", length(loo_studies), "\n\n")

# ==============================================================================
# LOO function (simplified — records model fit on reduced data, no predict)
# ==============================================================================

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
    
    # Check which held-out agents were removed entirely
    reduced_agents <- unique(reduced_mbnma$agent[reduced_mbnma$agent != "placebo"])
    held_agents <- held_out %>%
      filter(agent != "placebo") %>%
      pull(agent) %>%
      unique()
    removed_agents <- setdiff(held_agents, reduced_agents)
    
    # Fit model on reduced data
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
    converged <- max(rhat_vals) < 1.05
    totresdev <- model_loo$BUGSoutput$summary["totresdev", "50%"]
    tau_med <- model_loo$BUGSoutput$summary["sd", "50%"]
    
    return(list(
      study          = held_out_study,
      model          = model_name,
      dic            = model_loo$BUGSoutput$DIC,
      pD             = model_loo$BUGSoutput$pD,
      totresdev      = totresdev,
      tau            = tau_med,
      max_rhat       = max(rhat_vals),
      converged      = converged,
      removed_agents = removed_agents,
      n_arms_removed = nrow(held_out),
      error          = NULL
    ))
    
  }, error = function(e) {
    return(list(
      study          = held_out_study,
      model          = model_name,
      dic            = NA,
      pD             = NA,
      totresdev      = NA,
      tau            = NA,
      max_rhat       = NA,
      converged      = FALSE,
      removed_agents = NULL,
      n_arms_removed = NA,
      error          = e$message
    ))
  })
}

# ==============================================================================
# Run LOO loop (Models A and Ref only)
# ==============================================================================

cat("Starting LOO cross-validation (Models A and Ref)...\n")
cat("Estimated time: ~10 min per study x 47 studies = ~8 hours\n\n")

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

# ==============================================================================
# Summarize LOO results
# ==============================================================================

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

cat("\n=== LOO Summary ===\n")
cat("Studies completed:", length(loo_results), "of", length(loo_studies), "\n")
cat("Errors (Model A):", sum(!is.na(loo_summary$error[loo_summary$model == "Model A"])), "\n")
cat("Errors (Ref):", sum(!is.na(loo_summary$error[loo_summary$model == "Ref"])), "\n")
cat("Studies with agents removed:", sum(loo_summary$removed != "" & loo_summary$model == "Model A"), "\n")

# Compare DIC across held-out studies
loo_compare <- loo_summary %>%
  filter(!is.na(dic)) %>%
  select(study, model, dic) %>%
  pivot_wider(names_from = model, values_from = dic) %>%
  mutate(diff = `Model A` - Ref,
         A_better = diff < 0)

cat("\nModel A better in", sum(loo_compare$A_better, na.rm = TRUE), 
    "of", nrow(loo_compare), "studies\n")
cat("Ref better in", sum(!loo_compare$A_better, na.rm = TRUE), 
    "of", nrow(loo_compare), "studies\n")
cat("Mean DIC difference (A - Ref):", round(mean(loo_compare$diff, na.rm = TRUE), 2), "\n")

cat("\n=== All results saved ===\n")

# ==============================================================================
# 10. Summary Output
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  RESULTS SUMMARY — Chapter 4 GRISELDA Empirical Analysis\n")
cat("================================================================\n\n")

cat("--- Model Comparison ---\n")
print(comparison_table_g)

cat("\n--- Agent Rankings (SUCRA at 20mg FE) ---\n")
print(sucra_wide_g, n = Inf, width = Inf)

cat("\n--- Prior Sensitivity ---\n")
print(prior_summary_g)

cat("\n--- SUCRA Range Across Priors ---\n")
print(sucra_range_g, n = Inf)

cat("\n--- LOO Summary ---\n")
cat("Total LOO studies:", length(loo_results), "\n")
n_errors <- sum(sapply(loo_results, function(x) !is.null(x$A$error)))
cat("Studies with errors:", n_errors, "\n")

cat("\n================================================================\n")
cat("  Analysis complete.\n")
cat("================================================================\n")

# ==============================================================================
# 11. Save Results
# ==============================================================================

# ---- Model comparison table (DIC, pD, LPML for all 4 models) ----
write.csv(comparison_table_g, "results/griselda/model_comparison.csv", row.names = FALSE)

# ---- Convergence diagnostics ----
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

# ---- SUCRA comparison at max observed dose ----
write.csv(sucra_wide_g, "results/griselda/sucra_comparison_max_dose.csv", row.names = FALSE)

# ---- Full ranking results ----
write.csv(ranking_results_g, "results/griselda/ranking_results_full.csv", row.names = FALSE)

# ---- Max observed FE doses used for ranking ----
write.csv(max_rank_doses_g, "results/griselda/max_equiv_doses.csv", row.names = FALSE)

# ---- Lumped NMA treatment effects and ranks ----
lumped_results_g <- data.frame(
  agent = lumped_agent_names,
  d_median = lumped_median_d_g,
  d_lower = apply(d_agents_g, 2, quantile, 0.025),
  d_upper = apply(d_agents_g, 2, quantile, 0.975),
  mean_rank = lumped_mean_ranks_g,
  sucra = lumped_sucra_g
)
write.csv(lumped_results_g, "results/griselda/lumped_nma_results.csv", row.names = FALSE)

# ---- Heterogeneity (tau) posterior from each model ----
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

# ---- Prior sensitivity results ----
write.csv(prior_summary_g, "results/griselda/prior_sensitivity_summary.csv", row.names = FALSE)
write.csv(sucra_by_prior_g, "results/griselda/prior_sensitivity_sucra.csv", row.names = FALSE)

# ---- LOO cross-validation results ----
write.csv(loo_summary, "results/griselda/loo_summary.csv", row.names = FALSE)
write.csv(loo_compare, "results/griselda/loo_dic_comparison.csv", row.names = FALSE)
saveRDS(loo_results, "results/griselda/GRISELDA_LOO_final.rds")

# ---- Full model objects ----
save(model_A_g, model_B_g, model_Ref_g, model_Lumped_g,
     comparison_table_g, ranking_results_g, sucra_wide_g,
     tau_summary_g, prior_results_g, loo_results,
     file = "results/griselda/griselda_all_models.RData")

cat("\n=== All results saved to results/griselda/ ===\n")
cat("  model_comparison.csv          — DIC, pD, LPML for all 4 models\n")
cat("  convergence_diagnostics.csv   — Max Rhat and min n.eff per model\n")
cat("  sucra_comparison_max_dose.csv — SUCRA at max observed dose, wide format\n")
cat("  ranking_results_full.csv      — Full ranking details (all models)\n")
cat("  max_equiv_doses.csv           — Max FE dose per agent used for ranking\n")
cat("  lumped_nma_results.csv        — Lumped NMA treatment effects and ranks\n")
cat("  heterogeneity_tau.csv         — Tau posterior summaries per model\n")
cat("  prior_sensitivity_summary.csv — Tau and fit across prior specs\n")
cat("  prior_sensitivity_sucra.csv   — SUCRA stability across priors\n")
cat("  loo_summary.csv               — LOO study-level results\n")
cat("  loo_dic_comparison.csv        — LOO DIC comparison (A vs Ref)\n")
cat("  GRISELDA_LOO_final.rds        — Full LOO results object\n")
cat("  griselda_all_models.RData     — Full model objects for reproducibility\n")

# ==============================================================================
# Additional Rankings at 10mg and 40mg FE (run after loading saved results)
# ==============================================================================
# 
# Prerequisites: load("results/griselda/griselda_all_models.RData")
# This gives you: model_B_g, placebo_rate_g (via Griselda_clean), active_agents, etc.

# Reconstruct placebo rate and active agents if not in environment
if (!exists("placebo_rate_g")) {
  placebo_rate_g <- sum(Griselda_clean$r[Griselda_clean$agent == "placebo"]) /
    sum(Griselda_clean$n[Griselda_clean$agent == "placebo"])
}

if (!exists("active_agents")) {
  active_agents <- sort(unique(Griselda_clean$agent[Griselda_clean$agent != "placebo"]))
}

# ---- 10mg FE rankings ----
rank_doses_10 <- setNames(as.list(rep(10, length(active_agents))), active_agents)
pred_B_10 <- predict(model_B_g, E0 = placebo_rate_g, exact.doses = rank_doses_10)
rank_B_10 <- rank(pred_B_10, lower_better = FALSE)

cat("\n=== Rankings from Model B at 10mg FE ===\n")
print(rank_B_10)

# ---- 40mg FE rankings ----
rank_doses_40 <- setNames(as.list(rep(40, length(active_agents))), active_agents)
pred_B_40 <- predict(model_B_g, E0 = placebo_rate_g, exact.doses = rank_doses_40)
rank_B_40 <- rank(pred_B_40, lower_better = FALSE)

cat("\n=== Rankings from Model B at 40mg FE ===\n")
print(rank_B_40)

# ---- Compile mean ranks across all dose levels ----
# Extract mean ranks from summary objects
extract_mean_ranks <- function(rank_obj) {
  summ <- summary(rank_obj)$Predictions
  ranks <- setNames(summ$mean, as.character(summ$rank.param))
  # Remove placebo
  ranks <- ranks[!grepl("[Pp]lacebo", names(ranks))]
  # Clean agent names (remove _dose suffix)
  names(ranks) <- gsub("_[0-9.]+$", "", names(ranks))
  return(ranks)
}

ranks_10  <- extract_mean_ranks(rank_B_10)
ranks_20  <- extract_mean_ranks(rank_B_g)      # already computed in main script
ranks_40  <- extract_mean_ranks(rank_B_40)
ranks_max <- extract_mean_ranks(rank_B_max_g)   # already computed in main script

# Build comparison table
rank_dose_comparison <- data.frame(
  Agent   = names(ranks_10),
  `10mg`  = round(ranks_10, 2),
  `20mg`  = round(ranks_20[names(ranks_10)], 2),
  `40mg`  = round(ranks_40[names(ranks_10)], 2),
  Max     = round(ranks_max[names(ranks_10)], 2),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Sort by 20mg rank for readability
rank_dose_comparison <- rank_dose_comparison[order(rank_dose_comparison$`20mg`), ]

cat("\n=== Model B Mean Ranks Across Dose Levels ===\n")
print(rank_dose_comparison, row.names = FALSE)

# Save
write.csv(rank_dose_comparison, "results/griselda/rank_dose_comparison_modelB.csv", row.names = FALSE)
cat("\nSaved to results/griselda/rank_dose_comparison_modelB.csv\n")