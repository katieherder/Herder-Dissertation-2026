################################################################################
# Chapter 4: Simulation Study
# Functional Form Misspecification in MBNMA for Pharmacologic Networks
#
# This script implements:
#   1. Data-generating function calibrated to SSRI empirical results
#   2. Model-fitting wrapper (Models A, B, Ref, Lumped)
#   3. Metric extraction (bias, MSE, coverage, ranking accuracy, SUCRA,
#      per-agent rank match rate, DIC, LPML)
#   4. Replication loop (resumable from checkpoints)
#   5. Summary tables and figures
#
# Dependencies: MBNMAdose (>= 0.5.0), R2jags, ggplot2, dplyr, tidyr,
#               parallel, coda, MASS
################################################################################

# ==============================================================================
# 0. Setup
# ==============================================================================
library(MBNMAdose)
library(R2jags)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(coda)

set.seed(54321)
setwd("C:/Users/katie/Desktop/Depression-NMA/Chapter 4")

# ggplot theme (consistent with empirical analysis)
theme_diss <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )
theme_set(theme_diss)

# MCMC settings (consistent with empirical analysis)
n_iter   <- 20000
n_chains <- 3
n_burnin <- floor(n_iter / 2)   # 10,000
n_thin   <- max(1, floor((n_iter - n_burnin) / 1000))  # 10

cat("Simulation MCMC settings:\n",
    "  Iterations:", n_iter, "\n",
    "  Burn-in:", n_burnin, "\n",
    "  Thinning:", n_thin, "\n",
    "  Chains:", n_chains, "\n",
    "  Posterior samples:", n_chains * ((n_iter - n_burnin) / n_thin), "\n")


# ==============================================================================
# 1. True Parameters (Calibrated to SSRI Model B Posteriors)
# ==============================================================================

dose_equiv <- data.frame(
  agent    = c("citalopram", "escitalopram", "fluoxetine", "paroxetine", "sertraline"),
  equiv_mg = c(20, 9, 20, 17, 49.3),
  stringsAsFactors = FALSE
)

true_params <- data.frame(
  agent  = c("citalopram", "escitalopram", "fluoxetine", "paroxetine", "sertraline"),
  emax   = c(1.96, 1.13, 0.51, 2.40, 0.82),     # log-odds scale
  ed50   = c(91.0, 70.8, 77.1, 104.1, 58.4),     # fluoxetine-equivalent mg/day
  stringsAsFactors = FALSE
)

true_tau <- 0.27

true_mu_mean <- -0.44   # logit(0.393)
true_mu_sd   <- 0.30    # empirical variability across 60 trials

true_params$max_equiv_dose <- NA
true_params$d_true <- NA
true_params$true_rank <- NA

cat("\n=== True Parameters ===\n")
print(true_params, digits = 3)
cat("\nTrue tau:", true_tau, "\n")


# ==============================================================================
# 2. Extract Empirical Network Structure
# ==============================================================================

ssri_equiv <- ssri %>%
  left_join(dose_equiv, by = "agent") %>%
  mutate(dose_equiv = ifelse(agent == "placebo", 0, dose / equiv_mg * 20)) %>%
  select(-equiv_mg)

network_template <- ssri_equiv %>%
  group_by(studyID) %>%
  mutate(arm = row_number()) %>%
  ungroup() %>%
  select(studyID, arm, agent, dose_equiv, n, r)

empirical_arm_sizes <- network_template$n

cat("\n=== Network Template ===\n")
cat("  Studies:", length(unique(network_template$studyID)), "\n")
cat("  Total arms:", nrow(network_template), "\n")
cat("  Agents:", paste(sort(unique(network_template$agent[network_template$agent != "placebo"])),
                       collapse = ", "), "\n")
cat("  Arm size range:", min(empirical_arm_sizes), "--", max(empirical_arm_sizes), "\n")

max_equiv_doses <- network_template %>%
  filter(agent != "placebo") %>%
  group_by(agent) %>%
  summarise(max_dose_equiv = max(dose_equiv), .groups = "drop")

cat("\n=== Max Observed Equivalent Doses ===\n")
print(max_equiv_doses)

true_params <- true_params %>%
  dplyr::select(-max_equiv_dose, -d_true, -true_rank) %>%
  left_join(max_equiv_doses, by = "agent") %>%
  rename(max_equiv_dose = max_dose_equiv)

true_params$d_true <- with(true_params, emax * max_equiv_dose / (ed50 + max_equiv_dose))
true_params$true_rank <- base::rank(-true_params$d_true)

cat("\n=== True Parameters (at max observed dose) ===\n")
print(true_params, digits = 3)
cat("\nTrue ranking order:",
    paste(true_params$agent[order(true_params$true_rank)], collapse = " > "), "\n")


# ==============================================================================
# 3. Data-Generating Function
# ==============================================================================

generate_ssri_dataset <- function(template, true_params, true_tau,
                                  true_mu_mean, true_mu_sd,
                                  resample_n = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sim_data <- template
  studies <- unique(sim_data$studyID)
  n_studies <- length(studies)
  
  if (resample_n) {
    study_medians <- template %>%
      group_by(studyID) %>%
      summarise(med_n = median(n), .groups = "drop")
    
    drawn_medians <- sample(study_medians$med_n, n_studies, replace = TRUE)
    names(drawn_medians) <- studies
    
    sim_data <- sim_data %>%
      group_by(studyID) %>%
      mutate(
        orig_median = median(n),
        multiplier  = drawn_medians[as.character(studyID[1])] / orig_median,
        n           = pmax(10L, as.integer(round(n * multiplier)))
      ) %>%
      ungroup() %>%
      select(-orig_median, -multiplier)
  }
  # NOTE: One empirical study (citalopram) has no placebo arm; arm 1 is an active
  # dose. The contrast parameterization (d_model = f_dose - f_ref) handles this
  # correctly — mu_i absorbs the reference arm's dose-response contribution.
  # All evaluation metrics depend on contrasts, so drawing mu_i from the
  # placebo-centered distribution does not affect simulation results.
  mu <- rnorm(n_studies, mean = true_mu_mean, sd = true_mu_sd)
  names(mu) <- studies
  
  sim_data <- sim_data %>%
    left_join(true_params %>% select(agent, emax, ed50), by = "agent") %>%
    mutate(
      f_dose = ifelse(agent == "placebo", 0, emax * dose_equiv / (ed50 + dose_equiv))
    )
  
  sim_data <- sim_data %>%
    group_by(studyID) %>%
    mutate(
      f_ref = f_dose[arm == 1],
      d_model = f_dose - f_ref
    ) %>%
    ungroup()
  
  sim_data$epsilon <- 0
  sim_data$mu_study <- mu[as.character(sim_data$studyID)]
  
  for (s in studies) {
    idx <- which(sim_data$studyID == s & sim_data$arm != 1)
    k <- length(idx)
    
    if (k == 0) next
    
    if (k == 1) {
      sim_data$epsilon[idx] <- rnorm(1, 0, true_tau)
    } else {
      Sigma <- matrix(0.5 * true_tau^2, nrow = k, ncol = k)
      diag(Sigma) <- true_tau^2
      sim_data$epsilon[idx] <- MASS::mvrnorm(1, mu = rep(0, k), Sigma = Sigma)
    }
  }
  
  sim_data <- sim_data %>%
    mutate(
      logit_p = mu_study + d_model + epsilon,
      p = plogis(logit_p),
      r = rbinom(n(), size = n, prob = p)
    )
  
  out <- sim_data %>%
    select(studyID, agent, dose = dose_equiv, r, n) %>%
    as.data.frame()
  
  return(out)
}


# ==============================================================================
# 4. LPML Computation Function (same as empirical analysis)
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
  return(list(lpml = lpml, cpo = log_cpo))
}


# ==============================================================================
# 5. Metric Extraction Function
# ==============================================================================
#
# RANKING APPROACH (uniform across all models):
#   For every model (A, B, Ref, Lumped), we first construct d_hat — a matrix
#   of posterior samples of the RELATIVE treatment effect (log-odds vs placebo)
#   at each agent's maximum observed equivalent dose.
#
#   - Models A & B: d_hat[,a] = emax_samples * dose_a / (ed50_samples + dose_a)
#   - Model Ref:    d_hat[,a] = beta.1_samples * dose_a
#   - Lumped NMA:   d_hat[,a] = d_samples (direct NMA treatment effects)
#
#   Rankings are then computed IDENTICALLY for all models:
#     For each posterior draw, rank agents 1–5 (rank 1 = largest effect = best).
#     Report the mean rank across all posterior draws.
#
#   This ensures all models are compared on the same scale (relative treatment
#   effects) using the same ranking mechanism.
#
# METRICS RETURNED PER MODEL:
#   - bias, sq_error, coverage (per agent)
#   - d_median, d_lower, d_upper (per agent)
#   - est_ranks (mean posterior rank per agent)
#   - sucra (per agent)
#   - spearman (correlation of mean est_ranks with true_rank)
# ==============================================================================

extract_metrics <- function(fitted_models, true_params, placebo_rate = 0.393,
                            param_map = NULL) {
  if (is.null(fitted_models)) return(NULL)
  
  model_names <- c("A", "B", "Ref", "Lumped")
  metrics <- list()
  
  for (m in model_names) {
    mod <- fitted_models[[m]]
    if (is.null(mod)) {
      metrics[[m]] <- list(failed = TRUE)
      next
    }
    
    # --- Access BUGS output (different structure for nma vs mbnma objects) ---
    if (m == "Lumped") {
      bugs <- mod$jagsresult$BUGSoutput
    } else {
      bugs <- mod$BUGSoutput
    }
    
    # --- Convergence ---
    rhat_vals <- bugs$summary[, "Rhat"]
    rhat_vals <- rhat_vals[!is.na(rhat_vals)]
    max_rhat <- max(rhat_vals)
    converged <- max_rhat < 1.05
    
    # --- DIC ---
    dic <- bugs$DIC
    pD  <- bugs$pD
    
    # --- LPML ---
    if (m == "Lumped") {
      lpml <- NA
    } else {
      lpml_out <- compute_lpml(mod, paste("Model", m))
      lpml <- lpml_out$lpml
    }
    
    # --- Extract predicted treatment effects at each agent's max dose ---
    agent_names <- true_params$agent
    n_agents <- nrow(true_params)
    sims <- bugs$sims.list
    n_sims <- bugs$n.sims
    
    d_hat <- matrix(NA, nrow = n_sims, ncol = n_agents)
    colnames(d_hat) <- agent_names
    
    extraction_failed <- FALSE
    
    if (m %in% c("A", "B")) {
      for (a in seq_along(agent_names)) {
        col_idx <- a
        if (is.matrix(sims$emax) && col_idx <= ncol(sims$emax)) {
          emax_samps <- sims$emax[, col_idx]
        } else if (is.numeric(sims$emax) && n_agents == 1) {
          emax_samps <- sims$emax
        } else {
          cat(sprintf("  WARNING: Cannot extract emax for %s in Model %s\n",
                      agent_names[a], m))
          extraction_failed <- TRUE
          break
        }
        if (m == "A") {
          if (is.matrix(sims$ed50)) {
            ed50_samps <- sims$ed50[, 1]
          } else {
            ed50_samps <- as.numeric(sims$ed50)
          }
        } else {
          if (is.matrix(sims$ed50) && col_idx <= ncol(sims$ed50)) {
            ed50_samps <- sims$ed50[, col_idx]
          } else {
            cat(sprintf("  WARNING: Cannot extract ed50 for %s in Model %s\n",
                        agent_names[a], m))
            extraction_failed <- TRUE
            break
          }
        }
        d_a <- true_params$max_equiv_dose[a]
        d_hat[, a] <- emax_samps * d_a / (ed50_samps + d_a)
      }
      
    } else if (m == "Ref") {
      for (a in seq_along(agent_names)) {
        col_idx <- a
        if (is.matrix(sims$beta.1) && col_idx <= ncol(sims$beta.1)) {
          beta_samps <- sims$beta.1[, col_idx]
        } else {
          cat(sprintf("  WARNING: Cannot extract beta.1 for %s in Ref\n",
                      agent_names[a]))
          extraction_failed <- TRUE
          break
        }
        d_hat[, a] <- beta_samps * true_params$max_equiv_dose[a]
      }
      
    } else if (m == "Lumped") {
      for (a in seq_along(agent_names)) {
        col_idx <- a + 1   # column 1 is placebo
        if (is.matrix(sims$d) && col_idx <= ncol(sims$d)) {
          d_hat[, a] <- sims$d[, col_idx]
        } else {
          cat(sprintf("  WARNING: Cannot extract d for %s in Lumped\n",
                      agent_names[a]))
          extraction_failed <- TRUE
          break
        }
      }
    }
    
    if (extraction_failed) {
      metrics[[m]] <- list(failed = TRUE, dic = dic, lpml = lpml,
                           converged = converged, max_rhat = max_rhat,
                           note = "Parameter extraction failed")
      next
    }
    
    # --- Compute metrics per agent ---
    d_median <- apply(d_hat, 2, median)
    d_lower  <- apply(d_hat, 2, quantile, probs = 0.025)
    d_upper  <- apply(d_hat, 2, quantile, probs = 0.975)
    
    bias     <- d_median - true_params$d_true
    sq_error <- (d_median - true_params$d_true)^2
    covered  <- (true_params$d_true >= d_lower) & (true_params$d_true <= d_upper)
    
    # --- Rankings (UNIFORM across all models) ---
    # For each posterior draw, rank agents by relative treatment effect.
    # rank(-row) means: largest effect gets rank 1 (best).
    rank_mat <- t(apply(d_hat, 1, function(row) base::rank(-row)))
    est_mean_ranks <- colMeans(rank_mat)
    names(est_mean_ranks) <- agent_names
    
    # SUCRA: mean of (n_agents - rank) / (n_agents - 1) across posterior draws
    sucra <- colMeans((n_agents - rank_mat) / (n_agents - 1))
    names(sucra) <- agent_names
    
    # Spearman correlation with true ranks
    spearman <- cor(est_mean_ranks, true_params$true_rank, method = "spearman")
    
    metrics[[m]] <- list(
      failed       = FALSE,
      converged    = converged,
      max_rhat     = max_rhat,
      dic          = dic,
      pD           = pD,
      lpml         = lpml,
      bias         = setNames(bias, agent_names),
      sq_error     = setNames(sq_error, agent_names),
      coverage     = setNames(as.integer(covered), agent_names),
      d_median     = setNames(d_median, agent_names),
      d_lower      = setNames(d_lower, agent_names),
      d_upper      = setNames(d_upper, agent_names),
      est_ranks    = setNames(est_mean_ranks, agent_names),
      sucra        = sucra,
      spearman     = spearman
    )
  }
  
  return(metrics)
}


# ==============================================================================
# 6. Model-Fitting Wrapper
# ==============================================================================

fit_all_models <- function(sim_data, rep_seed) {
  net <- tryCatch(
    mbnma.network(sim_data),
    error = function(e) {
      cat("  Network creation failed:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(net)) return(NULL)
  
  results <- list()
  
  # ---- Model A: Shared Emax ----
  tryCatch({
    results$A <- mbnma.run(
      network  = net,
      fun      = demax(emax = "rel", ed50 = "common"),
      method   = "random",
      likelihood = "binomial",
      link     = "logit",
      parameters.to.save = c("resdev", "emax", "ed50", "sd", "totresdev"),
      n.iter   = n_iter,
      n.chains = n_chains,
      n.burnin = n_burnin,
      n.thin   = n_thin,
      jags.seed = rep_seed
    )
  }, error = function(e) {
    cat("  Model A failed:", e$message, "\n")
    results$A <<- NULL
  })
  
  # ---- Model B: Agent-specific Emax ----
  tryCatch({
    results$B <- mbnma.run(
      network  = net,
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
      parameters.to.save = c("resdev", "emax", "ed50", "sd", "totresdev"),
      n.iter   = n_iter,
      n.chains = n_chains,
      n.burnin = n_burnin,
      n.thin   = n_thin,
      jags.seed = rep_seed
    )
  }, error = function(e) {
    cat("  Model B failed:", e$message, "\n")
    results$B <<- NULL
  })
  
  # ---- Reference Model: Shared Linear ----
  tryCatch({
    results$Ref <- mbnma.run(
      network  = net,
      fun      = dpoly(degree = 1),
      method   = "random",
      likelihood = "binomial",
      link     = "logit",
      parameters.to.save = c("resdev", "beta.1", "sd", "totresdev"),
      n.iter   = n_iter,
      n.chains = n_chains,
      n.burnin = n_burnin,
      n.thin   = n_thin,
      jags.seed = rep_seed
    )
  }, error = function(e) {
    cat("  Ref model failed:", e$message, "\n")
    results$Ref <<- NULL
  })
  
  # ---- Lumped NMA: collapse all non-zero doses to a single node per agent ----
  tryCatch({
    sim_data_lumped <- sim_data
    sim_data_lumped$dose[sim_data_lumped$dose > 0] <- 1
    
    sim_data_lumped <- sim_data_lumped %>%
      group_by(studyID, agent, dose) %>%
      mutate(dup = n() > 1) %>%
      group_by(studyID) %>%
      filter(!any(dup & dose > 0)) %>%
      ungroup() %>%
      dplyr::select(-dup)
    
    sim_data_lumped <- sim_data_lumped %>%
      group_by(studyID) %>%
      filter(n() >= 2) %>%
      ungroup() %>%
      as.data.frame()
    
    net_lumped <- mbnma.network(sim_data_lumped)
    
    results$Lumped <- nma.run(
      network    = net_lumped,
      method     = "random",
      likelihood = "binomial",
      link       = "logit",
      n.iter     = n_iter,
      n.chains   = n_chains,
      n.burnin   = n_burnin,
      n.thin     = n_thin,
      jags.seed  = rep_seed
    )
  }, error = function(e) {
    cat("  Lumped NMA failed:", e$message, "\n")
    results$Lumped <<- NULL
  })
  
  results$network <- net
  return(results)
}


# ==============================================================================
# 6b. Parameter Structure Diagnostic (run ONCE before full simulation)
# ==============================================================================

cat("\n================================================================\n")
cat("  PARAMETER STRUCTURE DIAGNOSTIC\n")
cat("================================================================\n\n")

test_data <- generate_ssri_dataset(
  template = network_template, true_params = true_params,
  true_tau = true_tau, true_mu_mean = true_mu_mean,
  true_mu_sd = true_mu_sd, resample_n = TRUE, seed = 999999
)

test_fitted <- fit_all_models(test_data, rep_seed = 999999)

for (m in c("A", "B", "Ref", "Lumped")) {
  
  if (m == "Lumped") {
    cat(sprintf("\n--- Model %s: sims.list names ---\n", m))
    if (!is.null(test_fitted[[m]]) && !is.null(test_fitted[[m]]$jagsresult)) {
      sims <- test_fitted[[m]]$jagsresult$BUGSoutput$sims.list
      print(names(sims))
      cat(sprintf("\n--- Model %s: d dimensions ---\n", m))
      if ("d" %in% names(sims)) {
        cat("  d: ", paste(dim(sims$d), collapse = " x "), "\n")
      }
    } else {
      cat("  Lumped model is NULL\n")
    }
    next
  }
  
  cat(sprintf("\n--- Model %s: sims.list names ---\n", m))
  print(names(test_fitted[[m]]$BUGSoutput$sims.list))
  
  cat(sprintf("\n--- Model %s: emax/ed50/beta dimensions ---\n", m))
  sims <- test_fitted[[m]]$BUGSoutput$sims.list
  if ("emax" %in% names(sims)) {
    cat("  emax: ", paste(dim(sims$emax), collapse = " x "), "\n")
    emax_params <- grep("^emax\\[", rownames(test_fitted[[m]]$BUGSoutput$summary), value = TRUE)
    cat("  emax param names:", paste(emax_params, collapse = ", "), "\n")
  }
  if ("ed50" %in% names(sims)) {
    ed50_dim <- if (is.null(dim(sims$ed50))) length(sims$ed50) else dim(sims$ed50)
    cat("  ed50: ", paste(ed50_dim, collapse = " x "), "\n")
    ed50_params <- grep("^ed50", rownames(test_fitted[[m]]$BUGSoutput$summary), value = TRUE)
    cat("  ed50 param names:", paste(ed50_params, collapse = ", "), "\n")
  }
  if ("beta.1" %in% names(sims)) {
    cat("  beta.1: ", paste(dim(sims$beta.1), collapse = " x "), "\n")
    beta_params <- grep("^beta\\.1\\[", rownames(test_fitted[[m]]$BUGSoutput$summary), value = TRUE)
    cat("  beta.1 param names:", paste(beta_params, collapse = ", "), "\n")
  }
}

# --- Verify ranking with d_hat approach on test data ---
cat("\n--- Ranking verification (d_hat approach) ---\n")
test_metrics <- extract_metrics(test_fitted, true_params)
for (m in c("A", "B", "Ref", "Lumped")) {
  if (!is.null(test_metrics[[m]]) && !test_metrics[[m]]$failed) {
    cat(sprintf("\nModel %s (rho = %.3f):\n", m, test_metrics[[m]]$spearman))
    diag_df <- data.frame(
      agent       = true_params$agent,
      true_effect = round(true_params$d_true, 3),
      est_effect  = round(test_metrics[[m]]$d_median, 3),
      true_rank   = true_params$true_rank,
      est_rank    = round(test_metrics[[m]]$est_ranks, 2),
      sucra       = round(test_metrics[[m]]$sucra, 3)
    )
    print(diag_df[order(-diag_df$est_effect), ])
  }
}

cat("\n================================================================\n")
cat("  REVIEW OUTPUT ABOVE before running full simulation.\n")
cat("  Verify: emax[,1] = citalopram, emax[,2] = escitalopram, etc.\n")
cat("  Verify: rankings are sensible (positive Spearman expected).\n")
cat("================================================================\n")

DIAGNOSTIC_PASSED <- TRUE

if (!DIAGNOSTIC_PASSED) {
  stop("Review diagnostic output above and set DIAGNOSTIC_PASSED <- TRUE to proceed.")
}


# ==============================================================================
# 7. Run Simulation (resumable from checkpoint)
# ==============================================================================

n_reps <- 100

cat("\n================================================================\n")
cat("  SIMULATION: Starting", n_reps, "replications\n")
cat("  Each replication fits 4 models (A, B, Ref, Lumped)\n")
cat("  Estimated time: ~", round(n_reps * 15 / 60, 1), "hours\n")
cat("================================================================\n\n")

# Initialize or resume from most recent checkpoint
checkpoint_files <- list.files("results/simulation", pattern = "^sim_results_checkpoint_\\d+\\.RData$")
if (length(checkpoint_files) > 0) {
  latest_rep <- max(as.integer(gsub("\\D", "", checkpoint_files)))
  latest_file <- sprintf("results/simulation/sim_results_checkpoint_%d.RData", latest_rep)
  cat("Resuming from checkpoint:", latest_file, "\n")
  load(latest_file)
  if (length(all_metrics) < n_reps) {
    all_metrics <- c(all_metrics, vector("list", n_reps - length(all_metrics)))
  }
} else {
  all_metrics <- vector("list", n_reps)
}

for (rep in 1:n_reps) {
  
  # Skip already-completed replications (resume where left off)
  if (!is.null(all_metrics[[rep]])) {
    cat(sprintf("  Replication %d/%d: already complete, skipping\n", rep, n_reps))
    next
  }
  
  rep_seed <- 100000 + rep
  
  cat(sprintf("\n--- Replication %d/%d (seed: %d) ---\n", rep, n_reps, rep_seed))
  
  sim_data <- generate_ssri_dataset(
    template     = network_template,
    true_params  = true_params,
    true_tau     = true_tau,
    true_mu_mean = true_mu_mean,
    true_mu_sd   = true_mu_sd,
    resample_n   = TRUE,
    seed         = rep_seed
  )
  
  fitted <- fit_all_models(sim_data, rep_seed)
  
  all_metrics[[rep]] <- extract_metrics(fitted, true_params)
  
  # Progress summary
  if (!is.null(all_metrics[[rep]])) {
    for (m in c("A", "B", "Ref", "Lumped")) {
      if (!is.null(all_metrics[[rep]][[m]]) && !all_metrics[[rep]][[m]]$failed) {
        cat(sprintf("  %s: DIC=%.1f LPML=%s rho=%.3f Rhat_max=%.3f\n",
                    m,
                    all_metrics[[rep]][[m]]$dic,
                    ifelse(is.na(all_metrics[[rep]][[m]]$lpml), "N/A",
                           sprintf("%.1f", all_metrics[[rep]][[m]]$lpml)),
                    all_metrics[[rep]][[m]]$spearman,
                    all_metrics[[rep]][[m]]$max_rhat))
      } else {
        cat(sprintf("  %s: FAILED\n", m))
      }
    }
  }
  
  # Checkpoint every 5 replications (more frequent for 100-rep run)
  if (rep %% 5 == 0) {
    save(all_metrics, file = sprintf("results\simulation\sim_results_checkpoint_%d.RData", rep))
    cat(sprintf("  [Checkpoint saved: %d replications]\n", rep))
  }
}

save(all_metrics, true_params, true_tau, file = "results\simulation\sim_results_final.RData")
cat("\n=== Simulation complete. Results saved. ===\n")


# ==============================================================================
# 8. Compile Results
# ==============================================================================

cat("\n================================================================\n")
cat("  COMPILING SIMULATION RESULTS\n")
cat("================================================================\n\n")

agent_names <- true_params$agent
model_names <- c("A", "B", "Ref", "Lumped")

# ---- 8a. Convergence summary ----
convergence_summary <- expand.grid(model = model_names, stringsAsFactors = FALSE)
convergence_summary$n_converged <- 0
convergence_summary$n_failed <- 0
convergence_summary$n_total <- 0

for (m in model_names) {
  conv_flags <- sapply(all_metrics, function(rep) {
    if (is.null(rep) || is.null(rep[[m]])) return(NA)
    if (rep[[m]]$failed) return(NA)
    return(rep[[m]]$converged)
  })
  
  idx <- which(convergence_summary$model == m)
  convergence_summary$n_total[idx]     <- sum(!is.na(conv_flags))
  convergence_summary$n_converged[idx] <- sum(conv_flags, na.rm = TRUE)
  convergence_summary$n_failed[idx]    <- sum(is.na(conv_flags))
}

convergence_summary$pct_converged <- round(
  100 * convergence_summary$n_converged / convergence_summary$n_total, 1)

cat("--- Convergence Summary ---\n")
print(convergence_summary)


# ---- 8b. Bias and MSE by agent and model ----
bias_mse_results <- list()

for (m in model_names) {
  bias_mat <- do.call(rbind, lapply(all_metrics, function(rep) {
    if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
    if (!rep[[m]]$converged) return(NULL)
    return(t(rep[[m]]$bias))
  }))
  
  mse_mat <- do.call(rbind, lapply(all_metrics, function(rep) {
    if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
    if (!rep[[m]]$converged) return(NULL)
    return(t(rep[[m]]$sq_error))
  }))
  
  if (!is.null(bias_mat) && nrow(bias_mat) > 0) {
    bias_mse_results[[m]] <- data.frame(
      model    = m,
      agent    = agent_names,
      bias     = round(colMeans(bias_mat), 4),
      bias_se  = round(apply(bias_mat, 2, sd) / sqrt(nrow(bias_mat)), 4),
      mse      = round(colMeans(mse_mat), 4),
      mse_se   = round(apply(mse_mat, 2, sd) / sqrt(nrow(mse_mat)), 4),
      n_reps   = nrow(bias_mat),
      stringsAsFactors = FALSE
    )
  }
}

bias_mse_table <- bind_rows(bias_mse_results)

cat("\n--- Bias and MSE by Agent and Model ---\n")
print(bias_mse_table, digits = 4)


# ---- 8c. Coverage by agent and model ----
coverage_results <- list()

for (m in model_names) {
  cov_mat <- do.call(rbind, lapply(all_metrics, function(rep) {
    if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
    if (!rep[[m]]$converged) return(NULL)
    return(t(rep[[m]]$coverage))
  }))
  
  if (!is.null(cov_mat) && nrow(cov_mat) > 0) {
    coverage_results[[m]] <- data.frame(
      model    = m,
      agent    = agent_names,
      coverage = round(colMeans(cov_mat), 3),
      n_reps   = nrow(cov_mat),
      stringsAsFactors = FALSE
    )
  }
}

coverage_table <- bind_rows(coverage_results)

cat("\n--- Coverage by Agent and Model ---\n")
print(coverage_table)


# ---- 8d. Ranking accuracy (Spearman) ----
spearman_results <- list()

for (m in model_names) {
  rho_vals <- sapply(all_metrics, function(rep) {
    if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NA)
    if (!rep[[m]]$converged) return(NA)
    return(rep[[m]]$spearman)
  })
  rho_vals <- rho_vals[!is.na(rho_vals)]
  
  spearman_results[[m]] <- data.frame(
    model  = m,
    median = round(median(rho_vals), 3),
    mean   = round(mean(rho_vals), 3),
    q25    = round(quantile(rho_vals, 0.25), 3),
    q75    = round(quantile(rho_vals, 0.75), 3),
    n_reps = length(rho_vals),
    stringsAsFactors = FALSE
  )
}

spearman_table <- bind_rows(spearman_results)

cat("\n--- Ranking Accuracy (Spearman Correlation) ---\n")
print(spearman_table)


# ---- 8d2. Mean Absolute Rank Error by Agent ----
rank_error_results <- list()

for (m in model_names) {
  est_rank_mat <- do.call(rbind, lapply(all_metrics, function(rep) {
    if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
    if (!rep[[m]]$converged) return(NULL)
    return(t(rep[[m]]$est_ranks))
  }))
  
  if (!is.null(est_rank_mat) && nrow(est_rank_mat) > 0) {
    rank_err_mat <- sweep(est_rank_mat, 2, true_params$true_rank, "-")
    rank_error_results[[m]] <- data.frame(
      model          = m,
      agent          = agent_names,
      true_rank      = true_params$true_rank,
      mean_est_rank  = round(colMeans(est_rank_mat), 2),
      mean_abs_error = round(colMeans(abs(rank_err_mat)), 2),
      n_reps         = nrow(rank_err_mat),
      stringsAsFactors = FALSE
    )
  }
}

rank_error_table <- bind_rows(rank_error_results)

cat("\n--- Mean Absolute Rank Error by Agent ---\n")
print(rank_error_table)

cat("\n--- Mean Estimated Rank vs True Rank ---\n")
rank_error_wide <- rank_error_table %>%
  select(model, agent, true_rank, mean_est_rank) %>%
  pivot_wider(names_from = model, values_from = mean_est_rank,
              names_prefix = "est_rank_")
print(rank_error_wide)


# ---- 8d3. Per-agent rank match rate ----
# For each replication, the estimated integer rank of each agent is derived
# from est_ranks (rank of the mean posterior ranks). "Match" means the
# integer rank equals the true rank exactly.
rank_match_results <- list()

for (m in model_names) {
  match_mat <- do.call(rbind, lapply(all_metrics, function(rep) {
    if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
    if (!rep[[m]]$converged) return(NULL)
    est_int_rank <- base::rank(rep[[m]]$est_ranks, ties.method = "min")
    return(t(as.integer(est_int_rank == true_params$true_rank)))
  }))
  
  if (!is.null(match_mat) && nrow(match_mat) > 0) {
    rank_match_results[[m]] <- data.frame(
      model      = m,
      agent      = agent_names,
      true_rank  = true_params$true_rank,
      match_rate = round(colMeans(match_mat), 3),
      n_reps     = nrow(match_mat),
      stringsAsFactors = FALSE
    )
  }
}

rank_match_table <- bind_rows(rank_match_results)

cat("\n--- Per-Agent Rank Match Rate ---\n")
cat("(Proportion of replications where estimated rank matches true rank exactly)\n")
print(rank_match_table)


# ---- 8d4. SUCRA by agent and model ----
sucra_results <- list()

for (m in model_names) {
  sucra_mat <- do.call(rbind, lapply(all_metrics, function(rep) {
    if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
    if (!rep[[m]]$converged) return(NULL)
    return(t(rep[[m]]$sucra))
  }))
  
  if (!is.null(sucra_mat) && nrow(sucra_mat) > 0) {
    sucra_results[[m]] <- data.frame(
      model     = m,
      agent     = agent_names,
      true_rank = true_params$true_rank,
      sucra     = round(colMeans(sucra_mat), 3),
      sucra_sd  = round(apply(sucra_mat, 2, sd), 3),
      n_reps    = nrow(sucra_mat),
      stringsAsFactors = FALSE
    )
  }
}

sucra_table <- bind_rows(sucra_results)

cat("\n--- SUCRA by Agent and Model (averaged across replications) ---\n")
print(sucra_table)


# ---- 8e. Model selection distribution (DIC and LPML) ----
dic_selected  <- character(0)
lpml_selected <- character(0)

for (rep in seq_along(all_metrics)) {
  if (is.null(all_metrics[[rep]])) next
  
  all_ok <- all(sapply(model_names, function(m) {
    !is.null(all_metrics[[rep]][[m]]) &&
      !all_metrics[[rep]][[m]]$failed &&
      all_metrics[[rep]][[m]]$converged
  }))
  
  if (!all_ok) next
  
  # DIC: compare only dose-response models (Lumped fits different data)
  dic_models <- c("A", "B", "Ref")
  dics <- sapply(dic_models, function(m) all_metrics[[rep]][[m]]$dic)
  dic_selected <- c(dic_selected, dic_models[which.min(dics)])
  
  # LPML: compare only models that have it (exclude Lumped)
  lpml_models <- c("A", "B", "Ref")
  lpmls <- sapply(lpml_models, function(m) all_metrics[[rep]][[m]]$lpml)
  if (!any(is.na(lpmls))) {
    lpml_selected <- c(lpml_selected, lpml_models[which.max(lpmls)])
  }
}

cat("\n--- Model Selection Distribution ---\n")
cat("(Model B matches DGP but is most parameterized; selection is not guaranteed)\n\n")

cat("DIC selection:\n")
print(table(dic_selected))
cat(sprintf("  B selected: %d/%d (%.1f%%)\n",
            sum(dic_selected == "B"), length(dic_selected),
            100 * mean(dic_selected == "B")))

cat("\nLPML selection:\n")
print(table(lpml_selected))
cat(sprintf("  B selected: %d/%d (%.1f%%)\n",
            sum(lpml_selected == "B"), length(lpml_selected),
            100 * mean(lpml_selected == "B")))

concordance <- mean(dic_selected == lpml_selected)
cat(sprintf("\nDIC-LPML concordance: %.1f%%\n", 100 * concordance))


# ==============================================================================
# 9. Figures
# ==============================================================================

bias_mse_table$model <- factor(bias_mse_table$model, levels = c("A", "B", "Ref", "Lumped"))

# ---- 9a. Bias by agent and model ----
p_bias <- ggplot(bias_mse_table, aes(x = agent, y = bias, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbar(aes(ymin = bias - 1.96 * bias_se, ymax = bias + 1.96 * bias_se),
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c("A" = "#E41A1C", "B" = "#377EB8",
                               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
                    labels = c("A" = "Model A (Shared Emax)",
                               "B" = "Model B (Agent-specific Emax)",
                               "Ref" = "Ref (Linear)",
                               "Lumped" = "Lumped NMA")) +
  scale_x_discrete(labels = toupper) +
  labs(title = "Bias in Treatment Effect Estimation at Max Observed Dose",
       subtitle = paste0("Across ", n_reps, " simulation replications"),
       x = NULL, y = "Bias (log-odds scale)", fill = "Model") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_bias)
ggsave("figures/Figure 4.2.jpeg", p_bias, width = 9, height = 5, dpi = 300)


# ---- 9b. Coverage by agent and model ----
p_coverage <- ggplot(coverage_table, aes(x = agent, y = coverage, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("A" = "#E41A1C", "B" = "#377EB8",
                               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
                    labels = c("A" = "Model A (Shared Emax)",
                               "B" = "Model B (Agent-specific Emax)",
                               "Ref" = "Ref (Linear)",
                               "Lumped" = "Lumped NMA")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "95% Credible Interval Coverage by Agent",
       subtitle = paste0("Across ", n_reps, " simulation replications (nominal = 0.95)"),
       x = NULL, y = "Coverage", fill = "Model") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_coverage)

# ---- 9c. Spearman rank correlation boxplot ----
spearman_all <- do.call(rbind, lapply(model_names, function(m) {
  rho_vals <- sapply(all_metrics, function(rep) {
    if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NA)
    if (!rep[[m]]$converged) return(NA)
    return(rep[[m]]$spearman)
  })
  data.frame(model = m, spearman = rho_vals[!is.na(rho_vals)],
             stringsAsFactors = FALSE)
}))

p_spearman <- ggplot(spearman_all, aes(x = model, y = spearman, fill = model)) +
  geom_boxplot(width = 0.5, outlier.shape = 21) +
  scale_fill_manual(values = c("A" = "#E41A1C", "B" = "#377EB8",
                               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
                    labels = c("A" = "Model A (Shared Emax)",
                               "B" = "Model B (Agent-specific Emax)",
                               "Ref" = "Ref (Linear)",
                               "Lumped" = "Lumped NMA")) +
  scale_x_discrete(labels = c("A" = "Model A\n(Shared Emax)",
                              "B" = "Model B\n(Agent-specific)",
                              "Ref" = "Ref\n(Linear)",
                              "Lumped" = "Lumped\nNMA")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(title = "Ranking Accuracy Across Simulation Replications",
       subtitle = "Spearman correlation between estimated and true agent ranks",
       x = NULL, y = "Spearman \u03C1") +
  guides(fill = "none")

print(p_spearman)


# ---- 9d. Model selection barplot ----
dic_model_names <- c("A", "B", "Ref")

selection_df <- data.frame(
  criterion = c(rep("DIC", length(dic_model_names)), rep("LPML", length(dic_model_names))),
  selected  = rep(dic_model_names, 2),
  count     = c(
    sapply(dic_model_names, function(m) sum(dic_selected == m)),
    sapply(dic_model_names, function(m) sum(lpml_selected == m))
  ),
  stringsAsFactors = FALSE
)

p_selection <- ggplot(selection_df, aes(x = selected, y = count, fill = selected)) +
  geom_col(width = 0.6) +
  facet_wrap(~criterion) +
  scale_fill_manual(values = c("A" = "#E41A1C", "B" = "#377EB8", "Ref" = "#4DAF4A")) +
  labs(title = "Model Selection Across Replications",
       subtitle = "Number of replications selecting each model as preferred",
       x = "Model Selected", y = "Count") +
  guides(fill = "none")

print(p_selection)
ggsave("figures/Figure 4.3.jpeg", p_selection, width = 8, height = 5, dpi = 300)


# ==============================================================================
# 10. Summary Output for Results Section
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  SIMULATION RESULTS SUMMARY\n")
cat("================================================================\n\n")

cat("--- Convergence ---\n")
print(convergence_summary)

cat("\n--- Bias and MSE ---\n")
print(bias_mse_table, digits = 4)

cat("\n--- Coverage ---\n")
print(coverage_table)

cat("\n--- Ranking Accuracy (Spearman) ---\n")
print(spearman_table)

cat("\n--- Per-Agent Rank Match Rate ---\n")
print(rank_match_table)

cat("\n--- SUCRA ---\n")
print(sucra_table)

cat("\n--- Model Selection ---\n")
cat("DIC:  B selected", round(100 * mean(dic_selected == "B"), 1), "%\n")
cat("LPML: B selected", round(100 * mean(lpml_selected == "B"), 1), "%\n")
cat("Concordance:", round(100 * concordance, 1), "%\n")

cat("\n================================================================\n")
cat("  Simulation analysis complete.\n")
cat("  Results saved to sim_results_final.RData\n")
cat("  Figures saved to figures/\n")
cat("================================================================\n")


# ==============================================================================
# 11. Save Results
# ==============================================================================

dir.create("results/simulation", recursive = TRUE, showWarnings = FALSE)

write.csv(convergence_summary, "results/simulation/convergence_summary.csv", row.names = FALSE)
write.csv(bias_mse_table,      "results/simulation/bias_mse_by_agent.csv",   row.names = FALSE)
write.csv(coverage_table,      "results/simulation/coverage_by_agent.csv",   row.names = FALSE)
write.csv(spearman_table,      "results/simulation/spearman_summary.csv",    row.names = FALSE)
write.csv(rank_error_table,    "results/simulation/rank_error_by_agent.csv", row.names = FALSE)
write.csv(rank_error_wide,     "results/simulation/rank_comparison_wide.csv", row.names = FALSE)
write.csv(rank_match_table,    "results/simulation/rank_match_by_agent.csv", row.names = FALSE)
write.csv(sucra_table,         "results/simulation/sucra_by_agent.csv",      row.names = FALSE)

selection_summary <- data.frame(
  criterion = c(rep("DIC", length(table(dic_selected))),
                rep("LPML", length(table(lpml_selected)))),
  model     = c(names(table(dic_selected)), names(table(lpml_selected))),
  count     = c(as.integer(table(dic_selected)), as.integer(table(lpml_selected)))
)
write.csv(selection_summary, "results/simulation/model_selection_distribution.csv", row.names = FALSE)

write.csv(true_params, "results/simulation/simulation_true_params.csv", row.names = FALSE)

save(all_metrics, true_params, true_tau,
     bias_mse_table, coverage_table, spearman_table,
     rank_error_table, rank_match_table, sucra_table,
     convergence_summary,
     dic_selected, lpml_selected,
     file = "results/simulation/all_results.RData")

cat("\n=== All results saved to results/simulation/ ===\n")

writeLines(capture.output(sessionInfo()), "results/simulation/ssri_simulation.txt")