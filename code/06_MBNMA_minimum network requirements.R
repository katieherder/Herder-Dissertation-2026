################################################################################
# Chapter 5: Exploratory Simulations — Boundary Conditions for Reliable Rankings
#
# This script extends the Chapter 4 simulation framework to ask:
#   "Under what conditions can MBNMA rankings be trusted?"
#
# Two dimensions are varied ONE AT A TIME, holding others at empirical values:
#   5.4.2  Heterogeneity magnitude:  τ from 0.20 → 0.01 (0.27 is in Ch4)
#   5.4.3  Evidence density:         2× studies, 3× studies, and dose-enriched
#
# BACKBONE: Same SSRI network structure as Chapter 4.
# TRUTH:    Same Model B posterior medians as Chapter 4 DGP.
# MODELS:   Model A, Model B, Ref (Linear), Lumped NMA — same as Chapter 4.
# METRICS:  Spearman ρ, bias, coverage — same as Chapter 4.
#
# PREREQUISITES:
#   - The Chapter 4 simulation script must have been sourced or run first,
#     so that the following objects exist in the environment:
#       network_template, true_params, true_tau, true_mu_mean, true_mu_sd,
#       dose_equiv, generate_ssri_dataset(), fit_all_models(),
#       extract_metrics(), compute_lpml()
#   - If running standalone, source the Chapter 4 script through Section 6
#     (stopping before the replication loop in Section 7).
#
# OUTPUT:
#   results/adhoc/   — CSVs and .RData per simulation dimension
#   results/figures/ — Key figures for each dimension
#
################################################################################

# ==============================================================================
# 0. Setup and Verification
# ==============================================================================

stopifnot(
  exists("network_template"),
  exists("true_params"),
  exists("generate_ssri_dataset"),
  exists("fit_all_models"),
  exists("extract_metrics")
)

dir.create("results/adhoc", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

model_names <- c("A", "B", "Ref", "Lumped")
agent_names <- true_params$agent

cat("\n================================================================\n")
cat("  CHAPTER 5: EXPLORATORY SIMULATIONS\n")
cat("  Backbone: SSRI network (Chapter 4)\n")
cat("================================================================\n")


# ==============================================================================
# HELPER: Run one simulation scenario and compile results
# ==============================================================================

run_scenario <- function(scenario_label,
                         template_scenario   = network_template,
                         true_params_scenario = true_params,
                         true_tau_scenario    = true_tau,
                         n_reps              = 50,
                         seed_offset         = 0) {
  
  cat(sprintf("\n\n============================================================\n"))
  cat(sprintf("  SCENARIO: %s\n", scenario_label))
  cat(sprintf("  tau = %.3f | agents = %d | studies = %d | arms = %d\n",
              true_tau_scenario,
              length(unique(template_scenario$agent[template_scenario$agent != "placebo"])),
              length(unique(template_scenario$studyID)),
              nrow(template_scenario)))
  cat(sprintf("  Replications: %d\n", n_reps))
  cat(sprintf("============================================================\n"))
  
  all_metrics_scen <- vector("list", n_reps)
  
  for (rep in 1:n_reps) {
    rep_seed <- 200000 + seed_offset + rep
    
    cat(sprintf("\n--- %s | Rep %d/%d (seed: %d) ---\n",
                scenario_label, rep, n_reps, rep_seed))
    
    sim_data <- generate_ssri_dataset(
      template     = template_scenario,
      true_params  = true_params_scenario,
      true_tau     = true_tau_scenario,
      true_mu_mean = true_mu_mean,
      true_mu_sd   = true_mu_sd,
      resample_n   = TRUE,
      seed         = rep_seed
    )
    
    fitted <- fit_all_models(sim_data, rep_seed)
    all_metrics_scen[[rep]] <- extract_metrics(fitted, true_params_scenario)
    
    if (!is.null(all_metrics_scen[[rep]])) {
      for (m in model_names) {
        if (!is.null(all_metrics_scen[[rep]][[m]]) &&
            !all_metrics_scen[[rep]][[m]]$failed) {
          cat(sprintf("  %s: rho=%.3f Rhat_max=%.3f\n",
                      m,
                      all_metrics_scen[[rep]][[m]]$spearman,
                      all_metrics_scen[[rep]][[m]]$max_rhat))
        } else {
          cat(sprintf("  %s: FAILED\n", m))
        }
      }
    }
    
    if (rep %% 10 == 0) {
      save(all_metrics_scen,
           file = sprintf("results/adhoc/checkpoint_%s_rep%d.RData",
                          gsub("[^a-zA-Z0-9]", "_", scenario_label), rep))
    }
  }
  
  # --- Compile summary metrics ---
  
  spearman_by_model <- lapply(model_names, function(m) {
    rho_vals <- sapply(all_metrics_scen, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NA)
      if (!rep[[m]]$converged) return(NA)
      return(rep[[m]]$spearman)
    })
    rho_vals <- rho_vals[!is.na(rho_vals)]
    data.frame(
      model  = m,
      median_rho = ifelse(length(rho_vals) > 0, median(rho_vals), NA),
      mean_rho   = ifelse(length(rho_vals) > 0, mean(rho_vals), NA),
      q25        = ifelse(length(rho_vals) > 0, quantile(rho_vals, 0.25), NA),
      q75        = ifelse(length(rho_vals) > 0, quantile(rho_vals, 0.75), NA),
      n_reps     = length(rho_vals),
      stringsAsFactors = FALSE
    )
  })
  spearman_summary <- bind_rows(spearman_by_model)
  
  bias_by_model <- lapply(model_names, function(m) {
    bias_mat <- do.call(rbind, lapply(all_metrics_scen, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (!rep[[m]]$converged) return(NULL)
      return(t(rep[[m]]$bias))
    }))
    if (!is.null(bias_mat) && nrow(bias_mat) > 0) {
      data.frame(
        model  = m,
        agent  = agent_names,
        bias   = round(colMeans(bias_mat), 4),
        n_reps = nrow(bias_mat),
        stringsAsFactors = FALSE
      )
    } else NULL
  })
  bias_summary <- bind_rows(bias_by_model)
  
  cov_by_model <- lapply(model_names, function(m) {
    cov_mat <- do.call(rbind, lapply(all_metrics_scen, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (!rep[[m]]$converged) return(NULL)
      return(t(rep[[m]]$coverage))
    }))
    if (!is.null(cov_mat) && nrow(cov_mat) > 0) {
      data.frame(
        model    = m,
        agent    = agent_names,
        coverage = round(colMeans(cov_mat), 3),
        n_reps   = nrow(cov_mat),
        stringsAsFactors = FALSE
      )
    } else NULL
  })
  coverage_summary <- bind_rows(cov_by_model)
  
  conv_summary <- lapply(model_names, function(m) {
    conv_flags <- sapply(all_metrics_scen, function(rep) {
      if (is.null(rep) || is.null(rep[[m]])) return(NA)
      if (rep[[m]]$failed) return(NA)
      return(rep[[m]]$converged)
    })
    data.frame(
      model       = m,
      n_converged = sum(conv_flags, na.rm = TRUE),
      n_total     = sum(!is.na(conv_flags)),
      stringsAsFactors = FALSE
    )
  })
  conv_summary <- bind_rows(conv_summary)
  
  return(list(
    scenario_label   = scenario_label,
    all_metrics      = all_metrics_scen,
    spearman_summary = spearman_summary,
    bias_summary     = bias_summary,
    coverage_summary = coverage_summary,
    conv_summary     = conv_summary
  ))
}


################################################################################
#
# SIMULATION 5.4.2: EFFECT OF HETEROGENEITY MAGNITUDE
#
# NOTE: τ = 0.27 results are from Chapter 4 (already completed).
# This simulation covers τ = 0.20, 0.15, 0.10, 0.05, 0.01.
# 50 reps per scenario × 4 models × 5 τ levels = 1,000 model fits.
#
################################################################################

cat("\n\n################################################################\n")
cat("  SIMULATION 5.4.2: HETEROGENEITY MAGNITUDE\n")
cat("################################################################\n")
run_tau_sim <- FALSE   # Already completed — set TRUE to rerun

if (run_tau_sim) {

tau_values <- c(0.20, 0.15, 0.10, 0.05, 0.01)
n_reps_tau <- 50

tau_results <- list()

for (i in seq_along(tau_values)) {
  tau_val <- tau_values[i]
  label <- sprintf("tau_%.2f", tau_val)
  
  tau_results[[label]] <- run_scenario(
    scenario_label    = label,
    template_scenario = network_template,
    true_params_scenario = true_params,
    true_tau_scenario = tau_val,
    n_reps            = n_reps_tau,
    seed_offset       = i * 1000
  )
}

# --- Compile across τ levels ---
tau_spearman_compiled <- bind_rows(lapply(seq_along(tau_values), function(i) {
  label <- sprintf("tau_%.2f", tau_values[i])
  res <- tau_results[[label]]$spearman_summary
  res$tau <- tau_values[i]
  return(res)
}))

tau_bias_compiled <- bind_rows(lapply(seq_along(tau_values), function(i) {
  label <- sprintf("tau_%.2f", tau_values[i])
  res <- tau_results[[label]]$bias_summary
  res$tau <- tau_values[i]
  return(res)
}))

tau_coverage_compiled <- bind_rows(lapply(seq_along(tau_values), function(i) {
  label <- sprintf("tau_%.2f", tau_values[i])
  res <- tau_results[[label]]$coverage_summary
  res$tau <- tau_values[i]
  return(res)
}))

# Save
write.csv(tau_spearman_compiled,
          "results/adhoc/tau_spearman.csv", row.names = FALSE)
write.csv(tau_bias_compiled,
          "results/adhoc/tau_bias.csv", row.names = FALSE)
write.csv(tau_coverage_compiled,
          "results/adhoc/tau_coverage.csv", row.names = FALSE)

save(tau_results, tau_spearman_compiled, tau_bias_compiled, tau_coverage_compiled,
     file = "results/adhoc/sim_542_heterogeneity.RData")

# --- Figure 5.4.2: Spearman ρ vs τ ---
p_tau <- ggplot(tau_spearman_compiled,
                aes(x = tau, y = median_rho, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = model),
              alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.7, linetype = "dotted", color = "grey30") +
  annotate("text", x = max(tau_values), y = 0.73,
           label = "ρ = 0.7 threshold", hjust = 1, size = 3, color = "grey30") +
  scale_x_reverse(breaks = tau_values) +
  scale_color_manual(
    values = c("A" = "#E41A1C", "B" = "#377EB8",
               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
    labels = c("A" = "Model A (Shared Emax)",
               "B" = "Model B (Agent-specific)",
               "Ref" = "Linear", "Lumped" = "Lumped NMA")) +
  scale_fill_manual(
    values = c("A" = "#E41A1C", "B" = "#377EB8",
               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
    guide = "none") +
  labs(
    title = "Ranking Accuracy by Between-Study Heterogeneity",
    subtitle = "SSRI network structure; all other parameters at empirical values",
    x = expression(paste("Between-study heterogeneity (", tau, ")")),
    y = expression(paste("Median Spearman ", rho)),
    color = "Model"
  ) +
  theme_diss

print(p_tau)
ggsave("results/figures/fig_542_tau_vs_rho.jpeg", p_tau,
       width = 9, height = 5.5, dpi = 300)

cat("\n=== Simulation 5.4.2 complete ===\n")

} else {
  cat("\n=== Skipping tau simulation (already completed) ===\n")
  cat("  Loading saved results...\n")
  load("results/adhoc/sim_542_heterogeneity.RData")
}


################################################################################
#
# SIMULATION 5.4.3: EFFECT OF EVIDENCE DENSITY
#
# Three scenarios, all at τ = 0.27 (empirical):
#   (a) 2× studies: duplicate all studies at same dose levels
#   (b) 3× studies: triplicate all studies at same dose levels
#   (c) Dose-enriched: add midpoint doses between each existing pair
#
# Scenarios (a) and (b) test whether more data at the SAME doses helps.
# Scenario (c) tests whether more dose LEVELS (better curve coverage) helps.
#
# NOTE: Running 2 reps as placeholder. Increase n_reps_density for full run.
#
################################################################################

cat("\n\n################################################################\n")
cat("  SIMULATION 5.4.3: EVIDENCE DENSITY\n")
cat("################################################################\n")

n_reps_density <- 50   # PLACEHOLDER — increase to 50 for final run


# --- Helper: expand template by duplicating studies ---

expand_template <- function(template, multiplier) {
  if (multiplier == 1) return(template)
  
  original_studies <- unique(template$studyID)
  expanded <- template
  
  for (mult in 2:multiplier) {
    new_block <- template
    study_map <- setNames(
      paste0(original_studies, "_x", mult),
      original_studies
    )
    new_block$studyID <- study_map[as.character(new_block$studyID)]
    expanded <- rbind(expanded, new_block)
  }
  
  return(expanded)
}


# --- Helper: add midpoint doses between each existing pair per agent ---

add_midpoint_doses <- function(template) {
  agents <- unique(template$agent[template$agent != "placebo"])
  new_rows <- list()
  study_counter <- 9000
  
  for (ag in agents) {
    agent_doses <- sort(unique(template$dose_equiv[template$agent == ag]))
    agent_n <- median(template$n[template$agent == ag])
    
    for (j in 1:(length(agent_doses) - 1)) {
      midpoint <- (agent_doses[j] + agent_doses[j + 1]) / 2
      study_counter <- study_counter + 1
      
      new_rows[[length(new_rows) + 1]] <- data.frame(
        studyID = study_counter, agent = "placebo", dose_equiv = 0,
        r = NA_integer_, n = as.integer(round(agent_n)),
        stringsAsFactors = FALSE
      )
      new_rows[[length(new_rows) + 1]] <- data.frame(
        studyID = study_counter, agent = ag, dose_equiv = midpoint,
        r = NA_integer_, n = as.integer(round(agent_n)),
        stringsAsFactors = FALSE
      )
    }
  }
  
  new_df <- bind_rows(new_rows)
  expanded <- bind_rows(template, new_df)
  
  # Fix arm numbers and placeholder r values for new rows
  expanded <- expanded %>%
    group_by(studyID) %>%
    mutate(arm = row_number()) %>%
    ungroup() %>%
    mutate(r = ifelse(is.na(r), 0L, as.integer(r)))
  
  return(expanded)
}

# --- Run all three density scenarios ---

density_results <- list()

# (a) 2× studies
cat("\n--- Density: 2× studies ---\n")
template_2x <- expand_template(network_template, 2)
cat(sprintf("  %d studies, %d arms\n",
            length(unique(template_2x$studyID)), nrow(template_2x)))

density_results[["density_2x"]] <- run_scenario(
  scenario_label    = "density_2x",
  template_scenario = template_2x,
  true_params_scenario = true_params,
  true_tau_scenario = true_tau,
  n_reps            = n_reps_density,
  seed_offset       = 11000
)

# (b) 3× studies
cat("\n--- Density: 3× studies ---\n")
template_3x <- expand_template(network_template, 3)
cat(sprintf("  %d studies, %d arms\n",
            length(unique(template_3x$studyID)), nrow(template_3x)))

density_results[["density_3x"]] <- run_scenario(
  scenario_label    = "density_3x",
  template_scenario = template_3x,
  true_params_scenario = true_params,
  true_tau_scenario = true_tau,
  n_reps            = n_reps_density,
  seed_offset       = 12000
)

# (c) Dose-enriched (midpoint doses added)
cat("\n--- Density: dose-enriched ---\n")
template_enriched <- add_midpoint_doses(network_template)
cat(sprintf("  %d studies, %d arms\n",
            length(unique(template_enriched$studyID)), nrow(template_enriched)))

# Print new dose structure for verification
cat("\n  Dose structure after enrichment:\n")
template_enriched %>%
  filter(agent != "placebo") %>%
  group_by(agent) %>%
  summarise(
    n_doses = n_distinct(dose_equiv),
    n_studies = n_distinct(studyID),
    doses = paste(round(sort(unique(dose_equiv)), 1), collapse = ", "),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  print()

density_results[["density_enriched"]] <- run_scenario(
  scenario_label    = "density_enriched",
  template_scenario = template_enriched,
  true_params_scenario = true_params,
  true_tau_scenario = true_tau,
  n_reps            = n_reps_density,
  seed_offset       = 50000
)


# --- Compile across density scenarios ---

density_labels <- c("density_2x", "density_3x", "density_enriched")
density_display <- c("2× studies", "3× studies", "Dose-enriched")

density_spearman_compiled <- bind_rows(lapply(seq_along(density_labels), function(i) {
  res <- density_results[[density_labels[i]]]$spearman_summary
  res$scenario <- density_display[i]
  res$scenario_id <- density_labels[i]
  return(res)
}))

density_bias_compiled <- bind_rows(lapply(seq_along(density_labels), function(i) {
  res <- density_results[[density_labels[i]]]$bias_summary
  res$scenario <- density_display[i]
  return(res)
}))

density_coverage_compiled <- bind_rows(lapply(seq_along(density_labels), function(i) {
  res <- density_results[[density_labels[i]]]$coverage_summary
  res$scenario <- density_display[i]
  return(res)
}))

# Save
write.csv(density_spearman_compiled,
          "results/adhoc/density_spearman.csv", row.names = FALSE)
write.csv(density_bias_compiled,
          "results/adhoc/density_bias.csv", row.names = FALSE)
write.csv(density_coverage_compiled,
          "results/adhoc/density_coverage.csv", row.names = FALSE)

save(density_results, density_spearman_compiled,
     density_bias_compiled, density_coverage_compiled,
     template_enriched,
     file = "results/adhoc/sim_543_density.RData")

# --- Figure 5.4.3: Spearman ρ by density scenario ---
density_spearman_compiled$scenario <- factor(
  density_spearman_compiled$scenario,
  levels = density_display
)

p_density <- ggplot(density_spearman_compiled,
                    aes(x = scenario, y = median_rho, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.7, linetype = "dotted", color = "grey30") +
  scale_fill_manual(
    values = c("A" = "#E41A1C", "B" = "#377EB8",
               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
    labels = c("A" = "Model A (Shared Emax)",
               "B" = "Model B (Agent-specific)",
               "Ref" = "Linear", "Lumped" = "Lumped NMA")) +
  labs(
    title = "Ranking Accuracy by Evidence Density",
    subtitle = expression(paste(tau, " = 0.27 (empirical); baseline = 60 studies / 145 arms")),
    x = "Density scenario",
    y = expression(paste("Median Spearman ", rho)),
    fill = "Model"
  ) +
  theme_diss

print(p_density)
ggsave("results/figures/fig_543_density_vs_rho.jpeg", p_density,
       width = 9, height = 5.5, dpi = 300)

cat("\n=== Simulation 5.4.3 complete ===\n")


################################################################################
#
# COMBINED SUMMARY AND COMPILATION
#
################################################################################

cat("\n\n################################################################\n")
cat("  CHAPTER 5: COMPILING ALL RESULTS\n")
cat("################################################################\n\n")

# Load Chapter 4 base results for tau = 0.27 baseline
load("results/simulation/simulation_all_results.RData")
all_metrics_ch4 <- all_metrics
rm(all_metrics)

all_tau_values <- c(0.27, 0.20, 0.15, 0.10, 0.05, 0.01)

# ==============================================================================
# 1. SPEARMAN BY TAU (R-HAT < 1.10, INCLUDING CH4 BASELINE)
# ==============================================================================

cat("  1. Compiling Spearman by tau...\n")

tau_spearman_full <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) {
    res_list <- all_metrics_ch4
  } else {
    label <- sprintf("tau_%.2f", tv)
    res_list <- tau_results[[label]]$all_metrics
  }
  
  bind_rows(lapply(model_names, function(m) {
    rho_vals <- sapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NA)
      if (rep[[m]]$max_rhat >= 1.10) return(NA)
      return(rep[[m]]$spearman)
    })
    rho_vals <- rho_vals[!is.na(rho_vals)]
    data.frame(
      tau = tv, model = m,
      median_rho = ifelse(length(rho_vals) > 0, median(rho_vals), NA),
      mean_rho = ifelse(length(rho_vals) > 0, mean(rho_vals), NA),
      q25 = ifelse(length(rho_vals) > 0, quantile(rho_vals, 0.25), NA),
      q75 = ifelse(length(rho_vals) > 0, quantile(rho_vals, 0.75), NA),
      n_reps = length(rho_vals),
      stringsAsFactors = FALSE
    )
  }))
}))

# ==============================================================================
# 2. CONVERGENCE AT BOTH THRESHOLDS
# ==============================================================================

cat("  2. Compiling convergence...\n")

tau_convergence <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) {
    res_list <- all_metrics_ch4
  } else {
    label <- sprintf("tau_%.2f", tv)
    res_list <- tau_results[[label]]$all_metrics
  }
  
  bind_rows(lapply(model_names, function(m) {
    flags <- sapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(c(NA, NA))
      rhat <- rep[[m]]$max_rhat
      return(c(rhat < 1.05, rhat < 1.10))
    })
    if (is.matrix(flags)) {
      data.frame(
        tau = tv, model = m,
        n_converged_05 = sum(flags[1,], na.rm = TRUE),
        n_converged_10 = sum(flags[2,], na.rm = TRUE),
        n_total = sum(!is.na(flags[1,])),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(tau = tv, model = m, n_converged_05 = NA,
                 n_converged_10 = NA, n_total = NA, stringsAsFactors = FALSE)
    }
  }))
}))

# ==============================================================================
# 3. BIAS BY TAU, AGENT, AND MODEL (R-HAT < 1.10)
# ==============================================================================

cat("  3. Compiling bias...\n")

tau_bias_full <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) {
    res_list <- all_metrics_ch4
  } else {
    label <- sprintf("tau_%.2f", tv)
    res_list <- tau_results[[label]]$all_metrics
  }
  
  bind_rows(lapply(model_names, function(m) {
    bias_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.10) return(NULL)
      return(t(rep[[m]]$bias))
    }))
    if (!is.null(bias_mat) && nrow(bias_mat) > 0) {
      data.frame(
        tau = tv, model = m, agent = agent_names,
        bias = round(colMeans(bias_mat), 4),
        n_reps = nrow(bias_mat),
        stringsAsFactors = FALSE
      )
    } else NULL
  }))
}))

# ==============================================================================
# 4. COVERAGE BY TAU, AGENT, AND MODEL (R-HAT < 1.10)
# ==============================================================================

cat("  4. Compiling coverage...\n")

tau_coverage_full <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) {
    res_list <- all_metrics_ch4
  } else {
    label <- sprintf("tau_%.2f", tv)
    res_list <- tau_results[[label]]$all_metrics
  }
  
  bind_rows(lapply(model_names, function(m) {
    cov_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.10) return(NULL)
      return(t(rep[[m]]$coverage))
    }))
    if (!is.null(cov_mat) && nrow(cov_mat) > 0) {
      data.frame(
        tau = tv, model = m, agent = agent_names,
        coverage = round(colMeans(cov_mat), 3),
        n_reps = nrow(cov_mat),
        stringsAsFactors = FALSE
      )
    } else NULL
  }))
}))

# ==============================================================================
# 5. CRI WIDTHS BY TAU, AGENT, AND MODEL (R-HAT < 1.10)
# ==============================================================================

cat("  5. Compiling CRI widths...\n")

tau_cri_full <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) {
    res_list <- all_metrics_ch4
  } else {
    label <- sprintf("tau_%.2f", tv)
    res_list <- tau_results[[label]]$all_metrics
  }
  
  bind_rows(lapply(model_names, function(m) {
    width_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.10) return(NULL)
      return(t(rep[[m]]$d_upper - rep[[m]]$d_lower))
    }))
    if (!is.null(width_mat) && nrow(width_mat) > 0) {
      data.frame(
        tau = tv, model = m, agent = agent_names,
        median_width = round(apply(width_mat, 2, median), 3),
        mean_width = round(colMeans(width_mat), 3),
        n_reps = nrow(width_mat),
        stringsAsFactors = FALSE
      )
    } else NULL
  }))
}))

# ==============================================================================
# 6. MEAN ESTIMATED RANKS BY TAU, AGENT, AND MODEL (R-HAT < 1.10)
# ==============================================================================

cat("  6. Compiling ranks...\n")

tau_ranks_full <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) {
    res_list <- all_metrics_ch4
  } else {
    label <- sprintf("tau_%.2f", tv)
    res_list <- tau_results[[label]]$all_metrics
  }
  
  bind_rows(lapply(model_names, function(m) {
    rank_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.10) return(NULL)
      return(t(rep[[m]]$est_ranks))
    }))
    if (!is.null(rank_mat) && nrow(rank_mat) > 0) {
      data.frame(
        tau = tv, model = m, agent = agent_names,
        true_rank = true_params$true_rank,
        mean_est_rank = round(colMeans(rank_mat), 2),
        n_reps = nrow(rank_mat),
        stringsAsFactors = FALSE
      )
    } else NULL
  }))
}))

# ==============================================================================
# 7. DENSITY: SPEARMAN (R-HAT < 1.10)
# ==============================================================================

cat("  7. Compiling density Spearman...\n")

density_labels_all <- names(density_results)

density_spearman_full <- bind_rows(lapply(density_labels_all, function(label) {
  res_list <- density_results[[label]]$all_metrics
  
  bind_rows(lapply(model_names, function(m) {
    rho_vals <- sapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NA)
      if (rep[[m]]$max_rhat >= 1.10) return(NA)
      return(rep[[m]]$spearman)
    })
    rho_vals <- rho_vals[!is.na(rho_vals)]
    data.frame(
      scenario = label, model = m,
      median_rho = ifelse(length(rho_vals) > 0, median(rho_vals), NA),
      mean_rho = ifelse(length(rho_vals) > 0, mean(rho_vals), NA),
      n_reps = length(rho_vals),
      stringsAsFactors = FALSE
    )
  }))
}))

# ==============================================================================
# 8. DENSITY: CRI WIDTHS (R-HAT < 1.10)
# ==============================================================================

cat("  8. Compiling density CRI widths...\n")

density_cri_full <- bind_rows(lapply(density_labels_all, function(label) {
  res_list <- density_results[[label]]$all_metrics
  
  bind_rows(lapply(model_names, function(m) {
    width_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.10) return(NULL)
      return(t(rep[[m]]$d_upper - rep[[m]]$d_lower))
    }))
    if (!is.null(width_mat) && nrow(width_mat) > 0) {
      data.frame(
        scenario = label, model = m, agent = agent_names,
        median_width = round(apply(width_mat, 2, median), 3),
        mean_width = round(colMeans(width_mat), 3),
        n_reps = nrow(width_mat),
        stringsAsFactors = FALSE
      )
    } else NULL
  }))
}))

# ==============================================================================
# 9. PRINT ALL SUMMARIES
# ==============================================================================

cat("\n================================================================\n")
cat("  RESULTS SUMMARY\n")
cat("================================================================\n\n")

cat("--- Spearman by tau (R-hat < 1.10) ---\n")
print(as.data.frame(tau_spearman_full))

cat("\n--- Convergence by tau ---\n")
print(as.data.frame(tau_convergence))

cat("\n--- Bias by tau ---\n")
print(as.data.frame(tau_bias_full))

cat("\n--- Coverage by tau ---\n")
print(as.data.frame(tau_coverage_full))

cat("\n--- CRI widths by tau ---\n")
print(as.data.frame(tau_cri_full))

cat("\n--- Ranks by tau ---\n")
print(as.data.frame(tau_ranks_full))

cat("\n--- Density Spearman ---\n")
print(as.data.frame(density_spearman_full))

cat("\n--- Density CRI widths ---\n")
print(as.data.frame(density_cri_full))

cat("\n--- Density convergence ---\n")
for (label in density_labels_all) {
  cat(label, "\n")
  print(as.data.frame(density_results[[label]]$conv_summary))
  cat("\n")
}

cat("\n--- True parameters ---\n")
cat("True tau:", true_tau, "\n\n")
print(as.data.frame(true_params))

# ==============================================================================
# 10. SAVE COMPILED OUTPUT
# ==============================================================================

save(tau_spearman_full, tau_convergence, tau_bias_full,
     tau_coverage_full, tau_cri_full, tau_ranks_full,
     density_spearman_full, density_cri_full,
     true_params,
     file = "results/adhoc/ch5_compiled_for_writeup.RData")

# Also save CSVs for easy reference
write.csv(as.data.frame(tau_spearman_full),
          "results/adhoc/tau_spearman_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_convergence),
          "results/adhoc/tau_convergence_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_bias_full),
          "results/adhoc/tau_bias_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_coverage_full),
          "results/adhoc/tau_coverage_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_cri_full),
          "results/adhoc/tau_cri_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_ranks_full),
          "results/adhoc/tau_ranks_full.csv", row.names = FALSE)
write.csv(as.data.frame(density_spearman_full),
          "results/adhoc/density_spearman_full.csv", row.names = FALSE)
write.csv(as.data.frame(density_cri_full),
          "results/adhoc/density_cri_full.csv", row.names = FALSE)

save.image(file = "results/adhoc/ch5_full_workspace.RData")

cat("\n================================================================\n")
cat("  Chapter 5 simulations and compilation complete.\n")
cat("  Results: results/adhoc/\n")
cat("  Figures: results/figures/\n")
cat("================================================================\n")