################################################################################
# Chapter 5: Exploratory Simulations — Boundary Conditions for Reliable Rankings
#
# Three dimensions are varied ONE AT A TIME:
#   5.4.2  Heterogeneity magnitude:  τ from 0.20 → 0.01 (0.27 is in Ch4)
#   5.4.3  Evidence density:         2× studies, 3× studies, dose-enriched
#   5.4.4  Network connectivity:     0 → 5 → 10 → 20 → 40 head-to-head studies
#
# PREREQUISITES:
#   - Run Chapter 4 simulation script through Section 6, so that these exist in the environment:
#       network_template, true_params, true_tau, true_mu_mean, true_mu_sd,
#       dose_equiv, generate_ssri_dataset(), fit_all_models(),
#       extract_metrics(), compute_lpml()
#
# OUTPUT:
#   results/   — CSVs and .RData per simulation dimension
#   figures/ — Key figures for each dimension
################################################################################

setwd("C:/Users/katie/Desktop/Depression-NMA/Chapter 5")

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
model_names <- c("A", "B", "Ref", "Lumped")
agent_names <- true_params$agent

cat("\n================================================================\n")
cat("  CHAPTER 5: EXPLORATORY SIMULATIONS (CORRECTED)\n")
cat("  Backbone: SSRI network (Chapter 4)\n")
cat("  Rankings: unified d_hat approach\n")
cat("================================================================\n")


# ==============================================================================
# HELPER: Run one simulation scenario and compile results
# ==============================================================================

run_scenario <- function(scenario_label,
                         template_scenario   = network_template,
                         true_params_scenario = true_params,
                         true_tau_scenario    = true_tau,
                         n_reps              = 100,
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
  
  # --- Resume from checkpoint if one exists ---
  scen_safe <- gsub("[^a-zA-Z0-9]", "_", scenario_label)
  checkpoint_files <- list.files(
    path = "results",
    pattern = sprintf("^checkpoint_%s_rep\\d+\\.RData$", scen_safe),
    full.names = TRUE
  )
  
  if (length(checkpoint_files) > 0) {
    rep_nums <- as.integer(gsub(".*_rep(\\d+)\\.RData$", "\\1", checkpoint_files))
    latest_file <- checkpoint_files[which.max(rep_nums)]
    cat("  Resuming from:", basename(latest_file), "\n")
    load(latest_file)
    if (length(all_metrics_scen) < n_reps) {
      all_metrics_scen <- c(all_metrics_scen,
                            vector("list", n_reps - length(all_metrics_scen)))
    }
  } else {
    all_metrics_scen <- vector("list", n_reps)
  }
  
  for (rep in 1:n_reps) {
    
    # Skip already-completed replications
    if (!is.null(all_metrics_scen[[rep]])) {
      cat(sprintf("  Rep %d/%d: already complete, skipping\n", rep, n_reps))
      next
    }
    
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
    
    # Checkpoint every 5 reps (more frequent for 100-rep run)
    if (rep %% 5 == 0) {
      save(all_metrics_scen,
           file = sprintf("results/checkpoint_%s_rep%d.RData", scen_safe, rep))
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
      data.frame(model = m, agent = agent_names,
                 bias = round(colMeans(bias_mat), 4),
                 n_reps = nrow(bias_mat), stringsAsFactors = FALSE)
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
      data.frame(model = m, agent = agent_names,
                 coverage = round(colMeans(cov_mat), 3),
                 n_reps = nrow(cov_mat), stringsAsFactors = FALSE)
    } else NULL
  })
  coverage_summary <- bind_rows(cov_by_model)
  
  # Per-agent rank match rate
  rank_match_by_model <- lapply(model_names, function(m) {
    match_mat <- do.call(rbind, lapply(all_metrics_scen, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (!rep[[m]]$converged) return(NULL)
      est_int_rank <- base::rank(rep[[m]]$est_ranks, ties.method = "min")
      return(t(as.integer(est_int_rank == true_params_scenario$true_rank)))
    }))
    if (!is.null(match_mat) && nrow(match_mat) > 0) {
      data.frame(model = m, agent = agent_names,
                 true_rank = true_params_scenario$true_rank,
                 match_rate = round(colMeans(match_mat), 3),
                 n_reps = nrow(match_mat), stringsAsFactors = FALSE)
    } else NULL
  })
  rank_match_summary <- bind_rows(rank_match_by_model)
  
  # SUCRA averaged across replications
  sucra_by_model <- lapply(model_names, function(m) {
    sucra_mat <- do.call(rbind, lapply(all_metrics_scen, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (!rep[[m]]$converged) return(NULL)
      return(t(rep[[m]]$sucra))
    }))
    if (!is.null(sucra_mat) && nrow(sucra_mat) > 0) {
      data.frame(model = m, agent = agent_names,
                 true_rank = true_params_scenario$true_rank,
                 sucra = round(colMeans(sucra_mat), 3),
                 sucra_sd = round(apply(sucra_mat, 2, sd), 3),
                 n_reps = nrow(sucra_mat), stringsAsFactors = FALSE)
    } else NULL
  })
  sucra_summary <- bind_rows(sucra_by_model)
  
  conv_summary <- lapply(model_names, function(m) {
    conv_flags <- sapply(all_metrics_scen, function(rep) {
      if (is.null(rep) || is.null(rep[[m]])) return(NA)
      if (rep[[m]]$failed) return(NA)
      return(rep[[m]]$converged)
    })
    data.frame(model = m,
               n_converged = sum(conv_flags, na.rm = TRUE),
               n_total = sum(!is.na(conv_flags)),
               stringsAsFactors = FALSE)
  })
  conv_summary <- bind_rows(conv_summary)
  
  return(list(
    scenario_label     = scenario_label,
    all_metrics        = all_metrics_scen,
    spearman_summary   = spearman_summary,
    bias_summary       = bias_summary,
    coverage_summary   = coverage_summary,
    rank_match_summary = rank_match_summary,
    sucra_summary      = sucra_summary,
    conv_summary       = conv_summary
  ))
}


################################################################################
#
# SIMULATION 5.4.2: EFFECT OF HETEROGENEITY MAGNITUDE
#
################################################################################

cat("\n\n################################################################\n")
cat("  SIMULATION 5.4.2: HETEROGENEITY MAGNITUDE\n")
cat("################################################################\n")

tau_values <- c(0.20, 0.15, 0.10, 0.05)
n_reps_tau <- 100

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

# Compile across tau levels
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

write.csv(tau_spearman_compiled, "results/tau_spearman.csv", row.names = FALSE)
write.csv(tau_bias_compiled, "results/tau_bias.csv", row.names = FALSE)
write.csv(tau_coverage_compiled, "results/tau_coverage.csv", row.names = FALSE)

save(tau_results, tau_spearman_compiled, tau_bias_compiled, tau_coverage_compiled,
     file = "results/sim_542_heterogeneity.RData")

# Figure
p_tau <- ggplot(tau_spearman_compiled,
                aes(x = tau, y = median_rho, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = model),
              alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.7, linetype = "dotted", color = "grey30") +
  annotate("text", x = max(tau_values), y = 0.73,
           label = "\u03C1 = 0.7 threshold", hjust = 1, size = 3, color = "grey30") +
  scale_x_reverse(breaks = tau_values) +
  scale_color_manual(
    values = c("A" = "#E41A1C", "B" = "#377EB8",
               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
    labels = c("A" = "Model A (Shared Emax)", "B" = "Model B (Agent-specific)",
               "Ref" = "Linear", "Lumped" = "Lumped NMA")) +
  scale_fill_manual(
    values = c("A" = "#E41A1C", "B" = "#377EB8",
               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
    guide = "none") +
  labs(title = "Ranking Accuracy by Between-Study Heterogeneity",
       subtitle = "SSRI network structure; all other parameters at empirical values",
       x = expression(paste("Between-study heterogeneity (", tau, ")")),
       y = expression(paste("Median Spearman ", rho)),
       color = "Model") +
  theme_diss

print(p_tau)
ggsave("figures/fig_542_tau_vs_rho.jpeg", p_tau, width = 9, height = 5.5, dpi = 300)

cat("\n=== Simulation 5.4.2 complete ===\n")


################################################################################
#
# SIMULATION 5.4.3: EFFECT OF EVIDENCE DENSITY
#
################################################################################

cat("\n\n################################################################\n")
cat("  SIMULATION 5.4.3: EVIDENCE DENSITY\n")
cat("################################################################\n")

n_reps_density <- 100

# --- Helper: expand template by duplicating studies ---
expand_template <- function(template, multiplier) {
  if (multiplier == 1) return(template)
  original_studies <- unique(template$studyID)
  expanded <- template
  for (mult in 2:multiplier) {
    new_block <- template
    study_map <- setNames(paste0(original_studies, "_x", mult), original_studies)
    new_block$studyID <- study_map[as.character(new_block$studyID)]
    expanded <- rbind(expanded, new_block)
  }
  return(expanded)
}

# --- Helper: add midpoint doses ---
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
        r = NA_integer_, n = as.integer(round(agent_n)), stringsAsFactors = FALSE)
      new_rows[[length(new_rows) + 1]] <- data.frame(
        studyID = study_counter, agent = ag, dose_equiv = midpoint,
        r = NA_integer_, n = as.integer(round(agent_n)), stringsAsFactors = FALSE)
    }
  }
  new_df <- bind_rows(new_rows)
  expanded <- bind_rows(template, new_df)
  expanded <- expanded %>%
    group_by(studyID) %>% mutate(arm = row_number()) %>% ungroup() %>%
    mutate(r = ifelse(is.na(r), 0L, as.integer(r)))
  return(expanded)
}

density_results <- list()

# (a) 2x studies
cat("\n--- Density: 2x studies ---\n")
template_2x <- expand_template(network_template, 2)
cat(sprintf("  %d studies, %d arms\n",
            length(unique(template_2x$studyID)), nrow(template_2x)))

density_results[["density_2x"]] <- run_scenario(
  scenario_label = "density_2x", template_scenario = template_2x,
  true_params_scenario = true_params, true_tau_scenario = true_tau,
  n_reps = n_reps_density, seed_offset = 11000)

# (b) 3x studies
cat("\n--- Density: 3x studies ---\n")
template_3x <- expand_template(network_template, 3)
cat(sprintf("  %d studies, %d arms\n",
            length(unique(template_3x$studyID)), nrow(template_3x)))

density_results[["density_3x"]] <- run_scenario(
  scenario_label = "density_3x", template_scenario = template_3x,
  true_params_scenario = true_params, true_tau_scenario = true_tau,
  n_reps = n_reps_density, seed_offset = 12000)

# (c) Dose-enriched
cat("\n--- Density: dose-enriched ---\n")
template_enriched <- add_midpoint_doses(network_template)
cat(sprintf("  %d studies, %d arms\n",
            length(unique(template_enriched$studyID)), nrow(template_enriched)))

density_results[["density_enriched"]] <- run_scenario(
  scenario_label = "density_enriched", template_scenario = template_enriched,
  true_params_scenario = true_params, true_tau_scenario = true_tau,
  n_reps = n_reps_density, seed_offset = 50000)

# Compile
density_labels <- c("density_2x", "density_3x", "density_enriched")
density_display <- c("2\u00d7 studies", "3\u00d7 studies", "Dose-enriched")

density_spearman_compiled <- bind_rows(lapply(seq_along(density_labels), function(i) {
  res <- density_results[[density_labels[i]]]$spearman_summary
  res$scenario <- density_display[i]; res$scenario_id <- density_labels[i]
  return(res)
}))

density_bias_compiled <- bind_rows(lapply(seq_along(density_labels), function(i) {
  res <- density_results[[density_labels[i]]]$bias_summary
  res$scenario <- density_display[i]; return(res)
}))

density_coverage_compiled <- bind_rows(lapply(seq_along(density_labels), function(i) {
  res <- density_results[[density_labels[i]]]$coverage_summary
  res$scenario <- density_display[i]; return(res)
}))

write.csv(density_spearman_compiled, "results/density_spearman.csv", row.names = FALSE)
write.csv(density_bias_compiled, "results/density_bias.csv", row.names = FALSE)
write.csv(density_coverage_compiled, "results/density_coverage.csv", row.names = FALSE)

save(density_results, density_spearman_compiled,
     density_bias_compiled, density_coverage_compiled, template_enriched,
     file = "results/sim_543_density.RData")

# Figure
density_spearman_compiled$scenario <- factor(density_spearman_compiled$scenario,
                                             levels = density_display)
p_density <- ggplot(density_spearman_compiled,
                    aes(x = scenario, y = median_rho, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.7, linetype = "dotted", color = "grey30") +
  scale_fill_manual(
    values = c("A" = "#E41A1C", "B" = "#377EB8",
               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
    labels = c("A" = "Model A (Shared Emax)", "B" = "Model B (Agent-specific)",
               "Ref" = "Linear", "Lumped" = "Lumped NMA")) +
  labs(title = "Ranking Accuracy by Evidence Density",
       subtitle = expression(paste(tau, " = 0.27 (empirical); baseline = 60 studies / 145 arms")),
       x = "Density scenario",
       y = expression(paste("Median Spearman ", rho)),
       fill = "Model") +
  theme_diss

print(p_density)
ggsave("figures/fig_543_density_vs_rho.jpeg", p_density,
       width = 9, height = 5.5, dpi = 300)

cat("\n=== Simulation 5.4.3 complete ===\n")


################################################################################
#
# SIMULATION 5.4.4: EFFECT OF NETWORK CONNECTIVITY
#
################################################################################

cat("\n\n################################################################\n")
cat("  SIMULATION 5.4.4: NETWORK CONNECTIVITY\n")
cat("  Adding direct head-to-head comparisons\n")
cat("################################################################\n")

n_reps_connect <- 100

# All 10 possible pairwise comparisons among 5 SSRI agents
agent_pairs <- combn(true_params$agent, 2, simplify = FALSE)
cat("\nAll possible agent pairs (", length(agent_pairs), "):\n")
for (p in agent_pairs) cat("  ", p[1], "vs", p[2], "\n")

dose_1x_equiv <- 20
median_n <- median(network_template$n)
cat("Median arm size for new studies:", median_n, "\n")

create_h2h_study <- function(agent1, agent2, study_id, arm_n = median_n) {
  data.frame(
    studyID = rep(study_id, 2), arm = c(1, 2),
    agent = c(agent1, agent2), dose_equiv = rep(dose_1x_equiv, 2),
    n = rep(as.integer(arm_n), 2), r = rep(0L, 2),
    stringsAsFactors = FALSE)
}

# Level 1: cycle through all 5 agents
cycle_pairs <- list(
  c("citalopram", "escitalopram"), c("escitalopram", "fluoxetine"),
  c("fluoxetine", "paroxetine"), c("paroxetine", "sertraline"),
  c("sertraline", "citalopram"))

# Level 2: all 10 pairs (complete graph x1)
remaining_pairs <- setdiff(
  lapply(agent_pairs, function(p) paste(sort(p), collapse = "|")),
  lapply(cycle_pairs, function(p) paste(sort(p), collapse = "|")))
remaining_pairs <- agent_pairs[
  lapply(agent_pairs, function(p) paste(sort(p), collapse = "|")) %in% remaining_pairs]
level2_pairs <- c(cycle_pairs, remaining_pairs)

# Level 3: all 10 pairs x2
level3_pairs <- c(level2_pairs, level2_pairs)

# Level 4: all 10 pairs x4
level4_pairs <- rep(level2_pairs, 4)

connectivity_levels <- list(
  "h2h_05" = cycle_pairs,
  "h2h_10" = level2_pairs,
  "h2h_20" = level3_pairs,
  "h2h_40" = level4_pairs)

cat("\n=== Connectivity levels ===\n")
for (lvl in names(connectivity_levels)) {
  cat(sprintf("  %s: %d head-to-head studies\n", lvl, length(connectivity_levels[[lvl]])))
}

build_connectivity_template <- function(base_template, h2h_pairs, start_id = 80000) {
  if (length(h2h_pairs) == 0) return(base_template)
  new_studies <- list()
  for (i in seq_along(h2h_pairs)) {
    pair <- h2h_pairs[[i]]
    new_studies[[i]] <- create_h2h_study(pair[1], pair[2], start_id + i)
  }
  new_df <- bind_rows(new_studies)
  expanded <- bind_rows(base_template, new_df)
  expanded <- expanded %>%
    group_by(studyID) %>% mutate(arm = row_number()) %>% ungroup()
  return(expanded)
}

connectivity_templates <- lapply(names(connectivity_levels), function(lvl) {
  build_connectivity_template(network_template, connectivity_levels[[lvl]])
})
names(connectivity_templates) <- names(connectivity_levels)

cat("\n=== Template sizes ===\n")
for (lvl in names(connectivity_templates)) {
  tmpl <- connectivity_templates[[lvl]]
  cat(sprintf("  %s: %d studies, %d arms\n",
              lvl, length(unique(tmpl$studyID)), nrow(tmpl)))
}

# Run scenarios
connectivity_results <- list()

for (lvl in names(connectivity_levels)) {
  cat(sprintf("\n\n--- Connectivity: %s (%d h2h studies) ---\n",
              lvl, length(connectivity_levels[[lvl]])))
  connectivity_results[[lvl]] <- run_scenario(
    scenario_label = lvl,
    template_scenario = connectivity_templates[[lvl]],
    true_params_scenario = true_params,
    true_tau_scenario = true_tau,
    n_reps = n_reps_connect,
    seed_offset = as.integer(gsub("h2h_", "", lvl)) * 1000 + 60000)
}

# Compile
connect_labels <- names(connectivity_levels)
connect_display <- c("5 (cycle)", "10 (complete\u00d71)",
                     "20 (complete\u00d72)", "40 (complete\u00d74)")
n_h2h_studies <- c(5, 10, 20, 40)

connect_spearman <- bind_rows(lapply(seq_along(connect_labels), function(i) {
  res <- connectivity_results[[connect_labels[i]]]$spearman_summary
  res$n_h2h <- n_h2h_studies[i]; res$scenario <- connect_display[i]
  return(res)
}))

connect_bias <- bind_rows(lapply(seq_along(connect_labels), function(i) {
  res <- connectivity_results[[connect_labels[i]]]$bias_summary
  res$n_h2h <- n_h2h_studies[i]; res$scenario <- connect_display[i]
  return(res)
}))

connect_coverage <- bind_rows(lapply(seq_along(connect_labels), function(i) {
  res <- connectivity_results[[connect_labels[i]]]$coverage_summary
  res$n_h2h <- n_h2h_studies[i]; res$scenario <- connect_display[i]
  return(res)
}))

connect_convergence <- bind_rows(lapply(seq_along(connect_labels), function(i) {
  res <- connectivity_results[[connect_labels[i]]]$conv_summary
  res$n_h2h <- n_h2h_studies[i]
  return(res)
}))

connect_cri <- bind_rows(lapply(seq_along(connect_labels), function(i) {
  res_list <- connectivity_results[[connect_labels[i]]]$all_metrics
  bind_rows(lapply(model_names, function(m) {
    width_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.05) return(NULL)
      return(t(rep[[m]]$d_upper - rep[[m]]$d_lower))
    }))
    if (!is.null(width_mat) && nrow(width_mat) > 0) {
      data.frame(
        n_h2h = n_h2h_studies[i], scenario = connect_display[i],
        model = m, agent = agent_names,
        median_width = round(apply(width_mat, 2, median), 3),
        mean_width = round(colMeans(width_mat), 3),
        n_reps = nrow(width_mat), stringsAsFactors = FALSE)
    } else NULL
  }))
}))

write.csv(connect_spearman, "results/connectivity_spearman.csv", row.names = FALSE)
write.csv(connect_bias, "results/connectivity_bias.csv", row.names = FALSE)
write.csv(connect_coverage, "results/connectivity_coverage.csv", row.names = FALSE)
write.csv(connect_cri, "results/connectivity_cri.csv", row.names = FALSE)
write.csv(connect_convergence, "results/connectivity_convergence.csv", row.names = FALSE)

save(connectivity_results, connectivity_templates, connectivity_levels,
     connect_spearman, connect_bias, connect_coverage,
     connect_cri, connect_convergence,
     file = "results/sim_544_connectivity.RData")

# Figure
p_connect <- ggplot(connect_spearman,
                    aes(x = n_h2h, y = median_rho, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = model),
              alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.7, linetype = "dotted", color = "grey30") +
  annotate("text", x = 35, y = 0.73,
           label = expression(rho == 0.7), hjust = 0, size = 3, color = "grey30") +
  scale_x_continuous(breaks = n_h2h_studies) +
  scale_color_manual(
    values = c("A" = "#E41A1C", "B" = "#377EB8",
               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
    labels = c("A" = "Model A (Shared Emax)", "B" = "Model B (Agent-specific)",
               "Ref" = "Linear", "Lumped" = "Lumped NMA")) +
  scale_fill_manual(
    values = c("A" = "#E41A1C", "B" = "#377EB8",
               "Ref" = "#4DAF4A", "Lumped" = "#984EA3"),
    guide = "none") +
  labs(title = "Ranking Accuracy by Network Connectivity",
       subtitle = expression(paste(
         "Progressive addition of head-to-head studies; ",
         tau, " = 0.27, baseline = 60 placebo-controlled studies")),
       x = "Number of added head-to-head studies",
       y = expression(paste("Median Spearman ", rho)),
       color = "Model") +
  theme_diss

print(p_connect)
ggsave("figures/fig_544_connectivity_vs_rho.jpeg", p_connect,
       width = 10, height = 5.5, dpi = 300)

cat("\n=== Simulation 5.4.4 complete ===\n")


################################################################################
#
# COMBINED SUMMARY AND COMPILATION
#
################################################################################

cat("\n\n################################################################\n")
cat("  CHAPTER 5: COMPILING ALL RESULTS\n")
cat("################################################################\n\n")

# Load Chapter 4 base results for tau = 0.27 baseline
load("C:/Users/katie/Desktop/Depression-NMA/Chapter 4/results/simulation/all_results.RData")
all_metrics_ch4 <- all_metrics
rm(all_metrics)

# Add Ch4 baseline as h2h = 0 for connectivity compilation
ch4_connect_baseline <- bind_rows(lapply(model_names, function(m) {
  rho_vals <- sapply(all_metrics_ch4, function(rep) {
    if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NA)
    if (rep[[m]]$max_rhat >= 1.05) return(NA)
    return(rep[[m]]$spearman)
  })
  rho_vals <- rho_vals[!is.na(rho_vals)]
  data.frame(n_h2h = 0, scenario = "0 (star)", model = m,
             median_rho = ifelse(length(rho_vals) > 0, median(rho_vals), NA),
             mean_rho = ifelse(length(rho_vals) > 0, mean(rho_vals), NA),
             q25 = ifelse(length(rho_vals) > 0, quantile(rho_vals, 0.25), NA),
             q75 = ifelse(length(rho_vals) > 0, quantile(rho_vals, 0.75), NA),
             n_reps = length(rho_vals), stringsAsFactors = FALSE)
}))

all_tau_values <- c(0.27, 0.20, 0.15, 0.10, 0.05)

# ==============================================================================
# 1. SPEARMAN BY TAU (R-HAT < 1.05, INCLUDING CH4 BASELINE)
# ==============================================================================

cat("  1. Compiling Spearman by tau...\n")

tau_spearman_full <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) { res_list <- all_metrics_ch4
  } else { res_list <- tau_results[[sprintf("tau_%.2f", tv)]]$all_metrics }
  bind_rows(lapply(model_names, function(m) {
    rho_vals <- sapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NA)
      if (rep[[m]]$max_rhat >= 1.05) return(NA)
      return(rep[[m]]$spearman)
    })
    rho_vals <- rho_vals[!is.na(rho_vals)]
    data.frame(tau = tv, model = m,
               median_rho = ifelse(length(rho_vals) > 0, median(rho_vals), NA),
               mean_rho = ifelse(length(rho_vals) > 0, mean(rho_vals), NA),
               q25 = ifelse(length(rho_vals) > 0, quantile(rho_vals, 0.25), NA),
               q75 = ifelse(length(rho_vals) > 0, quantile(rho_vals, 0.75), NA),
               n_reps = length(rho_vals), stringsAsFactors = FALSE)
  }))
}))

# ==============================================================================
# 2. CONVERGENCE AT BOTH THRESHOLDS
# ==============================================================================

cat("  2. Compiling convergence...\n")

tau_convergence <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) { res_list <- all_metrics_ch4
  } else { res_list <- tau_results[[sprintf("tau_%.2f", tv)]]$all_metrics }
  bind_rows(lapply(model_names, function(m) {
    flags <- sapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(c(NA, NA))
      rhat <- rep[[m]]$max_rhat
      return(rhat < 1.05)
    })
    data.frame(tau = tv, model = m,
               n_converged = sum(flags, na.rm = TRUE),
               n_total = sum(!is.na(flags)),
               stringsAsFactors = FALSE)
  }))
}))

# ==============================================================================
# 3. BIAS BY TAU (R-HAT < 1.05)
# ==============================================================================

cat("  3. Compiling bias...\n")

tau_bias_full <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) { res_list <- all_metrics_ch4
  } else { res_list <- tau_results[[sprintf("tau_%.2f", tv)]]$all_metrics }
  bind_rows(lapply(model_names, function(m) {
    bias_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.05) return(NULL)
      return(t(rep[[m]]$bias))
    }))
    if (!is.null(bias_mat) && nrow(bias_mat) > 0) {
      data.frame(tau = tv, model = m, agent = agent_names,
                 bias = round(colMeans(bias_mat), 4),
                 n_reps = nrow(bias_mat), stringsAsFactors = FALSE)
    } else NULL
  }))
}))

# ==============================================================================
# 4. COVERAGE BY TAU (R-HAT < 1.05)
# ==============================================================================

cat("  4. Compiling coverage...\n")

tau_coverage_full <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) { res_list <- all_metrics_ch4
  } else { res_list <- tau_results[[sprintf("tau_%.2f", tv)]]$all_metrics }
  bind_rows(lapply(model_names, function(m) {
    cov_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.05) return(NULL)
      return(t(rep[[m]]$coverage))
    }))
    if (!is.null(cov_mat) && nrow(cov_mat) > 0) {
      data.frame(tau = tv, model = m, agent = agent_names,
                 coverage = round(colMeans(cov_mat), 3),
                 n_reps = nrow(cov_mat), stringsAsFactors = FALSE)
    } else NULL
  }))
}))

# ==============================================================================
# 5. CRI WIDTHS BY TAU (R-HAT < 1.05)
# ==============================================================================

cat("  5. Compiling CRI widths...\n")

tau_cri_full <- bind_rows(lapply(all_tau_values, function(tv) {
  if (tv == 0.27) { res_list <- all_metrics_ch4
  } else { res_list <- tau_results[[sprintf("tau_%.2f", tv)]]$all_metrics }
  bind_rows(lapply(model_names, function(m) {
    width_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.05) return(NULL)
      return(t(rep[[m]]$d_upper - rep[[m]]$d_lower))
    }))
    if (!is.null(width_mat) && nrow(width_mat) > 0) {
      data.frame(tau = tv, model = m, agent = agent_names,
                 median_width = round(apply(width_mat, 2, median), 3),
                 mean_width = round(colMeans(width_mat), 3),
                 n_reps = nrow(width_mat), stringsAsFactors = FALSE)
    } else NULL
  }))
}))

# ==============================================================================
# 6. RANKS BY TAU (R-HAT < 1.05)
# ==============================================================================

cat("  6. Compiling ranks...\n")

tau_ranks_full <- bind_rows(lapply(all_tau_values, function(tv) {
  # ... existing code ...
}))

# ==============================================================================
# 6b. RANK MATCH AND SUCRA BY TAU (R-HAT < 1.05)
# ==============================================================================

cat("  6b. Compiling rank match and SUCRA by tau...\n")

tau_rank_match_full <- bind_rows(lapply(all_tau_values, function(tv) {
  # ... new block ...
}))

tau_sucra_full <- bind_rows(lapply(all_tau_values, function(tv) {
  # ... new block ...
}))

write.csv(tau_rank_match_full, "results/tau_rank_match_full.csv", row.names = FALSE)
write.csv(tau_sucra_full, "results/tau_sucra_full.csv", row.names = FALSE)

# ==============================================================================
# 7-8. DENSITY COMPILATIONS (R-HAT < 1.05)
# ==============================================================================

cat("  7-8. Compiling density results...\n")

density_labels_all <- names(density_results)

density_spearman_full <- bind_rows(lapply(density_labels_all, function(label) {
  res_list <- density_results[[label]]$all_metrics
  bind_rows(lapply(model_names, function(m) {
    rho_vals <- sapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NA)
      if (rep[[m]]$max_rhat >= 1.05) return(NA)
      return(rep[[m]]$spearman)
    })
    rho_vals <- rho_vals[!is.na(rho_vals)]
    data.frame(scenario = label, model = m,
               median_rho = ifelse(length(rho_vals) > 0, median(rho_vals), NA),
               mean_rho = ifelse(length(rho_vals) > 0, mean(rho_vals), NA),
               n_reps = length(rho_vals), stringsAsFactors = FALSE)
  }))
}))

density_cri_full <- bind_rows(lapply(density_labels_all, function(label) {
  res_list <- density_results[[label]]$all_metrics
  bind_rows(lapply(model_names, function(m) {
    width_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.05) return(NULL)
      return(t(rep[[m]]$d_upper - rep[[m]]$d_lower))
    }))
    if (!is.null(width_mat) && nrow(width_mat) > 0) {
      data.frame(scenario = label, model = m, agent = agent_names,
                 median_width = round(apply(width_mat, 2, median), 3),
                 mean_width = round(colMeans(width_mat), 3),
                 n_reps = nrow(width_mat), stringsAsFactors = FALSE)
    } else NULL
  }))
}))

density_bias_full <- bind_rows(lapply(density_labels_all, function(label) {
  res_list <- density_results[[label]]$all_metrics
  bind_rows(lapply(model_names, function(m) {
    bias_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.05) return(NULL)
      return(t(rep[[m]]$bias))
    }))
    if (!is.null(bias_mat) && nrow(bias_mat) > 0) {
      data.frame(scenario = label, model = m, agent = agent_names,
                 bias = round(colMeans(bias_mat), 4),
                 n_reps = nrow(bias_mat), stringsAsFactors = FALSE)
    } else NULL
  }))
}))

density_coverage_full <- bind_rows(lapply(density_labels_all, function(label) {
  res_list <- density_results[[label]]$all_metrics
  bind_rows(lapply(model_names, function(m) {
    cov_mat <- do.call(rbind, lapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NULL)
      if (rep[[m]]$max_rhat >= 1.05) return(NULL)
      return(t(rep[[m]]$coverage))
    }))
    if (!is.null(cov_mat) && nrow(cov_mat) > 0) {
      data.frame(scenario = label, model = m, agent = agent_names,
                 coverage = round(colMeans(cov_mat), 3),
                 n_reps = nrow(cov_mat), stringsAsFactors = FALSE)
    } else NULL
  }))
}))

# ==============================================================================
# 9. CONNECTIVITY COMPILATIONS (R-HAT < 1.05)
# ==============================================================================

cat("  9. Compiling connectivity results (Rhat < 1.05)...\n")

connect_spearman_full <- bind_rows(lapply(seq_along(connect_labels), function(i) {
  res_list <- connectivity_results[[connect_labels[i]]]$all_metrics
  bind_rows(lapply(model_names, function(m) {
    rho_vals <- sapply(res_list, function(rep) {
      if (is.null(rep) || is.null(rep[[m]]) || rep[[m]]$failed) return(NA)
      if (rep[[m]]$max_rhat >= 1.05) return(NA)
      return(rep[[m]]$spearman)
    })
    rho_vals <- rho_vals[!is.na(rho_vals)]
    data.frame(n_h2h = n_h2h_studies[i], scenario = connect_display[i], model = m,
               median_rho = ifelse(length(rho_vals) > 0, median(rho_vals), NA),
               mean_rho = ifelse(length(rho_vals) > 0, mean(rho_vals), NA),
               n_reps = length(rho_vals), stringsAsFactors = FALSE)
  }))
}))
# Prepend Ch4 baseline (h2h = 0)
connect_spearman_full <- bind_rows(ch4_connect_baseline, connect_spearman_full)

# ==============================================================================
# 10. PRINT ALL SUMMARIES
# ==============================================================================

cat("\n================================================================\n")
cat("  RESULTS SUMMARY\n")
cat("================================================================\n\n")

cat("--- Spearman by tau (R-hat < 1.05) ---\n")
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

cat("\n--- Density bias ---\n")
print(as.data.frame(density_bias_full))

cat("\n--- Density coverage ---\n")
print(as.data.frame(density_coverage_full))

cat("\n--- Connectivity Spearman ---\n")
print(as.data.frame(connect_spearman))

cat("\n--- Connectivity Spearman (Rhat < 1.05) ---\n")
print(as.data.frame(connect_spearman_full))

cat("\n--- Connectivity bias ---\n")
print(as.data.frame(connect_bias))

cat("\n--- Connectivity coverage ---\n")
print(as.data.frame(connect_coverage))

cat("\n--- Connectivity CRI widths ---\n")
print(as.data.frame(connect_cri))

cat("\n--- Connectivity convergence ---\n")
print(as.data.frame(connect_convergence))

cat("\n--- True parameters ---\n")
cat("True tau:", true_tau, "\n\n")
print(as.data.frame(true_params))

cat("\n--- Rank match by tau ---\n")
print(as.data.frame(tau_rank_match_full))

cat("\n--- SUCRA by tau ---\n")
print(as.data.frame(tau_sucra_full))

# ==============================================================================
# 11. SAVE COMPILED OUTPUT
# ==============================================================================

save(tau_spearman_full, tau_convergence, tau_bias_full,
     tau_coverage_full, tau_cri_full, tau_ranks_full,
     tau_rank_match_full, tau_sucra_full,
     density_spearman_full, density_cri_full,
     density_bias_full, density_coverage_full,
     connect_spearman, connect_spearman_full,
     connect_bias, connect_coverage, connect_cri, connect_convergence,
     true_params,
     file = "results/ch5_compiled_for_writeup.RData")

write.csv(as.data.frame(tau_spearman_full), "results/tau_spearman_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_convergence), "results/tau_convergence_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_bias_full), "results/tau_bias_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_coverage_full), "results/tau_coverage_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_cri_full), "results/tau_cri_full.csv", row.names = FALSE)
write.csv(as.data.frame(tau_ranks_full), "results/tau_ranks_full.csv", row.names = FALSE)
write.csv(as.data.frame(density_spearman_full), "results/density_spearman_full.csv", row.names = FALSE)
write.csv(as.data.frame(density_cri_full), "results/density_cri_full.csv", row.names = FALSE)
write.csv(as.data.frame(density_bias_full), "results/density_bias_full.csv", row.names = FALSE)
write.csv(as.data.frame(density_coverage_full), "results/density_coverage_full.csv", row.names = FALSE)
write.csv(as.data.frame(connect_spearman_full), "results/connectivity_spearman_full.csv", row.names = FALSE)

save.image(file = "results/ch5_full_workspace.RData")

cat("\n================================================================\n")
cat("  Chapter 5 simulations and compilation complete.\n")
cat("  Results: results/\n")
cat("  Figures: figures/\n")
cat("================================================================\n")
################################################################################
# Chapter 5: Network Plots
#
# Figure 5.1: Evidence density scenarios (2x2)
#   (a) Baseline (60 studies, 145 arms)
#   (b) 2x studies (120 studies, 290 arms)
#   (c) 3x studies (180 studies, 435 arms)
#   (d) Dose-enriched (76 studies, 177 arms)
#
# Figure 5.2: Network connectivity scenarios (2x2)
#   (a) Baseline (60 studies, star-shaped)
#   (b) +5 head-to-head (cycle)
#   (c) +10 head-to-head (complete x1)
#   (d) +40 head-to-head (complete x4)
#
# PREREQUISITES: MBNMAdose loaded, ssri dataset available,
#                level2_pairs and median_n defined from Ch5 script
################################################################################

library(MBNMAdose)
library(dplyr)
setwd("C:/Users/katie/Desktop/Depression-NMA")

# ==============================================================================
# Helper: add midpoint doses (on original dose scale)
# ==============================================================================

add_midpoint_doses_orig <- function(df) {
  agents <- unique(df$agent[df$agent != "placebo"])
  new_rows <- list()
  study_counter <- 9000
  for (ag in agents) {
    agent_doses <- sort(unique(df$dose[df$agent == ag]))
    agent_n <- as.integer(median(df$n[df$agent == ag]))
    for (j in 1:(length(agent_doses) - 1)) {
      midpoint <- (agent_doses[j] + agent_doses[j + 1]) / 2
      study_counter <- study_counter + 1
      new_rows[[length(new_rows) + 1]] <- data.frame(
        studyID = study_counter, agent = "placebo",
        dose = 0, r = 0L, n = agent_n, stringsAsFactors = FALSE)
      new_rows[[length(new_rows) + 1]] <- data.frame(
        studyID = study_counter, agent = ag,
        dose = midpoint, r = 0L, n = agent_n, stringsAsFactors = FALSE)
    }
  }
  bind_rows(df %>% mutate(studyID = as.character(studyID)),
            bind_rows(new_rows) %>% mutate(studyID = as.character(studyID)))
}

# ==============================================================================
# Build density datasets
# ==============================================================================

# Baseline
net_baseline <- mbnma.network(data.ab = ssri)

# 2x
ssri_2x <- bind_rows(
  ssri %>% mutate(studyID = as.character(studyID)),
  ssri %>% mutate(studyID = paste0(studyID, "_x2")))
net_2x <- mbnma.network(data.ab = ssri_2x)

# 3x
ssri_3x <- bind_rows(
  ssri %>% mutate(studyID = as.character(studyID)),
  ssri %>% mutate(studyID = paste0(studyID, "_x2")),
  ssri %>% mutate(studyID = paste0(studyID, "_x3")))
net_3x <- mbnma.network(data.ab = ssri_3x)

# Dose-enriched
ssri_enriched <- add_midpoint_doses_orig(ssri)
net_enriched <- mbnma.network(data.ab = ssri_enriched)

# ==============================================================================
# FIGURE 5.1: Evidence Density (2x2)
# ==============================================================================

jpeg("figures/Figure 5.1.jpeg", width = 3600, height = 3000, res = 300)
par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))

plot(net_baseline, main = "")
title("(a) Baseline (60 studies, 145 arms)", adj = 0, cex.main = 0.9)

plot(net_2x, main = "")
title("(b) 2\u00d7 studies (120 studies, 290 arms)", adj = 0, cex.main = 0.9)

plot(net_3x, main = "")
title("(c) 3\u00d7 studies (180 studies, 435 arms)", adj = 0, cex.main = 0.9)

plot(net_enriched, main = "")
title("(d) Dose-enriched (76 studies, 177 arms)", adj = 0, cex.main = 0.9)

dev.off()
cat("Figure 5.1 saved to figures/Figure 5.1.jpeg\n")

# ==============================================================================
# Build connectivity datasets
# ==============================================================================

# All 10 possible pairwise combinations among 5 SSRI agents
agent_names_ssri <- sort(unique(ssri$agent[ssri$agent != "placebo"]))
agent_pairs <- combn(agent_names_ssri, 2, simplify = FALSE)

# Cycle: 5 pairs forming a single cycle
cycle_pairs <- list(
  c("citalopram", "escitalopram"),
  c("escitalopram", "fluoxetine"),
  c("fluoxetine", "paroxetine"),
  c("paroxetine", "sertraline"),
  c("sertraline", "citalopram"))

# Complete graph: all 10 pairs
remaining_pairs <- setdiff(
  lapply(agent_pairs, function(p) paste(sort(p), collapse = "|")),
  lapply(cycle_pairs, function(p) paste(sort(p), collapse = "|")))
remaining_pairs <- agent_pairs[
  lapply(agent_pairs, function(p) paste(sort(p), collapse = "|")) %in% remaining_pairs]
level2_pairs <- c(cycle_pairs, remaining_pairs)

# Complete x4: all 10 pairs repeated 4 times
level4_pairs <- rep(level2_pairs, 4)

median_n <- as.integer(median(ssri$n))

# Helper to add h2h studies to the ssri dataset
add_h2h_studies <- function(df, h2h_pairs, start_id = 90000) {
  ssri_h2h <- df %>% mutate(studyID = as.character(studyID))
  counter <- start_id
  for (pair in h2h_pairs) {
    counter <- counter + 1
    agent1_doses <- df$dose[df$agent == pair[1]]
    agent2_doses <- df$dose[df$agent == pair[2]]
    ssri_h2h <- bind_rows(ssri_h2h, data.frame(
      studyID = as.character(counter),
      agent = c(pair[1], pair[2]),
      dose = c(median(agent1_doses), median(agent2_doses)),
      r = c(0L, 0L),
      n = c(median_n, median_n),
      stringsAsFactors = FALSE))
  }
  return(ssri_h2h)
}

# +5 h2h (cycle)
ssri_h2h_5 <- add_h2h_studies(ssri, cycle_pairs)
net_h2h_5 <- mbnma.network(data.ab = ssri_h2h_5)

# +10 h2h (complete x1)
ssri_h2h_10 <- add_h2h_studies(ssri, level2_pairs)
net_h2h_10 <- mbnma.network(data.ab = ssri_h2h_10)

# +20 h2h (complete x2)
level3_pairs <- rep(level2_pairs, 2)
ssri_h2h_20 <- add_h2h_studies(ssri, level3_pairs)
net_h2h_20 <- mbnma.network(data.ab = ssri_h2h_20)

# +40 h2h (complete x4)
ssri_h2h_40 <- add_h2h_studies(ssri, level4_pairs)
net_h2h_40 <- mbnma.network(data.ab = ssri_h2h_40)

# ==============================================================================
# FIGURE 5.2: Network Connectivity (2x2)
# ==============================================================================

jpeg("figures/Figure 5.2.jpeg", width = 3600, height = 3000, res = 300)
par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))

plot(net_h2h_5, main = "")
title("(a) +5 head-to-head (cycle)", adj = 0, cex.main = 0.9)

plot(net_h2h_10, main = "")
title("(b) +10 head-to-head (complete \u00d71)", adj = 0, cex.main = 0.9)

plot(net_h2h_20, main = "")
title("(c) +20 head-to-head (complete \u00d72)", adj = 0, cex.main = 0.9)

plot(net_h2h_40, main = "")
title("(d) +40 head-to-head (complete \u00d74)", adj = 0, cex.main = 0.9)

dev.off()