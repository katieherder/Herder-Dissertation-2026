###############################################################################
# Chapter 4 Simulation: Data-Generating Function for Pharmacologic MBNMA
# 
# Purpose: Generate synthetic dose-response NMA datasets that mimic the
#          structure of the GRISELDA antidepressant network, with known
#          true parameters so we can evaluate bias, coverage, and model
#          selection across different MBNMA specifications.
#
# Author: Katherine Herder
# Date: February 2026
###############################################################################

library(tidyverse)

# =============================================================================
# PART 1: Define the "true" pharmacologic network structure
# =============================================================================

#' Define the true network structure
#'
#' This creates a simplified version of GRISELDA with:
#'   - 3 classes (SSRIs, SNRIs, Other) 
#'   - 8 agents total
#'   - Realistic dose levels per agent
#'   - Placebo as reference
#'
#' You can scale this up later to match full GRISELDA (21 agents, 4+ classes)

define_network_structure <- function() {
  
  agents <- tibble(
    agent_id   = 1:8,
    agent_name = c("Placebo", "fluoxetine", "sertraline", "escitalopram",
                   "venlafaxine", "duloxetine",
                   "mirtazapine", "bupropion"),
    class_name = c("Placebo", "SSRI", "SSRI", "SSRI",
                   "SNRI", "SNRI",
                   "Other", "Other"),
    class_id   = c(0, 1, 1, 1, 2, 2, 3, 3),
    # Typical dose levels tested in RCTs (mg/day)
    # These mimic what's actually in GRISELDA
    doses = list(
      c(0),                          # Placebo
      c(10, 20, 40, 60),             # fluoxetine
      c(50, 100, 150, 200),          # sertraline
      c(5, 10, 20),                  # escitalopram
      c(75, 150, 225, 375),          # venlafaxine
      c(30, 60, 120),               # duloxetine
      c(15, 30, 45),                # mirtazapine
      c(150, 300, 450)              # bupropion
    )
  )
  
  return(agents)
}


# =============================================================================
# PART 2: True dose-response functions
# =============================================================================

#' Emax dose-response function
#' @param dose Dose in mg
#' @param emax Maximum effect (asymptote)
#' @param ed50 Dose producing 50% of max effect
#' @return Effect size (SMD scale, negative = better for depression)
true_emax <- function(dose, emax, ed50) {
  emax * dose / (ed50 + dose)
}

#' Linear dose-response function
true_linear <- function(dose, beta) {
  beta * dose
}

#' Log-linear dose-response function  
true_loglinear <- function(dose, beta) {
  beta * log(dose + 1)
}


#' Define true dose-response parameters for each scenario
#'
#' Returns a list with:
#'   - func_by_class: which function each class truly follows
#'   - params: true parameter values per agent
#'   - tau: true between-study heterogeneity SD
#'
#' @param scenario Character: "all_emax", "mixed_forms", "sparse"

define_true_params <- function(scenario = "all_emax") {
  
  agents <- define_network_structure()
  
  if (scenario == "all_emax") {
    # -----------------------------------------------------------------
    # Scenario 1: All classes follow Emax (same form, different params)
    # Tests: Does class-effect pooling help? Does form selection matter
    #        when the form is actually the same?
    # -----------------------------------------------------------------
    
    true_params <- tibble(
      agent_name = agents$agent_name,
      class_name = agents$class_name,
      func       = "emax",  # everyone truly Emax
      # True Emax values (SMD scale; negative = improvement)
      # Within-class agents are similar but not identical
      emax = c(0,           # placebo
               -0.35, -0.30, -0.40,   # SSRIs (class mean ~ -0.35)
               -0.45, -0.50,          # SNRIs (class mean ~ -0.475)
               -0.25, -0.20),         # Other (class mean ~ -0.225)
      # True ED50 values (mg) 
      ed50 = c(NA,          # placebo
               15, 60, 8,            # SSRIs 
               100, 40,              # SNRIs
               20, 200),             # Other
      beta = NA  # not used in this scenario
    )
    
    tau <- 0.15  # between-study heterogeneity SD
    
  } else if (scenario == "mixed_forms") {
    # -----------------------------------------------------------------
    # Scenario 2: Classes follow DIFFERENT functional forms
    # Tests: Can model selection detect this? Does forcing shared Emax
    #        produce bias? This motivates the class-specific extension.
    # -----------------------------------------------------------------
    
    true_params <- tibble(
      agent_name = agents$agent_name,
      class_name = agents$class_name,
      func       = c("none",
                     "emax", "emax", "emax",       # SSRIs: truly Emax
                     "emax", "emax",                # SNRIs: truly Emax (different params)
                     "loglinear", "loglinear"),      # Other: truly log-linear
      emax = c(0, -0.35, -0.30, -0.40, -0.45, -0.50, NA, NA),
      ed50 = c(NA, 15, 60, 8, 100, 40, NA, NA),
      beta = c(NA, NA, NA, NA, NA, NA, -0.08, -0.05)  # log-linear slope for "Other"
    )
    
    tau <- 0.15
    
  } else if (scenario == "sparse") {
    # -----------------------------------------------------------------
    # Scenario 3: Same as all_emax but with reduced data
    # Tests: At what point does MBNMA break down? When do priors matter?
    # Fewer studies per comparison, fewer dose levels observed
    # -----------------------------------------------------------------
    
    # Same params as all_emax
    true_params <- tibble(
      agent_name = agents$agent_name,
      class_name = agents$class_name,
      func       = "emax",
      emax = c(0, -0.35, -0.30, -0.40, -0.45, -0.50, -0.25, -0.20),
      ed50 = c(NA, 15, 60, 8, 100, 40, 20, 200),
      beta = NA
    )
    
    tau <- 0.20  # slightly higher heterogeneity (more realistic for small networks)
  }
  
  return(list(true_params = true_params, tau = tau, scenario = scenario))
}


# =============================================================================
# PART 3: Generate a single simulated dataset
# =============================================================================

#' Compute the true effect for a given agent at a given dose
#'
#' @param agent_name Character
#' @param dose Numeric (mg)
#' @param true_params Tibble from define_true_params()
compute_true_effect <- function(agent_name, dose, true_params) {
  
  p <- true_params %>% filter(agent_name == !!agent_name)
  
  if (agent_name == "Placebo" | dose == 0) return(0)
  
  if (p$func == "emax") {
    return(true_emax(dose, p$emax, p$ed50))
  } else if (p$func == "loglinear") {
    return(true_loglinear(dose, p$beta))
  } else if (p$func == "linear") {
    return(true_linear(dose, p$beta))
  } else {
    stop(paste("Unknown function:", p$func))
  }
}


#' Generate one simulated NMA dataset
#'
#' Mimics the structure of a pharmacologic dose-response NMA:
#'   - Each study compares 2-3 arms (usually drug doses vs placebo)
#'   - Studies are mostly placebo-controlled (like GRISELDA)
#'   - Some head-to-head trials between active drugs
#'   - Continuous outcome (mean change in depression score, SMD scale)
#'
#' @param params Output from define_true_params()
#' @param n_studies_per_agent Number of studies per agent (approximate)
#' @param n_per_arm Sample size per arm (drawn from realistic range)
#' @param prop_headtohead Proportion of studies that are head-to-head
#' @param seed Random seed
#' @return A data frame in MBNMAdose long format

generate_dataset <- function(params,
                             n_studies_per_agent = 8,
                             n_per_arm = c(30, 150),  # range for Uniform draw
                             prop_headtohead = 0.15,
                             seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  true_params <- params$true_params
  tau <- params$tau
  agents <- define_network_structure()
  
  active_agents <- agents %>% filter(agent_name != "Placebo")
  
  studies <- list()
  study_counter <- 0
  
  for (i in 1:nrow(active_agents)) {
    ag <- active_agents$agent_name[i]
    ag_doses <- active_agents$doses[[i]]
    n_studies <- n_studies_per_agent
    
    # For sparse scenario, reduce studies for some agents
    if (params$scenario == "sparse") {
      n_studies <- sample(2:5, 1)  # 2-5 studies per agent
      # Also reduce dose levels observed
      if (length(ag_doses) > 2) {
        n_doses_observed <- sample(2:min(3, length(ag_doses)), 1)
        ag_doses <- sort(sample(ag_doses, n_doses_observed))
      }
    }
    
    for (j in 1:n_studies) {
      study_counter <- study_counter + 1
      
      # Decide if placebo-controlled or head-to-head
      is_h2h <- runif(1) < prop_headtohead
      
      # --- Study-level random effects ---
      # Study baseline (placebo response on SMD scale)
      mu_i <- rnorm(1, mean = -1.5, sd = 0.5)  # typical placebo response
      
      # Pick dose(s) for this study
      # Most studies test 1-2 doses of the same drug vs placebo
      n_doses_in_study <- sample(1:min(2, length(ag_doses)), 1, 
                                 prob = c(0.6, rep(0.4/(min(2, length(ag_doses))-1), 
                                                   min(2, length(ag_doses))-1)))
      study_doses <- sort(sample(ag_doses, n_doses_in_study))
      
      # Build arms
      arms <- list()
      
      if (!is_h2h) {
        # --- Placebo-controlled trial ---
        # Arm 1: Placebo
        n_k <- round(runif(1, n_per_arm[1], n_per_arm[2]))
        arms[[1]] <- tibble(
          studyID = study_counter,
          agent   = "Placebo",
          dose    = 0,
          class   = "Placebo",
          N       = n_k,
          true_effect = 0
        )
        
        # Active arms
        for (d in seq_along(study_doses)) {
          n_k <- round(runif(1, n_per_arm[1], n_per_arm[2]))
          te <- compute_true_effect(ag, study_doses[d], true_params)
          arms[[d + 1]] <- tibble(
            studyID = study_counter,
            agent   = ag,
            dose    = study_doses[d],
            class   = active_agents$class_name[i],
            N       = n_k,
            true_effect = te
          )
        }
        
      } else {
        # --- Head-to-head trial ---
        # Pick a second agent (different from first)
        other_agents <- active_agents %>% filter(agent_name != ag)
        other <- other_agents %>% slice_sample(n = 1)
        other_dose <- sample(other$doses[[1]], 1)
        
        n_k1 <- round(runif(1, n_per_arm[1], n_per_arm[2]))
        n_k2 <- round(runif(1, n_per_arm[1], n_per_arm[2]))
        
        te1 <- compute_true_effect(ag, study_doses[1], true_params)
        te2 <- compute_true_effect(other$agent_name, other_dose, true_params)
        
        arms[[1]] <- tibble(
          studyID = study_counter,
          agent   = ag,
          dose    = study_doses[1],
          class   = active_agents$class_name[i],
          N       = n_k1,
          true_effect = te1
        )
        arms[[2]] <- tibble(
          studyID = study_counter,
          agent   = other$agent_name,
          dose    = other_dose,
          class   = other$class_name,
          N       = n_k2,
          true_effect = te2
        )
      }
      
      study_df <- bind_rows(arms)
      
      # --- Generate observed data ---
      # For each arm, generate y (mean change) and se
      # y_{i,k} = mu_i + delta_{i,k} + noise
      # where delta_{i,k} ~ N(true_effect_k - true_effect_1, tau^2)
      
      # Reference arm effect
      ref_effect <- study_df$true_effect[1]
      
      study_df <- study_df %>%
        mutate(
          # Within-study SE (depends on sample size)
          se = sqrt(1 / N) * runif(n(), 1.5, 2.5),  # realistic SE range
          
          # Study-specific relative effect with heterogeneity
          delta = if_else(row_number() == 1, 0,
                          rnorm(n(), mean = true_effect - ref_effect, sd = tau)),
          
          # Observed outcome: baseline + relative effect + sampling error
          y = mu_i + delta + rnorm(n(), 0, se)
        ) %>%
        select(studyID, agent, dose, class, y, se, N)
      
      studies[[study_counter]] <- study_df
    }
  }
  
  dataset <- bind_rows(studies)
  
  # Ensure studyID is sequential
  dataset <- dataset %>%
    mutate(studyID = as.numeric(factor(studyID)))
  
  return(dataset)
}


# =============================================================================
# PART 4: Fit models and extract results
# =============================================================================

#' Fit all candidate MBNMA models to a simulated dataset
#'
#' @param data A simulated dataset from generate_dataset()
#' @param n_iter MCMC iterations (use fewer for simulation, more for real data)
#' @return A list of fitted model objects and summary statistics

fit_models <- function(data, n_iter = 15000, n_burnin = 5000) {
  
  library(MBNMAdose)
  
  # Create MBNMAdose network object
  network <- mbnma.network(data)
  
  results <- list()
  
  # --- Model 1: Shared Emax (all agents same parameters) ---
  tryCatch({
    mod_shared <- mbnma.run(network, 
                            fun = demax(), 
                            method = "random",
                            n.iter = n_iter, 
                            n.burnin = n_burnin,
                            n.chains = 3)
    results$shared_emax <- mod_shared
  }, error = function(e) {
    results$shared_emax <<- list(error = e$message)
  })
  
  # --- Model 2: Shared Emax + class effects on Emax parameter ---
  tryCatch({
    mod_class <- mbnma.run(network, 
                           fun = demax(),
                           class.effect = list(emax = "random"),
                           method = "random",
                           n.iter = n_iter, 
                           n.burnin = n_burnin,
                           n.chains = 3)
    results$class_emax <- mod_class
  }, error = function(e) {
    results$class_emax <<- list(error = e$message)
  })
  
  # --- Model 3: Shared Linear (wrong form if truth is Emax) ---
  tryCatch({
    mod_linear <- mbnma.run(network, 
                            fun = dpoly(degree = 1), 
                            method = "random",
                            n.iter = n_iter, 
                            n.burnin = n_burnin,
                            n.chains = 3)
    results$shared_linear <- mod_linear
  }, error = function(e) {
    results$shared_linear <<- list(error = e$message)
  })
  
  # --- Model 4: Shared Log-linear ---
  tryCatch({
    mod_loglin <- mbnma.run(network, 
                            fun = dloglin(), 
                            method = "random",
                            n.iter = n_iter, 
                            n.burnin = n_burnin,
                            n.chains = 3)
    results$shared_loglin <- mod_loglin
  }, error = function(e) {
    results$shared_loglin <<- list(error = e$message)
  })
  
  return(results)
}


#' Extract key metrics from fitted models
#'
#' @param fitted_models Output from fit_models()
#' @param true_params Output from define_true_params()
#' @return A tibble with one row per model containing DIC, convergence info, etc.

extract_metrics <- function(fitted_models, true_params) {
  
  metrics <- tibble()
  
  for (mod_name in names(fitted_models)) {
    mod <- fitted_models[[mod_name]]
    
    if ("error" %in% names(mod)) {
      metrics <- bind_rows(metrics, tibble(
        model = mod_name, 
        DIC = NA, pD = NA, totresdev = NA,
        converged = FALSE, 
        max_Rhat = NA
      ))
      next
    }
    
    # Extract DIC
    dic <- mod$BUGSoutput$DIC
    pd <- mod$BUGSoutput$pD
    
    # Total residual deviance
    totresdev <- mod$BUGSoutput$median$totresdev
    
    # Convergence: max Rhat across all parameters
    rhat_vals <- mod$BUGSoutput$summary[, "Rhat"]
    max_rhat <- max(rhat_vals, na.rm = TRUE)
    converged <- max_rhat < 1.05
    
    # Between-study heterogeneity estimate (if random effects)
    tau_est <- tryCatch(
      mod$BUGSoutput$median$sd,
      error = function(e) NA
    )
    
    metrics <- bind_rows(metrics, tibble(
      model = mod_name,
      DIC = dic,
      pD = pd,
      totresdev = totresdev,
      converged = converged,
      max_Rhat = max_rhat,
      tau_est = tau_est,
      tau_true = true_params$tau
    ))
  }
  
  return(metrics)
}


# =============================================================================
# PART 5: Run the simulation
# =============================================================================

#' Run a single simulation replicate
#'
#' @param rep_id Replicate number (used as seed)
#' @param scenario Scenario name
#' @param n_studies_per_agent Studies per agent
#' @param n_iter MCMC iterations
#' @return A tibble of metrics for this replicate

run_one_replicate <- function(rep_id, scenario, 
                              n_studies_per_agent = 8,
                              n_iter = 15000) {
  
  cat("  Rep", rep_id, "... ")
  
  # Generate true parameters
  params <- define_true_params(scenario)
  
  # Generate dataset
  data <- generate_dataset(params, 
                           n_studies_per_agent = n_studies_per_agent,
                           seed = rep_id * 1000 + as.numeric(Sys.time()) %% 1000)
  
  # Fit models
  fitted <- fit_models(data, n_iter = n_iter, n_burnin = floor(n_iter / 2))
  
  # Extract metrics
  metrics <- extract_metrics(fitted, params)
  metrics$replicate <- rep_id
  metrics$scenario <- scenario
  
  cat("done\n")
  return(metrics)
}


#' Run full simulation for one scenario
#'
#' @param scenario Scenario name
#' @param n_reps Number of replicates
#' @param n_studies_per_agent Studies per agent
#' @param n_iter MCMC iterations
#' @param parallel Use parallel processing?
#' @param n_cores Number of cores

run_simulation <- function(scenario, 
                           n_reps = 50,
                           n_studies_per_agent = 8,
                           n_iter = 15000,
                           parallel = FALSE,
                           n_cores = 4) {
  
  cat("=== Running scenario:", scenario, "with", n_reps, "replicates ===\n")
  
  if (parallel) {
    library(parallel)
    results <- mclapply(1:n_reps, function(i) {
      run_one_replicate(i, scenario, n_studies_per_agent, n_iter)
    }, mc.cores = n_cores)
  } else {
    results <- lapply(1:n_reps, function(i) {
      run_one_replicate(i, scenario, n_studies_per_agent, n_iter)
    })
  }
  
  all_results <- bind_rows(results)
  return(all_results)
}


# =============================================================================
# PART 6: Summarize simulation results
# =============================================================================

#' Summarize simulation results across replicates
#'
#' @param sim_results Output from run_simulation()

summarize_simulation <- function(sim_results) {
  
  summary <- sim_results %>%
    group_by(scenario, model) %>%
    summarize(
      n_reps        = n(),
      n_converged   = sum(converged, na.rm = TRUE),
      pct_converged = mean(converged, na.rm = TRUE) * 100,
      mean_DIC      = mean(DIC, na.rm = TRUE),
      sd_DIC        = sd(DIC, na.rm = TRUE),
      mean_resdev   = mean(totresdev, na.rm = TRUE),
      mean_tau_est  = mean(tau_est, na.rm = TRUE),
      tau_true      = first(tau_true),
      tau_bias      = mean(tau_est - tau_true, na.rm = TRUE),
      # How often is this model selected by DIC?
      .groups = "drop"
    )
  
  # Add DIC selection rate
  dic_selection <- sim_results %>%
    filter(converged) %>%
    group_by(scenario, replicate) %>%
    slice_min(DIC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    count(scenario, model, name = "n_selected")
  
  summary <- summary %>%
    left_join(dic_selection, by = c("scenario", "model")) %>%
    mutate(
      n_selected = replace_na(n_selected, 0),
      pct_selected = n_selected / n_reps * 100
    )
  
  return(summary)
}


# =============================================================================
# PART 7: Example usage
# =============================================================================

# --- Quick test (run this first to make sure everything works!) ---
# Takes ~5 minutes

if (FALSE) {  # Set to TRUE to run
  
  # Test data generation
  params <- define_true_params("all_emax")
  test_data <- generate_dataset(params, n_studies_per_agent = 5, seed = 42)
  
  cat("Generated dataset with", nrow(test_data), "rows,",
      length(unique(test_data$studyID)), "studies\n")
  cat("Agents:", paste(unique(test_data$agent), collapse = ", "), "\n")
  cat("Dose range per agent:\n")
  test_data %>% 
    group_by(agent) %>% 
    summarize(doses = paste(unique(dose), collapse = ", "),
              n_studies = length(unique(studyID))) %>%
    print()
  
  # Test fitting (just one dataset, ~3-5 min)
  fitted <- fit_models(test_data, n_iter = 10000)
  metrics <- extract_metrics(fitted, params)
  print(metrics)
}


# --- Preliminary simulation (for March 13 draft) ---
# ~50 reps x 3 scenarios x 4 models x ~3 min each = ~30 hours
# Run overnight with parallel processing

if (FALSE) {  # Set to TRUE to run
  
  # Start with just 10 reps to verify pipeline
  results_pilot <- run_simulation("all_emax", n_reps = 10, n_iter = 10000)
  print(summarize_simulation(results_pilot))
  
  # If that works, scale up
  results_s1 <- run_simulation("all_emax",     n_reps = 50, parallel = TRUE, n_cores = 4)
  results_s2 <- run_simulation("mixed_forms",  n_reps = 50, parallel = TRUE, n_cores = 4)
  results_s3 <- run_simulation("sparse",       n_reps = 50, parallel = TRUE, n_cores = 4)
  
  all_results <- bind_rows(results_s1, results_s2, results_s3)
  
  # Save results
  saveRDS(all_results, "ch4_simulation_results_preliminary.rds")
  
  # Summarize
  summary <- summarize_simulation(all_results)
  print(summary)
  write_csv(summary, "ch4_simulation_summary.csv")
}


# --- Full simulation (for April 13 defense) ---
# Scale up to 200 reps, add prior sensitivity

if (FALSE) {  # Set to TRUE when ready
  
  results_full_s1 <- run_simulation("all_emax",    n_reps = 200, parallel = TRUE, n_cores = 8)
  results_full_s2 <- run_simulation("mixed_forms", n_reps = 200, parallel = TRUE, n_cores = 8)
  results_full_s3 <- run_simulation("sparse",      n_reps = 200, parallel = TRUE, n_cores = 8)
  
  all_results_full <- bind_rows(results_full_s1, results_full_s2, results_full_s3)
  saveRDS(all_results_full, "ch4_simulation_results_full.rds")
}

###############################################################################
# Chapter 4: Prior Sensitivity Analysis for Heterogeneity in MBNMA
#
# Purpose: Systematically vary priors on between-study heterogeneity (tau)
#          and dose-response parameters to assess their impact on:
#          - Point estimates (Emax, ED50)
#          - Credible interval widths
#          - Treatment rankings
#          - Model fit (DIC, residual deviance)
#
# Extends Rosenberger et al. (2021) into the dose-response NMA context.
#
# Author: Katherine Herder
# Date: February 2026
###############################################################################

library(tidyverse)
library(MBNMAdose)
library(R2jags)
library(coda)

# =============================================================================
# PART 1: Define prior specifications
# =============================================================================

#' Get JAGS code from a fitted MBNMAdose model and modify priors
#'
#' MBNMAdose defaults:
#'   - Heterogeneity: sd ~ dunif(0, [user-specified or large])
#'   - Dose-response params: dnorm(0, 0.0001)  [variance = 10,000]
#'   - Class SDs: dnorm(0, 0.0025) T(0,)       [variance = 400]
#'
#' We test alternatives from the literature:

prior_specs <- list(
  
  # --- Prior Set 1: MBNMAdose defaults (vague) ---
  vague = list(
    name = "Vague (default)",
    description = "MBNMAdose defaults: Uniform(0,5) on tau",
    # This is the default - no modification needed
    tau_prior = "sd ~ dunif(0, 5)",
    modify = FALSE
  ),
  
  # --- Prior Set 2: Moderately informative half-normal ---
  # Rosenberger et al. (2021) showed this performed well
  halfnormal = list(
    name = "Half-Normal(0, 1)",
    description = "Rosenberger et al. recommended; tau ~ HN(0,1)",
    tau_prior = "sd ~ dnorm(0, 1) T(0,)",
    modify = TRUE,
    # In JAGS: replace the sd prior line
    old_prior_pattern = "sd ~ dunif\\(0,.*\\)",
    new_prior = "sd ~ dnorm(0, 1) T(0,)"
  ),
  
  # --- Prior Set 3: Turner et al. informative log-normal ---
  # For pharmacological interventions, subjective outcomes (like depression)
  # Turner et al. (2012): LN(-2.34, 1.62^2) for SMD
  # This gives: median tau ~ 0.096, 95% interval ~ [0.005, 1.98]
  turner = list(
    name = "Turner Informative",
    description = "Turner et al. (2012) LN(-2.34, 1.62^2) for pharma/subjective",
    tau_prior = "prec.tau <- pow(sd, -2)\nlog(sd) ~ dnorm(-2.34, pow(1.62, -2))",
    modify = TRUE,
    # Need to replace both the sd prior and how it enters the model
    # This is trickier - may need to replace the precision specification
    old_prior_pattern = "sd ~ dunif\\(0,.*\\)",
    new_prior = "sd ~ dlnorm(-2.34, 0.381)"  # 1/1.62^2 = 0.381
  ),
  
  # --- Prior Set 4: Weakly informative half-Cauchy ---
  # Recommended by Gelman (2006) for variance components
  # Approximated in JAGS via t-distribution
  halfcauchy = list(
    name = "Half-Cauchy(0, 0.5)",
    description = "Gelman (2006) recommendation, scale=0.5",
    tau_prior = "sd ~ dt(0, 4, 1) T(0,)",  # dt with precision=4 gives scale~0.5
    modify = TRUE,
    old_prior_pattern = "sd ~ dunif\\(0,.*\\)",
    new_prior = "sd ~ dt(0, 4, 1) T(0,)"
  )
)


# =============================================================================
# PART 2: Modify JAGS code with different priors
# =============================================================================

#' Extract and modify JAGS model code from a fitted MBNMAdose model
#'
#' Strategy: 
#'   1. Fit a model with MBNMAdose to get the JAGS code
#'   2. Extract it via result$model.arg$jagscode
#'   3. Use regex to swap out the prior specification
#'   4. Re-run via R2jags with the modified code
#'
#' @param base_model A fitted mbnma.run() object
#' @param prior_spec A prior specification from prior_specs list
#' @return Modified JAGS code as a string

modify_jags_priors <- function(base_model, prior_spec) {
  
  # Extract original JAGS code
  jags_code <- base_model$model.arg$jagscode
  
  if (!prior_spec$modify) {
    return(jags_code)  # defaults, no change needed
  }
  
  # Replace the heterogeneity prior
  # MBNMAdose typically writes: sd ~ dunif(0, <value>)
  # We need to find and replace this
  
  modified <- gsub(
    "sd ~ dunif\\(0,\\s*[0-9.]+\\)",
    prior_spec$new_prior,
    jags_code
  )
  
  # Verify the substitution happened
  if (modified == jags_code) {
    # Try alternative patterns MBNMAdose might use
    modified <- gsub(
      "sd\\s*~\\s*dunif\\(0\\s*,\\s*[^)]+\\)",
      prior_spec$new_prior,
      jags_code
    )
  }
  
  if (modified == jags_code) {
    warning("Could not find heterogeneity prior to replace. ",
            "Check the JAGS code manually.\n",
            "Pattern sought: sd ~ dunif(0, ...)")
  }
  
  return(modified)
}


# =============================================================================
# PART 3: Run prior sensitivity analysis
# =============================================================================

#' Run the full prior sensitivity analysis on a dataset
#'
#' @param data Dataset in MBNMAdose format (studyID, agent, dose, class, y, se, N)
#' @param base_fun Dose-response function (e.g., demax())
#' @param use_class_effect Whether to include class effects
#' @param n_iter MCMC iterations
#' @param n_burnin Burn-in iterations
#' @param priors_to_test Names of priors from prior_specs (default: all)
#' @return List with fitted models and comparison metrics

run_prior_sensitivity <- function(data, 
                                  base_fun = demax(),
                                  use_class_effect = FALSE,
                                  n_iter = 30000,
                                  n_burnin = 10000,
                                  n_thin = 2,
                                  priors_to_test = names(prior_specs)) {
  
  # Create network
  network <- mbnma.network(data)
  
  # --- Step 1: Fit baseline model with defaults ---
  cat("Fitting baseline model with default priors...\n")
  
  if (use_class_effect) {
    base_model <- mbnma.run(network, fun = base_fun,
                            class.effect = list(emax = "random"),
                            method = "random",
                            n.iter = n_iter, n.burnin = n_burnin, 
                            n.thin = n_thin, n.chains = 3)
  } else {
    base_model <- mbnma.run(network, fun = base_fun,
                            method = "random",
                            n.iter = n_iter, n.burnin = n_burnin,
                            n.thin = n_thin, n.chains = 3)
  }
  
  # Print the JAGS code for inspection
  cat("\n--- Base JAGS code (inspect for prior locations) ---\n")
  cat(base_model$model.arg$jagscode)
  cat("\n--- End JAGS code ---\n\n")
  
  results <- list()
  results[["vague"]] <- base_model
  
  # --- Step 2: Fit with each alternative prior ---
  for (prior_name in priors_to_test) {
    if (prior_name == "vague") next  # already fitted
    
    spec <- prior_specs[[prior_name]]
    cat("Fitting model with", spec$name, "prior...\n")
    
    # Modify JAGS code
    modified_jags <- modify_jags_priors(base_model, spec)
    
    # Write to temp file
    temp_file <- tempfile(fileext = ".jags")
    writeLines(modified_jags, temp_file)
    
    # Re-run with modified priors
    # Use mbnma.run with model.file argument if supported,
    # otherwise fall back to R2jags directly
    tryCatch({
      if (use_class_effect) {
        mod <- mbnma.run(network, fun = base_fun,
                         class.effect = list(emax = "random"),
                         method = "random",
                         model.file = temp_file,
                         n.iter = n_iter, n.burnin = n_burnin,
                         n.thin = n_thin, n.chains = 3)
      } else {
        mod <- mbnma.run(network, fun = base_fun,
                         method = "random",
                         model.file = temp_file,
                         n.iter = n_iter, n.burnin = n_burnin,
                         n.thin = n_thin, n.chains = 3)
      }
      results[[prior_name]] <- mod
    }, error = function(e) {
      cat("  Error with", spec$name, ":", e$message, "\n")
      cat("  Attempting direct R2jags fallback...\n")
      
      # Fallback: run via R2jags directly
      # You'll need to extract the data list from the base model
      tryCatch({
        jags_data <- base_model$model.arg$jagsdata
        jags_params <- base_model$parameters.to.save
        
        mod_jags <- jags(
          data = jags_data,
          parameters.to.save = jags_params,
          model.file = temp_file,
          n.iter = n_iter,
          n.burnin = n_burnin,
          n.thin = n_thin,
          n.chains = 3
        )
        results[[prior_name]] <- mod_jags
      }, error = function(e2) {
        cat("  R2jags fallback also failed:", e2$message, "\n")
        results[[prior_name]] <<- list(error = e2$message)
      })
    })
    
    # Clean up
    unlink(temp_file)
  }
  
  return(results)
}


# =============================================================================
# PART 4: Compare results across priors
# =============================================================================

#' Extract comparable metrics across prior specifications
#'
#' @param prior_results Output from run_prior_sensitivity()
#' @return A tibble comparing key metrics

compare_priors <- function(prior_results) {
  
  comparison <- tibble()
  
  for (prior_name in names(prior_results)) {
    mod <- prior_results[[prior_name]]
    
    # Skip errors
    if ("error" %in% names(mod)) {
      comparison <- bind_rows(comparison, tibble(
        prior = prior_name,
        prior_label = prior_specs[[prior_name]]$name,
        DIC = NA, pD = NA, totresdev = NA,
        tau_median = NA, tau_lower = NA, tau_upper = NA,
        max_Rhat = NA
      ))
      next
    }
    
    # Handle both mbnma objects and raw jags objects
    bugs <- if ("BUGSoutput" %in% names(mod)) mod$BUGSoutput else mod$BUGSoutput
    
    # DIC
    dic <- bugs$DIC
    pd <- bugs$pD
    resdev <- tryCatch(bugs$median$totresdev, error = function(e) NA)
    
    # Heterogeneity estimates
    # MBNMAdose stores this as "sd" 
    sd_summary <- tryCatch({
      idx <- grep("^sd$|^sd\\[", rownames(bugs$summary))
      if (length(idx) > 0) {
        list(
          median = bugs$summary[idx[1], "50%"],
          lower  = bugs$summary[idx[1], "2.5%"],
          upper  = bugs$summary[idx[1], "97.5%"]
        )
      } else {
        list(median = NA, lower = NA, upper = NA)
      }
    }, error = function(e) list(median = NA, lower = NA, upper = NA))
    
    # Convergence
    rhat_vals <- bugs$summary[, "Rhat"]
    max_rhat <- max(rhat_vals, na.rm = TRUE)
    
    comparison <- bind_rows(comparison, tibble(
      prior       = prior_name,
      prior_label = prior_specs[[prior_name]]$name,
      DIC         = dic,
      pD          = pd,
      totresdev   = resdev,
      tau_median  = sd_summary$median,
      tau_lower   = sd_summary$lower,
      tau_upper   = sd_summary$upper,
      tau_CrI_width = sd_summary$upper - sd_summary$lower,
      max_Rhat    = max_rhat
    ))
  }
  
  return(comparison)
}


#' Compare treatment rankings across prior specifications
#'
#' @param prior_results Output from run_prior_sensitivity()
#' @param agents Character vector of agent names to rank
#' @return A tibble with rank probabilities per prior

compare_rankings <- function(prior_results, agents = NULL) {
  
  rankings <- list()
  
  for (prior_name in names(prior_results)) {
    mod <- prior_results[[prior_name]]
    if ("error" %in% names(mod)) next
    
    # Try to extract dose-response parameter estimates
    # For Emax models, we care about emax (max effect) per agent
    bugs <- mod$BUGSoutput
    
    # Find emax parameters
    emax_idx <- grep("^d\\.1\\[|^emax\\[|^beta\\.1\\[", rownames(bugs$summary))
    
    if (length(emax_idx) > 0) {
      emax_df <- tibble(
        prior = prior_name,
        param_name = rownames(bugs$summary)[emax_idx],
        median = bugs$summary[emax_idx, "50%"],
        lower = bugs$summary[emax_idx, "2.5%"],
        upper = bugs$summary[emax_idx, "97.5%"],
        CrI_width = upper - lower
      )
      rankings[[prior_name]] <- emax_df
    }
  }
  
  if (length(rankings) > 0) {
    return(bind_rows(rankings))
  } else {
    return(tibble())
  }
}


# =============================================================================
# PART 5: Visualization functions
# =============================================================================

#' Plot heterogeneity posterior distributions across priors
#'
#' @param prior_results Output from run_prior_sensitivity()

plot_tau_posteriors <- function(prior_results) {
  
  posteriors <- list()
  
  for (prior_name in names(prior_results)) {
    mod <- prior_results[[prior_name]]
    if ("error" %in% names(mod)) next
    
    bugs <- mod$BUGSoutput
    
    # Extract MCMC samples for sd (tau)
    sd_idx <- grep("^sd$", colnames(bugs$sims.matrix))
    if (length(sd_idx) > 0) {
      posteriors[[prior_name]] <- tibble(
        prior = prior_specs[[prior_name]]$name,
        tau = bugs$sims.matrix[, sd_idx[1]]
      )
    }
  }
  
  if (length(posteriors) == 0) {
    cat("No tau posteriors found to plot.\n")
    return(invisible(NULL))
  }
  
  plot_data <- bind_rows(posteriors)
  
  p <- ggplot(plot_data, aes(x = tau, fill = prior, color = prior)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    labs(
      title = "Posterior Distribution of Between-Study Heterogeneity (τ)",
      subtitle = "Across prior specifications",
      x = expression(tau ~ "(between-study SD)"),
      y = "Posterior density",
      fill = "Prior", color = "Prior"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom") +
    coord_cartesian(xlim = c(0, 1.5))
  
  return(p)
}


#' Plot CrI widths for treatment effects across priors
#'
#' Shows how prior choice affects precision of estimates

plot_cri_comparison <- function(ranking_comparison) {
  
  if (nrow(ranking_comparison) == 0) return(invisible(NULL))
  
  p <- ggplot(ranking_comparison, 
              aes(x = param_name, y = CrI_width, fill = prior)) +
    geom_col(position = "dodge") +
    labs(
      title = "95% Credible Interval Width by Prior Specification",
      x = "Parameter", y = "CrI Width (SMD)",
      fill = "Prior"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


# =============================================================================
# PART 6: Example usage
# =============================================================================

if (FALSE) {  # Set to TRUE to run
  
  # ---- Using the built-in SSRI dataset as a test ----
  # (Use this to verify the pipeline before running on GRISELDA)
  
  data(ssri)
  head(ssri)
  
  # Note: ssri dataset is binary (r, N), not continuous
  # For testing the pipeline, this works fine
  # Your GRISELDA analysis will use continuous outcomes (y, se)
  
  # Quick test with SSRI data
  ssri_network <- mbnma.network(ssri)
  plot(ssri_network)
  
  # Fit base model
  base <- mbnma.run(ssri_network, fun = demax(),
                    method = "random",
                    n.iter = 20000, n.burnin = 10000, n.chains = 3)
  
  # Inspect the JAGS code to find prior locations
  cat(base$model.arg$jagscode)
  
  # Run prior sensitivity
  prior_results <- run_prior_sensitivity(
    data = ssri,
    base_fun = demax(),
    n_iter = 20000,
    n_burnin = 10000,
    priors_to_test = c("vague", "halfnormal", "turner")
  )
  
  # Compare
  comparison <- compare_priors(prior_results)
  print(comparison)
  
  # Rankings
  rank_comp <- compare_rankings(prior_results)
  print(rank_comp)
  
  # Plot
  p_tau <- plot_tau_posteriors(prior_results)
  print(p_tau)
  ggsave("ch4_tau_posterior_comparison.png", p_tau, width = 10, height = 6)
  
  
  # ---- For GRISELDA (once you have the data loaded) ----
  # 
  # griselda <- read.csv("griselda_dose_response.csv")
  # 
  # # You'll need columns: studyID, agent, dose, class, y, se, N
  # # (or r, N for binary outcome)
  # 
  # # Define classes
  # griselda <- griselda %>%
  #   mutate(class = case_when(
  #     agent %in% c("fluoxetine", "sertraline", "paroxetine", 
  #                  "citalopram", "escitalopram", "fluvoxamine") ~ "SSRI",
  #     agent %in% c("venlafaxine", "duloxetine", 
  #                  "desvenlafaxine", "milnacipran") ~ "SNRI",
  #     agent %in% c("amitriptyline", "clomipramine", "imipramine") ~ "TCA",
  #     TRUE ~ "Other"
  #   ))
  # 
  # prior_results <- run_prior_sensitivity(
  #   data = griselda,
  #   base_fun = demax(),
  #   use_class_effect = TRUE,
  #   n_iter = 50000,
  #   n_burnin = 20000,
  #   priors_to_test = c("vague", "halfnormal", "turner", "halfcauchy")
  # )
  # 
  # comparison <- compare_priors(prior_results)
  # write_csv(comparison, "ch4_prior_sensitivity_results.csv")
}