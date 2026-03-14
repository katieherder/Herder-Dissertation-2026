# ============================================================================
# DOSE-RESPONSE NETWORK META-ANALYSIS
# Testing different functions per agent to find optimal model
# ============================================================================

# Load packages ----
library(MBNMAdose)
library(tidyverse)
library(readxl)

# 1. DATA PREPARATION ----
cat("Loading and preparing data...\n")

tian_data <- read_excel("C://Users//katie//Box//Herder_Dissertation_Dose-Response NMA//Tian_Dataset.xlsx")

tian_data_clean <- tian_data %>%
  rename(
    studyID = Study, 
    agent = Agent, 
    se = SE, 
    dose = Dose,  # Using rounded dose (or use `Exact dose` for exact values)
    n = N
  ) %>%
  filter(studyID != "Gao et al, 2016") %>%  # Remove study with SE = 0
  filter(!is.na(agent) & !is.na(dose) & !is.na(y) & !is.na(se) & !is.na(n))


cat("Data prepared:", nrow(tian_data_clean), "observations from", 
    length(unique(tian_data_clean$studyID)), "studies\n\n")

# 2. CREATE NETWORK ----
network <- mbnma.network(tian_data_clean)
summary(network)

# 3. TEST ALL FUNCTIONS ON FULL NETWORK (FAIR COMPARISON) ----
cat("\n========================================\n")
cat("TESTING ALL FUNCTIONS ON FULL NETWORK\n")
cat("========================================\n")

test_all_functions <- function(network) {
  
  results <- data.frame(
    Function = character(),
    DIC = numeric(),
    Residual_Dev = numeric(),
    SD = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Linear
  cat("Testing linear (all agents)...\n")
  linear_fun <- ~ beta.1 * dose
  linear <- try(mbnma.run(network, 
                          link = "smd",
                          fun = duser(fun = linear_fun, beta.1 = "rel"),
                          method = "random"), silent = TRUE)
  if(inherits(linear, "mbnma")) {
    results <- rbind(results, data.frame(
      Function = "Linear",
      DIC = linear$BUGSoutput$DIC,
      Residual_Dev = linear$BUGSoutput$summary["totresdev", "mean"],
      SD = linear$BUGSoutput$summary["sd", "mean"]
    ))
  }
  
  # Quadratic
  cat("Testing quadratic (all agents)...\n")
  quad_fun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
  quadratic <- try(mbnma.run(network,
                             link = "smd",
                             fun = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
                             method = "random"), silent = TRUE)
  if(inherits(quadratic, "mbnma")) {
    results <- rbind(results, data.frame(
      Function = "Quadratic",
      DIC = quadratic$BUGSoutput$DIC,
      Residual_Dev = quadratic$BUGSoutput$summary["totresdev", "mean"],
      SD = quadratic$BUGSoutput$summary["sd", "mean"]
    ))
  }
  
  # Exponential
  cat("Testing exponential (all agents)...\n")
  exp_fun <- ~ (1 - exp(-beta.1 * dose))
  exponential <- try(mbnma.run(network,
                               link = "smd",
                               fun = duser(fun = exp_fun, beta.1 = "rel"),
                               method = "random"), silent = TRUE)
  if(inherits(exponential, "mbnma")) {
    results <- rbind(results, data.frame(
      Function = "Exponential",
      DIC = exponential$BUGSoutput$DIC,
      Residual_Dev = exponential$BUGSoutput$summary["totresdev", "mean"],
      SD = exponential$BUGSoutput$summary["sd", "mean"]
    ))
  }
  
  # Emax
  cat("Testing Emax (all agents)...\n")
  emax_fun <- ~ emax * dose / (ed50 + dose)
  emax <- try(mbnma.run(network,
                        link = "smd",
                        fun = duser(fun = emax_fun, emax = "rel", ed50 = "rel"),
                        method = "random"), silent = TRUE)
  if(inherits(emax, "mbnma")) {
    results <- rbind(results, data.frame(
      Function = "Emax",
      DIC = emax$BUGSoutput$DIC,
      Residual_Dev = emax$BUGSoutput$summary["totresdev", "mean"],
      SD = emax$BUGSoutput$summary["sd", "mean"]
    ))
  }
  
  results <- results[order(results$DIC), ]
  return(results)
}

overall_results <- test_all_functions(network)
print(overall_results)

# 4. TEST EACH AGENT SEPARATELY ----
cat("\n========================================\n")
cat("TESTING FUNCTIONS FOR EACH AGENT\n")
cat("========================================\n")

test_agent_functions <- function(network_data, agent_name) {
  
  cat("\nAgent:", agent_name, "\n")
  cat("-------------------\n")
  
  # Find studies that have BOTH this agent AND control
  studies_with_agent <- network_data %>%
    filter(agent == agent_name) %>%
    pull(studyID) %>%
    unique()
  
  studies_with_control <- network_data %>%
    filter(agent == "CON") %>%
    pull(studyID) %>%
    unique()
  
  # Keep only studies that have both
  valid_studies <- intersect(studies_with_agent, studies_with_control)
  
  cat("Studies with", agent_name, "AND control:", length(valid_studies), "\n")
  
  if(length(valid_studies) < 3) {
    cat("Insufficient studies for", agent_name, "\n")
    return(NULL)
  }
  
  # Subset to valid studies only
  agent_subset <- network_data %>%
    filter(studyID %in% valid_studies) %>%
    filter(agent %in% c("CON", agent_name))
  
  cat("Data points:", nrow(agent_subset), "\n")
  
  # Create network
  agent_network <- try(mbnma.network(agent_subset), silent = FALSE)
  
  if(!inherits(agent_network, "mbnma.network")) {
    cat("Network creation failed for", agent_name, "\n")
    return(NULL)
  }
  
  cat("Network created successfully!\n")
  
  results <- data.frame(
    Agent = character(),
    Function = character(),
    DIC = numeric(),
    Residual_Dev = numeric(),
    SD = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Test Linear
  cat("  Linear... ")
  linear_fun <- ~ beta.1 * dose
  linear <- try(mbnma.run(agent_network, 
                          link = "smd",
                          fun = duser(fun = linear_fun, beta.1 = "rel"),
                          method = "random"), silent = TRUE)
  if(inherits(linear, "mbnma")) {
    results <- rbind(results, data.frame(
      Agent = agent_name,
      Function = "Linear",
      DIC = linear$BUGSoutput$DIC,
      Residual_Dev = linear$BUGSoutput$summary["totresdev", "mean"],
      SD = linear$BUGSoutput$summary["sd", "mean"]
    ))
    cat("DIC =", round(linear$BUGSoutput$DIC, 1), "\n")
  } else {
    cat("Failed\n")
  }
  
  # Test Quadratic
  cat("  Quadratic... ")
  quad_fun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
  quadratic <- try(mbnma.run(agent_network,
                             link = "smd",
                             fun = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
                             method = "random"), silent = TRUE)
  if(inherits(quadratic, "mbnma")) {
    results <- rbind(results, data.frame(
      Agent = agent_name,
      Function = "Quadratic",
      DIC = quadratic$BUGSoutput$DIC,
      Residual_Dev = quadratic$BUGSoutput$summary["totresdev", "mean"],
      SD = quadratic$BUGSoutput$summary["sd", "mean"]
    ))
    cat("DIC =", round(quadratic$BUGSoutput$DIC, 1), "\n")
  } else {
    cat("Failed\n")
  }
  
  # Test Exponential
  cat("  Exponential... ")
  exp_fun <- ~ (1 - exp(-beta.1 * dose))
  exponential <- try(mbnma.run(agent_network,
                               link = "smd",
                               fun = duser(fun = exp_fun, beta.1 = "rel"),
                               method = "random"), silent = TRUE)
  if(inherits(exponential, "mbnma")) {
    results <- rbind(results, data.frame(
      Agent = agent_name,
      Function = "Exponential",
      DIC = exponential$BUGSoutput$DIC,
      Residual_Dev = exponential$BUGSoutput$summary["totresdev", "mean"],
      SD = exponential$BUGSoutput$summary["sd", "mean"]
    ))
    cat("DIC =", round(exponential$BUGSoutput$DIC, 1), "\n")
  } else {
    cat("Failed\n")
  }
  
  if(nrow(results) > 0) {
    results <- results[order(results$DIC), ]
  }
  return(results)
}

# Test all agents
agents <- c("AE", "MBE", "ME", "RE")
agent_results_list <- list()

for (agent in agents) {
  agent_results_list[[agent]] <- test_agent_functions(tian_data_clean, agent)
}

# Combine results
agent_results <- do.call(rbind, agent_results_list)

if(!is.null(agent_results) && nrow(agent_results) > 0) {
  print(agent_results)
  
  # Show best function for each agent
  cat("\n========================================\n")
  cat("BEST FUNCTION FOR EACH AGENT:\n")
  cat("========================================\n")
  best_per_agent <- agent_results %>%
    group_by(Agent) %>%
    slice_min(DIC, n = 1) %>%
    select(Agent, Function, DIC, SD)
  print(best_per_agent)
} else {
  cat("No valid agent-specific results obtained.\n")
  cat("Recommendation: Use the overall network results instead.\n")
}

# 5. BUILD OPTIMIZED MULTI-FUNCTION MODEL ----
cat("\n========================================\n")
cat("BUILDING OPTIMIZED MULTI-FUNCTION MODEL\n")
cat("========================================\n")

# YOU NEED TO MANUALLY SET THESE BASED ON RESULTS ABOVE!
# Example structure - modify based on your actual results:

# Define the functions
linear_fun <- ~ beta.1 * dose
quad_fun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
exp_fun <- ~ (1 - exp(-beta.1 * dose))
emax_fun <- ~ emax * dose / (ed50 + dose)

# MODIFY THIS BASED ON YOUR RESULTS FROM AGENT-SPECIFIC TESTING
# Example: if AE=linear, MBE=quad, ME=linear, RE=exp
cat("Building multi-function model...\n")
cat("(Modify function assignments based on results above!)\n\n")

network$agents

multi_fun <- dmulti(list(
  Placebo = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),  # placebo, no dose effect
  AE  = duser(fun = linear_fun, beta.1 = "rel"),
  MBE = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
  ME  = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
  RE  = duser(fun = linear_fun, beta.1 = "rel")
))
set.seed(1234)

# Fit optimized model with custom prior

nma_optimized <- mbnma.run(
  network,
  link = "smd",
  fun = multi_fun,
  method = "random",
  n.iter = 50000,   # increase iterations
  n.burnin = 10000,
  n.chains = 3
)

summary(nma_optimized)

# 6. COMPARE ALL MODELS ----
cat("\n========================================\n")
cat("MODEL COMPARISON SUMMARY\n")
cat("========================================\n")

comparison <- data.frame(
  Model = c("Linear (all)", "Quadratic (all)", "Exponential (all)", 
            "Emax (all)", "Optimized (mixed)"),
  DIC = c(
    if(nrow(overall_results[overall_results$Function == "Linear",]) > 0) 
      overall_results[overall_results$Function == "Linear", "DIC"] else NA,
    if(nrow(overall_results[overall_results$Function == "Quadratic",]) > 0) 
      overall_results[overall_results$Function == "Quadratic", "DIC"] else NA,
    if(nrow(overall_results[overall_results$Function == "Exponential",]) > 0) 
      overall_results[overall_results$Function == "Exponential", "DIC"] else NA,
    if(nrow(overall_results[overall_results$Function == "Emax",]) > 0) 
      overall_results[overall_results$Function == "Emax", "DIC"] else NA,
    nma_optimized$BUGSoutput$DIC
  )
)

comparison <- comparison[order(comparison$DIC), ]
print(comparison)

cat("\n*** BEST MODEL:", comparison$Model[1], "with DIC =", 
    round(comparison$DIC[1], 2), "***\n")


## OPTIMIZE FURTHER
print(nma_optimized)

set.seed(12345)
multi_fun <- dmulti(list(
  Placebo = duser(fun = linear_fun, beta.1 = "rel"),  # placebo, no dose effect
  AE  = duser(fun = linear_fun, beta.1 = "rel"),
  MBE = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
  ME  = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
  RE  = duser(fun = linear_fun, beta.1 = "rel")
))

custom_priors <- list(
  beta.1 = "dnorm(0, 0.1)",
  beta.2 = "dnorm(0, 0.05)",
  sd     = "dunif(0, 60)",
  D      = "dnorm(0, 0.05)"
)

# Fit optimized model with custom prior

nma_optimized_2 <- mbnma.run(
  network,
  link = "smd",
  fun = multi_fun,
  method = "random",
  n.iter = 50000,   # increase iterations
  n.burnin = 10000,
  n.chains = 3,
  priors = custom_priors
)
summary(nma_optimized)
summary(nma_optimized_2)
summary(overall_results)

# For nma_optimized
if("BUGSoutput" %in% names(nma_optimized)) {
  pd_optimized <- nma_optimized$BUGSoutput$pD
  dic_optimized <- nma_optimized$BUGSoutput$DIC
  
  cat("nma_optimized:\n")
  cat(sprintf("  pD: %.1f\n", pd_optimized))
  cat(sprintf("  DIC: %.1f\n", dic_optimized))
}

# For nma_optimized_2
if("BUGSoutput" %in% names(nma_optimized_2)) {
  pd_optimized_2 <- nma_optimized_2$BUGSoutput$pD
  dic_optimized_2 <- nma_optimized_2$BUGSoutput$DIC
  
  cat("\nnma_optimized_2:\n")
  cat(sprintf("  pD: %.1f\n", pd_optimized_2))
  cat(sprintf("  DIC: %.1f\n", dic_optimized_2))
}

# Re-run just quadratic to extract pD
cat("Re-running quadratic model to extract pD...\n")
quad_fun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
quadratic <- mbnma.run(network,
                       link = "smd",
                       fun = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
                       method = "random")

# Extract all metrics
quad_metrics <- data.frame(
  Function = "Quadratic",
  DIC = quadratic$BUGSoutput$DIC,
  Deviance = quadratic$BUGSoutput$summary["deviance", "mean"],
  Residual_Dev = quadratic$BUGSoutput$summary["totresdev", "mean"],
  pD = quadratic$BUGSoutput$pD,
  SD = quadratic$BUGSoutput$summary["sd", "mean"]
)

print(quad_metrics)

# For your slide - this is the row from the original paper
cat("\n=== Paper's reported values (from your slide) ===\n")
cat("2 MBE Quadratic 107. 3.58\n")
cat("\n=== Your replication ===\n")
cat(sprintf("2 MBE Quadratic %.1f %.3f (pD: %.1f)\n", 
            quad_metrics$DIC, 
            quad_metrics$SD,
            quad_metrics$pD))
