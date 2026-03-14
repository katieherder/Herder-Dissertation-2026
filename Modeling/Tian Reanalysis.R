# Load the packages.
library(MBNMAdose)
library(tidyverse)
library(cowplot)
library(janitor)
library(patchwork)
library(readxl)
library(mcmcplots)
library(rjags)
library(multinma)
library(overlapping)
library(ggthemes)
# Read the Excel file
tian_data <- read_excel("C://Users//katie//Box//Herder_Dissertation_Dose-Response NMA//Tian_Dataset.xlsx")

# Use Exact dose
tian_data_clean <- tian_data %>%
  rename(studyID = Study, agent = Agent, se = SE, dose=Dose, n=N) %>%
  filter(studyID != "Gao et al, 2016") %>%
  filter(!is.na(agent) & !is.na(dose) & !is.na(y) & !is.na(se))

# Create network
network <- mbnma.network(tian_data_clean)

# Define quadratic function THEIR WAY
dpoly_fun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))

# Fit using their syntax
nma_quad <- mbnma.run(network, 
                      link = "smd",  # or try "smd"
                      fun = duser(fun = dpoly_fun, beta.1 = "rel", beta.2 = "rel"), 
                      method = "random")
summary(nma_quad)


# More comprehensive function testing
test_functions <- function(network, agent_name = NULL) {
  
  results <- data.frame(
    Function = character(),
    DIC = numeric(),
    Residual_Dev = numeric(),
    SD = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 1. Linear
  cat("Testing linear...\n")
  linear <- try(mbnma.run(network, fun = dpoly(degree = 1), method = "random"), silent = TRUE)
  if(class(linear)[1] != "try-error") {
    results <- rbind(results, data.frame(
      Function = "Linear",
      DIC = linear$BUGSoutput$DIC,
      Residual_Dev = linear$BUGSoutput$summary["totresdev", "mean"],
      SD = linear$BUGSoutput$summary["sd", "mean"]
    ))
  }
  
  # 2. Quadratic with duser
  cat("Testing quadratic...\n")
  quad_fun <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
  quadratic <- try(mbnma.run(network, link = "smd",
                             fun = duser(fun = quad_fun, beta.1 = "rel", beta.2 = "rel"),
                             method = "random"), silent = TRUE)
  if(class(quadratic)[1] != "try-error") {
    results <- rbind(results, data.frame(
      Function = "Quadratic",
      DIC = quadratic$BUGSoutput$DIC,
      Residual_Dev = quadratic$BUGSoutput$summary["totresdev", "mean"],
      SD = quadratic$BUGSoutput$summary["sd", "mean"]
    ))
  }
  
  # 3. Emax
  cat("Testing Emax...\n")
  emax <- try(mbnma.emax(network, emax = "rel", ed50 = "rel", method = "random"), silent = TRUE)
  if(class(emax)[1] != "try-error") {
    results <- rbind(results, data.frame(
      Function = "Emax",
      DIC = emax$BUGSoutput$DIC,
      Residual_Dev = emax$BUGSoutput$summary["totresdev", "mean"],
      SD = emax$BUGSoutput$summary["sd", "mean"]
    ))
  }
  
  # 4. Exponential
  cat("Testing exponential...\n")
  exponential <- try(mbnma.run(network, fun = dexp(), method = "random"), silent = TRUE)
  if(class(exponential)[1] != "try-error") {
    results <- rbind(results, data.frame(
      Function = "Exponential",
      DIC = exponential$BUGSoutput$DIC,
      Residual_Dev = exponential$BUGSoutput$summary["totresdev", "mean"],
      SD = exponential$BUGSoutput$summary["sd", "mean"]
    ))
  }
  
  # Sort by DIC
  results <- results[order(results$DIC), ]
  return(results)
}

# Test all functions on full network
all_results <- test_functions(network)
print(all_results)