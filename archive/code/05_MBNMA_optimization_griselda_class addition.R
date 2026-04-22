################################################################################
# Chapter 4: GRISELDA Class-Level ED50 Analysis
#
# CORRECTED VERSION — March 2026
# FIX: Unified ranking approach using d_hat for all models including C1, C2.
#
# PREREQUISITE: Run the corrected main GRISELDA analysis first, then:
#   load("results/griselda/griselda_all_models.RData")
#   Griselda_clean must be available
#
# Dependencies: MBNMAdose (>= 0.5.0), R2jags, ggplot2, dplyr, tidyr,
#               patchwork, coda, stringr
################################################################################

# ==============================================================================
# 0. Setup & Reload
# ==============================================================================

library(MBNMAdose)
library(R2jags)
library(ggplot2)
library(dplyr)
library(tidyr)
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

n_iter   <- 20000
n_chains <- 3
n_burnin <- floor(n_iter / 2)
n_thin   <- max(1, floor((n_iter - n_burnin) / 1000))

load("results/griselda/griselda_all_models.RData")
cat("Loaded saved model objects from main GRISELDA analysis.\n")

placebo_rate_g <- sum(Griselda_clean$r[Griselda_clean$agent == "placebo"]) /
  sum(Griselda_clean$n[Griselda_clean$agent == "placebo"])
cat("Pooled placebo response rate:", round(placebo_rate_g, 3), "\n")

active_agents <- sort(unique(Griselda_clean$agent[Griselda_clean$agent != "placebo"]))
cat("Active agents (", length(active_agents), "):", paste(active_agents, collapse = ", "), "\n")


# ==============================================================================
# 1. Rebuild Network with 'class' Variable
# ==============================================================================

griselda_mbnma_class <- Griselda_clean %>%
  mutate(
    dose  = ifelse(agent == "placebo", 0, dose_FE),
    class = ifelse(agent == "placebo", "placebo", Class)
  ) %>%
  select(studyID, agent, dose, n, r, class) %>%
  as.data.frame()

cat("\n=== Class assignments ===\n")
class_check <- griselda_mbnma_class %>%
  filter(agent != "placebo") %>%
  distinct(agent, class) %>%
  arrange(class, agent)
print(as.data.frame(class_check))
cat("\nAgents per class:\n")
print(table(class_check$class))

griselda_net_class <- mbnma.network(griselda_mbnma_class)
cat("\nNetwork created with class variable.\n")


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
  if (is_nma) { bugs <- model$jagsresult$BUGSoutput
  } else { bugs <- model$BUGSoutput }
  if (!is.null(bugs)) {
    rhat_vals <- bugs$summary[, "Rhat"]
    rhat_vals <- rhat_vals[!is.na(rhat_vals)]
    cat("  Max Rhat:", round(max(rhat_vals), 4), "\n")
    n_above <- sum(rhat_vals > 1.05)
    cat("  Parameters with Rhat > 1.05:", n_above, "of", length(rhat_vals), "\n")
    if (n_above > 0) {
      flagged <- rhat_vals[rhat_vals > 1.05]
      cat("  Flagged parameters:\n")
      print(round(sort(flagged, decreasing = TRUE), 3))
    }
    neff_vals <- bugs$summary[, "n.eff"]
    neff_vals <- neff_vals[!is.na(neff_vals) & neff_vals > 0]
    cat("  Min n.eff:", min(neff_vals), "\n")
    cat("  Median n.eff:", median(neff_vals), "\n")
  }
  invisible(NULL)
}

# If extract_d_hat and compute_ranks_from_d_hat were saved in RData, they're
# already loaded. Otherwise define them here:
if (!exists("extract_d_hat")) {
  extract_d_hat <- function(model, model_type, agent_names, doses, is_nma = FALSE) {
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
    }
    return(d_hat)
  }
}

if (!exists("compute_ranks_from_d_hat")) {
  compute_ranks_from_d_hat <- function(d_hat, model_name) {
    agent_names <- colnames(d_hat)
    n_agents <- ncol(d_hat)
    rank_mat <- t(apply(d_hat, 1, function(row) base::rank(-row)))
    mean_ranks   <- colMeans(rank_mat)
    median_ranks <- apply(rank_mat, 2, median)
    sucra        <- (n_agents - mean_ranks) / (n_agents - 1)
    data.frame(
      agent = agent_names, median_rank = median_ranks, mean_rank = mean_ranks,
      sucra = sucra, d_median = apply(d_hat, 2, median),
      d_lower = apply(d_hat, 2, quantile, 0.025),
      d_upper = apply(d_hat, 2, quantile, 0.975),
      model = model_name, stringsAsFactors = FALSE
    )
  }
}


# ==============================================================================
# 3. Fit Class-Effect Models
# ==============================================================================

cat("\n================================================================\n")
cat("  Fitting Model C1: Common Class Effect on ED50\n")
cat("================================================================\n")

model_C1_g <- mbnma.run(
  network    = griselda_net_class,
  fun        = demax(emax = "rel", ed50 = "rel"),
  class.effect = list(ed50 = "common"),
  method     = "random",
  likelihood = "binomial",
  link       = "logit",
  parameters.to.save = c("resdev", "emax", "ed50", "sd", "totresdev", "theta", "ED50"),
  n.iter     = n_iter,
  n.chains   = n_chains,
  n.burnin   = n_burnin,
  n.thin     = n_thin,
  jags.seed  = 123456
)
print(model_C1_g)
check_convergence(model_C1_g, "Model C1 (Class-common ED50)")

cat("\n================================================================\n")
cat("  Fitting Model C2: Random Class Effect on ED50\n")
cat("================================================================\n")

model_C2_g <- mbnma.run(
  network    = griselda_net_class,
  fun        = demax(emax = "rel", ed50 = "rel"),
  class.effect = list(ed50 = "random"),
  method     = "random",
  likelihood = "binomial",
  link       = "logit",
  parameters.to.save = c("resdev", "emax", "ed50", "sd", "totresdev", "theta", "ED50", "sd.ED50"),
  n.iter     = n_iter,
  n.chains   = n_chains,
  n.burnin   = n_burnin,
  n.thin     = n_thin,
  jags.seed  = 123456
)
print(model_C2_g)
check_convergence(model_C2_g, "Model C2 (Class-random ED50)")


# ==============================================================================
# 4. Diagnostic: Check ed50 dimensions for C1 and C2
# ==============================================================================
# This determines whether to use "emax_shared" or "emax_agent" for ranking

cat("\n=== ED50 posterior structure diagnostic ===\n")
cat("C1 ed50 dimensions:", paste(dim(model_C1_g$BUGSoutput$sims.list$ed50), collapse = " x "), "\n")
cat("C2 ed50 dimensions:", paste(dim(model_C2_g$BUGSoutput$sims.list$ed50), collapse = " x "), "\n")

# If C1 has per-agent columns (class-common means same within class, different
# across classes), use "emax_agent". If single column, use "emax_shared".
c1_ed50_dim <- dim(model_C1_g$BUGSoutput$sims.list$ed50)
c1_model_type <- if (!is.null(c1_ed50_dim) && length(c1_ed50_dim) == 2 && c1_ed50_dim[2] > 1) {
  "emax_agent"
} else {
  "emax_shared"
}
cat("C1 ranking model_type:", c1_model_type, "\n")

# C2 should always have per-agent columns
cat("C2 ranking model_type: emax_agent\n")


# ==============================================================================
# 5. Model Comparison: DIC and LPML
# ==============================================================================

cat("\n================================================================\n")
cat("  Model Comparison: Class-Effect Models vs. Original Models\n")
cat("================================================================\n")

lpml_C1_g <- compute_lpml(model_C1_g, "Model C1")
lpml_C2_g <- compute_lpml(model_C2_g, "Model C2")

comparison_class <- data.frame(
  Model = c("A: Shared Emax (network ED50)",
            "B: Agent-specific Emax",
            "C1: Class-common ED50",
            "C2: Class-random ED50",
            "Ref: Shared Linear"),
  DIC   = round(c(model_A_g$BUGSoutput$DIC, model_B_g$BUGSoutput$DIC,
                  model_C1_g$BUGSoutput$DIC, model_C2_g$BUGSoutput$DIC,
                  model_Ref_g$BUGSoutput$DIC), 1),
  pD    = round(c(model_A_g$BUGSoutput$pD, model_B_g$BUGSoutput$pD,
                  model_C1_g$BUGSoutput$pD, model_C2_g$BUGSoutput$pD,
                  model_Ref_g$BUGSoutput$pD), 1),
  LPML  = round(c(lpml_A_g$lpml, lpml_B_g$lpml,
                  lpml_C1_g$lpml, lpml_C2_g$lpml, lpml_Ref_g$lpml), 1),
  stringsAsFactors = FALSE
)
comparison_class$Dbar <- comparison_class$DIC - comparison_class$pD

cat("\n=== Full Model Comparison (GRISELDA) ===\n")
print(comparison_class)


# ==============================================================================
# 6. Extract Class-Level ED50 Estimates
# ==============================================================================

cat("\n================================================================\n")
cat("  Class-Level ED50 Parameter Estimates\n")
cat("================================================================\n")

cat("\n--- Model C1: Class-common ED50 estimates ---\n")
ed50_params_C1 <- model_C1_g$BUGSoutput$summary
ed50_rows_C1 <- ed50_params_C1[grep("^ED50", rownames(ed50_params_C1)), ]
print(round(ed50_rows_C1[, c("50%", "2.5%", "97.5%", "Rhat", "n.eff")], 2))

cat("\n--- Model C2: Class-random ED50 estimates ---\n")
ed50_params_C2 <- model_C2_g$BUGSoutput$summary
ed50_class_C2 <- ed50_params_C2[grep("^ED50", rownames(ed50_params_C2)), ]

cat("\nWithin-class SD of ED50:\n")
if ("sd.ED50" %in% rownames(ed50_params_C2)) {
  print(round(ed50_params_C2["sd.ED50", c("50%", "2.5%", "97.5%", "Rhat", "n.eff")], 2))
}

cat("\n--- Agent-level ED50 estimates ---\n")
agent_ed50_C1 <- ed50_params_C1[grep("^ed50\\[", rownames(ed50_params_C1)), ]
agent_ed50_C2 <- ed50_params_C2[grep("^ed50\\[", rownames(ed50_params_C2)), ]

cat("\nModel C1 agent-level ED50:\n")
print(round(agent_ed50_C1[, c("50%", "2.5%", "97.5%")], 2))

cat("\nModel C2 agent-level ED50:\n")
print(round(agent_ed50_C2[, c("50%", "2.5%", "97.5%", "Rhat", "n.eff")], 2))

cat("\n--- For reference: Model A network-wide shared ED50 ---\n")
ed50_A <- model_A_g$BUGSoutput$summary["ed50", c("50%", "2.5%", "97.5%")]
print(round(ed50_A, 2))


# ==============================================================================
# 7. Heterogeneity (tau) Comparison
# ==============================================================================

cat("\n================================================================\n")
cat("  Between-Study Heterogeneity (tau) Across All Models\n")
cat("================================================================\n")

get_tau_summary <- function(model, model_name, is_nma = FALSE) {
  if (is_nma) { tau_samps <- model$jagsresult$BUGSoutput$sims.list$sd
  } else { tau_samps <- model$BUGSoutput$sims.list$sd }
  if (is.matrix(tau_samps)) tau_samps <- tau_samps[, 1]
  data.frame(Model = model_name,
             tau_median = round(median(tau_samps), 3),
             tau_lower = round(quantile(tau_samps, 0.025), 3),
             tau_upper = round(quantile(tau_samps, 0.975), 3),
             stringsAsFactors = FALSE)
}

tau_all <- bind_rows(
  get_tau_summary(model_A_g,      "A: Shared Emax"),
  get_tau_summary(model_B_g,      "B: Agent-specific Emax"),
  get_tau_summary(model_C1_g,     "C1: Class-common ED50"),
  get_tau_summary(model_C2_g,     "C2: Class-random ED50"),
  get_tau_summary(model_Ref_g,    "Ref: Shared Linear"),
  get_tau_summary(model_Lumped_g, "Lumped NMA", is_nma = TRUE)
)
print(tau_all, row.names = FALSE)


# ==============================================================================
# 8. Dose-Response Curves
# ==============================================================================

pred_C1_g <- predict(model_C1_g, E0 = placebo_rate_g, n.doses = 30)
pred_C2_g <- predict(model_C2_g, E0 = placebo_rate_g, n.doses = 30)

jpeg("figures/griselda/GRISELDA_DR_ModelC1.jpeg", width = 14, height = 10, units = "in", res = 300)
plot(pred_C1_g) + ggtitle("GRISELDA Model C1: Class-Common ED50")
dev.off()

jpeg("figures/griselda/GRISELDA_DR_ModelC2.jpeg", width = 14, height = 10, units = "in", res = 300)
plot(pred_C2_g) + ggtitle("GRISELDA Model C2: Class-Random ED50")
dev.off()

jpeg("figures/griselda/GRISELDA_Fitplot_C1.jpeg", width = 14, height = 12, units = "in", res = 300)
fitplot(model_C1_g)
dev.off()

jpeg("figures/griselda/GRISELDA_Fitplot_C2.jpeg", width = 14, height = 12, units = "in", res = 300)
fitplot(model_C2_g)
dev.off()


# ==============================================================================
# 9. Agent Rankings — Unified d_hat Approach
# ==============================================================================

cat("\n================================================================\n")
cat("  Agent Rankings: All Models (unified d_hat)\n")
cat("================================================================\n")

max_rank_doses_g <- griselda_mbnma_class %>%
  filter(agent != "placebo") %>%
  group_by(agent) %>%
  summarise(max_dose = max(dose), .groups = "drop")
max_rank_doses_g_list <- as.list(setNames(max_rank_doses_g$max_dose, max_rank_doses_g$agent))

rank_doses_g    <- setNames(as.list(rep(20, length(active_agents))), active_agents)
rank_doses_10_g <- setNames(as.list(rep(10, length(active_agents))), active_agents)
rank_doses_40_g <- setNames(as.list(rep(40, length(active_agents))), active_agents)

# ---- Rankings at max observed dose (all 6 models) ----
d_hat_A_max   <- extract_d_hat(model_A_g,   "emax_shared", active_agents, max_rank_doses_g_list)
d_hat_B_max   <- extract_d_hat(model_B_g,   "emax_agent",  active_agents, max_rank_doses_g_list)
d_hat_C1_max  <- extract_d_hat(model_C1_g,  c1_model_type, active_agents, max_rank_doses_g_list)
d_hat_C2_max  <- extract_d_hat(model_C2_g,  "emax_agent",  active_agents, max_rank_doses_g_list)
d_hat_Ref_max <- extract_d_hat(model_Ref_g, "linear",      active_agents, max_rank_doses_g_list)

lumped_trt_labs <- model_Lumped_g$trt.labs
lumped_active_labs <- lumped_trt_labs[!grepl("Placebo|placebo", lumped_trt_labs)]
lumped_agent_names <- gsub("_1$", "", lumped_active_labs)
d_hat_Lmp_max <- extract_d_hat(model_Lumped_g, "lumped", lumped_agent_names,
                               max_rank_doses_g_list, is_nma = TRUE)

ranks_A_max   <- compute_ranks_from_d_hat(d_hat_A_max,   "A: Shared Emax")
ranks_B_max   <- compute_ranks_from_d_hat(d_hat_B_max,   "B: Agent-specific Emax")
ranks_C1_max  <- compute_ranks_from_d_hat(d_hat_C1_max,  "C1: Class-common ED50")
ranks_C2_max  <- compute_ranks_from_d_hat(d_hat_C2_max,  "C2: Class-random ED50")
ranks_Ref_max <- compute_ranks_from_d_hat(d_hat_Ref_max, "Ref: Linear")
ranks_Lmp_max <- compute_ranks_from_d_hat(d_hat_Lmp_max, "Lumped NMA")

ranking_all <- bind_rows(ranks_A_max, ranks_B_max, ranks_C1_max,
                         ranks_C2_max, ranks_Ref_max, ranks_Lmp_max)

# Add class information
class_lookup <- Griselda_clean %>%
  filter(agent != "placebo") %>%
  distinct(agent, Class)

ranking_all <- ranking_all %>%
  left_join(class_lookup, by = "agent")

sucra_wide_all <- ranking_all %>%
  select(agent, Class, model, sucra) %>%
  pivot_wider(names_from = model, values_from = sucra) %>%
  arrange(Class, desc(`B: Agent-specific Emax`))

cat("\n=== SUCRA Comparison Across All Models (max observed dose) ===\n")
print(sucra_wide_all, n = Inf, width = Inf)

# ---- Rankings at 20mg FE ----
d_hat_C1_20 <- extract_d_hat(model_C1_g, c1_model_type, active_agents, rank_doses_g)
d_hat_C2_20 <- extract_d_hat(model_C2_g, "emax_agent",  active_agents, rank_doses_g)
ranks_C1_20 <- compute_ranks_from_d_hat(d_hat_C1_20, "C1 20mg")
ranks_C2_20 <- compute_ranks_from_d_hat(d_hat_C2_20, "C2 20mg")


# ==============================================================================
# 10. SUCRA Plot Including Class-Effect Models
# ==============================================================================

ranking_all$agent_clean <- ranking_all$agent %>%
  str_replace_all("_", " ") %>%
  str_to_title()

c2_order <- ranking_all %>%
  filter(model == "C2: Class-random ED50") %>%
  arrange(sucra) %>%
  pull(agent_clean)

ranking_all$agent_clean <- factor(ranking_all$agent_clean, levels = unique(c2_order))

ranking_dr <- ranking_all %>% filter(model != "Lumped NMA")

p_sucra_all <- ggplot(ranking_dr,
                      aes(x = agent_clean, y = sucra, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), name = "SUCRA") +
  scale_fill_manual(
    values = c("A: Shared Emax"        = "#E41A1C",
               "B: Agent-specific Emax" = "#377EB8",
               "C1: Class-common ED50"  = "#FF7F00",
               "C2: Class-random ED50"  = "#984EA3",
               "Ref: Linear"            = "#4DAF4A")
  ) +
  labs(title = "GRISELDA: Agent Rankings at Max Observed FE Dose",
       subtitle = "SUCRA values across model specifications",
       x = NULL, fill = "Model") +
  theme_minimal(base_size = 12)

print(p_sucra_all)
ggsave("figures/griselda/GRISELDA_SUCRA_with_class_effects.jpeg",
       p_sucra_all, width = 12, height = 12, dpi = 300)

# Version including lumped
p_sucra_all_lumped <- ggplot(ranking_all,
                             aes(x = agent_clean, y = sucra, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), name = "SUCRA") +
  scale_fill_manual(
    values = c("A: Shared Emax"        = "#E41A1C",
               "B: Agent-specific Emax" = "#377EB8",
               "C1: Class-common ED50"  = "#FF7F00",
               "C2: Class-random ED50"  = "#984EA3",
               "Ref: Linear"            = "#4DAF4A",
               "Lumped NMA"             = "#999999")
  ) +
  labs(title = "GRISELDA: Agent Rankings at Max Observed FE Dose",
       subtitle = "SUCRA values across all model specifications (unified d_hat ranking)",
       x = NULL, fill = "Model") +
  theme_minimal(base_size = 12)

ggsave("figures/griselda/GRISELDA_SUCRA_all_models.jpeg",
       p_sucra_all_lumped, width = 12, height = 14, dpi = 300)


# ==============================================================================
# 11. Model C2 Rank Stability Across Dose Levels
# ==============================================================================

cat("\n================================================================\n")
cat("  Model C2 Rank Stability Across Dose Levels\n")
cat("================================================================\n")

d_hat_C2_10 <- extract_d_hat(model_C2_g, "emax_agent", active_agents, rank_doses_10_g)
d_hat_C2_40 <- extract_d_hat(model_C2_g, "emax_agent", active_agents, rank_doses_40_g)

ranks_C2_10  <- compute_ranks_from_d_hat(d_hat_C2_10,  "C2 10mg")
ranks_C2_40  <- compute_ranks_from_d_hat(d_hat_C2_40,  "C2 40mg")
ranks_C2_max_tbl <- compute_ranks_from_d_hat(d_hat_C2_max, "C2 max")

rank_dose_C2 <- data.frame(
  Agent  = active_agents,
  `10mg` = round(ranks_C2_10$mean_rank, 2),
  `20mg` = round(ranks_C2_20$mean_rank, 2),
  `40mg` = round(ranks_C2_40$mean_rank, 2),
  Max    = round(ranks_C2_max_tbl$mean_rank, 2),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rank_dose_C2 <- rank_dose_C2[order(rank_dose_C2$`20mg`), ]

cat("\n=== Model C2 Mean Ranks Across Dose Levels ===\n")
print(rank_dose_C2, row.names = FALSE)


# ==============================================================================
# 12. Deviance-Deviance Comparison Plots
# ==============================================================================

cat("\n=== Generating deviance-deviance plots ===\n")

jpeg("figures/griselda/GRISELDA_DevDev_C1_vs_A.jpeg", width = 8, height = 6, units = "in", res = 300)
devdev(model_C1_g, model_A_g, dev.type = "resdev")
dev.off()

jpeg("figures/griselda/GRISELDA_DevDev_C2_vs_A.jpeg", width = 8, height = 6, units = "in", res = 300)
devdev(model_C2_g, model_A_g, dev.type = "resdev")
dev.off()

jpeg("figures/griselda/GRISELDA_DevDev_C2_vs_C1.jpeg", width = 8, height = 6, units = "in", res = 300)
devdev(model_C2_g, model_C1_g, dev.type = "resdev")
dev.off()

jpeg("figures/griselda/GRISELDA_DevDev_C2_vs_B.jpeg", width = 8, height = 6, units = "in", res = 300)
devdev(model_C2_g, model_B_g, dev.type = "resdev")
dev.off()


# ==============================================================================
# 13. Save All Results
# ==============================================================================

dir.create("results/griselda/class_effects", showWarnings = FALSE, recursive = TRUE)

write.csv(comparison_class, "results/griselda/class_effects/model_comparison_with_class.csv",
          row.names = FALSE)
write.csv(sucra_wide_all, "results/griselda/class_effects/sucra_comparison_all_models.csv",
          row.names = FALSE)
write.csv(ranking_all, "results/griselda/class_effects/ranking_results_all_models.csv",
          row.names = FALSE)
write.csv(rank_dose_C2, "results/griselda/class_effects/rank_dose_comparison_C2.csv",
          row.names = FALSE)
write.csv(tau_all, "results/griselda/class_effects/heterogeneity_tau_all_models.csv",
          row.names = FALSE)

save(model_C1_g, model_C2_g,
     lpml_C1_g, lpml_C2_g,
     comparison_class, ranking_all, sucra_wide_all,
     tau_all, rank_dose_C2,
     file = "results/griselda/class_effects/class_effect_models.RData")

cat("\n=== All class-effect results saved to results/griselda/class_effects/ ===\n")


# ==============================================================================
# 14. Summary
# ==============================================================================

cat("\n================================================================\n")
cat("  RESULTS SUMMARY \u2014 Class-Effect Analysis\n")
cat("================================================================\n\n")

cat("--- Model Comparison (DIC / LPML) ---\n")
print(comparison_class[, c("Model", "DIC", "pD", "LPML")], row.names = FALSE)

cat("\n--- Between-Study Heterogeneity ---\n")
print(tau_all, row.names = FALSE)

cat("\n--- Class-Level ED50 Estimates (Model C1) ---\n")
print(round(ed50_rows_C1[, c("50%", "2.5%", "97.5%")], 2))

cat("\n--- Class-Level ED50 Estimates (Model C2: class means) ---\n")
print(round(ed50_class_C2[, c("50%", "2.5%", "97.5%")], 2))

cat("\n--- Top 5 Agents by SUCRA (Model C2 at max dose) ---\n")
top5_C2 <- ranks_C2_max %>% arrange(desc(sucra)) %>% head(5)
print(top5_C2[, c("agent", "sucra", "mean_rank")], row.names = FALSE)

cat("\n--- Bottom 5 Agents by SUCRA (Model C2 at max dose) ---\n")
bot5_C2 <- ranks_C2_max %>% arrange(sucra) %>% head(5)
print(bot5_C2[, c("agent", "sucra", "mean_rank")], row.names = FALSE)

cat("\n================================================================\n")
cat("  Class-effect analysis complete.\n")
cat("================================================================\n")

cat("\n=== All results saved to results/simulation/ ===\n")

cat("\n=== Session Info ===\n")
writeLines(capture.output(sessionInfo()), "results/griselda/session_info_class.txt")

