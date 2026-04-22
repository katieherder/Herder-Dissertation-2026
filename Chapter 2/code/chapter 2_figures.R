################################################################################
# Chapter 2: Empirical Review — Figures
#
# Produces:
#   Figure 2: Dose specification across intervention types (heatmap)
#   Figure 3: Primary analysis approach by intervention category (bar chart)
#
# Data: Manually coded from 22 depression NMAs (2020–2025)
################################################################################

library(dplyr)
library(ggplot2)

setwd("C:/Users/katie/Desktop/Depression-NMA/Chapter 2")

# ==============================================================================
# Data (shared across both figures)
# ==============================================================================

nma_data <- data.frame(
  network = 1:22,
  intervention_type = c(
    "Single Drug Class", "Single Drug Class", "Multiple Drug Classes",
    "Multiple Drug Classes", "Single Drug Class", "Single Drug Class",
    "Single Drug", "Single Drug", "Single Drug", "Single Drug Class",
    "Multiple Drug Classes", "Multiple Drug Classes", "Single Exercise Type",
    "Single Supplement", "Multiple Exercise Types", "Multiple Exercise Types",
    "Multiple Exercise Types", "Multiple Exercise Types",
    "Multiple Exercise Types", "Multiple Exercise Types",
    "Multiple Drug Classes and Therapy Types",
    "Single Drug Classe and Multiple Therapy Types"
  ),
  dose_handling = c(
    "Fixed", "Categorized", "Fixed", "Mixed",
    "Mixed", "Categorized", "Fixed", "Fixed",
    "Categorized", "Fixed", "Fixed", "Fixed",
    "Continuous", "Categorized", "Continuous", "Categorized",
    "Categorized", "Categorized", "Continuous", "Fixed",
    "Categorized", "Continuous"
  ),
  Primary_analysis = c(
    "Split", "Split", "Model-Based", "Split",
    "Split", "Split", "Split", "Split",
    "Split", "Model-Based", "Lumped", "Lumped",
    "Split", "Split", "Lumped", "Lumped",
    "Lumped", "Lumped", "Lumped", "Lumped",
    "Split", "Model-Based"
  )
)


# ==============================================================================
# Figure 2: Dose Specification Across Intervention Types (Heatmap)
# ==============================================================================

heatmap_data <- nma_data %>%
  mutate(
    intervention_group = factor(case_when(
      grepl("Therapy Types", intervention_type) ~
        "Multiple Therapy Types\nCombined with Drugs",
      intervention_type %in% c("Single Drug", "Single Supplement") ~
        "Single Drug/Supplement",
      intervention_type %in% c("Single Drug Class", "Multiple Drug Classes") ~
        "Multiple Drugs",
      grepl("Exercise", intervention_type) ~
        "Exercise/Mindfulness\nInterventions"
    ), levels = c("Multiple Therapy Types\nCombined with Drugs",
                  "Exercise/Mindfulness\nInterventions",
                  "Multiple Drugs",
                  "Single Drug/Supplement")),
    dose_handling = factor(dose_handling,
                           levels = c("Fixed", "Categorized", "Mixed", "Continuous"))
  ) %>%
  group_by(intervention_group, dose_handling) %>%
  summarise(n = n(), .groups = "drop")

fig2 <- ggplot(heatmap_data,
               aes(x = dose_handling, y = intervention_group, fill = n)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = n), size = 4.5, fontface = "bold", color = "black") +
  scale_fill_gradient(low = "#FFF5F0", high = "#CB181D",
                      name = "Number of\nNetworks",
                      breaks = 1:6, limits = c(0, 6)) +
  labs(x = "Dose Handling Method", y = "Intervention Type") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 10, color = "black"),
    axis.text.y  = element_text(size = 10, color = "black"),
    axis.title   = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    panel.grid   = element_blank()
  )

ggsave("figures/Figure 2.2.jpeg", fig2,
       width = 170, height = 100, units = "mm", dpi = 300)

cat("Figure 2 saved.\n")


# ==============================================================================
# Figure 3: Primary Analysis Approach by Intervention Category (Bar Chart)
# ==============================================================================

analysis_data <- nma_data %>%
  mutate(
    intervention_category = factor(case_when(
      grepl("Drug", intervention_type) & !grepl("Therapy", intervention_type) ~
        "Pharmacologic\nOnly",
      grepl("Exercise|Supplement", intervention_type) ~
        "Non-Pharmacologic\nOnly",
      grepl("and", intervention_type) ~
        "Mixed\nInterventions"
    )),
    Primary_analysis = factor(Primary_analysis,
                              levels = c("Model-Based", "Lumped", "Split"))
  ) %>%
  group_by(intervention_category, Primary_analysis) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(intervention_category) %>%
  mutate(total = sum(n)) %>%
  ungroup()

fig3 <- ggplot(analysis_data,
               aes(x = reorder(intervention_category, -total),
                   y = n, fill = Primary_analysis)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.4, width = 0.7) +
  geom_text(data = analysis_data %>%
              group_by(intervention_category) %>%
              summarise(total = sum(n), .groups = "drop"),
            aes(x = intervention_category, y = total,
                label = paste0("n=", total), fill = NULL),
            vjust = -0.5, fontface = "bold", size = 3.5) +
  scale_fill_manual(
    values = c("Model-Based" = "#66C2A5",
               "Lumped"       = "#FC8D62",
               "Split"        = "#8DA0CB"),
    name = "Primary Analysis\nApproach"
  ) +
  labs(x = "Intervention Category", y = "Number of Networks") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     breaks = seq(0, 12, 3)) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x  = element_text(size = 9, color = "black"),
    axis.text.y  = element_text(size = 9, color = "black"),
    axis.title   = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 9, face = "bold", lineheight = 1.2),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(4, "mm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25)
  )

ggsave("figures/Figure 2.3.jpeg", fig3,
       width = 170, height = 90, units = "mm", dpi = 300)

cat("Figure 3 saved.\n")