# Check these:
length(unique(Li$Study))           # how many studies?
table(Li$Intervention_Type)        # what are the agents?
range(Li$METs_min_per_week)        # dose range?
table(Li$Intensity_Code)           # what intensities?

# Critical: do any studies have multiple arms of the same type?
Li %>% group_by(Study, Intervention_Type) %>% 
  filter(n() > 1) %>% 
  select(Study, Intervention_Type, METs_min_per_week)