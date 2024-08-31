library(dplyr)
library(tidyr)
library(ggplot2)
library(OlinkAnalyze)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

# Plasma <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/Proteomics_Plasma_INF_Transposed.csv")
df <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/Proteomics_pheno_plasma.csv")
# Plasma <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/PlasmaNPX_Transposed.csv")
df <- na.omit(df)

#Group the prodromal 
df <- df%>%
  mutate(subgroup_simplified = case_when(
    subgroup_simplified == "RBD +/- Hyposomia" ~ "RBD +/- Hyposomia", 
    subgroup_simplified == "Hyposmia" ~ "RBD +/- Hyposomia",
    TRUE ~ subgroup_simplified
  ))

# Rename columns for consistency
names(df)[names(df) == 'OLINKID'] <- 'OlinkID'
names(df)[names(df) == 'ASSAY'] <- 'Assay'
names(df)[names(df) == 'PANEL'] <- 'Panel'
names(df)[names(df) == 'UNIPROT'] <- 'UniProt'
names(df)[names(df) == 'PATNO'] <- 'SampleID'
event_ids <- unique(df$EVENT_ID)

####RUN STATS AND PLOTS#####
# T-test and z-score function
perform_ttest <- function(df, cohort1_label, cohort2_label, event_id) {
  filtered_data <- df %>%
    filter(subgroup_simplified %in% c(cohort1_label, cohort2_label), EVENT_ID == event_id)
      
  if (nlevels(factor(filtered_data$subgroup_simplified)) == 2 && nrow(filtered_data) > 1) {
    tryCatch({
      # Perform t-test
      ttest_results <- olink_ttest(df = filtered_data, variable = "subgroup_simplified")

      ttest_results <- ttest_results %>%
        mutate(Comparison = paste(cohort1_label, "vs", cohort2_label),
               EVENT_ID = event_id)
      return(ttest_results)
    }, error = function(e) {
      message(paste("Error performing t-test for cohorts", cohort1_label, "vs", cohort2_label, "at", event_id, ":", e$message))
      return(NULL)
    })
  } else {
    message(paste("Skipping t-test for cohorts", cohort1_label, "vs", cohort2_label, "at", event_id, ": insufficient data or missing cohorts"))
    return(NULL)
  }
}
# Cohort comparisons
cohort_comparisons <- list(
  list("Healthy Control", "Sporadic"),
  list("Healthy Control", "RBD +/- Hyposomia"),
  list("Sporadic", "RBD +/- Hyposomia")
  # list("RBD +/- Hyposomia", "Sporadic")  # Reverse
)

# Loop through all combinations and calculate z-scores
ttest_results_list <- list()
for (event_id in event_ids) {
  for (comparison in cohort_comparisons) {
    cohort1 <- comparison[[1]]
    cohort2 <- comparison[[2]]
    result <- perform_ttest(df, cohort1, cohort2, event_id)
    if (!is.null(result)) {
      comparison_label <- paste(cohort1, "vs", cohort2, sep = "_")
      ttest_results_list[[paste(comparison_label, event_id, sep = "_")]] <- result
    }
  }
}

# Combine all results
ttest_results_all <- bind_rows(ttest_results_list)

# Custom color palette ensuring 0 = white
heatmap_colors <- colorRampPalette(c("#ea9139", "white", "#2f618c"))(100)

# Define color breaks to ensure consistency
breaks <- seq(-3, 3, length.out = 101)

# Function to generate heatmap data for each comparison
prepare_heatmap_data <- function(ttest_results, comparison) {
  heatmap_data <- ttest_results %>%
    filter(Comparison == comparison, !EVENT_ID %in% c("V12")) %>%
    select(Assay, EVENT_ID, estimate) %>%
    spread(EVENT_ID, estimate)
  
  # Conditionally invert the estimate values for specific comparison
  if (comparison == "Sporadic vs RBD +/- Hyposomia") {
    heatmap_data[,-1] <- -heatmap_data[,-1]
  }
  
  # Ensure EVENT_ID order with safe subsetting of existing columns only
  expected_cols <- c("Assay", "BL", "V02", "V04", "V06", "V08", "V10")
  existing_cols <- expected_cols[expected_cols %in% colnames(heatmap_data)]
  heatmap_data <- heatmap_data[, existing_cols, drop = FALSE]
  
  # Convert to matrix
  heatmap_matrix <- as.matrix(heatmap_data[,-1])
  rownames(heatmap_matrix) <- heatmap_data$Assay
  return(heatmap_matrix)
}
heatmap_data_sporadic_vs_prodromal <- prepare_heatmap_data(ttest_results_all, "Sporadic vs RBD +/- Hyposomia")

pheatmap(heatmap_data_sporadic_vs_prodromal,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = heatmap_colors,
         breaks = breaks,
         main = "Sporadic vs RBD +/- Hyposomia",
         fontsize = 6)

# Generate heatmap data
heatmap_data_healthy_vs_sporadic <- prepare_heatmap_data(ttest_results_all, "Healthy Control vs Sporadic")
heatmap_data_healthy_vs_prodromal <- prepare_heatmap_data(ttest_results_all, "Healthy Control vs RBD +/- Hyposomia")
heatmap_data_sporadic_vs_prodromal <- prepare_heatmap_data(ttest_results_all, "Sporadic vs RBD +/- Hyposomia")

# Plot each heatmap
pheatmap(heatmap_data_healthy_vs_sporadic,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = heatmap_colors,
         breaks = breaks,
         main = "Healthy Control vs Sporadic",
         fontsize = 6)

pheatmap(heatmap_data_healthy_vs_prodromal,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = heatmap_colors,
         breaks = breaks,
         main = "Healthy Control vs RBD +/- Hyposomia",
         fontsize = 6)

pheatmap(heatmap_data_sporadic_vs_prodromal,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = heatmap_colors,
         breaks = breaks,
         main = "Sporadic vs RBD +/- Hyposomia",
         fontsize = 6)

######COMBINED PLOT#######
#invert sporadic and prodromal
ttest_results_all <- ttest_results_all %>%
  mutate(estimate = ifelse(Comparison == "Sporadic vs RBD +/- Hyposomia", -estimate, estimate))
# Spread data to wide format
ttest_results_long <- ttest_results_all %>%
  select(Assay, EVENT_ID, estimate, Comparison) %>%
  filter(EVENT_ID %in% c("BL", "V02", "V04", "V06", "V08", "V10")) %>%
  spread(Comparison, estimate)  # Spread based on estimate

# Melt the dataframe for ggplot2
ttest_results_melt <- melt(ttest_results_long, id.vars = c("Assay", "EVENT_ID"), variable.name = "Comparison", value.name = "estimate")

# Ensure EVENT_ID order
ttest_results_melt$EVENT_ID <- factor(ttest_results_melt$EVENT_ID, levels = c("BL", "V02", "V04", "V06", "V08", "V10"))

# Order Assay by the mean estimate value
ordered_assays <- ttest_results_melt %>%
  group_by(Assay) %>%
  summarize(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
  arrange(mean_estimate) %>%
  pull(Assay)

# Convert Assay to factor with ordered levels
ttest_results_melt$Assay <- factor(ttest_results_melt$Assay, levels = ordered_assays)

# Plot using ggplot2 with facets
ggplot(ttest_results_melt, aes(x = EVENT_ID, y = Assay, fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "#ce6713", mid = "white", high = "#285577", midpoint = 0,
                       limits = c(-2, 2),  # Extending limits to emphasize extreme values
                       oob = scales::squish) +  # Squish values beyond the limits to the max colors
  facet_wrap(~ Comparison, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 1)) +  # Adjust y-axis text size
  labs(title = "Heatmaps for Cohort Comparisons", x = "EVENT_ID", y = "Assay", fill = "Estimate")


#SAVE DATAFRAME STATS
write.csv(heatmap_data_healthy_vs_prodromal, "/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/pvalues_Plasma_healthy_vs_prodromal.csv", row.names=TRUE, quote=TRUE) 
write.csv(heatmap_data_healthy_vs_sporadic, "/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/pvalues_Plasma_healthy_vs_sporadic.csv", row.names=TRUE, quote=TRUE) 
write.csv(heatmap_data_sporadic_vs_prodromal, "/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/pvalues_Plasma_sporadic_vs_prodromal.csv", row.names=TRUE, quote=TRUE) 