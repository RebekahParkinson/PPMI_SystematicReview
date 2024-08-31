library(ggplot2) # for plots
library(dplyr) #for filtering and
# library(car)  # for Levene's test
# library(stats)  # for Shapiro-Wilk test
library(tidyverse)
library(rstatix)

################# PREPROCESSING #####################
df <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/FinalConcat.csv")
df[df == '.' | df == '-' | df == ' - ' | df == '#DIV/0!' | df == ''| df == ' '] <- NA 

# Step 1: Group prodromal categories
df_cleaned <- df %>%
  mutate(subgroup_simplified = case_when(
    subgroup_simplified == "RBD +/- Hyposomia" ~ "RBD +/- Hyposomia", 
    subgroup_simplified == "Hyposmia" ~ "RBD +/- Hyposomia",
    TRUE ~ subgroup_simplified
  ))

# Step 2: Filter and group the data by PATNO and EVENT_ID
filtered_df <- df %>%
  group_by(PATNO, EVENT_ID) %>%
  filter(n() == 1) 
# %>%
# filter(PRIMDIAG == 1)
# &!is.na(QUALRep1) 
# & QUALRep1 == "Positive") ####If only want positive SAAs
filtered_df$mean_putamen <- as.numeric(filtered_df$mean_putamen)

#including grouping prodromal
filtered_prodgroup_df <- df_cleaned %>%
  group_by(PATNO, EVENT_ID) %>%
  filter(n() == 1) 
# %>%
  # filter(PRIMDIAG == 1)
         # &!is.na(QUALRep1) 
         # & QUALRep1 == "Positive") ####If only want positive SAAs
filtered_prodgroup_df$mean_putamen <- as.numeric(filtered_prodgroup_df$mean_putamen)


################# NLR vs DATSCAN #####################
# Step 3: Calculate average metrics for each patient (DATSCAN only)
datscan_patient_averages <- filtered_prodgroup_df%>%
  group_by(PATNO) %>%
  summarise(
    mean_NLR = mean(NLR, na.rm = TRUE), 
    mean_SEX  = mean(SEX_x, na.rm = TRUE), 
    mean_age = mean(age, na.rm = TRUE),
    mean_mean_putamen = mean(mean_putamen, na.rm = TRUE)
  )
#Step 4: STATS
#ShapiroWilksTest
shapiro.test(datscan_patient_averages$mean_NLR)
shapiro.test(datscan_patient_averages$mean_mean_putamen)

compute_statistics <- function(datscan_patient_averages) {
  if (nrow(datscan_patient_averages) < 2 || is.na(var(datscan_patient_averages$mean_mean_putamen, na.rm = TRUE)) || var(datscan_patient_averages$mean_mean_putamen, na.rm = TRUE) == 0 || 
      is.na(var(datscan_patient_averages$mean_NLR, na.rm = TRUE)) || var(datscan_patient_averages$mean_NLR, na.rm = TRUE) == 0) {
    return(data.frame(correlation = NA, p_value = NA))
  }
  
  model <- lm(mean_NLR ~ mean_mean_putamen + mean_SEX + mean_age, data = datscan_patient_averages)
  summary_model <- summary(model)
  # correlation <- cor(datscan_patient_averages$mean_mean_putamen, datscan_patient_averages$mean_NLR, use = "pairwise.complete.obs")
  correlation_test <- cor.test(x = datscan_patient_averages$mean_mean_putamen, y = datscan_patient_averages$mean_NLR, method = 'spearman', exact = FALSE)
  correlation <- correlation_test$estimate
  p_value <- summary_model$coefficients[2, 4]
  
  return(data.frame(correlation = correlation, p_value = p_value))
}

# Step 5: Apply statistics function and merge results
statistics_results <- datscan_patient_averages %>%
  do(compute_statistics(.)) %>%
  ungroup()

datscan_patient_averages <- datscan_patient_averages %>%
  cross_join(statistics_results)

# Step 6: Plot correlation between mean_putamen and NLR
plot_mean_putamen_vs_NLR <- datscan_patient_averages %>%
  ggplot(aes(x = mean_mean_putamen, y = mean_NLR)) +
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) + 
  labs(x = "Mean Putamen", y = "Mean NLR", title = "Correlation Plot") +
  geom_text(aes(label = paste("r:", round(correlation, 2), "\n", "p:", ifelse(is.na(p_value), "NA", signif(p_value, 2)))), 
            x = Inf, y = Inf, hjust = 1, vjust = 1, 
            size = 4, color = "black", fontface = "italic")
print(plot_mean_putamen_vs_NLR)


################# NLR vs updrs #####################
#if testing upsit, not updrs
filtered_df$upsit <- as.numeric(filtered_df$upsit)
filtered_prodgroup_df$upsit <- as.numeric(filtered_prodgroup_df$upsit)
#begin
updrs_patient_averages <- filtered_df %>%
  group_by(PATNO) %>%
  summarise(
    mean_NLR = mean(NLR, na.rm = TRUE), 
    mean_SEX  = mean(SEX_x, na.rm = TRUE), 
    mean_age = mean(age, na.rm = TRUE),
    mean_updrs3_score_on = mean(updrs3_score, na.rm = TRUE), 
    subgroup_simplified = first(subgroup_simplified)
  )

# #ShapiroWilksTest
# shapiro_df <- updrs_patient_averages %>%
#   group_by(subgroup_simplified) %>%
#   shapiro_test(mean_NLR)
# shapiro_df_2 <- updrs_patient_averages %>%
#   group_by(subgroup_simplified) %>%
#   shapiro_test(mean_updrs3_score_on)

# Step 9: Compute statistics for updrs3_score_on
compute_statistics_updrs <- function(updrs_patient_averages) {
  if (nrow(updrs_patient_averages) < 2 || is.na(var(updrs_patient_averages$mean_updrs3_score_on, na.rm = TRUE)) || var(updrs_patient_averages$mean_updrs3_score_on, na.rm = TRUE) == 0 || 
      is.na(var(updrs_patient_averages$mean_NLR, na.rm = TRUE)) || var(updrs_patient_averages$mean_NLR, na.rm = TRUE) == 0) {
    return(data.frame(correlation = NA, p_value = NA))
  }
  
  model <- lm(mean_NLR ~ mean_updrs3_score_on + mean_SEX + mean_age, data = updrs_patient_averages)
  summary_model <- summary(model)
  correlation_test <- cor.test(x = updrs_patient_averages$mean_updrs3_score_on, y = updrs_patient_averages$mean_NLR, method = 'spearman', exact = FALSE)
  correlation <- correlation_test$estimate
  p_value <- summary_model$coefficients[2, 4]
  
  return(data.frame(correlation = correlation, p_value = p_value))
}

# Step 10: Apply to each subgroup and merge results
subgroup_statistics <- updrs_patient_averages %>%
  group_by(subgroup_simplified) %>%
  do(compute_statistics_updrs(.)) %>%
  ungroup()

updrs_patient_averages <- updrs_patient_averages %>%
  left_join(subgroup_statistics, by = "subgroup_simplified")

# Step 11: Plot correlation between updrs3_score_on and NLR by subgroup
plot_updrs3_vs_NLR <- updrs_patient_averages %>%
  ggplot(aes(x = mean_updrs3_score_on, y = mean_NLR)) +
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) + 
  labs(x = "Mean UPDRS3 Score", y = "Mean NLR", title = "Correlation Plot with Spearman's rank") +
  facet_wrap(~ subgroup_simplified, labeller = label_both) + 
  geom_text(aes(label = paste("r:", round(correlation, 2), "\n", "p:", ifelse(is.na(p_value), "NA", signif(p_value, 2)))), 
            x = Inf, y = Inf, hjust = 1, vjust = 1, 
            size = 4, color = "black", fontface = "italic")
print(plot_updrs3_vs_NLR)
# 
# # NOT NORMALLY DISTRIBUTED:
# compute_statistics <- function(updrs_patient_averages) {
#   if (nrow(updrs_patient_averages) < 2 || is.na(var(updrs_patient_averages$mean_updrs3_score_on, na.rm = TRUE)) || var(updrs_patient_averages$mean_updrs3_score_on, na.rm = TRUE) == 0 || 
#       is.na(var(updrs_patient_averages$mean_NLR, na.rm = TRUE)) || var(updrs_patient_averages$mean_NLR, na.rm = TRUE) == 0) {
#     return(data.frame(correlation = NA, p_value = NA))
#   }
#   model <- lm(mean_NLR ~ mean_updrs3_score_on + mean_SEX + mean_age, data = updrs_patient_averages)
#   summary_model <- summary(model)
#   correlation_test <- cor.test(x = updrs_patient_averages$mean_updrs3_score_on, y = updrs_patient_averages$mean_NLR, method = 'spearman', exact = FALSE)
#   correlation <- correlation_test$estimate
#   p_value <- summary_model$coefficients[2, 4]
#   
#   return(data.frame(correlation = correlation, p_value = p_value))
# }
# # Apply to each subgroup_simplified and merge results
# statistics_df <- updrs_patient_averages %>%
#   group_by(subgroup_simplified) %>%
#   do(compute_statistics(.)) %>%
#   ungroup()
# # Merge the statistics back into the original dataframe
# updrs_patient_averages <- updrs_patient_averages %>%
#   left_join(statistics_df, by = "subgroup_simplified")
# # Plot
# correlation_plots_NLRvsupdrs3_score_on_Bysubgroup_simplified <- updrs_patient_averages %>%
#   ggplot(aes(x = mean_updrs3_score_on, y = mean_NLR)) +
#   geom_point() + 
#   geom_smooth(method = "lm", se = TRUE) + 
#   labs(x = "mean_updrs3_score_on", y = "mean_NLR", title = "Correlation Plot with Spearman's rank") +
#   facet_wrap(~ subgroup_simplified, labeller = label_both) + 
#   geom_text(aes(label = paste("r:", round(correlation, 2), "\n", "p:", ifelse(is.na(p_value), "NA", signif(p_value, 2)))), 
#             x = Inf, y = Inf, hjust = 1, vjust = 1, 
#             size = 4, color = "black", fontface = "italic")
# print(correlation_plots_NLRvsupdrs3_score_on_Bysubgroup_simplified)
# 
# 
