library(ggplot2)
library(dplyr)

# Load the data
df <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/FinalConcat.csv")

# Filter and group the data
df[df == '.' | df == '-' | df == ' - ' | df == '#DIV/0!' | df == ''| df == ' '] <- NA 

filtered_df <- df %>%
  group_by(PATNO, EVENT_ID) %>%
  filter(n() == 1) %>%
  filter(PRIMDIAG == 1)
         # &!is.na(QUALRep1) )
         # & QUALRep1 == "Positive") ####If only want positive SAAs


################# NLR vs updrs #####################
# Compute averages for each subgroup_simplified
data_avg <- filtered_df %>%
  group_by(PATNO) %>%
  summarise(mean_NLR = mean(NLR, na.rm = TRUE), 
            mean_updrs3_score_on = mean(updrs3_score_on, na.rm = TRUE), 
            subgroup_simplified = first(subgroup_simplified))
# compute correlation and p-value
compute_statistics <- function(df) {
  # Check if there are enough data points and if the variances are not NA or zero
  if (nrow(df) < 2 || is.na(var(df$mean_updrs3_score_on, na.rm = TRUE)) || var(df$mean_updrs3_score_on, na.rm = TRUE) == 0 || 
      is.na(var(df$mean_NLR, na.rm = TRUE)) || var(df$mean_NLR, na.rm = TRUE) == 0) {
    return(data.frame(correlation = NA, p_value = NA))
  }
  # Fit a linear model
  model <- lm(mean_NLR ~ mean_updrs3_score_on, data = df)
  summary_model <- summary(model)
  correlation <- cor(df$mean_updrs3_score_on, df$mean_NLR, use = "pairwise.complete.obs")
  p_value <- summary_model$coefficients[2, 4]
  return(data.frame(correlation = correlation, p_value = p_value))
}
# Apply to each subgroup_simplified and merge results
statistics_df <- data_avg %>%
  group_by(subgroup_simplified) %>%
  do(compute_statistics(.)) %>%
  ungroup()
# Merge the statistics back into the original dataframe
data_avg <- data_avg %>%
  left_join(statistics_df, by = "subgroup_simplified")
# Plot with ggplot2
correlation_plots_NLRvsupdrs3_score_on_Bysubgroup_simplified <- data_avg %>%
  ggplot(aes(x = mean_updrs3_score_on, y = mean_NLR)) +
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) + 
  labs(x = "mean_updrs3_score_on", y = "mean_NLR", title = "Correlation Plot") +
  facet_wrap(~ subgroup_simplified, labeller = label_both) + 
  geom_text(aes(label = paste("r:", round(correlation, 2), "\n", "p:", ifelse(is.na(p_value), "NA", signif(p_value, 2)))), 
            x = Inf, y = Inf, hjust = 1, vjust = 1, 
            size = 4, color = "black", fontface = "italic")
print(correlation_plots_NLRvsupdrs3_score_on_Bysubgroup_simplified)

################# NLR vs DATSCAN #####################
# Compute averages for each subgroup_simplified
filtered_df <- df %>%
  group_by(PATNO, EVENT_ID) %>%
  filter(n() == 1)%>%
  filter(!is.na(mean_putamen) & !is.na(NLR) & !is.na(subgroup_simplified))


filtered_df$mean_putamen <- as.double(filtered_df$mean_putamen)

# Ensure there is data left after filtering
print(nrow(filtered_df))
print(summary(filtered_df))

data_avg <- filtered_df %>%
  group_by(PATNO) %>%
  summarise(mean_NLR = mean(NLR, na.rm = TRUE), 
            mean_mean_putamen = mean(mean_putamen, na.rm = TRUE), 
            subgroup_simplified = first(subgroup_simplified))

print(nrow(data_avg))
print(summary(data_avg))

compute_statistics <- function(df) {
  print(nrow(df))
  print(summary(df))
  
  if (nrow(df) < 2 || is.na(var(df$mean_mean_putamen, na.rm = TRUE)) || var(df$mean_mean_putamen, na.rm = TRUE) == 0 || 
      is.na(var(df$mean_NLR, na.rm = TRUE)) || var(df$mean_NLR, na.rm = TRUE) == 0) {
    return(data.frame(correlation = NA, p_value = NA))
  }
  
  model <- lm(mean_NLR ~ mean_mean_putamen, data = df)
  summary_model <- summary(model)
  
  correlation <- cor(df$mean_mean_putamen, df$mean_NLR, use = "pairwise.complete.obs")
  p_value <- summary_model$coefficients[2, 4]
  
  return(data.frame(correlation = correlation, p_value = p_value))
}

statistics_df <- data_avg %>%
  group_by(subgroup_simplified) %>%
  do(compute_statistics(.)) %>%
  ungroup()

data_avg <- data_avg %>%
  left_join(statistics_df, by = "subgroup_simplified")

print(summary(statistics_df))

correlation_plots_NLRvsmean_putamen_Bysubgroup_simplified <- data_avg %>%
  ggplot(aes(x = mean_mean_putamen, y = mean_NLR)) +
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) + 
  labs(x = "mean_putamen", y = "mean_NLR", title = "Correlation Plot") +
  facet_wrap(~ subgroup_simplified, labeller = label_both) + 
  geom_text(data = filter(data_avg, !is.na(correlation) & !is.na(p_value)),
            aes(label = paste("r:", round(correlation, 2), "\n", "p:", ifelse(is.na(p_value), "NA", signif(p_value, 2)))), 
            x = Inf, y = Inf, hjust = 1, vjust = 1, 
            size = 4, color = "black", fontface = "italic")

print(correlation_plots_NLRvsmean_putamen_Bysubgroup_simplified)

data_avg_filtered <- data_avg %>%
  filter(!is.na(mean_mean_putamen) & !is.na(mean_NLR))
table(data_avg_filtered$subgroup_simplified)


statistics_df <- data_avg_filtered %>%
  group_by(subgroup_simplified) %>%
  do(compute_statistics(.)) %>%
  ungroup()

print(statistics_df)




################# NLR vs fmax #####################
# Compute averages for each subgroup_simplified
data_avg <- filtered_df %>%
  group_by(PATNO) %>%
  summarise(mean_NLR = mean(NLR, na.rm = TRUE),
            mean_FmaxRep_av = mean(FmaxRep_av, na.rm = TRUE),
            subgroup_simplified = first(subgroup_simplified))
# compute correlation and p-value
compute_statistics <- function(df) {
  # Check if there are enough data points and if the variances are not NA or zero
  if (nrow(df) < 2 || is.na(var(df$mean_FmaxRep_av, na.rm = TRUE)) || var(df$mean_FmaxRep_av, na.rm = TRUE) == 0 || 
      is.na(var(df$mean_NLR, na.rm = TRUE)) || var(df$mean_NLR, na.rm = TRUE) == 0) {
    return(data.frame(correlation = NA, p_value = NA))
  }
  # Fit a linear model
  model <- lm(mean_NLR ~ mean_FmaxRep_av, data = df)
  summary_model <- summary(model)
  correlation <- cor(df$mean_FmaxRep_av, df$mean_NLR, use = "pairwise.complete.obs")
  p_value <- summary_model$coefficients[2, 4]
  return(data.frame(correlation = correlation, p_value = p_value))
}
# Apply to each subgroup_simplified and merge results
statistics_df <- data_avg %>%
  group_by(subgroup_simplified) %>%
  do(compute_statistics(.)) %>%
  ungroup()
# Merge the statistics back into the original dataframe
data_avg <- data_avg %>%
  left_join(statistics_df, by = "subgroup_simplified")
# Plot with ggplot2
correlation_plots_NLRvsFmaxRep_av_Bysubgroup_simplified <- data_avg %>%
  ggplot(aes(x = mean_FmaxRep_av, y = mean_NLR)) +
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) + 
  labs(x = "FmaxRep_av", y = "mean_NLR", title = "Correlation Plot") +
  facet_wrap(~ subgroup_simplified, labeller = label_both) + 
  geom_text(aes(label = paste("r:", round(correlation, 2), "\n", "p:", ifelse(is.na(p_value), "NA", signif(p_value, 2)))), 
            x = Inf, y = Inf, hjust = 1, vjust = 1, 
            size = 4, color = "black", fontface = "italic")
print(correlation_plots_NLRvsFmaxRep_av_Bysubgroup_simplified)


################# Neutrophils vs fmax #####################
# Compute averages for each subgroup_simplified
data_avg <- filtered_df %>%
  group_by(PATNO) %>%
  summarise(mean_Neutrophils = mean(Neutrophils, na.rm = TRUE),
            mean_FmaxRep_av = mean(FmaxRep_av, na.rm = TRUE),
            subgroup_simplified = first(subgroup_simplified))
# compute correlation and p-value
compute_statistics <- function(df) {
  # Check if there are enough data points and if the variances are not NA or zero
  if (nrow(df) < 2 || is.na(var(df$mean_FmaxRep_av, na.rm = TRUE)) || var(df$mean_FmaxRep_av, na.rm = TRUE) == 0 || 
      is.na(var(df$mean_Neutrophils, na.rm = TRUE)) || var(df$mean_Neutrophils, na.rm = TRUE) == 0) {
    return(data.frame(correlation = NA, p_value = NA))
  }
  # Fit a linear model
  model <- lm(mean_Neutrophils ~ mean_FmaxRep_av, data = df)
  summary_model <- summary(model)
  correlation <- cor(df$mean_FmaxRep_av, df$mean_Neutrophils, use = "pairwise.complete.obs")
  p_value <- summary_model$coefficients[2, 4]
  return(data.frame(correlation = correlation, p_value = p_value))
}
# Apply to each subgroup_simplified and merge results
statistics_df <- data_avg %>%
  group_by(subgroup_simplified) %>%
  do(compute_statistics(.)) %>%
  ungroup()
# Merge the statistics back into the original dataframe
data_avg <- data_avg %>%
  left_join(statistics_df, by = "subgroup_simplified")
# Plot with ggplot2
correlation_plots_NeutrophilsvsFmaxRep_av_Bysubgroup_simplified <- data_avg %>%
  ggplot(aes(x = mean_FmaxRep_av, y = mean_Neutrophils)) +
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) + 
  labs(x = "FmaxRep_av", y = "mean_Neutrophils", title = "Correlation Plot") +
  facet_wrap(~ subgroup_simplified, labeller = label_both) + 
  geom_text(aes(label = paste("r:", round(correlation, 2), "\n", "p:", ifelse(is.na(p_value), "NA", signif(p_value, 2)))), 
            x = Inf, y = Inf, hjust = 1, vjust = 1, 
            size = 4, color = "black", fontface = "italic")
print(correlation_plots_NeutrophilsvsFmaxRep_av_Bysubgroup_simplified)


################# DATSCAN_PUTAMEN_L vs fmax #####################
# Compute averages for each subgroup_simplified

data_avg <- filtered_df %>%
  group_by(PATNO) %>%
  summarise(mean_DATSCAN_PUTAMEN_L = mean(DATSCAN_PUTAMEN_L, na.rm = TRUE),
            mean_FmaxRep_av = mean(FmaxRep_av, na.rm = TRUE),
            subgroup_simplified = first(subgroup_simplified))
# compute correlation and p-value
compute_statistics <- function(df) {
  # Check if there are enough data points and if the variances are not NA or zero
  if (nrow(df) < 2 || is.na(var(df$mean_FmaxRep_av, na.rm = TRUE)) || var(df$mean_FmaxRep_av, na.rm = TRUE) == 0 || 
      is.na(var(df$mean_DATSCAN_PUTAMEN_L, na.rm = TRUE)) || var(df$mean_DATSCAN_PUTAMEN_L, na.rm = TRUE) == 0) {
    return(data.frame(correlation = NA, p_value = NA))
  }
  # Fit a linear model
  model <- lm(mean_DATSCAN_PUTAMEN_L ~ mean_FmaxRep_av, data = df)
  summary_model <- summary(model)
  correlation <- cor(df$mean_FmaxRep_av, df$mean_DATSCAN_PUTAMEN_L, use = "pairwise.complete.obs")
  p_value <- summary_model$coefficients[2, 4]
  return(data.frame(correlation = correlation, p_value = p_value))
}
# Apply to each subgroup_simplified and merge results
statistics_df <- data_avg %>%
  group_by(subgroup_simplified) %>%
  do(compute_statistics(.)) %>%
  ungroup()
# Merge the statistics back into the original dataframe
data_avg <- data_avg %>%
  left_join(statistics_df, by = "subgroup_simplified")
# Plot with ggplot2
correlation_plots_DATSCAN_PUTAMEN_LvsFmaxRep_av_Bysubgroup_simplified <- data_avg %>%
  ggplot(aes(x = mean_FmaxRep_av, y = mean_DATSCAN_PUTAMEN_L)) +
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) + 
  labs(x = "FmaxRep_av", y = "mean_DATSCAN_PUTAMEN_L", title = "Correlation Plot") +
  facet_wrap(~ subgroup_simplified, labeller = label_both) + 
  geom_text(aes(label = paste("r:", round(correlation, 2), "\n", "p:", ifelse(is.na(p_value), "NA", signif(p_value, 2)))), 
            x = Inf, y = Inf, hjust = 1, vjust = 1, 
            size = 4, color = "black", fontface = "italic")
print(correlation_plots_DATSCAN_PUTAMEN_LvsFmaxRep_av_Bysubgroup_simplified)




##EXPERIMENTING ENDS HERE
##odds ratio_upsit_pctl
filtered_df$upsit_pctl <- as.integer(filtered_df$upsit_pctl)

filtered_df2 <- subset(filtered_df, QUALRep1 == 'Positive' | QUALRep1 == 'Negative')
ggplot(filtered_df2 , aes(x = upsit_pctl, y = NLR)) +
  geom_point(aes(fill = QUALRep1, color = QUALRep1), shape = 21, size = 1) +
  # scale_color_manual(values = c("Positive" = "black", "Negative" = "black")) +
  # scale_fill_manual(values = c("Positive" = "white", "Negative" = "black")) +
  theme_minimal() +
  labs(
    x = "upsit_pctl",
    y = "NLR",
    title = "Odds Ratio Plot",
    fill = "QUALRep1",
    color = "QUALRep1"
  ) +
  facet_wrap(~subgroup_simplified) +
  geom_vline(xintercept = 15, linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(6, 9), linetype = "dashed", color = "black")

##odds ratio_"updrs3_score_on"
ggplot(filtered_df2 , aes(x = updrs3_score_on, y = NLR)) +
  geom_point(aes(fill = QUALRep1, color = QUALRep1), shape = 21, size = 1) +
  theme_minimal() +
  labs(
    x = "updrs3_score_on",
    y = "NLR",
    title = "Odds Ratio Plot",
    fill = "QUALRep1",
    color = "QUALRep1"
  ) +
  facet_wrap(~subgroup_simplified) +
  geom_vline(xintercept = 15, linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(6, 9), linetype = "dashed", color = "black")