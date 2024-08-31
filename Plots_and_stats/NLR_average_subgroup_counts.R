library(ggplot2)
library(dplyr)
library(ggsignif)
library(epitools)
library(ggmosaic)
library(ggpattern)
data <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/FinalConcat.csv")

#######PREP DATA#######
data[data == '.' | data == '-' | data == ' - ' | data == '#DIV/0!' | data == '' | data == ' '] <- NA 
data <- data %>%
  filter(!is.na(subgroup_simplified) &
           !is.na(NLR))
data$subgroup_simplified <- as.factor(data$subgroup_simplified) 
subgroup_simplified_order <- c("Healthy Control", "Genetic", "Sporadic", "RBD +/- Hyposomia","Hyposmia")

# data <- df %>%
#   mutate(subgroup_simplified = case_when(
#     subgroup_simplified == "Hyposmia" ~ "Prodromal", 
#     subgroup_simplified == "RBD +/- Hyposomia" ~ "Prodromal",
#     TRUE ~ subgroup_simplified  # Default case to keep existing value
#   ))

# Average by individual
averaged_data <- data %>%
  group_by(PATNO, subgroup_simplified) %>%
  summarise(NLR = mean(NLR, na.rm = TRUE)) 
# Categorize NLR values
averaged_data <- averaged_data %>%
  mutate(NLR_category = case_when(
    NLR < 0.7 ~ "Below 0.7",
    NLR > 3 ~ "Above 3",
    TRUE ~ "Normal"
  ))
averaged_data <- averaged_data %>%
  filter(!is.na(subgroup_simplified) &
           !is.na(NLR))

#######Central tendency calculations#######
# Count the number of patients in each category
category_counts <- averaged_data %>%
  group_by(NLR_category) %>%
  summarise(count = n())
print(category_counts)
# Calculate SEM
data_summary <- averaged_data %>%
  group_by(subgroup_simplified) %>%
  summarise(mean = mean(NLR), 
            sem = sd(NLR)/sqrt(n()))

#######STATS#######
# One-way ANOVA
anova_result <- aov(NLR ~ subgroup_simplified, data = averaged_data)
summary(anova_result)
# Tukey's HSD
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
str(tukey_result)
print(names(tukey_result$cohort))
# ANOVA for NLR
anova_result_category <- aov(NLR ~ NLR_category, data = averaged_data)
summary(anova_result_category)
# # Annotate the ggplot with p-values (example with one comparison)
# # Assuming comparison between first two cohorts for demonstration
# p_value <- tukey_result$subgroup_simplified["p adj"][5] # Replace with actual comparison index
# plot_label <- paste("p-value:", round(p_value, 3))

#######Plot#######
ggplot(averaged_data, aes(x = factor(subgroup_simplified, subgroup_simplified_order), y = NLR, fill = subgroup_simplified)) +
  geom_violin() +
  geom_point(aes(color = NLR_category), position = position_jitter(width = 0.1), alpha = 0.5) +
  geom_errorbar(data = data_summary, aes(y = mean, ymin = mean - sem, ymax = mean + sem), width = 0.2) +
  geom_signif(comparisons = list(c("Healthy Control", "Sporadic"), c("RBD+/- Hyposomia", "Healthy Control"), c("RBD+/- Hyposomia", "Sporadic")), 
              map_signif_level = TRUE) +
  scale_color_manual(values = c("Below 0.7" = "red", "Above 3" = "blue", "Normal" = "black")) +
  theme_minimal() +
  labs(title = "Violin Plot of NLR by subgroup_simplified", x = "subgroup_simplified", y = "NLR")


#######NLR Threshold Stats#######
#order the categories
averaged_data$subgroup_simplified <- factor(averaged_data$subgroup_simplified,
                                            levels = c("Healthy Control", "Genetic", "Sporadic", "RBD +/- Hyposomia", "Hyposmia"))
averaged_data$NLR_category <- factor(averaged_data$NLR_category,
                                            levels = c( "Above 3", "Normal","Below 0.7" ))
#make a chi table for chi/odds ratio tests
contingency_table <- table(averaged_data$subgroup_simplified, averaged_data$NLR_category)
print(contingency_table)
chisq_test <- chisq.test(contingency_table)
print(chisq_test)
odds_ratio <- oddsratio(contingency_table)
print(odds_ratio)
#Chi and mosaic plots
mosaicplot(contingency_table, main = "Mosaic Plot of NLR Categories by Subgroup",
           col = c("lightblue", "pink", "lightgreen"), ylab = "NLR")
# counts <- as.numeric(contingency_table)
# coords <- expand.grid(x = 1:nrow(contingency_table), y = 1:ncol(contingency_table))
# midpoints <- apply(coords, 1, function(x) {
#   list(x = x[1] - 0.5 + runif(1, 0, 1), y = x[2] - 0.5 + runif(1, 0, 1))
# })
# mapply(function(x, y, label) text(x, y, labels = label, cex = 0.8), 
#        sapply(midpoints, "[[", "x"), 
#        sapply(midpoints, "[[", "y"), 
#        counts)
contingency_df <- as.data.frame(as.table(contingency_table))
ggplot(contingency_df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Chi-square Analysis of NLR Categories by Subgroup",
       x = "Subgroup",
       y = "Count",
       fill = "NLR Category") +
  theme_minimal()
