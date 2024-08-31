df <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/FinalConcat.csv")
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpattern)

df[df %in% c('-', '.', ' - ', '#DIV/0!', '', ' ')] <- NA

# Filter the data
df_filtered <- df %>%
  filter(!is.na(QUALRep1) & 
           QUALRep1 != 'Inconclusive' & 
           !is.na(subgroup_simplified) & 
           !is.na(PRIMDIAG))

# Calculate frequency and proportion within each subgroup
df_summary <- df_filtered %>%
  group_by(subgroup_simplified, EVENT_NUM, PRIMDIAG) %>%
  tally() %>%
  group_by(subgroup_simplified, EVENT_NUM) %>%
  mutate(proportion = n / sum(n))

# Ensure PRIMDIAG is a factor with the specified order
df_summary$PRIMDIAG <- factor(df_summary$PRIMDIAG,
                              levels = c(1, 23, 24, 2, 3, 4, 5, 6, 7, 9, 11, 12, 15, 17, 25, 97))

# Define colors for specific categories
custom_colors <- c("1" = "green", 
                   "23" = "blue", 
                   "24" = "purple")

# Get remaining categories
remaining_categories <- setdiff(levels(df_summary$PRIMDIAG), names(custom_colors))

# Create a named vector of colors
color_scale <- c(custom_colors, 
                 setNames(viridis(length(remaining_categories)), remaining_categories))

# Create the stacked bar plot
ggplot(df_summary, aes(x = EVENT_NUM, y = proportion, fill = PRIMDIAG)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ subgroup_simplified) +
  scale_fill_manual(values = color_scale) +
  theme_minimal() +
  labs(
    x = "Event Number",
    y = "Proportion",
    fill = "Primary Diagnosis",
    title = "Stacked Bar Plot of PRIMDIAG Proportions by EVENT_NUM and Subgroup"
  )



###PLOT 2 (including SAA+)####
df_summary2 <- df_filtered %>%
  group_by(subgroup_simplified, EVENT_NUM, PRIMDIAG, QUALRep1) %>%
  tally() %>%
  group_by(subgroup_simplified, EVENT_NUM, QUALRep1) %>%
  mutate(proportion = n / sum(n))

# Ensure PRIMDIAG is a factor with the specified order
df_summary2$PRIMDIAG <- factor(df_summary2$PRIMDIAG,
                              levels = c(1, 23, 24, 2, 3, 4, 5, 6, 7, 9, 11, 12, 15, 17, 25, 97))

# Define colors for specific categories
custom_colors <- c("1" = "green", 
                   "23" = "blue", 
                   "24" = "purple")

# Get remaining categories
remaining_categories <- setdiff(levels(df_summary2$PRIMDIAG), names(custom_colors))

# Create a named vector of colors
color_scale <- c(custom_colors, 
                 setNames(viridis(length(remaining_categories)), remaining_categories))

# Create the stacked bar plot with patterns
ggplot(df_summary2, aes(x = EVENT_NUM, y = proportion, fill = PRIMDIAG, alpha = QUALRep1)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ subgroup_simplified) +
  scale_fill_manual(values = color_scale) +
  scale_alpha_manual(values = c("Positive" = 1, "Negative" = 0.5)) +
  theme_minimal() +
  labs(
    x = "Event Number",
    y = "Proportion",
    fill = "Primary Diagnosis",
    alpha = "QUALRep1",
    title = "Stacked Bar Plot of PRIMDIAG and QUALRep1 Proportions by EVENT_NUM and Subgroup"
  )


###PLOT 3 (including SAA+) - RAW FREQ####
df_summary3 <- df_filtered %>%
  group_by(subgroup_simplified, EVENT_NUM, PRIMDIAG, QUALRep1) %>%
  tally()

# Ensure PRIMDIAG is a factor with the specified order
df_summary3$PRIMDIAG <- factor(df_summary3$PRIMDIAG,
                              levels = c(1, 23, 24, 2, 3, 4, 5, 6, 7, 9, 11, 12, 15, 17, 25, 97))

# Define colors for specific categories
custom_colors <- c("1" = "green", 
                   "23" = "blue", 
                   "24" = "purple")

# Get remaining categories
remaining_categories <- setdiff(levels(df_summary3$PRIMDIAG), names(custom_colors))

# Create a named vector of colors
color_scale <- c(custom_colors, 
                 setNames(viridis(length(remaining_categories)), remaining_categories))

# Create the stacked bar plot with patterns
ggplot(df_summary3, aes(x = EVENT_NUM, y = n, fill = PRIMDIAG, alpha = QUALRep1)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ subgroup_simplified, scales = "free_y") +
  scale_fill_manual(values = color_scale) +
  scale_alpha_manual(values = c("Positive" = 1, "Negative" = 0.5)) +
  theme_minimal() +
  labs(
    x = "Event Number",
    y = "Frequency",
    fill = "Primary Diagnosis",
    alpha = "QUALRep1",
    title = "Stacked Bar Plot of PRIMDIAG Frequencies by EVENT_NUM and Subgroup"
  )