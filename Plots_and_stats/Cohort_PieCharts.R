df_FinalConcat <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/FinalConcat.csv")
# df_Phenoconverted <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/Phenoconverted.csv")

library(dplyr)
library(ggplot2)
library(scales)


##########By SAA+##########
# summarise
df_FinalConcat[df_FinalConcat == '-' | df_FinalConcat == ' - ' | df_FinalConcat == '#DIV/0!' | df_FinalConcat == ''| df_FinalConcat == ' '] <- NA 
df_FinalConcat_summary0 <- df_FinalConcat %>%
  group_by(PATNO) %>%
  summarise(mean_QUALRep1 = first(QUALRep1), #keep first instance of qualrep for now, but consider doing 'any instance' like with primdiag below
            subgroup_simplified = first(subgroup_simplified)) %>%
  ungroup() %>%
  filter(!is.na(mean_QUALRep1) & !is.na(subgroup_simplified))

# counts and percentages for the pie charts
df_FinalConcat_pie <- df_FinalConcat_summary0 %>%
  group_by(subgroup_simplified, mean_QUALRep1) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(subgroup_simplified) %>%
  mutate(percentage = count / sum(count) * 100)
#Colours+order
custom_colors <- c("Positive" = "#c2971f", "Negative" = "black", "Inconclusive" = "#939393")
df_FinalConcat_pie$subgroup_simplified <- factor(df_FinalConcat_pie$subgroup_simplified, 
                                                 levels = c("Healthy Control", "Genetic", "Sporadic", "Hyposmia", "RBD +/- Hyposomia"))
#sum to 100
df_FinalConcat_pie <- df_FinalConcat_pie %>%
  group_by(subgroup_simplified) %>%
  mutate(percentage = percentage / sum(percentage) * 100) %>%
  ungroup()
#Plot
p0 <- ggplot(df_FinalConcat_pie, aes(x = "", y = percentage, fill = factor(mean_QUALRep1))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  facet_wrap(~ subgroup_simplified) +
  geom_text(aes(label = paste0(count, " (", round(percentage, 1), "%)")), color = "white", 
            position = position_stack(vjust = 0.5)) +
  theme_void() +
  scale_fill_manual(values = custom_colors) +  # Use the custom colors
  theme(legend.position = "bottom") +
  labs(fill = "Mean QUALRep1")
print(p0)

##########By PRIMDIAG##########
# summarise
# df_FinalConcat[df_FinalConcat == '-' | df_FinalConcat == ' - ' | df_FinalConcat == '#DIV/0!' | df_FinalConcat == ''| df_FinalConcat == ' '] <- NA 
# df_FinalConcat_summary <- df_FinalConcat %>%
#   group_by(PATNO) %>%
#   summarise(PRIMDIAG = last(PRIMDIAG),
#             subgroup_simplified = last(subgroup_simplified)) %>%
#   ungroup() %>%
#   filter(!is.na(PRIMDIAG) & !is.na(subgroup_simplified)) 

#ALTERNATIVE!
df_FinalConcat[df_FinalConcat == '-' | df_FinalConcat == ' - ' | df_FinalConcat == '#DIV/0!' | df_FinalConcat == ''| df_FinalConcat == ' '] <- NA 
df_FinalConcat_summary1 <- df_FinalConcat %>% #filter NaNs
  filter(!is.na(PRIMDIAG) & !is.na(subgroup_simplified))
df_FinalConcat_summary1 <- df_FinalConcat_summary1 %>% #Select only 1 PRIMDIAG values
  group_by(PATNO) %>%
  summarise(PRIMDIAG_1 = any(PRIMDIAG == "1"),
            subgroup_simplified = first(subgroup_simplified)) %>%
  ungroup()
df_FinalConcat_summary1 <- df_FinalConcat_summary1 %>% #Sets true/false = 1/Other in PRIMDIAG
  mutate(PRIMDIAG = ifelse(PRIMDIAG_1, "1", "Other")) %>%
  select(-PRIMDIAG_1)

# counts and percentages for the pie charts
df_FinalConcat_pie <- df_FinalConcat_summary1 %>%
  group_by(subgroup_simplified, PRIMDIAG) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(subgroup_simplified) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
# Colour
custom_colors <- c("1" = "#c2971f", "Other" = "black")
#Order
df_FinalConcat_pie$subgroup_simplified <- factor(df_FinalConcat_pie$subgroup_simplified, 
                                                 levels = c("Healthy Control", "Genetic", "Sporadic", "Hyposmia", "RBD +/- Hyposomia"))
#sum to 100
df_FinalConcat_pie <- df_FinalConcat_pie %>%
  group_by(subgroup_simplified) %>%
  mutate(percentage = percentage / sum(percentage) * 100) %>%
  ungroup()
#plot
p1 <- ggplot(df_FinalConcat_pie, aes(x = "", y = percentage, fill = factor(PRIMDIAG))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  facet_wrap(~ subgroup_simplified) +
  geom_text(aes(label = paste0(count, " (", round(percentage, 1), "%)")), color = "white", 
            position = position_stack(vjust = 0.5)) +
  theme_void() +
  scale_fill_manual(values = custom_colors) +  # Use the custom colors
  theme(legend.position = "bottom") +
  labs(fill = "PRIMDIAG")
print(p1)
# 





##########By mean_putamen##########
df_FinalConcat <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/FinalConcat.csv")
df_FinalConcat[df_FinalConcat %in% c('-', '.', ' - ', '#DIV/0!', '', ' ')] <- NA

#QC: Check putamen values in specific subgroups
hyposmia_data <- df_FinalConcat %>%
  filter(subgroup_simplified == 'Hyposmia') %>%
  select(PATNO, mean_putamen, EVENT_ID, subgroup_simplified)
print(hyposmia_data)
hyposmia_by_patno <- hyposmia_data %>%
  group_by(PATNO) %>%
  summarise(all_mean_putamen_values = list(mean_putamen))
print(hyposmia_by_patno)

# QC1
nan_count_before <- df_FinalConcat %>%
  group_by(PATNO) %>%
  summarise(nan_count = sum(is.na(mean_putamen)))
print(nan_count_before)
# Convert the column mean_putamen to double from strings
df_FinalConcat$mean_putamen <- as.double(df_FinalConcat$mean_putamen)
# QC2
df_FinalConcat$mean_putamen <- as.numeric(df_FinalConcat$mean_putamen)
nan_count_after <- df_FinalConcat %>%
  group_by(PATNO) %>%
  summarise(nan_count = sum(is.na(mean_putamen)))
print(nan_count_after)


# Attempt to summarise the df by grouping mean_putamen
df_FinalConcat_summary <- df_FinalConcat %>%
  group_by(PATNO) %>%
  summarise(mean_putamen = mean(mean_putamen, na.rm = TRUE),
            #STILL NOT SURE IF THIS IS WORKING PROPERLY. COMPLETELY DIFFERENT N IF CHANGE THE MEAN FUNCTION (e.g., to min, or max, or first, etc)
            subgroup_simplified = first(subgroup_simplified)) %>%
  ungroup() %>%
  filter(!is.na(mean_putamen) & !is.na(subgroup_simplified))
# counts and percentages for pie charts
df_FinalConcat_pie <- df_FinalConcat_summary %>%
  mutate(putamen_group = ifelse(mean_putamen <= 0.65,"abnormal1", "normal1")) %>%
  group_by(subgroup_simplified, putamen_group) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(subgroup_simplified) %>%
  mutate(percentage = count / sum(count) * 100)

#colour
custom_colors <- c("abnormal1" = "#c2971f", "normal1" = "black")
#Order
df_FinalConcat_pie$subgroup_simplified <- factor(df_FinalConcat_pie$subgroup_simplified,
                                                 levels = c("Healthy Control", "Genetic", "Sporadic", "Hyposmia", "RBD +/- Hyposomia"))
# sum to 100
df_FinalConcat_pie <- df_FinalConcat_pie %>%
  group_by(subgroup_simplified) %>%
  mutate(percentage = percentage / sum(percentage) * 100) %>%
  ungroup()
# plot
p2 <- ggplot(df_FinalConcat_pie, aes(x = "", y = percentage, fill = factor(putamen_group))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  facet_wrap(~ subgroup_simplified) +
  geom_text(aes(label = paste0(count, " (", round(percentage, 1), "%)")), color = "white", 
            position = position_stack(vjust = 0.5)) +
  theme_void() +
  scale_fill_manual(values = custom_colors) +  # Use the custom colors
  theme(legend.position = "bottom") +
  labs(fill = "Mean Putamen")
print(p2)

##########By mean_putamen_experimental attempt##########
