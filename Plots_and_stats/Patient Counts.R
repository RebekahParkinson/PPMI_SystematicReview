df_FinalConcat <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/FinalConcat.csv")
# df_FinalMerge <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/Final_Merge_v2.csv")
# df_FinalMergeOld <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/Final_Merge.csv")
df_Phenoconverted <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/Phenoconverted.csv")

library(dplyr)

df_FinalConcat[df_FinalConcat %in% c('-', '.', ' - ', '#DIV/0!', '', ' ')] <- NA
df_FinalConcat$mean_putamen <- as.double(df_FinalConcat$mean_putamen)
df_FinalConcat$NLR <- as.numeric(df_FinalConcat$NLR)
df_FinalConcat$FmaxRep_av <- as.numeric(df_FinalConcat$FmaxRep_av)
# Filter for 'prodromal' subgroup and calculate metrics
unique_patients <- df_FinalConcat %>%
  filter(subgroup_simplified == 'Hyposmia' & 
           !is.na(NLR) & 
           !is.na(FmaxRep_av) & 
           # !is.na(normalised_mean_putamen) & 
           mean_putamen <= 0.65) %>%
  distinct(PATNO) %>%
  nrow()
print(unique_patients)

unique_patients <- df_Phenoconverted %>%
  filter(subgroup_simplified == 'RBD +/- Hyposomia' & 
           !is.na(NLR) & 
           !is.na(FmaxRep_av) & 
           # !is.na(normalised_mean_putamen) & 
           normalised_mean_putamen <= 0.65) %>%
  distinct(PATNO) %>%
  nrow()
print(unique_patients)


non_na_counts <- df_FinalConcat %>%
  group_by(subgroup_simplified) %>%
  summarise(non_na_count = sum(!is.na(mean_putamen)))
print(non_na_counts)

summary(df_FinalConcat$mean_putamen)








#Count of PATNO's 
count_unique <- df_FinalConcat %>%
  filter(subgroup== 'Prodromal') %>%
  # filter(ageonset) %>%
  summarise(count_unique_values = n_distinct(PATNO))

count_unique_values <- count_unique$count_unique_values
print(count_unique_values)


colnames(df_FinalConcat)
unique(df_FinalConcat$subgroup)


#Means and Ranges
summary_stats <- df %>%
  filter(condition_column == "A") %>%
  summarise(
    mean_numeric = mean(numeric_column),
    range_numeric = max(numeric_column) - min(numeric_column)
  )
mean_value <- summary_stats$mean_numeric
range_value <- summary_stats$range_numeric
print(mean_value)
print(range_value)


#Non NaNs
count_non_na <- sum(!is.na(df_FinalConcat$NLR[df_FinalMerge$subgroup == "Prodromal"]))
print(count_non_na)

