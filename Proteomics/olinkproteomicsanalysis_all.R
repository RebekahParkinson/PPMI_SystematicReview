library(dplyr)
options(dplyr.summarise.inform = TRUE)
library(tidyr)
library(ggplot2)
library(OlinkAnalyze)
library(rlang)

# CSF_INF <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/Proteomics_CSF_INF_Transposed.csv")
CSF_INF <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/Proteomics_pheno_CSF.csv")
# Plasma_INF <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/Proteomics_Plasma_INF_Transposed.csv")
Plasma_INF <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/Proteomics_pheno_plasma.csv")

# CSF_NEU <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/Proteomics_CSF_NEU_Transposed.csv")
# Plasma_NEU <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Code_Output/Proteomics_Plasma_NEU_Transposed.csv")

#Group the prodromal 
CSF_INF <- CSF_INF %>%
  mutate(subgroup_simplified = case_when(
    subgroup_simplified == "RBD +/- Hyposomia" ~ "RBD +/- Hyposomia", 
    subgroup_simplified == "Hyposmia" ~ "RBD +/- Hyposomia",
    TRUE ~ subgroup_simplified
  ))
Plasma_INF <- Plasma_INF %>%
  mutate(subgroup_simplified = case_when(
    subgroup_simplified == "RBD +/- Hyposomia" ~ "RBD +/- Hyposomia", 
    subgroup_simplified == "Hyposmia" ~ "RBD +/- Hyposomia",
    TRUE ~ subgroup_simplified
  ))

######## CSF_INF ######## 
CSF_INF <- na.omit(CSF_INF)
names(CSF_INF)[names(CSF_INF) == 'OLINKID'] <- 'OlinkID'
names(CSF_INF)[names(CSF_INF) == 'ASSAY'] <- 'Assay'
names(CSF_INF)[names(CSF_INF) == 'PANEL'] <- 'Panel'
names(CSF_INF)[names(CSF_INF) == 'UNIPROT'] <- 'UniProt'
names(CSF_INF)[names(CSF_INF) == 'PATNO'] <- 'SampleID'

filtered_CSF_INF <- CSF_INF %>%
  filter(subgroup_simplified %in% c("Healthy Control", "RBD +/- Hyposomia"))
        #& (QUALRep1 == 'Positive'))

averaged_CSF_INF <- filtered_CSF_INF %>%
  ungroup() %>%
  group_by(subgroup_simplified, OlinkID, Assay, Panel, UniProt, SampleID) %>%
  summarise(NPX = mean(NPX, na.rm = TRUE), .groups = 'drop')

if (!is.factor(averaged_CSF_INF$subgroup_simplified) && !is.character(averaged_CSF_INF$subgroup_simplified)) {
  averaged_CSF_INF$subgroup_simplified <- as.factor(averaged_CSF_INF$subgroup_simplified)
}
ttest_results <- olink_ttest(df = averaged_CSF_INF, variable = "subgroup_simplified")

# select names of the top #10 most significant proteins
ttest_significant <- ttest_results %>%
  arrange(p.value) %>%
  head(n = 20) %>%
  pull(OlinkID)

# volcano plot with annotated top #10 most significant proteins
volcano_plot <- olink_volcano_plot(p.val_tbl = ttest_results, olinkid_list = ttest_significant) +
  scale_color_manual(values = c('black', 'green')) +
  ggtitle("CSF_INF, Sporadic vs RBD +/- Hyposomia")
print(volcano_plot)

######## Plasma_INF######## 
Plasma_INF <- na.omit(Plasma_INF)

names(Plasma_INF)[names(Plasma_INF) == 'OLINKID'] <- 'OlinkID'
names(Plasma_INF)[names(Plasma_INF) == 'ASSAY'] <- 'Assay'
names(Plasma_INF)[names(Plasma_INF) == 'PANEL'] <- 'Panel'
names(Plasma_INF)[names(Plasma_INF) == 'UNIPROT'] <- 'UniProt'
names(Plasma_INF)[names(Plasma_INF) == 'PATNO'] <- 'SampleID'

filtered_Plasma_INF <- Plasma_INF %>%
  filter(subgroup_simplified %in% c("Healthy Control", "RBD +/- Hyposomia"))
           # & (QUALRep1 == 'Positive'))

averaged_Plasma_INF <- filtered_Plasma_INF %>%
  ungroup() %>%
  group_by(subgroup_simplified, OlinkID, Assay, Panel, UniProt, SampleID) %>%
  summarise(NPX = mean(NPX, na.rm = TRUE), .groups = 'drop')

if (!is.factor(averaged_Plasma_INF$subgroup_simplified) && !is.character(averaged_Plasma_INF$subgroup_simplified)) {
  averaged_Plasma_INF$subgroup_simplified <- as.factor(averaged_Plasma_INF$subgroup_simplified)
}
ttest_results <- olink_ttest(df = averaged_Plasma_INF, variable = "subgroup_simplified")

# select names of the top #10 most significant proteins
ttest_significant <- ttest_results %>%
  arrange(p.value) %>%
  head(n = 20) %>%
  pull(OlinkID)

# volcano plot with annotated top #10 most significant proteins
volcano_plot <- olink_volcano_plot(p.val_tbl = ttest_results, olinkid_list = ttest_significant) +
  scale_color_manual(values = c('black', 'green')) +
  ggtitle("Plasma_INF, Sporadic vs RBD +/- Hyposomia")
print(volcano_plot)

######## CSF_NEU ######## 
CSF_NEU <- na.omit(CSF_NEU)
names(CSF_NEU)[names(CSF_NEU) == 'OLINKID'] <- 'OlinkID'
names(CSF_NEU)[names(CSF_NEU) == 'ASSAY'] <- 'Assay'
names(CSF_NEU)[names(CSF_NEU) == 'PANEL'] <- 'Panel'
names(CSF_NEU)[names(CSF_NEU) == 'UNIPROT'] <- 'UniProt'
names(CSF_NEU)[names(CSF_NEU) == 'PATNO'] <- 'SampleID'

filtered_CSF_NEU <- CSF_NEU %>%
  filter(subgroup_simplified %in% c("Sporadic", "RBD +/- Hyposomia"))

averaged_CSF_NEU <- filtered_CSF_NEU %>%
  ungroup() %>%
  group_by(subgroup_simplified, OlinkID, Assay, Panel, UniProt, SampleID) %>%
  summarise(NPX = mean(NPX, na.rm = TRUE), .groups = 'drop')

if (!is.factor(averaged_CSF_NEU$subgroup_simplified) && !is.character(averaged_CSF_NEU$subgroup_simplified)) {
  averaged_CSF_NEUF$subgroup_simplified <- as.factor(averaged_CSF_NEU$subgroup_simplified)
}
ttest_results <- olink_ttest(df = averaged_CSF_NEU, variable = "subgroup_simplified")

# select names of the top #10 most significant proteins
ttest_significant <- ttest_results %>%
  arrange(p.value) %>%
  head(n = 20) %>%
  pull(OlinkID)

# volcano plot with annotated top #10 most significant proteins
volcano_plot <- olink_volcano_plot(p.val_tbl = ttest_results, olinkid_list = ttest_significant) +
  scale_color_manual(values = c('black', 'green')) +
  ggtitle("CSF_NEU, Sporadic vs RBD +/- Hyposomia")
print(volcano_plot)

######## Plasma_NEU ######
Plasma_NEU <- na.omit(Plasma_NEU)

names(Plasma_NEU)[names(Plasma_NEU) == 'OLINKID'] <- 'OlinkID'
names(Plasma_NEU)[names(Plasma_NEU) == 'ASSAY'] <- 'Assay'
names(Plasma_NEU)[names(Plasma_NEU) == 'PANEL'] <- 'Panel'
names(Plasma_NEU)[names(Plasma_NEU) == 'UNIPROT'] <- 'UniProt'
names(Plasma_NEU)[names(Plasma_NEU) == 'PATNO'] <- 'SampleID'

filtered_Plasma_NEU <- Plasma_NEU %>%
  filter(subgroup_simplified %in% c("Sporadic", "RBD +/- Hyposomia"))

averaged_Plasma_NEU <- filtered_Plasma_NEU %>%
  ungroup() %>%
  group_by(subgroup_simplified, OlinkID, Assay, Panel, UniProt, SampleID) %>%
  summarise(NPX = mean(NPX, na.rm = TRUE), .groups = 'drop')

if (!is.factor(averaged_Plasma_NEU$subgroup_simplified) && !is.character(averaged_Plasma_NEU$subgroup_simplified)) {
  averaged_Plasma_NEUF$subgroup_simplified <- as.factor(averaged_Plasma_NEU$subgroup_simplified)
}
ttest_results <- olink_ttest(df = averaged_Plasma_NEU, variable = "subgroup_simplified")

# select names of the top #10 most significant proteins
ttest_significant <- ttest_results %>%
  arrange(p.value) %>%
  head(n = 20) %>%
  pull(OlinkID)

# volcano plot with annotated top #10 most significant proteins
volcano_plot <- olink_volcano_plot(p.val_tbl = ttest_results, olinkid_list = ttest_significant) +
  scale_color_manual(values = c('black', 'green')) +
  ggtitle("Plasma_NEU, Sporadic vs RBD +/- Hyposomia")
print(volcano_plot)