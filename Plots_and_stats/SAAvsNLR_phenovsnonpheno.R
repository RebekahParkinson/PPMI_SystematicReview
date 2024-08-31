library(ggplot2)
library(dplyr)
library(ggsignif)
library(broom) 
library(ggpubr) # need this for stat_pvalue_manual

# Load the dataa
df <- read.csv("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/FinalConcat.csv")
df[df == '.' | df == '-' | df == ' - ' | df == '#DIV/0!' | df == '' | df == ' '] <- NA 

#######PREP DATA#######
filtered_df <- df %>%
  group_by(PATNO, EVENT_ID) %>%
  filter(n() == 1) %>%
  filter(!is.na(QUALRep1) 
         & !is.na(NLR) 
         & !is.na(subgroup_simplified) 
         & QUALRep1 != "Inconclusive")

#######Central tendency calculations#######
# exclude outliers of NLR based on the IQR
exclude_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

##All
#av by patient
# post_pheno <- filtered_df[filtered_df$PRIMDIAG == 1 & filtered_df$subgroup_simplified %in% c("Sporadic", "Genetic", "Hyposmia", "RBD +/- Hyposomia"), ]
post_pheno <- filtered_df[filtered_df$PRIMDIAG == 1 & filtered_df$subgroup_simplified %in% c("Healthy Control","Sporadic","Genetic", "Hyposmia", "RBD +/- Hyposomia"), ]
pre_pheno <- filtered_df[filtered_df$PRIMDIAG %in% c(23, 24) & filtered_df$subgroup_simplified %in% c("Healthy Control","Genetic", "Hyposmia", "RBD +/- Hyposomia"), ]
pre_pheno <- filtered_df[filtered_df$PRIMDIAG != 1 & filtered_df$subgroup_simplified %in% c( "Genetic","Hyposmia", "RBD +/- Hyposomia"), ]
data_avg <- filtered_df %>%
  group_by(PATNO) %>%
  summarise(mean_NLR = mean(NLR, na.rm = TRUE),
            mean_QUALRep1 = first(QUALRep1),
            mean_subgroup_simplified = first(subgroup_simplified))
data_avg <- data_avg %>% #apply outlier function
  mutate(mean_NLR = exclude_outliers(mean_NLR))
NLR_stats_by_SAA <- data_avg %>% # Mean and sd
  group_by(mean_QUALRep1) %>%
  summarise(SAA_NLR_mean = mean(mean_NLR, na.rm = TRUE),
            SAA_NLR_std = sd(mean_NLR, na.rm = TRUE),
            SAA_n = n())
NLR_stats_by_subgroup <- data_avg %>%
  group_by(mean_subgroup_simplified, mean_QUALRep1) %>%
  summarise(subgroup_NLR_mean = mean(mean_NLR, na.rm = TRUE),
            subgroup_NLR_std = sd(mean_NLR, na.rm = TRUE),
            subgroup_n = n())
SAAvNLR <- ggplot(NLR_stats_by_SAA, aes(x = mean_QUALRep1, y = SAA_NLR_mean)) + ###PLOT
  geom_col(aes(fill = mean_QUALRep1), width = 0.85) +  # Adjusted width and position
  geom_errorbar(aes(ymin = SAA_NLR_mean - SAA_NLR_std,
                    ymax = SAA_NLR_mean + SAA_NLR_std),
                width = 0.3, position = position_dodge(width = 0.6)) +  # Adjusted width and position
  geom_jitter(data = data_avg, aes(x = mean_QUALRep1, y = mean_NLR), 
              width = 0.2, height = 0, color = "white", fill = "black", 
              size = 4, shape = 21, stroke = 0.3, alpha = 0.5) +
  theme(
    axis.text = element_text(size = 8, family = "Arial"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    panel.background = element_blank()
  ) +
  scale_x_discrete(labels = c("Negative", "Positive")) +
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
  scale_fill_manual(values = c("black", "gold"), labels = c("Negative", "Positive"))
t_test <- t.test(mean_NLR ~ mean_QUALRep1, data = data_avg) ##STATS
p_value <- signif(t_test$p.value, digits = 4)
test_result <- data.frame(group1 = "Negative",
                          group2 = "Positive",
                          p.value = p_value)
SAAvNLR <- SAAvNLR + 
  stat_pvalue_manual(test_result, y.position = 5.5,
                     tip.length = 0.05, label.size = 4.5, bracket.size = 0.5)
SAAvNLR #visualise plot
ggsave("/Volumes/PARK2023-Q5758/SAA Systematic Review (Human Data)/Becky/Code_Output/SAAvsNLR_PD.PRIMDIAG1.svg", plot = SAAvNLR, width = 5, height = 10, dpi = 300)








#####BY PRIMDIAG
data_avg <- filtered_df %>%
  group_by(PATNO) %>%
  summarise(mean_NLR = mean(NLR, na.rm = TRUE),
            mean_PRIMDIAG = mean(PRIMDIAG),
            mean_QUALRep1 = first(QUALRep1),
            mean_subgroup_simplified = first(subgroup_simplified)) %>%
  mutate(mean_NLR = exclude_outliers(mean_NLR))
t_test_results <- data_avg %>%
  group_by(mean_PRIMDIAG) %>%
  filter(n_distinct(mean_QUALRep1) == 2) %>%  # Ensure there are exactly 2 levels
  do(tidy(t.test(mean_NLR ~ mean_QUALRep1, data = .))) %>%
  mutate(group1 = "Negative", group2 = "Positive")

NLR_stats_by_PRIMDIAG <- data_avg %>%
  group_by(mean_PRIMDIAG, mean_QUALRep1) %>%
  summarise(SAA_NLR_mean = mean(mean_NLR, na.rm = TRUE),
            SAA_NLR_std = sd(mean_NLR, na.rm = TRUE),
            SAA_n = n())
#PLOT
SAAvNLR <- ggplot(NLR_stats_by_PRIMDIAG, aes(x = mean_PRIMDIAG, y = SAA_NLR_mean, fill = mean_QUALRep1)) +
  geom_col(width = 0.85, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = SAA_NLR_mean - SAA_NLR_std, ymax = SAA_NLR_mean + SAA_NLR_std),
                width = 0.3, position = position_dodge(width = 0.7)) +
  geom_jitter(data = data_avg, aes(x = mean_PRIMDIAG, y = mean_NLR, color = mean_QUALRep1),
              width = 0.2, height = 0, size = 1.5, shape = 21, stroke = 0.3, alpha = 0.5) +  # Removed 'position'
  theme(
    axis.text = element_text(size = 8, family = "Arial"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    panel.background = element_blank()
  ) +
  scale_x_discrete(labels = levels(data_avg$mean_PRIMDIAG)) +
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
  scale_fill_manual(values = c("Negative" = "black", "Positive" = "gold")) +
  scale_color_manual(values = c("Negative" = "black", "Positive" = "gold")) +
  labs(y = "Mean NLR")

t_test_results <- t_test_results %>%
  mutate(x = mean_PRIMDIAG,
         y.position = 5.5,  # Set a constant y-position for the p-value annotations
         label = paste0("p = ", signif(p.value, digits = 4)))  # Create label for p-values
SAAvNLR <- SAAvNLR + 
  stat_pvalue_manual(data = t_test_results, 
                     aes(x = x, y.position = y.position, label = label),
                     tip.length = 0.05, 
                     label.size = 4.5, 
                     bracket.size = 0.5)
SAAvNLR
ggsave("SAAvNLRbyPRIMDIAG.svg", plot = SAAvNLR, width = 8, height = 6, dpi = 300)



##Pre phenoconverted
# pre_pheno <- filtered_df[filtered_df$PRIMDIAG != 1, ] #option 1
pre_pheno <- filtered_df[filtered_df$PRIMDIAG %in% c(23, 24), ]#option 2
pre_pheno_avg <- pre_pheno  %>% #av by patient
  group_by(PATNO) %>%
  summarise(mean_NLR = mean(NLR, na.rm = TRUE),
            mean_QUALRep1 = first(QUALRep1),
            mean_subgroup_simplified = first(subgroup_simplified))
pre_pheno_avg <- pre_pheno_avg %>% #apply outlier function
  mutate(mean_NLR = exclude_outliers(mean_NLR))
pre_pheno_NLR <- pre_pheno_avg %>% # Mean and sd
  group_by(mean_QUALRep1) %>%
  summarise(pre_pheno_SAA_NLR_mean = mean(mean_NLR, na.rm = TRUE),
            pre_pheno_SAA_NLR_std = sd(mean_NLR, na.rm = TRUE),
            SAA_n = n())
pre_pheno_NLR_by_subgroup<- pre_pheno_avg %>%
  group_by(mean_subgroup_simplified, mean_QUALRep1) %>%
  summarise(subgroup_NLR_mean = mean(mean_NLR, na.rm = TRUE),
            subgroup_NLR_std = sd(mean_NLR, na.rm = TRUE),
            subgroup_n = n())
pre_pheno_SAAvsNLR <- ggplot(pre_pheno_NLR, aes(x = mean_QUALRep1, y = pre_pheno_SAA_NLR_mean)) + ###PLOT
  geom_col(aes(fill = mean_QUALRep1), width = 0.85) +  # Adjusted width and position
  geom_errorbar(aes(ymin = pre_pheno_SAA_NLR_mean - pre_pheno_SAA_NLR_std,
                    ymax = pre_pheno_SAA_NLR_mean + pre_pheno_SAA_NLR_std),
                width = 0.3, position = position_dodge(width = 0.6)) +  # Adjusted width and position
  geom_jitter(data = data_avg, aes(x = mean_QUALRep1, y = mean_NLR), 
              width = 0.2, height = 0, color = "white", fill = "black", 
              size = 1.5, shape = 21, stroke = 0.3, alpha = 0.5) +
  theme(
    axis.text = element_text(size = 8, family = "Arial"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    panel.background = element_blank()
  ) +
  scale_x_discrete(labels = c("Negative", "Positive")) +
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
  scale_fill_manual(values = c("black", "gold"), labels = c("Negative", "Positive"))
t_test <- t.test(mean_NLR ~ mean_QUALRep1, data = data_avg) ##STATS
p_value <- signif(t_test$p.value, digits = 4)
test_result <- data.frame(group1 = "Negative",
                          group2 = "Positive",
                          p.value = p_value)
pre_pheno_SAAvsNLR <- pre_pheno_SAAvsNLR  + 
  stat_pvalue_manual(test_result, y.position = 5.5,
                     tip.length = 0.05, label.size = 4.5, bracket.size = 0.5)
pre_pheno_SAAvsNLR  #visualise plot


##Post phenoconverted
post_pheno <- filtered_df[filtered_df$PRIMDIAG == 1, ]
post_pheno_avg <- post_pheno  %> %#av by patient
  group_by(PATNO) %>%
  summarise(mean_NLR = mean(NLR, na.rm = TRUE),
            mean_QUALRep1 = first(QUALRep1),
            mean_subgroup_simplified = first(subgroup_simplified))
post_pheno_avg <- post_pheno_avg %>% #apply outlier function
  mutate(mean_NLR = exclude_outliers(mean_NLR))
post_pheno_NLR <- post_pheno_avg %>% # Mean and sd
  group_by(mean_QUALRep1) %>%
  summarise(post_pheno_SAA_NLR_mean = mean(mean_NLR, na.rm = TRUE),
            post_pheno_SAA_NLR_std = sd(mean_NLR, na.rm = TRUE),
            SAA_n = n())
post_pheno_NLR_by_subgroup<- post_pheno_avg %>%
  group_by(mean_subgroup_simplified, mean_QUALRep1) %>%
  summarise(subgroup_NLR_mean = mean(mean_NLR, na.rm = TRUE),
            subgroup_NLR_std = sd(mean_NLR, na.rm = TRUE),
            subgroup_n = n())
post_pheno_SAAvsNLR <- ggplot(post_pheno_NLR, aes(x = mean_QUALRep1, y = post_pheno_SAA_NLR_mean)) + ###PLOT
  geom_col(aes(fill = mean_QUALRep1), width = 0.85) +  # Adjusted width and position
  geom_errorbar(aes(ymin = post_pheno_SAA_NLR_mean - post_pheno_SAA_NLR_std,
                    ymax = post_pheno_SAA_NLR_mean + post_pheno_SAA_NLR_std),
                width = 0.3, position = position_dodge(width = 0.6)) +  # Adjusted width and position
  geom_jitter(data = data_avg, aes(x = mean_QUALRep1, y = mean_NLR), 
              width = 0.2, height = 0, color = "white", fill = "black", 
              size = 1.5, shape = 21, stroke = 0.3, alpha = 0.5) +
  theme(
    axis.text = element_text(size = 8, family = "Arial"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    panel.background = element_blank()
  ) +
  scale_x_discrete(labels = c("Negative", "Positive")) +
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
  scale_fill_manual(values = c("black", "gold"), labels = c("Negative", "Positive"))
t_test <- t.test(mean_NLR ~ mean_QUALRep1, data = data_avg) ##STATS
p_value <- signif(t_test$p.value, digits = 4)
test_result <- data.frame(group1 = "Negative",
                          group2 = "Positive",
                          p.value = p_value)
post_pheno_SAAvsNLR <- post_pheno_SAAvsNLR + 
  stat_pvalue_manual(test_result, y.position = 5.5,
                     tip.length = 0.05, label.size = 4.5, bracket.size = 0.5)
post_pheno_SAAvsNLR #visualise plot





ggsave("barplot_SAAvNLR_ALL.png", plot = SAAvNLR, width = 8, height = 6, dpi = 300)


######barplot######
SAAvNLR <- ggplot(NLR_stats_by_SAA, aes(x = mean_QUALRep1, y = SAA_NLR_mean)) +
  geom_col(aes(fill = mean_QUALRep1), width = 0.85) +  # Adjusted width and position
  geom_errorbar(aes(ymin = SAA_NLR_mean - SAA_NLR_std,
                    ymax = SAA_NLR_mean + SAA_NLR_std),
                width = 0.3, position = position_dodge(width = 0.6)) +  # Adjusted width and position
  geom_jitter(data = data_avg, aes(x = mean_QUALRep1, y = mean_NLR), 
              width = 0.2, height = 0, color = "white", fill = "black", 
              size = 1.5, shape = 21, stroke = 0.3, alpha = 0.5) +
  theme(
    axis.text = element_text(size = 8, family = "Arial"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    panel.background = element_blank()
  ) +
  scale_x_discrete(labels = c("Negative", "Positive")) +
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
  scale_fill_manual(values = c("black", "gold"), labels = c("Negative", "Positive"))
# scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 5)) 
# geom_text(aes(label = n, y = 5), size = 3, color = "black") 

# Add stats
t_test <- t.test(mean_NLR ~ mean_QUALRep1, data = data_avg)
p_value <- signif(t_test$p.value, digits = 4)
test_result <- data.frame(group1 = "Negative",
                          group2 = "Positive",
                          p.value = p_value)
SAAvNLR <- SAAvNLR + 
  stat_pvalue_manual(test_result, y.position = 5.5,
                     tip.length = 0.05, label.size = 4.5, bracket.size = 0.5)
SAAvNLR
ggsave("barplot_SAAvNLR_ALL.png", plot = SAAvNLR, width = 8, height = 6, dpi = 300)

######barplot_by subgroup######
SAAvNLR_subgroup <- ggplot(NLR_stats_by_subgroup, aes(x = mean_QUALRep1, y = subgroup_NLR_mean)) +
  geom_col(aes(fill = mean_QUALRep1), width = 0.85) +
  geom_errorbar(aes(ymin = subgroup_NLR_mean - subgroup_NLR_std, ymax = subgroup_NLR_mean + subgroup_NLR_std),
                width = 0.3, position = position_dodge(width = 0.6)) +
  facet_wrap(~ mean_subgroup_simplified) +
  geom_jitter(data = data_avg, aes(x = mean_QUALRep1, y = mean_NLR),
              width = 0.2, height = 0, color = "white", fill = "black",
              size = 1.5, shape = 21, stroke = 0.3, alpha = 0.5) +
  theme(
    axis.text = element_text(size = 8, family = "Arial"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10)
  ) +
  scale_x_discrete(labels = c("Negative", "Positive")) +
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
  scale_fill_manual(values = c("black", "gold"), labels = c("Negative", "Positive"))

# Perform t-test for subgroup
t_test_subgroup <- t.test(subgroup_NLR_mean ~ mean_QUALRep1, data = NLR_stats_by_subgroup)
p_value_subgroup <- signif(t_test_subgroup$p.value, digits = 4)
test_result_subgroup <- data.frame(group1 = "Negative",
                                   group2 = "Positive",
                                   p.value = p_value_subgroup)

# Add stats to the plot
SAAvNLR_subgroup <- SAAvNLR_subgroup + 
  stat_pvalue_manual(test_result_subgroup, y.position = max(NLR_stats_by_subgroup$subgroup_NLR_mean) + 3,
                     tip.length = 0.05, label.size = 4.5, bracket.size = 0.5)
SAAvNLR_subgroup

######violin plot######
SAAvNLR_violin <- ggplot(data_avg, aes(x = mean_QUALRep1, y = mean_NLR)) +
  geom_violin(aes(fill = mean_QUALRep1), width = 0.85) +
  geom_jitter(width = 0.2, height = 0, color = "white", fill = "black", 
              size = 1.5, shape = 21, stroke = 0.5, alpha = 0.5) +
  theme(
    axis.text = element_text(size = 8, family = "Arial"), 
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    panel.background = element_blank()
  ) + 
  scale_x_discrete(labels = c("Negative", "Positive")) +
  scale_fill_manual(values = c("black", "gold"), labels = c("Negative", "Positive")) 
# scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 5)) +
# geom_text(data = NLR_stats_by_SAA, aes(label = n, y = 5), size = 3, color = "black")
# Add stats
SAAvNLR_violin <- SAAvNLR_violin + 
  stat_pvalue_manual(test_result, y.position = 6, 
                     tip.length = 0.05, label.size = 4.5, bracket.size = 0.5)
SAAvNLR_violin