library(ggplot2)
library(dplyr)
library(ggsignif)
library(broom) 

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
#av by patient
data_avg <- filtered_df %>%
  group_by(PATNO) %>%
  summarise(mean_NLR = mean(NLR, na.rm = TRUE),
            mean_QUALRep1 = first(QUALRep1),
            mean_subgroup_simplified = first(subgroup_simplified))
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
data_avg <- data_avg %>%
  mutate(mean_NLR = exclude_outliers(mean_NLR))
# Calculate mean and standard deviation for error bars
NLR_stats_by_SAA <- data_avg %>%
  group_by(mean_QUALRep1) %>%
  summarise(SAA_NLR_mean = mean(mean_NLR, na.rm = TRUE),
            SAA_NLR_std = sd(mean_NLR, na.rm = TRUE),
            SAA_n = n())

NLR_stats_by_subgroup <- data_avg %>%
  group_by(mean_subgroup_simplified, mean_QUALRep1) %>%
  summarise(subgroup_NLR_mean = mean(mean_NLR, na.rm = TRUE),
            subgroup_NLR_std = sd(mean_NLR, na.rm = TRUE),
            subgroup_n = n())
# Merge
# data_avg <- left_join(data_avg, NLR_stats_by_SAA, by = "mean_QUALRep1")
# data_avg <- left_join(data_avg, NLR_stats_by_subgroup, by = "mean_subgroup_simplified")

######POST phenoconversion barplot######
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
ggsave("post-phenoconversion_barplot_SAAvNLR_ALL.png", plot = SAAvNLR, width = 8, height = 6, dpi = 300)


######PRE phenoconversion barplot######
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
ggsave("pre-phenoconversion_barplot_SAAvNLR_ALL.png", plot = SAAvNLR, width = 8, height = 6, dpi = 300)

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