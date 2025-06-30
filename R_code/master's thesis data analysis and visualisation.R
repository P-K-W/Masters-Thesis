## Load packages
library(readxl)
library(corrplot)
library(car)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)

## Read data
data <- read_excel("indecies_data.xlsx")
data$Group <- factor(data$Group, levels = c("HC","OCD"))



## Chi-square test for gender
gender_table <- table(data$Group, data$Gender)
print(gender_table)
chi2_result <- chisq.test(gender_table)
print(chi2_result)



## Levene's test between HC & OCD (bonferroni holm corrected)
# Run Levene’s test for each spatial variable
levene_res <- lapply(spatial_vars, function(var){
  fmla <- as.formula(paste(var, "~ Group"))
  tst  <- leveneTest(fmla, data = data, center = median)
  data.frame(
    Measure = var,
    F.value = tst$`F value`[1],
    p.raw   = tst$`Pr(>F)`[1],
    stringsAsFactors = FALSE
  )
}) %>% 
  do.call(rbind, .)

# Apply Bonferroni–Holm correction across all measures
levene_res$p.adj <- p.adjust(levene_res$p.raw, method = "holm")
print(levene_res, row.names = FALSE)



## Levene's test across four groups (bonferroni-holm corrected)

# Reshape data 
data <- read_excel("indecies_data.xlsx") %>%
  mutate(
    Group  = factor(Group,   levels = c("HC","OCD")),
    Gender = factor(Gender,  levels = c(0,1), labels = c("Female","Male"))
  )

spatial_vars <- c("GOF", "medGOF", "weightGOF",
                  "NumDisjfROIs", "isoQ", "FractalDim", "Lacunarity")

long_df <- data %>%
  pivot_longer(
    cols      = all_of(spatial_vars),
    names_to  = "Measure",
    values_to = "Value"
  ) %>%
  mutate(
    GroupGender = factor(
      paste(Group, Gender, sep = "-"),
      levels = c("HC-Female","HC-Male","OCD-Female","OCD-Male")
    )
  )

# Global Brown–Forsythe
bf_global <- long_df %>%
  group_by(Measure) %>%
  do({
    tst <- leveneTest(Value ~ GroupGender, data = ., center = median)
    data.frame(
      Df1     = tst$Df[1],
      Df2     = tst$Df[2],
      F.value = tst$`F value`[1],
      p.value = tst$`Pr(>F)`[1]
    )
  }) %>%
  ungroup() %>%
  mutate(
    p.holm = p.adjust(p.value, method = "holm")
  )

print(bf_global)


# Pairwise Levene tests on FractalDim & isoQ 
pair_measures <- c("FractalDim","isoQ")
groups       <- levels(long_df$GroupGender)
comparisons  <- combn(groups, 2, simplify = FALSE)

pairwise_results <- data.frame(
  Measure = character(),
  group1  = character(),
  group2  = character(),
  df1     = integer(),
  df2     = integer(),
  F.value = double(),
  p.raw   = double(),
  p.adj   = double(),
  stringsAsFactors = FALSE
)

for (m in pair_measures) {
  sub_m <- filter(long_df, Measure == m)
  
  # get a matrix (6 × 4) of stats: df1, df2, F.value, p.raw
  stats_mat <- t(sapply(comparisons, function(cmp) {
    df2 <- filter(sub_m, GroupGender %in% cmp)
    df2$GroupGender <- factor(df2$GroupGender, levels = cmp)
    tst <- leveneTest(Value ~ GroupGender, data = df2, center = median)
    c(
      df1     = tst$Df[1],
      df2     = tst$Df[2],
      F.value = tst$`F value`[1],
      p.raw   = tst$`Pr(>F)`[1]
    )
  }))
  
  stats_df <- as.data.frame(stats_mat, stringsAsFactors = FALSE)
  stats_df$p.adj <- p.adjust(stats_df$p.raw, method = "holm")
  stats_df$group1  <- sapply(comparisons, `[`, 1)
  stats_df$group2  <- sapply(comparisons, `[`, 2)
  stats_df$Measure <- m
  
  # reorder columns by name
  stats_df <- stats_df[, c("Measure","group1","group2",
                           "df1","df2","F.value","p.raw","p.adj")]
  
  pairwise_results <- bind_rows(pairwise_results, stats_df)
}

print(pairwise_results)


## plot (HC vs. OCD)
# Compute mean and variance for each group
summary_df <- data %>%
  group_by(Group) %>%
  summarise(across(all_of(spatial_vars),
                   list(mean = ~mean(.),
                        var  = ~var(.)),
                   .names = "{col}_{fn}")) %>%
  ungroup()

print(summary_df)

# Separate plots for each spatial heterogeneity measures
for (v in spatial_vars) {
  p <- ggplot(data, aes(x = Group, y = .data[[v]], fill = Group)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "black", fill = "white") +
    labs(title = paste0(v, " by Group"),
         y = v, x = NULL) +
    theme_minimal() +
    theme(legend.position = "none")
  print(p)
}

# Combined plot with t-test
ggplot(long_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, color = "black", fill = "white") +
  stat_compare_means(
    aes(group = Group), 
    method      = "t.test",
    comparisons = list(c("HC", "OCD")), 
    label       = "p.signif",
    hide.ns     = TRUE,
    label.size  = 3
  ) +
  facet_wrap(~ Measure, scales = "free_y") +
  labs(title = "Spatial Measures by Group",
       y = "Value", x = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Gender ratio bar plot
data$Gender <- factor(data$Gender,
                      levels = c(0, 1),
                      labels = c("Female", "Male"))

gender_df <- data %>%
  count(Group, Gender) %>%
  group_by(Group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(gender_df, aes(x = Group, y = prop, fill = Gender)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6) +
  # ggplot2 defaults
  # add labels on top of bars:
  geom_text(aes(label = percent(prop, accuracy = 1)),
            position = position_dodge(width = 0.6),
            vjust = -0.5, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(c(0, 0.1))) +
  labs(title = "Gender Composition by Group",
       y = "Proportion", x = NULL, fill = "Gender") +
  theme_minimal()


## Plot (gender+group)

# Pivot to long form 
long_df_gender <- data %>%
  pivot_longer(
    all_of(spatial_vars),
    names_to  = "Measure",
    values_to = "Value"
  )


# Boxplot + jitter + mean + HC vs OCD significance test
ggplot(long_df_gender, aes(x = Group, y = Value, fill = Gender)) +
  geom_boxplot(
    alpha = 0.3,
    outlier.shape = NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = Gender),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    ),
    alpha = 0.5,
    size = 0.8
  ) +
  stat_summary(
    fun      = mean,
    geom     = "point",
    position = position_dodge(width = 0.75),
    shape    = 23,
    size     = 2,
    color    = "black",
    fill     = "white"
  ) +
  stat_compare_means(
    aes(group = Group), 
    method      = "t.test",
    comparisons = list(c("HC", "OCD")), 
    label       = "p.signif",
    hide.ns     = TRUE,
    label.size  = 3
  ) +
  facet_wrap(~ Measure, scales = "free_y") +
  labs(
    title = "Spatial Measures by Group × Gender\n(* indicates HC vs OCD, p < 0.05)",
    y     = "Value",
    x     = NULL,
    fill  = "Gender"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")




## Plot (within group)
# HC vs OCD (uncorrected)
group_test <- compare_means(
  formula   = Value ~ Group,
  data      = long_df_gender,
  method    = "t.test",      
  group.by  = "Measure"      
)


# Within HC (uncorrected)
within_HC <- compare_means(
  formula   = Value ~ Gender,
  data      = filter(long_df_gender, Group == "HC"),
  method    = "t.test",
  group.by  = "Measure"
)

# Within OCD (uncorrected)
within_OCD <- compare_means(
  formula   = Value ~ Gender,
  data      = filter(long_df_gender, Group == "OCD"),
  method    = "t.test",
  group.by  = "Measure"
)

# Apply Bonferroni-Holm correction to group_test
raw_p_group <- group_test$p
holm_p_group <- p.adjust(raw_p_group, method = "holm")

group_test <- group_test %>%
  mutate(
    p_holm         = holm_p_group,
    p.signif_holm  = case_when(
      p_holm < 0.001 ~ "***",
      p_holm < 0.01  ~ "**",
      p_holm < 0.05  ~ "*",
      TRUE           ~ ""
    )
  )

# Apply Bonferroni-Holm correction to within_HC
raw_p_HC <- within_HC$p
holm_p_HC <- p.adjust(raw_p_HC, method = "holm")

within_HC <- within_HC %>%
  mutate(
    p_holm         = holm_p_HC,
    p.signif_holm  = case_when(
      p_holm < 0.001 ~ "***",
      p_holm < 0.01  ~ "**",
      p_holm < 0.05  ~ "*",
      TRUE           ~ ""
    )
  )

# Apply Bonferroni-Holm correction to within_OCD
raw_p_OCD <- within_OCD$p
holm_p_OCD <- p.adjust(raw_p_OCD, method = "holm")

within_OCD <- within_OCD %>%
  mutate(
    p_holm         = holm_p_OCD,
    p.signif_holm  = case_when(
      p_holm < 0.001 ~ "***",
      p_holm < 0.01  ~ "**",
      p_holm < 0.05  ~ "*",
      TRUE           ~ ""
    )
  )

# Compute maximum for each maximum to place the asterisk
ylim_df <- long_df_gender %>%
  group_by(Measure) %>%
  summarise(y_max = max(Value, na.rm = TRUE)) %>%
  ungroup()

# Annotation table for group_test
group_annot <- group_test %>%
  select(Measure, p.signif_holm) %>%
  left_join(ylim_df, by = "Measure") %>%
  mutate(
    x     = 1.5,                 
    y     = y_max * 1.05,         
    label = p.signif_holm
  )

# Annotation table for within_HC
hc_annot <- within_HC %>%
  select(Measure, p.signif_holm) %>%
  left_join(ylim_df, by = "Measure") %>%
  mutate(
    x     = 1,                   
    y     = y_max * 1.10,       
    label = p.signif_holm
  )

# Annotation table for within_OCD
ocd_annot <- within_OCD %>%
  select(Measure, p.signif_holm) %>%
  left_join(ylim_df, by = "Measure") %>%
  mutate(
    x     = 2,                   
    y     = y_max * 1.10,
    label = p.signif_holm
  )

# Stack all annotations 
all_annot <- bind_rows(
  group_annot %>% mutate(type = "Group"),
  hc_annot    %>% mutate(type = "HC"),
  ocd_annot   %>% mutate(type = "OCD")
)


# Final plot
ggplot(long_df_gender, aes(x = Group, y = Value, fill = Gender)) +
  geom_boxplot(
    alpha         = 0.3,
    outlier.shape = NA,
    position      = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = Gender),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    ),
    alpha = 0.5,
    size  = 0.8
  ) +
  stat_summary(
    fun       = mean,
    geom      = "point",
    position  = position_dodge(width = 0.75),
    shape     = 23,
    size      = 2,
    color     = "black",
    fill      = "white"
  ) +
  facet_wrap(~ Measure, scales = "free_y") +
  geom_text(
    data        = all_annot,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size        = 3
  ) +
  labs(
    title    = "Spatial Measures by Group × Gender\n(* p<0.05, ** p<0.01, *** p<0.001)",
    subtitle = "Top asterisk = HC vs OCD; Side asterisks (x=1 or x=2) = Male vs Female within each Group",
    x        = NULL,
    y        = "Value",
    fill     = "Gender",
    color    = "Gender"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")




## Plotting female comparison

long_female <- filter(long_df_gender, Gender == "Female")

female_group_test <- compare_means(
  Value ~ Group,
  data     = long_female,
  method   = "t.test",
  group.by = "Measure"
) %>%
  mutate(
    p_holm        = p.adjust(p, method = "holm"),
    p.signif_holm = case_when(
      p_holm < 0.001 ~ "***",
      p_holm < 0.01  ~ "**",
      p_holm < 0.05  ~ "*",
      TRUE           ~ ""
    )
  )

female_ylim <- long_female %>%
  group_by(Measure) %>%
  summarise(y_max = max(Value, na.rm = TRUE)) %>%
  ungroup()

annot_female <- female_group_test %>%
  filter(group1 == "HC", group2 == "OCD", p.signif_holm != "") %>%
  select(Measure, p.signif_holm) %>%
  left_join(female_ylim, by = "Measure") %>%
  mutate(
    x     = 1.5,
    y     = y_max * 1.05,
    label = p.signif_holm
  )

ggplot(long_female, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  stat_summary(
    fun   = mean,
    geom  = "point",
    shape = 23,
    size  = 2,
    color = "black",
    fill  = "white"
  ) +
  geom_text(
    data        = annot_female,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size        = 3
  ) +
  facet_wrap(~ Measure, scales = "free_y") +
  scale_fill_discrete(
    name   = NULL,
    labels = c("HC female", "OCD female")
  ) +
  labs(
    title = "HC ♀ vs OCD ♀  (asterisk if Holm-adjusted p < 0.05)",
    y     = "Value",
    x     = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_blank(),  
    axis.ticks.x     = element_blank()
  )


##Plotting male comparison
long_male <- long_df_gender %>%
  filter(Gender == "Male")

male_group_test <- compare_means(
  Value ~ Group,
  data     = long_male,
  method   = "t.test",
  group.by = "Measure"
)

male_ylim <- long_male %>%
  group_by(Measure) %>%
  summarise(y_max = max(Value, na.rm = TRUE)) %>%
  ungroup()

annot_male <- male_group_test %>%
  filter(
    group1   == "HC",
    group2   == "OCD",
    p.signif != "ns"       
  ) %>%
  select(Measure, p.signif) %>%
  left_join(male_ylim, by = "Measure") %>%
  mutate(
    x     = 1.5,            
    y     = y_max * 1.05,   
    label = p.signif
  )

ggplot(long_male, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  stat_summary(
    fun      = mean,
    geom     = "point",
    shape    = 23,
    size     = 2,
    color    = "black",
    fill     = "white"
  ) +
  geom_text(
    data        = annot_male,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size        = 3
  ) +
  facet_wrap(~ Measure, scales = "free_y") +
  scale_fill_discrete(
    name   = NULL,
    labels = c("HC male", "OCD male")
  ) +
  labs(
    title = "HC ♂ vs OCD ♂  (asterisk if p < 0.05)",
    y     = "Value",
    x     = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_blank(),  
    axis.ticks.x     = element_blank()
  )
