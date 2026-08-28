################################################################################
################################################################################
###################  MASTER FIGURE DOCUMENT ####################################
################################################################################

library(systemfonts)
library(ggplot2)
library(extrafont)
library("ggpubr")
library(tidyr)
library(dplyr)
library(stringr)
library("plotrix")
library(gridExtra)
library(magick)
library(showtext)
library(patchwork)
library(ggh4x)
library(forcats)
library(ggpattern)
library(rstatix)
library(forcats)
library(glmnet)
library(rsample)
library(purrr)
library(pROC)
library(tidyverse)
library(ggsci)
library(distillery)
library(ggpubr)
library(grid)
library(gridExtra)

##INITIALIZATION ----------------------------------------------------------------
# Times New Roman is usually installed as "Times New Roman"
system_fonts()
showtext_auto()

#Mac
font_add(
  "Times",
  regular    = "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
  bold       = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
  italic     = "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf",
  bolditalic = "/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf"
)

#PC
font_add(
  "Times",
  regular    = "C:/WINDOWS/Fonts/TIMES.TTF",
  bold       = "C:/WINDOWS/Fonts/TIMESBD.TTF",
  italic     = "C:/WINDOWS/Fonts/TIMESI.TTF",
  bolditalic = "C:/WINDOWS/Fonts/TIMESBI.TTF"
)

showtext_opts(dpi = 320)

base_theme <- theme_minimal() + 
  theme(
    text = element_text(family = "Times", color = "black"),
    plot.title = element_text(family = "Times", face = "bold",
                              hjust = 0.5, size = 14),
    legend.title = element_text(family = "Times", hjust = 0.5, face = 'bold', size = 12),
    legend.text = element_text(family = "Times", size = 12),
    axis.text.x = element_text(size = 12, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black'),
    axis.title.x = element_text(size = 12, color = 'black'),
    axis.title.y = element_text(size = 12, color = 'black'),
    strip.text = element_text(size = 12, face = 'bold')
    
  )




## STUDY ONE --------------------------------------------------------------

#Assign Model column as a factor under forcats 
results <- read.csv('lca.fit.csv', header = TRUE) %>%
  select(-X, -cAIC, -resid..df, -likelihood.ratio) %>%
  mutate(Model = case_when(
    Model == '2 class' ~ 2,
    Model == '3 class' ~ 3,
    Model == '4 class' ~ 4,
    Model == '5 class' ~ 5,
    Model == '6 class' ~ 6
  )) %>%
  rename(logLike = log.likelihood, entropy = Entropy, ng = Model)

#Convert to long format

results_long <- results %>%
  pivot_longer(
    cols = c(AIC, BIC, aBIC, logLike, entropy),
    names_to = "Criterion",
    values_to = "Value"
  )

#Plot aBIC, BIC, cAIC and likelihood ratio data for each class 
# Faceted elbow plot
results_long <- results_long %>%
  mutate(Criterion = factor(Criterion,
                            levels = c('aBIC', 'BIC', 'AIC', 'logLike',
                                       'entropy'),
                            labels = c('Sample Adjusted BIC',
                                       'BIC',
                                       'AIC',
                                       'Log-likelihood',
                                       'Entropy')))

aic_plot <- ggplot(results_long %>%
                     filter(Criterion == 'AIC' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'AIC'
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )

bic_plot <- ggplot(results_long %>%
                     filter(Criterion == 'BIC' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'BIC'
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )

abic_plot <- ggplot(results_long %>%
                      filter(Criterion == 'Sample Adjusted BIC' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'Sample Adjusted BIC'
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )

loglike_plot <- ggplot(results_long %>%
                         filter(Criterion == 'Log-likelihood' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'Log-likelihood'
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )

ent_plot <- ggplot(results_long %>%
                     filter(Criterion == 'Entropy' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'Entropy'
  ) +
  scale_y_continuous(limits = c(0.6,0.9)) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )



final_plot <- (aic_plot | loglike_plot) / (bic_plot | abic_plot) / ent_plot
final_plot
#End

ggsave('Thesis Figures/Output/Study 1/elbow_plot.png', final_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')



#CBCL Sub-scale Probability by LCA
summary.plot.data <- read.csv("Study 1 RmD/summary.plot.data.csv", header = TRUE) %>%
  mutate(Class = factor(Class))

CBCL.LCA <- ggplot(summary.plot.data, aes(x=measure, y=probability, group=Class, color = Class)) + 
  geom_point(aes(shape=Class), size = 6) + 
  geom_line(linewidth = 2) + 
  scale_shape_manual(labels=c('2' = 'Low Symptom (88.4%)', '3' = 'Predominantly\nInternalising (7.1%)', '4' = 'Predominantly\nExternalising (2.5%)', '1' = 'Universal\nDifficulties (2.0%)'), 
                     values = c(15, 16, 17, 18)) + 
  scale_color_manual(values = c("2" = "red", '3' = 'green', '4'= 'purple', '1' = 'black'),
                     labels=c('2' = 'Low Symptom (88.4%)', '3' = 'Predominantly\nInternalising (7.1%)', '4' = 'Predominantly\nExternalising (2.5%)', '1' = 'Universal\nDifficulties (2.0%)')) +
  scale_x_discrete(limits = c('Anxious','Withdrawn','Somatic','Social','Thought','Attention','Delinquent','Aggressive')) + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x = "CBCL Subscale", y = "Probability", title = "CBCL Subscale Probability by Latent Class") + 
  base_theme + 
  theme(plot.title = element_blank(), 
        plot.margin = unit(c(22.5,4,8.48,4), "pt"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position="bottom")
CBCL.LCA
#End

ggsave('Thesis Figures/Output/Study 1/cbcl_plot.png', CBCL.LCA, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')


#Expectancies by LCA
mean_valuesLCA <- read.csv('Study 1 RmD/mean_valuesLCA2.csv', header = TRUE)


mean_valuesLCA <- mean_valuesLCA %>%
  filter(!is.na(Class))


mean_valuesLCA$Time = factor(mean_valuesLCA$Time,
                             levels = c(1,3),
                             labels = c("Year 1", "Year 3"))

mean_valuesLCA <- mean_valuesLCA %>%
  mutate(
    Direction = case_when(
      grepl("Positive", Group) ~ "Positive",
      grepl("Negative", Group) ~ "Negative"
    ),
    Substance = case_when(
      grepl("Alcohol", Group)   ~ "Alcohol",
      grepl("Marijuana", Group) ~ "Marijuana",
      grepl("Tobacco", Group)   ~ "Tobacco"
    ),
    Direction = factor(Direction,
                       levels = c("Positive", "Negative")),
    Class = factor(Class, levels = c("Low Symptom", "Internalising",
                                     "Externalising", "Universal"))
  )

Exp.LCA.bar <- ggplot(
  mean_valuesLCA,
  aes(
    x = Class,
    y = Mean_Value,
    fill = Class,
    pattern = Time
  )
) +
  
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = .85),
    width = .75,
    colour = "black",
    
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = .08,
    pattern_spacing = .03,
    
    pattern_key_scale_factor = .6
  ) +
  
  geom_errorbar(
    aes(
      ymin = Mean_Value - SE_Value,
      ymax = Mean_Value + SE_Value
    ),
    position = position_dodge(width = .85),
    width = .2,
    linewidth = .5
  ) +
  
  scale_fill_manual(
    values = c(
      "Low Symptom" = "#D55E00",
      "Internalising" = "#009E73",
      "Externalising" = "#7B61FF",
      "Universal" = "#8a8a8a"
    )
  ) +
  
  scale_pattern_manual(
    values = c(
      "Year 1" = "none",
      "Year 3" = "stripe"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Average Score",
    fill = "Class",
    pattern = "Follow-up"
  ) +
  
  facet_nested(
    ~ Direction + Substance,
    scales = "fixed",
    strip = strip_nested(
      text_x = list(
        element_text(size = 16, face = "bold"),
        element_text(size = 16, face = "bold"),
        element_text(size = 12),
        element_text(size = 12),
        element_text(size = 12),
        element_text(size = 12),
        element_text(size = 12),
        element_text(size = 12)
      ),
      background_x = list(
        element_rect(fill = "grey90"),
        element_rect(fill = "grey90"),
        element_rect(fill = "grey90"),
        element_rect(fill = "grey90"),
        element_rect(fill = "grey90"),
        element_rect(fill = "grey90"),
        element_rect(fill = "grey90"),
        element_rect(fill = "grey90")
      )
    )
  ) +
  
  base_theme +
  
  theme(
    plot.margin = unit(c(1,1,0.5,1), "cm"),
    
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    
    legend.background = element_rect(fill = "white"),
    
    strip.background = element_blank(),
    
    ggh4x.facet.nestline = element_line(colour = "grey40"),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.grid.major.y =
      element_line(colour = "grey85", linewidth = 0.5),
    
    legend.position = "bottom",
    
    legend.title = element_text(face = "bold"),
    
    panel.spacing.x =
      unit(c(0.2,0.2,1.5,0.2,0.2), "cm"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()) +
  annotate(
    "text",
    x = 2.5, y = -Inf,
    label = "Class & Time",
    size = 5,
    vjust = -0.5,
    family = 'Times'
  )

Exp.LCA.bar




ggsave('Thesis Figures/Output/Study 1/expectancies_barplot_naFix.png', Exp.LCA.bar, 'png',
       path = NULL,
       scale = 1,
       width = 11,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')







## STUDY TWO -------------------------------------------------------------------

###### Figures #################################################################
load("output/study2/data_management/elastic_net_data.RData")

###### Sample Characteristics: Categorical Variables ############################
# split by DV
demos_DV0 <- data.2 %>% 
  subset(DV == 0) %>% 
  select(-src_subject_id, -DV)

demos_DV1 <- data.2 %>% 
  subset(DV == 1) %>% 
  select(-src_subject_id, -DV)

# summary stats
obs_num_DV0 <- get_summary_stats(demos_DV0, type = "common") %>% 
  select(variable, n, mean, sd, min, max) %>% 
  left_join(table_names_numeric, by = 'variable') %>% 
  relocate(domain_name, table_name, variable) %>% 
  rename(
    n_0 = n,
    mean_0 = mean,
    sd_0 = sd,
    min_0 = min,
    max_0 = max)

obs_num_DV1 <- get_summary_stats(demos_DV1, type = "common") %>% 
  select(variable, n, mean, sd, min, max) %>% 
  left_join(table_names_numeric, by = 'variable') %>% 
  relocate(domain_name, table_name, variable) %>% 
  rename(
    n_1 = n,
    mean_1 = mean,
    sd_1 = sd,
    min_1 = min,
    max_1 = max)

# merge DV0 and DV1 summary stats
obs_num_DV <- obs_num_DV0 %>% 
  full_join(by = c('domain_name', 'table_name', 'variable'), obs_num_DV1)
names(obs_num_DV)

factor_sum <- function(x) {
  n <- table(x) # frequency
  proportion <- round(prop.table(n), 4) # proportion of non-missing values
  percentage <- (round(prop.table(n), 4))*100
  OUT <- cbind(n, proportion, percentage)
}

# split by DV
p1_DV0 <- factor_sum(demos_DV0$race_4l) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = race_4l', .before = 1) 
p1_DV1 <- factor_sum(demos_DV1$race_4l) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = race_4l', .before = 1) 

p2_DV0 <- factor_sum(demos_DV0$eth_hisp) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = eth_hisp', .before = 1) 
p2_DV1 <- factor_sum(demos_DV1$eth_hisp) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = eth_hisp', .before = 1) 

p3_DV0 <- factor_sum(demos_DV0$income) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = income', .before = 1) 
p3_DV1 <- factor_sum(demos_DV1$income) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = income', .before = 1) 

p4_DV0 <- factor_sum(demos_DV0$religion) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = religion', .before = 1) 
p4_DV1 <- factor_sum(demos_DV1$religion) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = religion', .before = 1) 

p5_DV0 <- factor_sum(demos_DV0$p_edu) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = p_edu', .before = 1) 
p5_DV1 <- factor_sum(demos_DV1$p_edu) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = p_edu', .before = 1) 

p6_DV0 <- factor_sum(demos_DV0$sex_2l) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = sex_2l', .before = 1) 
p6_DV1 <- factor_sum(demos_DV1$sex_2l) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = sex_2l', .before = 1) 

p7_DV0 <- factor_sum(demos_DV0$rec_bin) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = rec_bin', .before = 1) 
p7_DV1 <- factor_sum(demos_DV1$rec_bin) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = rec_bin', .before = 1) 

p8_DV0 <- factor_sum(demos_DV0$exp_sub) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = exp_sub', .before = 1) 
p8_DV1 <- factor_sum(demos_DV1$exp_sub) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = exp_sub', .before = 1) 

p9_DV0 <- factor_sum(demos_DV0$tbi_injury) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = tbi_injury', .before = 1) 
p9_DV1 <- factor_sum(demos_DV1$tbi_injury) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = tbi_injury', .before = 1) 

p10_DV0 <- factor_sum(demos_DV0$kbi_p_grades_in_school) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = kbi_p_grades_in_school', .before = 1) 
p10_DV1 <- factor_sum(demos_DV1$kbi_p_grades_in_school) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = kbi_p_grades_in_school', .before = 1) 

p11_DV0 <- factor_sum(demos_DV0$det_susp) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = det_susp', .before = 1) 
p11_DV1 <- factor_sum(demos_DV1$det_susp) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = det_susp', .before = 1) 

p12_DV0 <- factor_sum(demos_DV0$se_services) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = se_services', .before = 1) 
p12_DV1 <- factor_sum(demos_DV1$se_services) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = se_services', .before = 1) 

p13_DV0 <- factor_sum(demos_DV0$cct) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 0; Variable = cct', .before = 1) 
p13_DV1 <- factor_sum(demos_DV1$cct) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  add_row(variable = 'DV = 1; Variable = cct', .before = 1) 

p1 <- cbind(p1_DV0, p1_DV1)
p2 <- cbind(p2_DV0, p2_DV1) 
p3 <- cbind(p3_DV0, p3_DV1)
p4 <- cbind(p4_DV0, p4_DV1)
p5 <- cbind(p5_DV0, p5_DV1) 
p6 <- cbind(p6_DV0, p6_DV1) 
p7 <- cbind(p7_DV0, p7_DV1) 
p8 <- cbind(p8_DV0, p8_DV1) 
p9 <- cbind(p9_DV0, p9_DV1) 
p10 <- cbind(p10_DV0, p10_DV1) 
p11 <- cbind(p11_DV0, p11_DV1) 
p12 <- cbind(p12_DV0, p12_DV1) 
p13 <- cbind(p13_DV0, p13_DV1) 

obs_factor_DV <- rbind(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13)

write.csv(obs_factor_DV, 'output.alc2/data_management/obs_factor_DV.csv') 

rm(p1_DV0, p1_DV1, p2_DV0, p2_DV1, p3_DV0, p3_DV1, p4_DV0, p4_DV1,
   p5_DV0, p5_DV1, p6_DV0, p6_DV1, p7_DV0, p7_DV1, p8_DV0, p8_DV1,
   p9_DV0, p9_DV1, p10_DV0, p10_DV1, p11_DV0, p11_DV1, p12_DV0, p12_DV1,
   p13_DV0, p13_DV1,
   p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13,
   obs_factor_DV)

# chi-square and Fisher's Test

chisq_race <- tidy(chisq.test(data.2$DV, data.2$race_4l)) %>% 
  mutate(variable = 'Race') %>% 
  relocate(variable) 
chisq_eth_hisp <- tidy(chisq.test(data.2$DV, data.2$eth_hisp)) %>% 
  mutate(variable = 'Ethnicity') %>% 
  relocate(variable) 
chisq_income <- tidy(chisq.test(data.2$DV, data.2$income)) %>% 
  mutate(variable = 'Income') %>% 
  relocate(variable)
chisq_religion <- tidy(chisq.test(data.2$DV, data.2$religion))  
# follow-up with Fisher's Test due to small cell sizes: 
fisher_religion <- as.matrix(
  table(data.2$DV, data.2$religion, useNA = 'ifany'))
fisher.test(fisher_religion, simulate.p.value=TRUE) 
rm(fisher_religion)

chisq_sex <- tidy(chisq.test(data.2$DV, data.2$sex_2l)) %>% 
  mutate(variable = 'Sex') %>% 
  relocate(variable)
chisq_p_edu <- tidy(chisq.test(data.2$DV, data.2$p_edu)) %>% 
  mutate(variable = 'P_edu') %>% 
  relocate(variable)
chisq_rec_bin <- tidy(chisq.test(data.2$DV, data.2$rec_bin)) %>% 
  mutate(variable = 'rec_bin') %>% 
  relocate(variable)
chisq_exp_sub <- tidy(chisq.test(data.2$DV, data.2$exp_sub)) %>% 
  mutate(variable = 'exp_sub') %>% 
  relocate(variable)
chisq_tbi_injury <- tidy(chisq.test(data.2$DV, data.2$tbi_injury)) %>% 
  mutate(variable = 'tbi_injury') %>% 
  relocate(variable)
chisq_grades <- tidy(chisq.test(data.2$DV, data.2$kbi_p_grades_in_school)) %>% 
  mutate(variable = 'grades') %>% 
  relocate(variable) 
chisq_det_susp <- tidy(chisq.test(data.2$DV, data.2$det_susp)) %>% 
  mutate(variable = 'det_susp') %>% 
  relocate(variable)
chisq_se_services <- tidy(chisq.test(data.2$DV, data.2$se_services)) %>% 
  mutate(variable = 'se_services') %>% 
  relocate(variable)
chisq_cct <- tidy(chisq.test(data.2$DV, data.2$cct)) %>% 
  mutate(variable = 'cct') %>% 
  relocate(variable)

obs_factor_diff <- rbind(
  chisq_race, chisq_eth_hisp, chisq_income, chisq_sex,
  chisq_p_edu, chisq_rec_bin, chisq_exp_sub, chisq_tbi_injury, chisq_grades,
  chisq_det_susp, chisq_se_services, chisq_cct) %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  mutate(p.value = ifelse(p.value == 0.000, '<0.001', p.value))

write.csv(obs_factor_diff, 'output.alc2/data_management/obs_factor_diff.csv')

rm(demos_DV0, demos_DV1, 
   chisq_race, chisq_eth_hisp, chisq_income, chisq_religion,
   chisq_sex, chisq_p_edu, chisq_rec_bin, chisq_exp_sub, chisq_tbi_injury, 
   chisq_grades, chisq_det_susp, chisq_se_services, chisq_cct,
   obs_factor_diff)

###### Horizontal Boxplot ######################################################

# update names for plot
m1_coef_top5 <- m1_coef_top5 %>%
  mutate(across(where(is.character), ~ .x %>%
                  str_replace_all("\\bBehaviors\\b", "Behaviours") %>%
                  str_replace_all("\\bBehavior\\b", "Behaviour"))) %>%
  mutate(Domain = ifelse(Domain == 'Culture & Environment',
         'Community & School', Domain)) %>%
  mutate(Domain = ifelse(Variable == 'CRPF', 'Community & School', Domain),
         Variable = ifelse(Variable == 'CRPF', 'Perceived Access to Substances', Variable),
         Domain = ifelse(Domain == 'Neurocognitive Factors', 'Cognitive Ability', Domain),
         Domain = ifelse(Domain == 'Physical Health', 'Health and Physical Activity', Domain),
         Domain = ifelse(Domain == 'Mental Health', 'Psychological Traits', Domain))

top_predictors <- m1_coef_top5 %>%
  mutate(abs_mean_coef = abs(`Average Across Splits`)) %>%
  slice_max(order_by = abs_mean_coef, n = 25) %>%
  mutate(Variable = ifelse(Variable == 'Religion - Mormon (Church of Jesus Christ of Latter Day Saints/LDS)', 'Religion - Church of Jesus Christ of Latter Day Saints', Variable))



Fig_1 <- ggplot(top_predictors, 
                aes(fill = Domain, y = `Average Across Splits`, 
                    x = fct_reorder(Variable, abs_mean_coef))) + 
  geom_bar(stat = 'identity') +
  geom_text(aes(label = round(`Average Across Splits`, 2)),
            hjust = ifelse(top_predictors$`Average Across Splits` > 0, -0.2, 1.2),  # text offset depends on sign
            size = 5,
            family = 'Times') +
  coord_flip() +
  labs(y = 'Coefficient', x = 'Predictor') +
  scale_y_continuous(limits = c(-2, 1)) +
  base_theme +
  theme(axis.title.x = element_text(
          color="black", size = 14, face = 'bold'),
        axis.title.y = element_text(
          color="black", size = 14, face = 'bold'),
        axis.text.y = element_text(
          color="black", size=14),
        axis.text.x = element_text(
          color="black", size=14),
        legend.text = element_text(
          color="black", size = 14),
        legend.title = element_text(
          color="black", size = 16))
Fig_1

ggsave('Thesis Figures/Output/Study 2/study2_horizontal_boxplot.png', Fig_1, 'png',
       path = NULL,
       scale = 1,
       width = 42,
       height = 29.7,
       units = 'cm',
       dpi = 320,
       limitsize = FALSE,
       bg = NULL,
       create.dir = FALSE)


## STUDY THREE -----------------------------------------------------------------
################################################################################


load("output/study3/study3_data.RData")

###### Visualize with Neuroconductor ###########################################
# Install if not already
install.package(devtools)
library(devtools)
Sys.setenv(GITHUB_PAT = "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")

install.packages("remotes")
install.packages("neurobase",
                 repos = c("https://muschellij2.r-universe.dev", "https://cloud.r-project.org")
)
install.packages("RNiftyReg")

library(neurobase)
library(oro.nifti)
library(RNiftyReg)
library(png)
library(grDevices)

yeo17 <- readNIfTI("BrainAtlas/Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz", reorient = FALSE)

# Load a structural image (e.g., MNI152 T1-weighted)
t1 <- readNIfTI("BrainAtlas/FSL_MNI152_FreeSurferConformed_1mm.nii.gz", reorient = FALSE)


# Network names and corresponding Yeo 17 labels
networks <- list(
  DMN = c(14, 15, 16, 17),
  SN = c(7, 8),
  FPN = c(11, 12, 13)
)

# Colors for visualization
colors <- c(DMN = "red", SN = "blue", FPN = "green")

# Loop through networks
for (net_name in names(networks)) {
  
  # Binary mask for this network
  mask_array <- as.numeric(yeo17[] %in% networks[[net_name]])
  mask <- yeo17
  mask[] <- mask_array
  
  if (sum(mask) == 0) {
    warning(paste("Network", net_name, "has no voxels. Skipping."))
    next
  }
  
  # Pick axial slice (orthogonal to z-axis)
  z <- floor(dim(t1)[2] / 2)  # [, z, ] indexing
  t1_slice <- t1[, z, ]
  mask_slice <- mask[, z, ]
  
  # Create brain mask (all voxels inside the brain)
  brain_mask <- t1_slice != 0
  
  # Normalize T1 slice to 0–1
  t1_norm <- t1_slice / max(t1_slice, na.rm = TRUE)
  
  # Initialize RGBA array
  img <- array(0, dim = c(nrow(t1_norm), ncol(t1_norm), 4))  # R,G,B,A
  img[,,1] <- t1_norm  # Red
  img[,,2] <- t1_norm  # Green
  img[,,3] <- t1_norm  # Blue
  img[,,4] <- brain_mask  # Alpha channel: transparent outside brain
  
  # Apply network overlay: color only where mask_slice==1
  overlay_color <- col2rgb(adjustcolor(colors[net_name], alpha.f = 0.5)) / 255
  overlay_alpha <- mask_slice != 0
  
  for (i in 1:3) {  # RGB channels
    img[,,i][overlay_alpha] <- overlay_color[i]
  }
  img[,,4][overlay_alpha] <- 0.5  # semi-transparent overlay
  
  img <- aperm(img, c(2,1,3))       # swap rows/cols
  img <- img[nrow(img):1,,] # flip vertically
  
  # Save PNG
  writePNG(img, target = paste0("output/study3/", net_name, "_axial.png"))
}

###### Create Triple Network Triangle Graph ####################################
library(igraph)
library(ggraph)
library(patchwork)
library(rlang)
library(ggimage)


# 1. Compute mean connectivity per DV

mean_conn <- data.2 %>%
  group_by(DV) %>%
  summarise(across(DMN_DMN:SN_SN, \(x) mean(x, na.rm = TRUE)))

uninitiated <- mean_conn %>% filter(DV == "Uninitiated") %>% select(-DV)
initiated   <- mean_conn %>% filter(DV == "Initiated") %>% select(-DV)
# 2. Convert to edge list

make_edges <- function(df_row, color = NULL) {
  df_long <- df_row %>%
    pivot_longer(cols = everything(), names_to = "connection", values_to = "weight") %>%
    separate(connection, into = c("from_temp","to_temp"), sep = "_") %>%
    rename(from = from_temp, to = to_temp) %>%
    mutate(weight = abs(as.numeric(weight)))
  
  if(!is.null(color)) df_long$color <- color
  return(df_long)
}

edges_uninit <- make_edges(uninitiated, color = "black")   # baseline black
edges_init   <- make_edges(initiated)                       # temporary, will add color


# 3. Add difference-based color to initiated edges

# Compute difference
diff_conn <- initiated - uninitiated
diff_conn <- diff_conn %>%
  pivot_longer(cols = everything(), names_to = "connection", values_to = "diff") %>%
  separate(connection, into = c("from","to"), sep="_") %>%
  mutate(diff = as.numeric(diff))

# Separate loops vs non-loops
edges_links <- diff_conn %>% filter(from != to)
edges_loops <- diff_conn %>% filter(from == to)

# Duplicate non-loop edges to get bidirectional edges
diff_conn <- bind_rows(edges_links, edges_links %>% rename(from_tmp = from) %>%
                         rename(from = to, to = from_tmp) %>% select(from, to, diff))

# Combine with loops
diff_conn <- bind_rows(diff_conn, edges_loops)

diff_conn <- diff_conn %>%
  mutate(
    sig_edge = (from == "DMN" & to == "SN") |
      (from == "FPN" & to == "SN") |
      (from == "SN" & to == "DMN") |
      (from == "SN" & to == "FPN")
  )

diff_conn <- diff_conn %>%
  mutate(
    diff_scaled = ifelse(sig_edge,
                         scales::rescale(diff, to = c(6, 7)),  # amplify significant edges
                         scales::rescale(diff, to = c(2,3)))   # keep other edges subtle
  )

edges_init <- edges_init %>%
  left_join(diff_conn %>% select(from, to, diff, diff_scaled), by = c("from","to"))

edges_init <- edges_init %>%
  mutate(
    # Identify significant edges
    sig_edge = (from == "DMN" & to == "SN") |
      (from == "FPN" & to == "SN") |
      (from == "SN" & to == "DMN") |
      (from == "SN" & to == "FPN"),
    edge_color = ifelse(sig_edge, 'red', 'black')
  )

# Make a named vector for loop angles
loop_angles <- c(DMN = -90, FPN = -90, SN = 90)

# Add a column to the edge dataframe before creating the graph
edges_init <- edges_init %>%
  mutate(loop_angle = ifelse(from == to, loop_angles[from], NA))

network_positions <- tibble::tibble(
  network = c("DMN", "FPN", "SN"),
  x = c(-1, 1, 0),
  y = c(1, 1, -1),
  img = paste0("output/study3/", network, "_axial.png")  # PNG files saved from previous step
)

# 4. Function to plot network

# Apply Times New Roman globally in the theme
base_theme <- theme_void() + 
  theme(
    text = element_text(family = "Times New Roman", color = "black"),
    plot.title = element_text(family = "Times New Roman", face = "bold", size = 24,
                              hjust = 0.5),
    legend.title = element_text(family = "Times New Roman", hjust = 0.5, size = 16),
    legend.title.position = 'top',
    legend.text = element_text(family = "Times New Roman", size = 14),
    legend.key.width = unit(2, "in"),
    plot.background   = element_rect(fill = NA, colour = NA),  # no background
    panel.background  = element_rect(fill = NA, colour = NA),  # no panel
    legend.background = element_rect(fill = NA, colour = NA),  # transparent legend
    legend.key        = element_rect(fill = NA, colour = NA)   # transparent legend keys
  )


# Initiated plot (edges colored by difference)
plot_network_init <- function(edges, title) {
  # Create graph
  g <- graph_from_data_frame(edges[, c("from","to","weight","diff","sig_edge","edge_color", 'loop_angle')], directed = TRUE)
  
  # Layout with scaled coordinates
  layout <- create_layout(g, layout = "circle")
  
  # Rotate coordinates by angle theta
  theta <- -5*pi/6  # 60 degrees clockwise, adjust as needed
  x_rot <- layout$x * cos(theta) - layout$y * sin(theta)
  y_rot <- layout$x * sin(theta) + layout$y * cos(theta)
  layout$x <- x_rot
  layout$y <- y_rot
  
  x_center <- mean(range(layout$x))
  y_center <- mean(range(layout$y))
  
  layout$x <- layout$x - x_center
  layout$y <- layout$y - y_center
  
  network_positions <- tibble::tibble(
    network = layout$name,
    x = layout$x,
    y = layout$y,
    img = paste0("output/study3/", layout$name, "_axial.png")
  )
  
  x_range <- range(layout$x)
  y_range <- range(layout$y)
  
  # Prepare edge positions and labels
  edge_positions <- igraph::as_data_frame(g, what = "edges") %>%
    left_join(layout %>% select(name, x, y), by = c("from" = "name")) %>%
    rename(x_from = x, y_from = y) %>%
    left_join(layout %>% select(name, x, y), by = c("to" = "name")) %>%
    rename(x_to = x, y_to = y) %>%
    mutate(weight = as.numeric(weight)) %>%
    mutate(angle = atan2(y_to - y_from, x_to - x_from) * 180 / pi,  # convert radians to degrees
           # Optionally flip upside-down labels
           angle = ifelse(angle < -90 | angle > 90, angle + 180, angle)) %>%
    mutate(
      x_mid = (x_from + x_to)/2,
      y_mid = (y_from + y_to)/2,
      # For loops, nudge label
      nudge_x = ifelse(from == to, case_when(
        from == 'FPN' ~ 0,
        from == 'DMN' ~ 0,
        from == 'SN'  ~ 0
      ), -0.1 * sin(angle*pi/180)),
      nudge_y = ifelse(from == to, case_when(
        from == 'FPN' ~ -0.6,
        from == 'DMN' ~ -0.6,
        from == 'SN'  ~ 0.6
      ), ifelse(from == 'DMN' & to == 'FPN', -(0.1 * cos(angle*pi/180)), 0.1 * cos(angle*pi/180))),
      # Conditional bold for diff if significant
      diff_label = ifelse(sig_edge,
                          paste0("<b>", ifelse(round(diff,3)>0, paste0("+", round(diff,3)), round(diff,3)), "</b>"),
                          ifelse(round(diff,3)>0, paste0("+", round(diff,3)), round(diff,3))),
      label = ifelse(from == to,
                     paste0(
                       "<i>r</i> = ", round(weight,3),
                       "; \u0394<i>r</i> = ", diff_label), paste0(
                         "<i>r</i> = ", round(weight,3),
                         "; \u0394<i>r</i> = ", diff_label, ifelse(sig_edge,
                                                                   '*',''))
      )
    )
  
  # Extract numeric coordinates
  x_coords <- as.numeric(layout$x)
  y_coords <- as.numeric(layout$y)
  
  margin <- 0.5
  
  # Plot
  ggraph(graph = g, layout = layout) +
    geom_edge_fan(aes(width = diff, color = sig_edge),
                  arrow = arrow(length = unit(2, "mm"), type='closed', ends = 'both'),
                  start_cap = circle(14,'mm'),
                  end_cap = circle(14, 'mm'),
                  strength = 0) +
    geom_richtext(data = edge_positions,
                  aes(x = x_mid, y = y_mid, label = label, angle = angle),
                  nudge_x = edge_positions$nudge_x,
                  nudge_y = edge_positions$nudge_y,
                  fill = "white",
                  label.color = 'black',
                  label.size = 0, 
                  label.padding = unit(0.1, "lines"),
                  size = 4,
                  family = "Times New Roman",
                  label.r = unit(0, "lines"),
                  vjust = 0.5,
                  hjust = 0.5,
                  inherit.aes	= FALSE) +
    geom_image(
      data = network_positions,
      aes(x = x, y = y, image = img),
      size = 0.26,
      by = "width",
      asp = 1       # maintain aspect ratio
    ) +
    geom_edge_loop(aes(width = weight, color = sig_edge, direction = loop_angle),
                   arrow = arrow(length = unit(2, "mm"), type='closed', ends= 'both'),
                   start_cap = circle(14,'mm'),
                   end_cap = circle(14, 'mm')) +
    geom_node_label(aes(label = name),
                    fill = "white",
                    label.padding = unit(0.05, "lines"),
                    size = 5,
                    family = "Times New Roman",
                    label.size = 0) +
    
    scale_edge_width(name = 'Difference from Uninitiated', range = c(1,3)) +
    expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) +
    base_theme +
    theme(panel.border = element_blank(),
          text = element_text(family = "Times New Roman", color = 'black')) +
    coord_fixed(
      xlim = c(min(layout$x) - margin, max(layout$x) + margin),
      ylim = c(min(layout$y) - margin, max(layout$y) + margin),
      clip = 'off'
    ) + 
    guides(
      edge_width = 'none',
      edge_color = 'none'
    ) +
    scale_edge_color_manual(values = c("FALSE" = "black", "TRUE" = "red"))
}

# Generate plots
plot_init   <- plot_network_init(edges_init,   "Drug-Use Initiated Adolescent\nTriple Network Dynamics")

ggsave('Thesis Figures/Output/Study 3/study 3 triple network plot init_only.png', plot_init, 'png',
       path = NULL,
       scale = 1,
       width = 6,
       height = 6,
       units = 'in',
       dpi = 320,
       limitsize = TRUE,
       bg = 'transparent',
       create.dir = FALSE)



## STUDY FOUR ------------------------------------------------------------------

###### MPLUS OUTPUT FOR ELBOW PLOT #############################################
###### ELBOW PLOT ##############################################################
##Information Criteria
mplus_ic <- data.frame(
  ng = 1:10,
  AIC = NA,
  BIC = NA,
  aBIC = NA,
  logLike = NA,
  entropy = NA
)

# Combine all ng vectors into a list, info manually retrieved from mplus output
ng_list <- list(
  c(118816.961, 118886.129, 118854.351, -59398.481, NA),  # ng1
  c(84893.908, 85018.405, 84961.205, -42428.954, 0.880),  # ng2
  c(80820.450, 81014.112, 80925.134, -40382.225, 0.878),  # ng3
  c(79076.220, 79339.047, 79218.291, -39500.110, 0.827),  # ng4
  c(78092.239, 78424.231, 78271.697, -38998.119, 0.799),  # ng5
  c(77486.450, 77887.607, 77703.295, -38685.225, 0.809),  # ng6
  c(77136.650, 77606.972, 77390.883, -38500.325, 0.799),  # ng7
  c(76889.399, 77428.887, 77181.020, -38366.700, 0.805),  # ng8
  c(76792.801, 77401.453, 77121.808, -38308.400, 0.835),  # ng9
  c(76656.283, 77334.101, 77022.678, -38230.141, 0.687)   # ng10
)

# Loop through each class and fill mplus_ic
for (i in 1:10) {
  mplus_ic[mplus_ic$ng == i, c("AIC","BIC",'aBIC',"logLike", 'entropy')] <- ng_list[[i]]
}

mplus_poly_ic_long <- mplus_ic %>%
  pivot_longer(
    cols = c(AIC, BIC, aBIC, logLike, entropy),
    names_to = "Criterion",
    values_to = "Value"
  )

mplus_ic <- mplus_ic %>%
  arrange(ng) %>%
  mutate(
    deltaAIC = ifelse( ng == 1, NA, AIC - lag(AIC)),
    deltaBIC = ifelse( ng == 1, NA, BIC - lag(BIC)),
    deltaABIC = ifelse(ng == 1, NA, aBIC - lag(aBIC)),
    deltalogLike = ifelse( ng == 1, NA, logLike - lag(logLike)),
    deltaEntropy = ifelse(ng == 1, NA, entropy - lag(entropy))
  ) %>%
  select(ng, AIC, deltaAIC, BIC, deltaBIC, aBIC, deltaABIC, logLike, deltalogLike,
         entropy, deltaEntropy)



# Define the file pattern and empty list for results
ng_files <- sprintf("output/study4/mPlus/poly_wBL/poly_ng%d_mplus_wMRI.dat", 1:10)
assignment_counts <- list()

assignment_counts <- lapply(ng_files, function(file) {
  message("Reading: ", file)
  dat <- read.table(file, header = TRUE, fill = TRUE)
  lst <- ncol(dat)
  last_col <- dat[[lst]]
  table(last_col)
})

names(assignment_counts) <- paste0("ng", 1:10)
assignment_counts

mplus_ic$class_dist <- assignment_counts

combos <- lapply(1:10, function(num) {
  choose(num, 2)
})

mplus_ic$unique_tests <- combos

mplus_ic <- mplus_ic %>%
  select(ng, unique_tests, class_dist, everything())

mplus_ic$class_dist <- sapply(mplus_ic$class_dist, 
                              function(x) paste(x, collapse = "; "))

mplus_ic$unique_tests <- as.character(mplus_ic$unique_tests)

# Faceted elbow plot
mplus_poly_ic_long <- mplus_poly_ic_long %>%
  mutate(Criterion = factor(Criterion,
                            levels = c('aBIC', 'BIC', 'AIC', 'logLike',
                                       'entropy'),
                            labels = c('Sample Adjusted BIC',
                                       'BIC',
                                       'AIC',
                                       'Log-likelihood',
                                       'Entropy')))

aic_plot <- ggplot(mplus_poly_ic_long %>%
                     filter(Criterion == 'AIC' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'AIC'
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )

bic_plot <- ggplot(mplus_poly_ic_long %>%
                     filter(Criterion == 'BIC' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'BIC'
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )

abic_plot <- ggplot(mplus_poly_ic_long %>%
                      filter(Criterion == 'Sample Adjusted BIC' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'Sample Adjusted BIC'
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )

loglike_plot <- ggplot(mplus_poly_ic_long %>%
                         filter(Criterion == 'Log-likelihood' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'Log-likelihood'
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )

ent_plot <- ggplot(mplus_poly_ic_long %>%
                     filter(Criterion == 'Entropy' & ng != 1), aes(x = factor(ng), y = Value)) +
  geom_line(aes(group=1),color = "steelblue", linewidth = 1) +
  geom_point(size = 2, color = "darkred") +
  base_theme +
  labs(
    title = 'Entropy'
  ) +
  scale_y_continuous(limits = c(0.65,0.9)) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = -40)
  )



final_plot <- (aic_plot | loglike_plot) / (bic_plot | abic_plot) / ent_plot



ggsave('Thesis Figures/Output/Study 4/elbow_plot.png', final_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')

###### Plot best model #########################################################
data.long <- readRDS('output/study4/data_long.rds')
id_lookup<-read.csv('output/study4/mPlus/id_lookup.csv', header=TRUE)

two.model <- read.table("output/study4/mPlus/poly_wBL/poly_ng2_mplus_wMRI.dat", header = TRUE)
comb_ng2_assignments <- data.frame(
  mplus_id = as.numeric(two.model$X1),
  ng2_class = factor(two.model[[ncol(two.model)]])
)

three.model <- read.table("output/study4/mPlus/poly_wBL/poly_ng3_mplus_wMRI.dat", header = TRUE)
comb_ng3_assignments <- data.frame(
  mplus_id = as.numeric(three.model$X1),
  ng3_class = factor(three.model[[ncol(three.model)]])
)

four.model <- read.table("output/study4/mPlus/poly_wBL/poly_ng4_mplus_wMRI.dat", header = TRUE)
comb_ng4_assignments <- data.frame(
  mplus_id = as.numeric(four.model$X1),
  ng4_class = factor(four.model[[ncol(four.model)]])
)

five.model <- read.table("output/study4/mPlus/poly_wBL/poly_ng5_mplus_wMRI.dat", header = TRUE)
comb_ng5_assignments <- data.frame(
  mplus_id = as.numeric(five.model$X1),
  ng5_class = factor(five.model[[ncol(five.model)]])
)

six.model <- read.table("output/study4/mPlus/poly_wBL/poly_ng6_mplus_wMRI.dat", header = TRUE)
comb_ng6_assignments <- data.frame(
  mplus_id = as.numeric(six.model$X1),
  ng6_class = factor(six.model[[ncol(six.model)]])
)

seven.model <- read.table("output/study4/mPlus/poly_wBL/poly_ng7_mplus_wMRI.dat", header = TRUE,
                          fill = TRUE)
comb_ng7_assignments <- data.frame(
  mplus_id = as.numeric(seven.model$X1),
  ng7_class = factor(seven.model[[ncol(seven.model)]])
)

eight.model <- read.table("output/study4/mPlus/poly_wBL/poly_ng8_mplus_wMRI.dat", header = TRUE)
comb_ng8_assignments <- data.frame(
  mplus_id = as.numeric(eight.model$X1),
  ng8_class = factor(eight.model[[ncol(eight.model)]])
)

nine.model <- read.table("output/study4/mPlus/poly_wBL/poly_ng9_mplus_wMRI.dat", header = TRUE,
                         fill = TRUE)
comb_ng9_assignments <- data.frame(
  mplus_id = as.numeric(nine.model$X1),
  ng9_class = factor(nine.model[[ncol(nine.model)]])
)

ten.model <- read.table("output/study4/mPlus/poly_wBL/poly_ng10_mplus_wMRI.dat", header = TRUE)
comb_ng10_assignments <- data.frame(
  mplus_id = as.numeric(ten.model$X1),
  ng10_class = factor(ten.model[[ncol(ten.model)]])
)


id_lookup <- id_lookup %>%
  left_join(comb_ng2_assignments, by = 'mplus_id')%>%
  left_join(comb_ng3_assignments, by = 'mplus_id') %>%
  left_join(comb_ng4_assignments, by = 'mplus_id')%>%
  left_join(comb_ng5_assignments, by = 'mplus_id')%>%
  left_join(comb_ng6_assignments, by = 'mplus_id')%>%
  left_join(comb_ng7_assignments, by = 'mplus_id')%>%
  left_join(comb_ng8_assignments, by = 'mplus_id')%>%
  left_join(comb_ng9_assignments, by = 'mplus_id')%>%
  left_join(comb_ng10_assignments, by = 'mplus_id')

data.long <- data.long %>%
  left_join(id_lookup, by = 'participant_id') %>%
  select(-mplus_id) %>%
  filter(participant_id != 'sub-003RTV85')

data.long <- data.long %>%
  mutate(
    ng7_class = ifelse(is.na(ng7_class), 7, ng7_class),
    ng9_class = ifelse(is.na(ng9_class), 9, ng9_class)
  )

rm(comb_ng3_assignments, three.model, four.model, comb_ng4_assignments, id_lookup)

# Compute means per drug per class per time
draw_key_rect <- function(data, params, size) {
  grid::rectGrob(
    width = unit(0.5, "cm"),
    height = unit(0.3, "cm"),  # smaller height
    gp = grid::gpar(
      col = NA,
      fill = alpha(data$colour, data$alpha %||% 1)
    )
  )
}




plot_ng2_data <- data.long %>%
  filter(time_mplus != 0 ) %>%
  pivot_longer(cols = c(alc_use, nic_use, mj_use, oth_use),
               names_to = "drug",
               values_to = "use") %>%
  group_by(ng2_class, drug, time_mplus) %>%
  summarise(mean_use = mean(use, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    drug = case_when(
      drug == 'alc_use' ~ 'Alcohol',
      drug == 'mj_use' ~ 'Marijuana',
      drug == 'nic_use' ~ 'Nicotine',
      drug == 'oth_use' ~ 'Other'
    ),
    time_mplus = time_mplus * 6,
    ng2_class = factor(ng2_class, levels = c(1,2))
  )

poly_ng2_plot <- ggplot(plot_ng2_data[!is.na(plot_ng2_data$ng2_class),], aes(x = factor(time_mplus), y = mean_use, 
                                                                             color = factor(ng2_class), group = ng2_class)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) +
  facet_wrap(~drug, ncol = 1, scales = "fixed") +  # one row per drug
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Time from Baseline (in Months)",
       y = "Probability to Use",
       color = "Profile") +
  base_theme

poly_ng2_plot

ggsave('Thesis Figures/Output/Study 4/ng2_traj_plot.png', poly_ng2_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')




plot_ng3_data <- data.long %>%
  filter(time_mplus != 0 ) %>%
  pivot_longer(cols = c(alc_use, nic_use, mj_use, oth_use),
               names_to = "drug",
               values_to = "use") %>%
  group_by(ng3_class, drug, time_mplus) %>%
  summarise(mean_use = mean(use, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    drug = case_when(
      drug == 'alc_use' ~ 'Alcohol',
      drug == 'mj_use' ~ 'Marijuana',
      drug == 'nic_use' ~ 'Nicotine',
      drug == 'oth_use' ~ 'Other'
    ),
    time_mplus = time_mplus * 6,
    ng3_class = factor(ng3_class, levels = c(1,2,3))
  )

poly_ng3_plot <- ggplot(plot_ng3_data[!is.na(plot_ng3_data$ng3_class),], aes(x = factor(time_mplus), y = mean_use, 
                                                                             color = factor(ng3_class), group = ng3_class)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) +
  facet_wrap(~drug, ncol = 1, scales = "fixed") +  # one row per drug
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Time from Baseline (in Months)",
       y = "Probability to Use",
       color = "Profile") +
  base_theme

poly_ng3_plot

ggsave('Thesis Figures/Output/Study 4/ng3_traj_plot.png', poly_ng3_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')


plot_ng4_data <- data.long %>%
  filter(time_mplus != 0 ) %>%
  pivot_longer(cols = c(alc_use, nic_use, mj_use, oth_use),
               names_to = "drug",
               values_to = "use") %>%
  group_by(ng4_class, drug, time_mplus) %>%
  summarise(mean_use = mean(use, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    drug = case_when(
      drug == 'alc_use' ~ 'Alcohol',
      drug == 'mj_use' ~ 'Marijuana',
      drug == 'nic_use' ~ 'Nicotine',
      drug == 'oth_use' ~ 'Other'
    ),
    time_mplus = time_mplus * 6,
    ng4_class = factor(ng4_class, levels = c(1,2,3,4), labels = c('No/Low Use',
                                                                  'Alcohol Low-increasing-to-Moderate',
                                                                  'Alcohol Moderate-increasing-to-High',
                                                                  'Polysubstance Low-increasing-to-High'))
  )

poly_ng4_plot <- ggplot(plot_ng4_data[!is.na(plot_ng4_data$ng4_class),], aes(x = factor(time_mplus), y = mean_use, 
                                                                             color = factor(ng4_class), group = ng4_class)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) +
  facet_wrap(~drug, ncol = 1, scales = "fixed") +  # one row per drug
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Time from Baseline (in Months)",
       y = "Probability of Use",
       color = "Profile") +
  base_theme +
  scale_color_npg()

poly_ng4_plot

ggsave('Thesis Figures/Output/Study 4/ng4_traj_plot.png', poly_ng4_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')


plot_ng5_data <- data.long %>%
  filter(time_mplus != 0 ) %>%
  pivot_longer(cols = c(alc_use, nic_use, mj_use, oth_use),
               names_to = "drug",
               values_to = "use") %>%
  group_by(ng5_class, drug, time_mplus) %>%
  summarise(mean_use = mean(use, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    drug = case_when(
      drug == 'alc_use' ~ 'Alcohol',
      drug == 'mj_use' ~ 'Marijuana',
      drug == 'nic_use' ~ 'Nicotine',
      drug == 'oth_use' ~ 'Other'
    ),
    time_mplus = time_mplus * 6,
    ng5_class = factor(ng5_class, levels = c(1,2,3,4,5))
  )

poly_ng5_plot <- ggplot(plot_ng5_data[!is.na(plot_ng5_data$ng5_class),], aes(x = factor(time_mplus), y = mean_use, 
                                                                             color = factor(ng5_class), group = ng5_class)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) +
  facet_wrap(~drug, ncol = 1, scales = "fixed") +  # one row per drug
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Time from Baseline (in Months)",
       y = "Probability to Use",
       color = "Profile") +
  base_theme

poly_ng5_plot

ggsave('Thesis Figures/Output/Study 4/ng5_traj_plot.png', poly_ng5_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')



plot_ng6_data <- data.long %>%
  filter(time_mplus != 0 ) %>%
  pivot_longer(cols = c(alc_use, nic_use, mj_use, oth_use),
               names_to = "drug",
               values_to = "use") %>%
  group_by(ng6_class, drug, time_mplus) %>%
  summarise(mean_use = mean(use, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    drug = case_when(
      drug == 'alc_use' ~ 'Alcohol',
      drug == 'mj_use' ~ 'Marijuana',
      drug == 'nic_use' ~ 'Nicotine',
      drug == 'oth_use' ~ 'Other'
    ),
    time_mplus = time_mplus * 6,
    ng6_class = factor(ng6_class, levels = c(1,2,3,4,5,6))
  )

poly_ng6_plot <- ggplot(plot_ng6_data[!is.na(plot_ng6_data$ng6_class),], aes(x = factor(time_mplus), y = mean_use, 
                                                                             color = factor(ng6_class), group = ng6_class)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) +
  facet_wrap(~drug, ncol = 1, scales = "fixed") +  # one row per drug
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Time from Baseline (in Months)",
       y = "Probability to Use",
       color = "Profile") +
  base_theme

poly_ng6_plot

ggsave('Thesis Figures/Output/Study 4/ng6_traj_plot.png', poly_ng6_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')



plot_ng7_data <- data.long %>%
  filter(time_mplus != 0 ) %>%
  pivot_longer(cols = c(alc_use, nic_use, mj_use, oth_use),
               names_to = "drug",
               values_to = "use") %>%
  group_by(ng7_class, drug, time_mplus) %>%
  summarise(mean_use = mean(use, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    drug = case_when(
      drug == 'alc_use' ~ 'Alcohol',
      drug == 'mj_use' ~ 'Marijuana',
      drug == 'nic_use' ~ 'Nicotine',
      drug == 'oth_use' ~ 'Other'
    ),
    time_mplus = time_mplus * 6,
    ng7_class = factor(ng7_class, levels = c(1,2,3,4,5,6,7),
                       labels = c(1,2,3,4,5,6,'7*'))
  )

poly_ng7_plot <- ggplot(plot_ng7_data[!is.na(plot_ng7_data$ng7_class),], aes(x = factor(time_mplus), y = mean_use, 
                                                                             color = factor(ng7_class), group = ng7_class)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) +
  facet_wrap(~drug, ncol = 1, scales = "fixed") +  # one row per drug
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Time from Baseline (in Months)",
       y = "Probability to Use",
       color = "Profile") +
  base_theme

poly_ng7_plot

ggsave('Thesis Figures/Output/Study 4/ng7_traj_plot.png', poly_ng7_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')


plot_ng8_data <- data.long %>%
  filter(time_mplus != 0 ) %>%
  pivot_longer(cols = c(alc_use, nic_use, mj_use, oth_use),
               names_to = "drug",
               values_to = "use") %>%
  group_by(ng8_class, drug, time_mplus) %>%
  summarise(mean_use = mean(use, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    drug = case_when(
      drug == 'alc_use' ~ 'Alcohol',
      drug == 'mj_use' ~ 'Marijuana',
      drug == 'nic_use' ~ 'Nicotine',
      drug == 'oth_use' ~ 'Other'
    ),
    time_mplus = time_mplus * 6,
    ng8_class = factor(ng8_class, levels = c(1,2,3,4,5,6,7,8))
  )

poly_ng8_plot <- ggplot(plot_ng8_data[!is.na(plot_ng8_data$ng8_class),], aes(x = factor(time_mplus), y = mean_use, 
                                                                             color = factor(ng8_class), group = ng8_class)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) +
  facet_wrap(~drug, ncol = 1, scales = "fixed") +  # one row per drug
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Time from Baseline (in Months)",
       y = "Probability to Use",
       color = "Profile") +
  base_theme

poly_ng8_plot

ggsave('Thesis Figures/Output/Study 4/ng8_traj_plot.png', poly_ng8_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')


plot_ng9_data <- data.long %>%
  filter(time_mplus != 0 ) %>%
  pivot_longer(cols = c(alc_use, nic_use, mj_use, oth_use),
               names_to = "drug",
               values_to = "use") %>%
  group_by(ng9_class, drug, time_mplus) %>%
  summarise(mean_use = mean(use, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    drug = case_when(
      drug == 'alc_use' ~ 'Alcohol',
      drug == 'mj_use' ~ 'Marijuana',
      drug == 'nic_use' ~ 'Nicotine',
      drug == 'oth_use' ~ 'Other'
    ),
    time_mplus = time_mplus * 6,
    ng9_class = factor(ng9_class, levels = c(1,2,3,4,5,6,7,8,9),
                       labels = c(1,2,3,4,5,6,7,8,'9*'))
  )

poly_ng9_plot <- ggplot(plot_ng9_data[!is.na(plot_ng9_data$ng9_class),], aes(x = factor(time_mplus), y = mean_use, 
                                                                             color = factor(ng9_class), group = ng9_class)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) +
  facet_wrap(~drug, ncol = 1, scales = "fixed") +  # one row per drug
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Time from Baseline (in Months)",
       y = "Probability to Use",
       color = "Profile") +
  base_theme

poly_ng9_plot

ggsave('Thesis Figures/Output/Study 4/ng9_traj_plot.png', poly_ng9_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')



plot_ng10_data <- data.long %>%
  filter(time_mplus != 0 ) %>%
  pivot_longer(cols = c(alc_use, nic_use, mj_use, oth_use),
               names_to = "drug",
               values_to = "use") %>%
  group_by(ng10_class, drug, time_mplus) %>%
  summarise(mean_use = mean(use, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    drug = case_when(
      drug == 'alc_use' ~ 'Alcohol',
      drug == 'mj_use' ~ 'Marijuana',
      drug == 'nic_use' ~ 'Nicotine',
      drug == 'oth_use' ~ 'Other'
    ),
    time_mplus = time_mplus * 6,
    ng10_class = factor(ng10_class, levels = c(1,2,3,4,5,6,7,8,9,10))
  )

poly_ng10_plot <- ggplot(plot_ng10_data[!is.na(plot_ng10_data$ng10_class),], aes(x = factor(time_mplus), y = mean_use, 
                                                                                 color = factor(ng10_class), group = ng10_class)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) +
  facet_wrap(~drug, ncol = 1, scales = "fixed") +  # one row per drug
  scale_y_continuous(limits = c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Time from Baseline (in Months)",
       y = "Probability to Use",
       color = "Profile") +
  base_theme

poly_ng10_plot

ggsave('Thesis Figures/Output/Study 4/ng10_traj_plot.png', poly_ng10_plot, 'png',
       path = NULL,
       scale = 1,
       width = 10,
       height = 8,
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')


##AVERAGE POSTERIOR PROBABILITIES

get_app <- function(df, model_name) {
  
  # Assigned class column (last column)
  class_col <- ncol(df)
  
  # Assigned classes
  assigned_class <- df[[class_col]]
  
  # Number of classes
  k <- length(unique(assigned_class))
  
  # Probability columns = k columns immediately before class column
  prob_cols <- (class_col - k):(class_col - 1)
  
  # Probability matrix
  probs <- as.matrix(df[, prob_cols])
  
  # Extract probability corresponding to assigned class
  assigned_prob <- probs[cbind(seq_len(nrow(probs)), assigned_class)]
  
  # Build output dataframe
  out <- data.frame(
    class = assigned_class,
    assigned_prob = assigned_prob
  ) %>%
    group_by(class) %>%
    summarise(
      APP = mean(assigned_prob, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(model = model_name)
  
  return(out)
}

models <- list(
  two = two.model,
  three = three.model,
  four = four.model,
  five = five.model,
  six = six.model,
  seven = seven.model,
  eight = eight.model,
  nine = nine.model,
  ten = ten.model
)

results <- imap_dfr(models, get_app)

write.csv(results, 'output/study4/mPlus/average_probabilities.csv')


###### Summary Statistics ######################################################
cat("\n\n=== ELASTIC NET SUMMARY ===\n")

load("output/study4/en_results/claude_elastic_net_results_wBL.RData")
data.3 <- readRDS('output/study4/data_management/wBL_fulldata.rds')
table_names_all <- readRDS('output/study4/data_management/table_names_wBL.rds')
table_names_all$table_name <- gsub("CRPF", "Access to substances (CRPF)", table_names_all$table_name)
table_names_all$domain_name <- ifelse(table_names_all$table_name == "Access to substances (CRPF)", "Culture & Environment",
                                      table_names_all$domain_name)
table_names_all$domain_name <- ifelse(table_names_all$table_name == "Family History of Substance Use", "Parenting Behaviors",
                                      table_names_all$domain_name)
table_names_all$domain_name <- gsub("Culture & Environment", "Community & School", table_names_all$domain_name)
table_names_all$table_name <- gsub("Family History of Substance Use", "Family history of substance use", table_names_all$table_name)


pairs <- combn(1:4, 2, simplify = FALSE)
# Extract AUCs for each comparison
summary_df <- data.frame()

for (comp_idx in seq_along(pairs)) {
  
  p <- pairs[[comp_idx]]
  level_order <- p[order(-table(data.3[data.3$ng4_class %in% p,]$ng4_class)[p])]
  level1 <- level_order[1]
  level2 <- level_order[2]
  
  aucs <- numeric(10)
  
  for (split_idx in 1:10) {
    result_name <- sprintf("comp%d_split%d", comp_idx, split_idx)
    if (result_name %in% names(enet_all_results)) {
      aucs[split_idx] <- enet_all_results[[result_name]]$metrics$auc
    }
  }
  
  summary_df <- rbind(summary_df, data.frame(
    comparison = sprintf("%d vs %d", level1, level2),
    mean_auc = mean(aucs),
    sd_auc = sd(aucs),
    min_auc = min(aucs),
    max_auc = max(aucs)
  ))
}

print(summary_df)


# Define pairwise comparisons
pairs <- combn(1:4, 2, simplify = FALSE)

# Initialize storage
all_results <- list()

cat("=== Extracting Elastic Net Results from Environment ===\n\n")

# Loop through all comparisons and splits
for (comp_idx in seq_along(pairs)) {
  
  p <- pairs[[comp_idx]]
  level_order <- p[order(-table(data.3[data.3$ng4_class %in% p,]$ng4_class)[p])]
  level1 <- level_order[1]
  level2 <- level_order[2]
  
  cat(sprintf("\nComparison %d: %d vs %d\n", comp_idx, level1, level2))
  
  # Storage for this comparison's splits
  comp_results <- list()
  
  for (split_idx in 1:10) {
    
    # Get object name
    enet_name <- sprintf("enet_%d_%dvs%d", split_idx, level1, level2)
    
    # Check if object exists
    if (exists(enet_name, envir = .GlobalEnv)) {
      
      enet_obj <- get(enet_name, envir = .GlobalEnv)
      
      # Extract key information
      comp_results[[paste0("split_", split_idx)]] <- list(
        split = split_idx,
        comparison = sprintf("%d vs %d", level1, level2),
        levels = c(level1, level2),
        auc = enet_obj$metrics$auc,
        accuracy = enet_obj$metrics$accuracy,
        sensitivity = enet_obj$metrics$sensitivity,
        specificity = enet_obj$metrics$specificity,
        ppv = enet_obj$metrics$ppv,
        npv = enet_obj$metrics$npv,
        threshold = enet_obj$metrics$threshold,
        lambda_min = enet_obj$enet_results$lambda_final,
        n_predictors = nrow(enet_obj$pooled_coef$results_df),
        coefficients = enet_obj$pooled_coef$results_df,
        roc_object = enet_obj$metrics$roc_obj,
        confusion_matrix = enet_obj$metrics$confusion_matrix
      )
      
      cat(sprintf("  Split %d: AUC = %.3f, Accuracy = %.3f, %d predictors\n",
                  split_idx, 
                  enet_obj$metrics$auc,
                  enet_obj$metrics$accuracy,
                  nrow(enet_obj$pooled_coef$results_df)))
      
    } else {
      cat(sprintf("  Split %d: NOT FOUND\n", split_idx))
    }
  }
  
  # Store results for this comparison
  if (length(comp_results) > 0) {
    all_results[[sprintf("comp_%d_%dvs%d", comp_idx, level1, level2)]] <- comp_results
  }
}

# Check if we have any results
if (length(all_results) == 0) {
  stop("No elastic net results found in environment!")
}

cat("\n\n=== Creating Summary Tables ===\n")

######  PERFORMANCE SUMMARY BY SPLIT ###########################################

performance_by_split <- list()

for (comp_name in names(all_results)) {
  comp_results <- all_results[[comp_name]]
  
  for (split_name in names(comp_results)) {
    split_result <- comp_results[[split_name]]
    
    performance_by_split[[paste(comp_name, split_name, sep = "_")]] <- data.frame(
      comparison = split_result$comparison,
      split = split_result$split,
      auc = split_result$auc,
      accuracy = split_result$accuracy,
      sensitivity = split_result$sensitivity,
      specificity = split_result$specificity,
      ppv = split_result$ppv,
      npv = split_result$npv,
      lambda_min = split_result$lambda_min,
      n_predictors = split_result$n_predictors,
      stringsAsFactors = FALSE
    ) %>%
      mutate(comparison = case_when(
        comparison == '1 vs 2' ~ 'No/Low Use (Reference) v L-M Alcohol',
        comparison == '1 vs 3' ~ 'No/Low Use (Reference) v M-H Alcohol',
        comparison == '2 vs 3' ~ 'L-M Alcohol (Reference) v M-H Alcohol',
        comparison == '1 vs 4' ~ 'No/Low Use (Reference) v L-H Poly',
        comparison == '2 vs 4' ~ 'L-M Alcohol (Reference) v L-H Poly',
        comparison == '3 vs 4' ~ 'M-H Alcohol (Reference) v L-H Poly'
      ))
  }
}

performance_df <- bind_rows(performance_by_split) %>%
  arrange(comparison, split)

cat("\n--- Performance by Split ---\n")
print(head(performance_df, 20))


######  AGGREGATED PERFORMANCE SUMMARY (MEAN ACROSS SPLITS) ####################

performance_summary <- performance_df %>%
  group_by(comparison) %>%
  summarise(
    n_splits = n(),
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE),
    min_auc = min(auc, na.rm = TRUE),
    max_auc = max(auc, na.rm = TRUE),
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    sd_accuracy = sd(accuracy, na.rm = TRUE),
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    mean_n_predictors = mean(n_predictors, na.rm = TRUE),
    sd_n_predictors = sd(n_predictors, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_auc))

cat("\n--- Aggregated Performance Summary ---\n")
print(performance_summary)

######  TOP PREDICTORS BY COMPARISON (ACROSS ALL SPLITS) #######################

cat("\n--- Extracting Coefficients ---\n")

all_coefficients <- list()

for (comp_name in names(all_results)) {
  comp_results <- all_results[[comp_name]]
  
  for (split_name in names(comp_results)) {
    split_result <- comp_results[[split_name]]
    
    if (nrow(split_result$coefficients) > 0) {
      coefs <- split_result$coefficients %>%
        mutate(
          comparison = split_result$comparison,
          split = split_result$split
        ) %>%
        select(comparison, split, variable, coefficient, stability_pct)
      
      all_coefficients[[paste(comp_name, split_name, sep = "_")]] <- coefs
    }
  }
}

all_coefs_df <- bind_rows(all_coefficients)

cat(sprintf("Total coefficient entries: %d\n", nrow(all_coefs_df)))

# Calculate variable importance: how often selected and average coefficient
variable_importance <- all_coefs_df %>%
  mutate(comparison = case_when(
    comparison == '1 vs 2' ~ 'No/Low Use (Reference) v Alcohol L-M',
    comparison == '1 vs 3' ~ 'No/Low Use (Reference) v Alcohol M-H',
    comparison == '2 vs 3' ~ 'Alcohol L-M (Reference) v Alcohol M-H',
    comparison == '1 vs 4' ~ 'No/Low Use (Reference) v Poly L-H',
    comparison == '2 vs 4' ~ 'Alcohol L-M (Reference) v Poly L-H',
    comparison == '3 vs 4' ~ 'Alcohol M-H (Reference) v Poly L-H'
  )) %>%
  group_by(comparison, variable) %>%
  summarise(
    n_splits_selected = n(),
    mean_coefficient = mean(coefficient, na.rm = TRUE),
    sd_coefficient = sd(coefficient, na.rm = TRUE),
    mean_stability = mean(stability_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(abs_mean_coef = abs(mean_coefficient)) %>%
  arrange(comparison, desc(n_splits_selected), desc(abs_mean_coef)) %>%
  left_join(table_names_all[,c('variable', 'table_name', 'domain_name')], by = 'variable') %>%
  select(comparison, domain_name, table_name, variable, everything()) %>%
  arrange(comparison, rev(abs_mean_coef))

supplement_predictors <- variable_importance %>%
  filter(n_splits_selected >= 5) %>%
  filter(variable != "(Intercept)") %>%
  group_by(comparison) %>%
  arrange(desc(abs_mean_coef), .by_group = TRUE) %>%
  mutate(rank = dense_rank(desc(abs_mean_coef))) %>%
  ungroup()

# Top predictors (selected in >=7 splits with high coefficients)
top_predictors <- variable_importance %>%
  filter(n_splits_selected >= 5) %>%
  filter(variable != '(Intercept)') %>%
  group_by(comparison) %>%
  slice_max(order_by = abs_mean_coef, n = 20) %>%
  ungroup() %>%
  arrange(comparison, rev(abs_mean_coef))

cat("\n--- Top Predictors (selected in >=7 splits) ---\n")
print(head(top_predictors, 30))
#write.csv(supplement_predictors, 'ABCD6_Study1/en_results_claude/supplementary_table.csv')


######  HORIZONTAL BOXPLOTS ####################################################




domains <- unique(top_predictors$domain_name)
chroma <- data.frame(
  domain = domains,
  colors = ggsci::pal_npg()(8))%>%
  mutate(
    domain = ifelse(domain == 'Self and Peer Involvement with Substance Use', 'Involvement with\nSubstance Use', domain),
    colors = 
      ifelse(domain == 'Parenting Behaviors', '#7E614899', colors)
  )

domain_palette <- setNames(chroma$colors, chroma$domain)

#Shrink label titles for figure
top_predictors$domain_name <- ifelse(top_predictors$domain_name == 'Self and Peer Involvement with Substance Use', 'Involvement with\nSubstance Use', top_predictors$domain_name)
top_predictors$table_name <- gsub(" of Desikan ROI:", ":", top_predictors$table_name)
top_predictors$table_name <- gsub(" in Desikan ROI:", ":", top_predictors$table_name)
top_predictors$table_name <- gsub("Average correlation between Gordon networks", "rs-FC", top_predictors$table_name)
top_predictors$table_name <- gsub("Average correlation between Gordon network", "rs-FC", top_predictors$table_name)
top_predictors$table_name <- gsub(" Subcortical ROI:", "", top_predictors$table_name)
top_predictors$table_name <- gsub(" \\(full shell DTI\\) in AtlasTrack fiber tract", "", top_predictors$table_name)
top_predictors$table_name <- gsub(" Test", "", top_predictors$table_name)
top_predictors$table_name <- gsub("Total v", "V", top_predictors$table_name)
top_predictors$table_name <- gsub("Total s", "S", top_predictors$table_name)
top_predictors$table_name <- gsub("Average s", "S", top_predictors$table_name)
top_predictors$table_name <- gsub("Average c", "C", top_predictors$table_name)
top_predictors$table_name <- gsub("Weighted average f", "F", top_predictors$table_name)
top_predictors$table_name <- gsub(" hemisphere)", ")", top_predictors$table_name)
top_predictors$table_name <- gsub("(right)", "R", top_predictors$table_name)
top_predictors$table_name <- gsub("(left)", "L", top_predictors$table_name)
top_predictors$table_name <- gsub("Primary Caregiver's Maximum Education - Some College or Associate's Degree",
                                  "Caregiver education - some college", top_predictors$table_name)
top_predictors$table_name <- gsub("NIH Toolbox - Dimensional Change Card Sort", "Card sort (NIH toolbox)", top_predictors$table_name)
top_predictors$table_name <- gsub("NIH Toolbox - Flanker Inhibitory Control and Attention", "Flanker task (NIH toolbox)", top_predictors$table_name)
top_predictors$table_name <- gsub("sensorimotor", "SM", top_predictors$table_name)
top_predictors$table_name <- gsub("SDS - Disorders of Initiating and Maintaining Sleep", "Disorders of maintaining sleep (SDS)", top_predictors$table_name)
top_predictors$table_name <- gsub("Baseline Use - Alcohol", "Baseline alcohol use", top_predictors$table_name)
top_predictors$table_name <- gsub("WISC-V Matrix Reasoning", "Matrix reasoning (WISC-V)", top_predictors$table_name)
top_predictors$table_name <- gsub("NIH Toolbox - Picture Vocabulary", "Picture vocabulary (NIH toolbox)", top_predictors$table_name)
top_predictors$table_name <- gsub("CBCL - Rule-breaking Behavior", "Rule-breaking behavior (CBCL)", top_predictors$table_name)
top_predictors$table_name <- gsub("CBCL - Externalizing Disorders", "Externalizing disorders (CBCL)", top_predictors$table_name)
top_predictors$table_name <- gsub("PDS - Puberty Development", "Puberty development (PDS)", top_predictors$table_name)
top_predictors$table_name <- gsub("Religious Affiliation -", "Religious affiliation -", top_predictors$table_name)
top_predictors$table_name <- gsub("Intent to Use - Alcohol", "Intent to use alcohol", top_predictors$table_name)
top_predictors$table_name <- gsub("UPPS-P - Sensation Seeking", "Sensation-seeking (UPPS-P)", top_predictors$table_name)
top_predictors$table_name <- gsub("RAVLT - Short Delay", "Short delay (RAVLT)", top_predictors$table_name)
top_predictors$table_name <- gsub("UPPS-P - Lack of Planning", "Lack of planning (UPPS-P)", top_predictors$table_name)
top_predictors$table_name <- gsub("Prenatal Exposure - Caffeine", "Prenatal exposure to caffeine", top_predictors$table_name)
top_predictors$table_name <- gsub("Prenatal Exposure - Substance Use - Yes", "Prenatal exposure to substances", top_predictors$table_name)
top_predictors$table_name <- gsub("Intent to Use - Nicotine", "Intent to use nicotine", top_predictors$table_name)
top_predictors$table_name <- gsub("Intent to Use - Cannabis", "Intent to use cannabis", top_predictors$table_name)
top_predictors$table_name <- gsub("PMQ - Parental Monitoring", "Parental monitoring (PMQ)", top_predictors$table_name)
top_predictors$table_name <- gsub("SRPF - School Involvement Subscale", "School involvement (SRPF)", top_predictors$table_name)
top_predictors$table_name <- gsub("SDS - Disorders of Excessive Somnolence", "Disorders of excessive somnolence (SDS)", top_predictors$table_name)
top_predictors$table_name <- gsub("Peer Substance Use - Nicotine", "Peer nicotine use", top_predictors$table_name)
top_predictors$table_name <- gsub("NIH Toolbox - Oral Reading Recognition", "Oral reading recognition (NIH toolbox)", top_predictors$table_name)
top_predictors$table_name <- gsub("Parent Rules", "Parent rules", top_predictors$table_name)
top_predictors$table_name <- gsub("Recreational Activities (continuous)", "Recreational activities (continuous)", top_predictors$table_name)
top_predictors$table_name <- gsub("none", "unassigned", top_predictors$table_name)
top_predictors$table_name <- gsub("Physical Activity - Days of 60-min.", "Physical activity - days of 60-min.", top_predictors$table_name)
top_predictors$table_name <- gsub("WISC-V - Matrix Reasoning", "Matrix reasoning (WISC-V)", top_predictors$table_name)
top_predictors$table_name <- gsub("CBCL - Withdraw / Depression", "Withdraw/depression (CBCL)", top_predictors$table_name)
top_predictors$table_name <- gsub("inferior parietal", "inferior parietal region", top_predictors$table_name)
top_predictors$table_name <- gsub("inferior temporal", "inferior temporal region", top_predictors$table_name)


for(i in 1:6) {
  
  comps <- unique(top_predictors$comparison)
  comp <- comps[i]
  box_data <- top_predictors %>%
    filter(comparison == comp) %>%
    mutate(table_name = if_else(
      nchar(table_name) > 49,
      sub("(.{49,}?)\\s", "\\1\n", table_name),
      table_name
    ),
    domain_name = ifelse(domain_name == 'Self and Peer Involvement with Substance Use', 'Self and Peer Involvement\nwith Substance Use', 
                         domain_name),
    mean_coefficient = ifelse(mean_coefficient > 2, 2, mean_coefficient)
    )
  if(i == 3 | i == 1 | i == 5){
    temp <- ggplot(box_data, 
                   aes(fill = domain_name, y = mean_coefficient, 
                       x = reorder(table_name, abs_mean_coef))) + 
      geom_bar(stat = 'identity') +
      geom_text(aes(label = round(mean_coefficient, 2)),
                hjust = ifelse(box_data$mean_coefficient > 0, -0.2, 1.2),  # text offset depends on sign
                size = 3,
                family = 'Times') +
      coord_flip() +
      labs(title = comp, y = 'Coefficient') +
      scale_y_continuous(limits = c(-1,2.1)) +
      scale_fill_manual(values = domain_palette) +
      scale_x_discrete(position = "top") +
      base_theme +
      theme(
        legend.position = 'none',
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0, margin = margin(t= 5, b = 5)),
        axis.text.y = element_text(size = 10,
                                   margin = margin(r = 10)),
        plot.margin = margin(0, 0, 0, 0),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 10), hjust = 0.4))
  } else{
    temp <- ggplot(box_data, 
                   aes(fill = domain_name, y = mean_coefficient, 
                       x = reorder(table_name, abs_mean_coef))) + 
      geom_bar(stat = 'identity') +
      geom_text(aes(label = round(mean_coefficient, 2)),
                hjust = ifelse(box_data$mean_coefficient > 0, -0.2, 1.2),  # text offset depends on sign
                size = 3,
                family = 'Times') +
      coord_flip() +
      labs(title = comp, y = 'Coefficient', 
           x = 'Predictors') +
      scale_y_continuous(limits = c(-1,2.1)) +
      scale_fill_manual(values = domain_palette) +
      base_theme +
      theme(
        legend.position = 'none',
        legend.title = element_blank(),
        axis.text.y = element_text(size = 10,
                                   margin = margin(r = 10)),
        plot.margin = margin(0, 0, 0, 0),
        plot.title = element_text(hjust = 1, margin = margin(t= 5, b = 5)),
        axis.title.x = element_text(size = 14, margin = margin(t = 10), hjust = 0.4), 
        axis.title.y = element_text(size = 14, margin = margin(b = 10)))
  }
  
  
  assign(sprintf("boxplot_comp_%d", i),
         temp, envir = .GlobalEnv)
  
}


temp

# Create manual legend

names(domain_palette)[1] <- 'Involvement with\nSubstance Use'
names(domain_palette)[5] <- 'Cognitive\nAbility'
names(domain_palette)[4] <- 'Health &\nPhysical Activity'
names(domain_palette)[2] <- 'Neuroimaging'
names(domain_palette)[3] <- 'Psychological\nTraits'
names(domain_palette)[6] <- 'Community &\nSchool'
names(domain_palette)[8] <- 'Parenting\nBehaviours'
legend_labels <- names(domain_palette)
legend_colors <- domain_palette

# Create a simple legend grob
make_legend <- function(labels, colors, spacing = 0.8) {
  n <- length(labels)
  
  # Calculate text widths to ensure enough space
  max_width <- max(strwidth(labels, units = "inches", family = "Times", cex = 0.8))
  
  # Create legend as a table
  legend_grobs <- lapply(1:n, function(i) {
    grobTree(
      rectGrob(x = 0, width = unit(0.4, "cm"), height = unit(0.4, "cm"), 
               gp = gpar(fill = colors[i], col = NA),
               just = "left"),
      textGrob(labels[i], x = unit(0.6, "cm"), just = "left",
               gp = gpar(fontsize = 10, fontfamily = "Times"))
    )
  })
  
  # Arrange horizontally with explicit spacing
  do.call(arrangeGrob, c(legend_grobs, list(
    ncol = n,
    widths = unit(rep(spacing, n), "inches")  # Control spacing between items
  )))
}

manual_legend <- make_legend(legend_labels, legend_colors, spacing = 1.3)  # Adjust spacing value

sep <- ggplot() +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  theme_void() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    panel.spacing = unit(0, "pt"),
    axis.ticks.length = unit(0, "pt")
  )



master_plot <- ((boxplot_comp_6 + theme(axis.title.x=element_blank(), axis.title.y = element_blank(),
                                        axis.text.x=element_blank())) + sep + (boxplot_comp_5 + theme(axis.title.x=element_blank(),
                                                                                                      axis.text.x=element_blank())) + 
                  plot_layout(widths = c(1, 0.25, 1)))  /
  ((boxplot_comp_4 + theme(axis.title.x=element_blank(),
                           axis.text.x=element_blank())) + sep + (boxplot_comp_1 + theme(axis.title.x=element_blank(),
                                                                                         axis.text.x=element_blank())) + 
     plot_layout(widths = c(1, 0.25, 1))) /
  ((boxplot_comp_2 + theme(axis.title.y = element_blank())) + sep + boxplot_comp_3 + plot_layout(widths = c(1, 0.25, 1))) / 
  wrap_elements(manual_legend) +
  plot_layout(heights = c(1,1,1, 0.08))

ggsave("Thesis Figures/Output/Study 4/master_plot.png", master_plot, 'png',
       path = NULL,
       scale = 1,
       width = 29.7,
       height = 42,
       unit = 'cm',
       dpi = 320,
       limitsize = FALSE,
       bg = 'white')


######  ROC CURVES VISUALIZATION ###############################################


cat("\n--- Creating ROC Curves ---\n")

roc_data_list <- list()

for (comp_name in names(all_results)) {
  comp_results <- all_results[[comp_name]]
  
  for (split_name in names(comp_results)) {
    split_result <- comp_results[[split_name]]
    
    if (!is.null(split_result$roc_object)) {
      roc_obj <- split_result$roc_object
      
      roc_data_list[[paste(comp_name, split_name, sep = "_")]] <- data.frame(
        comparison = split_result$comparison,
        split = split_result$split,
        sensitivity = roc_obj$sensitivities,
        specificity = roc_obj$specificities,
        auc = split_result$auc,
        stringsAsFactors = FALSE
      ) %>%
        mutate(comparison = case_when(
          comparison == '1 vs 2' ~ 'No/Low Use (Reference) v Alcohol L-M',
          comparison == '1 vs 3' ~ 'No/Low Use (Reference) v Alcohol M-H',
          comparison == '2 vs 3' ~ 'Alcohol L-M (Reference) v Alcohol M-H',
          comparison == '1 vs 4' ~ 'No/Low Use (Reference) v Poly L-H',
          comparison == '2 vs 4' ~ 'Alcohol L-M (Reference) v Poly L-H',
          comparison == '3 vs 4' ~ 'Alcohol M-H (Reference) v Poly L-H'
        ))
    }
  }
}

roc_data <- bind_rows(roc_data_list) %>%
  mutate(
    comparison_split = paste0(comparison, " (Split ", split, ")"),
    fpr = 1 - specificity
  )

# Plot: One facet per comparison, showing all splits
p_roc_facets <- ggplot(roc_data, 
                       aes(x = fpr, y = sensitivity, 
                           color = factor(split), group = comparison_split)) +
  geom_line(alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  facet_wrap(~ comparison, ncol = 2) +
  labs(
    title = "ROC Curves by Comparison (All Splits)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Split"
  ) +
  base_theme +
  theme(legend.position = 'none',
        text = element_text(family = "Times", color = "black")
  )

ggsave("Thesis Figures/Output/Study 4/roc_curves_by_comparison_all_splits_wMRI.png",
       p_roc_facets, width = 12, height = 10, dpi = 300)

# Plot: Mean ROC curve per comparison
mean_roc_data <- roc_data %>%
  group_by(comparison) %>%
  summarise(
    mean_auc = mean(auc),
    .groups = "drop"
  ) %>%
  distinct()

# For mean curves, we need to interpolate at common FPR points
fpr_seq <- seq(0, 1, by = 0.01)

# Interpolate each split's ROC curve to common FPR grid
mean_roc_curves <- roc_data %>%
  group_by(comparison, split) %>%
  arrange(fpr) %>%
  reframe(
    fpr_grid = fpr_seq,
    interp_sens = approx(fpr, sensitivity, xout = fpr_seq, rule = 2)$y
  ) %>%
  group_by(comparison, fpr_grid) %>%
  summarise(
    mean_sensitivity = mean(interp_sens, na.rm = TRUE),
    sd_sensitivity = sd(interp_sens, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(mean_roc_data, by = "comparison") %>%
  mutate(comparison_auc = sprintf("%s\n(Mean AUC = %.3f)", comparison, mean_auc))

draw_key_rect <- function(data, params, size) {
  grid::rectGrob(
    width = unit(0.5, "cm"),
    height = unit(0.3, "cm"),  # smaller height
    gp = grid::gpar(
      col = NA,
      fill = alpha(data$colour, data$alpha %||% 1)
    )
  )
}

p_roc_mean <- ggplot(mean_roc_curves, 
                     aes(x = fpr_grid, y = mean_sensitivity, color = comparison_auc)) +
  geom_line(linewidth = 1.2, key_glyph = draw_key_rect) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Mean ROC Curves Across Splits",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Comparison"
  ) +
  base_theme +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, margin = margin(t = 5, l = 5, b = 5, unit = "pt")),
    text = element_text(family = "Times", color = "black"),
    plot.title = element_blank()
  ) +
  scale_color_npg()

ggsave("Thesis Figures/Output/Study 4/roc_curves_mean_by_comparison_wMRI.png",
       p_roc_mean, width = 10, height = 8, dpi = 320)

cat("ROC curves saved\n")


