# Mediation Approach with Sensation Seeking -------------------------------


#                               Mediators                                         #


mh_y_upps <- read_csv('data/5.1/mh_y_upps.csv') %>% 
  subset (eventname == 'baseline_year_1_arm_1') %>%
  select(src_subject_id, upps_y_ss_negative_urgency,
         upps_y_ss_lack_of_perseverance, upps_y_ss_lack_of_planning,
         upps_y_ss_sensation_seeking, upps_y_ss_positive_urgency)


# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(mh_y_upps, by = 'src_subject_id')

rm(mh_y_upps)

mh_p_cbcl <- read_csv('data/5.1/mh_p_cbcl.csv') %>% 
  subset (eventname == 'baseline_year_1_arm_1') %>%
  select(src_subject_id, starts_with('cbcl_scr_syn_') & ends_with('_r'))

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(mh_p_cbcl, by = 'src_subject_id')

rm(mh_p_cbcl)


###DMN:SN Mediation -> Insignificant
# Model for mediator (linear regression)
model.M1 <- lm(upps_y_ss_sensation_seeking ~ DMN_SN, data = data.2)

# Model for outcome (logistic regression)
model.Y1 <- glm(DV ~ DMN_SN + upps_y_ss_sensation_seeking,
                family = binomial(link = "logit"), data = data.2)

# Mediation analysis
med.out.dmnsn <- mediation::mediate(model.M1, model.Y1,
                                    treat = "DMN_SN",
                                    mediator = "upps_y_ss_sensation_seeking",
                                    sims = 5000)

summary(med.out.dmnsn)

# --- 1. Extract coefficients + p-values ---
coef.M <- summary(model.M1)$coefficients
coef.Y <- summary(model.Y1)$coefficients

a <- coef.M["DMN_SN", "Estimate"]
a_p <- coef.M["DMN_SN", "Pr(>|t|)"]

b <- coef.Y["upps_y_ss_sensation_seeking", "Estimate"]
b_p <- coef.Y["upps_y_ss_sensation_seeking", "Pr(>|z|)"]

cprime <- coef.Y["DMN_SN", "Estimate"]
cprime_p <- coef.Y["DMN_SN", "Pr(>|z|)"]

indirect <- med.out.dmnsn$d0
total <- med.out.dmnsn$tau.coef

# helper function for stars
pstars <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}

# format labels
a_lab <- paste0("a = ", round(a, 2), " ", pstars(a_p))
b_lab <- paste0("b = ", round(b, 2), " ", pstars(b_p))
cprime_lab <- paste0("c' = ", round(cprime, 2), " ", pstars(cprime_p))

# --- 2. Define node positions (triangle layout) ---
nodes <- data.frame(
  name = c("DMN:SN", "Sensation Seeking", "Drug\nInitiation"),
  x = c(0, 0.5, 1),
  y = c(0, 0.5, 0)
)

# --- 3. Define edges ---
edges <- data.frame(
  from = c("DMN:SN", "Sensation Seeking", "DMN:SN"),
  to   = c("Sensation Seeking", "Drug\nInitiation", "Drug\nInitiation"),
  label = c(a_lab, b_lab, cprime_lab),
  label_dx = c(-0.03, 0.03, 0),    # horizontal offsets
  label_dy = c(0.03, 0.03, 0.03)         # vertical offsets
) %>%
  rowwise() %>%
  mutate(
    # original coordinates
    x0 = nodes$x[match(from, nodes$name)],
    y0 = nodes$y[match(from, nodes$name)],
    x1 = nodes$x[match(to, nodes$name)],
    y1 = nodes$y[match(to, nodes$name)],
    
    # shorten arrows by 10% at each end (adjust factor if needed)
    shorten = 0.15,
    x_start = x0 + shorten * (x1 - x0),
    y_start = y0 + shorten * (y1 - y0),
    x_end   = x1 - shorten * (x1 - x0),
    y_end   = y1 - shorten * (y1 - y0),
    
    # label positions (still midpoint, plus offsets you already have)
    label_x = (x0 + x1) / 2 + label_dx,
    label_y = (y0 + y1) / 2 + label_dy
  ) %>%
  ungroup()

edges$angle <- atan2(edges$y_end - edges$y_start,
                     edges$x_end - edges$x_start) * 180 / pi

# --- 4. Plot diagram ---
sem_plot1 <- ggplot() +
  # arrows
  geom_curve(
    data = edges,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    curvature = 0,
    arrow = arrow(length = unit(0.03, "npc")),
    lineend = "round", linewidth = 1
  ) +
  # edge labels
  geom_text(
    data = edges,
    aes(x = label_x, y = label_y, label = label, angle = angle),
    size = 4, family = 'Times New Roman'
  ) +
  # nodes
  geom_label(
    data = nodes,
    aes(x = x, y = y, label = name),
    size = 5, fill = "white", label.size = 0, color = "black", family = 'Times New Roman'
  ) +
  theme_void() +
  coord_fixed(xlim = c(-0.1, 1.1), ylim = c(-0.1, 0.6))+
  labs(title = "DMN:SN to Initiation via Sensation Seeking") +
  theme(plot.title = element_text(
    hjust = 0.5,        # center
    size = 16,          # adjust font size
    family = "Times New Roman"),
    panel.border = element_blank())
sem_plot1

ggsave('Thesis Figures/Output/Study 3/study 3 dmn_sn mediation plot.png', sem_plot1, 'png',
       path = NULL,
       scale = 1,
       width = 6,
       height = 4,
       units = 'in',
       dpi = 320,
       limitsize = FALSE,
       bg = 'transparent',
       create.dir = FALSE)

###FPN:SN Mediation -> Significant
# Model for mediator (linear regression)
model.M <- lm(upps_y_ss_sensation_seeking ~ FPN_SN, data = data.2)

# Model for outcome (logistic regression)
model.Y <- glm(DV ~ FPN_SN + upps_y_ss_sensation_seeking,
               family = binomial(link = "logit"), data = data.2)

# Mediation analysis
med.out.fpnsn <- mediation::mediate(model.M, model.Y,
                                    treat = "FPN_SN",
                                    mediator = "upps_y_ss_sensation_seeking",
                                    sims = 5000)

summary(med.out.fpnsn)




# --- 1. Extract coefficients + p-values ---
coef.M <- summary(model.M)$coefficients
coef.Y <- summary(model.Y)$coefficients

a <- coef.M["FPN_SN", "Estimate"]
a_p <- coef.M["FPN_SN", "Pr(>|t|)"]

b <- coef.Y["upps_y_ss_sensation_seeking", "Estimate"]
b_p <- coef.Y["upps_y_ss_sensation_seeking", "Pr(>|z|)"]

cprime <- coef.Y["FPN_SN", "Estimate"]
cprime_p <- coef.Y["FPN_SN", "Pr(>|z|)"]

indirect <- med.out.fpnsn$d0
total <- med.out.fpnsn$tau.coef

# helper function for stars
pstars <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}

# format labels
a_lab <- paste0("a = ", round(a, 2), " ", pstars(a_p))
b_lab <- paste0("b = ", round(b, 2), " ", pstars(b_p))
cprime_lab <- paste0("c' = ", round(cprime, 2), " ", pstars(cprime_p))

# --- 2. Define node positions (triangle layout) ---
nodes <- data.frame(
  name = c("FPN:SN", "Sensation Seeking", "Drug\nInitiation"),
  x = c(0, 0.5, 1),
  y = c(0, 0.5, 0)
)

# --- 3. Define edges ---
edges <- data.frame(
  from = c("FPN:SN", "Sensation Seeking", "FPN:SN"),
  to   = c("Sensation Seeking", "Drug\nInitiation", "Drug\nInitiation"),
  label = c(a_lab, b_lab, cprime_lab),
  label_dx = c(-0.03, 0.03, 0),    # horizontal offsets
  label_dy = c(0.03, 0.03, 0.03)         # vertical offsets
) %>%
  rowwise() %>%
  mutate(
    # original coordinates
    x0 = nodes$x[match(from, nodes$name)],
    y0 = nodes$y[match(from, nodes$name)],
    x1 = nodes$x[match(to, nodes$name)],
    y1 = nodes$y[match(to, nodes$name)],
    
    # shorten arrows by 10% at each end (adjust factor if needed)
    shorten = 0.15,
    x_start = x0 + shorten * (x1 - x0),
    y_start = y0 + shorten * (y1 - y0),
    x_end   = x1 - shorten * (x1 - x0),
    y_end   = y1 - shorten * (y1 - y0),
    
    # label positions (still midpoint, plus offsets you already have)
    label_x = (x0 + x1) / 2 + label_dx,
    label_y = (y0 + y1) / 2 + label_dy
  ) %>%
  ungroup()

edges$angle <- atan2(edges$y_end - edges$y_start,
                     edges$x_end - edges$x_start) * 180 / pi

# --- 4. Plot diagram ---
sem_plot2 <- ggplot() +
  # arrows
  geom_curve(
    data = edges,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    curvature = 0,
    arrow = arrow(length = unit(0.03, "npc")),
    lineend = "round", linewidth = 1
  ) +
  # edge labels
  geom_text(
    data = edges,
    aes(x = label_x, y = label_y, label = label, angle = angle),
    size = 4, family = 'Times New Roman'
  ) +
  # nodes
  geom_label(
    data = nodes,
    aes(x = x, y = y, label = name),
    size = 5, fill = "white", label.size = 0, color = "black", family = 'Times New Roman'
  ) +
  theme_void() +
  coord_fixed(xlim = c(-0.1, 1.1), ylim = c(-0.1, 0.6))+
  labs(title = "FPN:SN to Initiation via Sensation Seeking") +
  theme(plot.title = element_text(
    hjust = 0.5,        # center
    size = 16,          # adjust font size
    family = "Times New Roman"),
    panel.border = element_blank())

sem_plot2

ggsave('Thesis Figures/Output/Study 3/study 3 fpn_sn mediation plot.png', sem_plot2, 'png',
       path = NULL,
       scale = 1,
       width = 6,
       height = 4,
       units = 'in',
       dpi = 320,
       limitsize = FALSE,
       bg = 'transparent',
       create.dir = FALSE)


comb_plot <- sem_plot1 / sem_plot2
comb_plot <- comb_plot + theme(panel.border = element_blank())

ggsave('Thesis Figures/Output/Study 3/study 3 dual mediation plot.png', comb_plot, 'png',
       path = NULL,
       scale = 1,
       width = 6,
       height = 9,
       units = 'in',
       dpi = 320,
       limitsize = FALSE,
       bg = 'transparent',
       create.dir = FALSE)

rm(coef.M, coef.Y, edges, med.out.dmnsn, med.out.fpnsn, med.out.interact, model.M,
   model.M1, model.M2, model.Y, model.Y1, model.Y2, nodes)

# Function to convert mediation summary to dataframe
med_summary_to_df <- function(med_obj) {
  s <- summary(med_obj)
  data.frame(
    ACME_Estimate = s$d0,
    ACME_CI_Lower = s$d0.ci[1],
    ACME_CI_Upper = s$d0.ci[2],
    ACME_p = s$d0.p,
    
    ADE_Estimate = s$z0,
    ADE_CI_Lower = s$z0.ci[1],
    ADE_CI_Upper = s$z0.ci[2],
    ADE_p = s$z0.p,
    
    Total_Effect = s$tau.coef,
    Total_CI_Lower = s$tau.ci[1],
    Total_CI_Upper = s$tau.ci[2],
    Total_p = s$tau.p,
    
    Prop_Mediated = s$n0,
    Prop_CI_Lower = s$n0.ci[1],
    Prop_CI_Upper = s$n0.ci[2],
    Prop_p = s$n0.p
  )
}

# Convert both mediation results
df_dmnsn <- med_summary_to_df(med.out.dmnsn)
df_fpnsn <- med_summary_to_df(med.out.fpnsn)
df_dmfpnsn <- med_summary_to_df(med.out.interact)

# Add model name column for clarity
df_dmnsn$Model <- "DMN:SN"
df_fpnsn$Model <- "FPN:SN"

# Combine
final_df <- rbind(df_dmnsn, df_fpnsn, df_dmfpnsn)
final_df <- t(final_df)

# Write to Excel
write.xlsx(final_df, file = "output/study3/study3_mediation_results.xlsx", rowNames = TRUE)

save.image("output/study3/study3_data.RData")
