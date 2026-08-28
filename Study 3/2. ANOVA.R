# ANOVA ---------------------------------------------------------------------
##RAW ANOVA
data.1 <- data.1 %>%
  rename(
    DMN_DMN = rsfmri_c_ngd_dt_ngd_dt,
    FPN_FPN = rsfmri_c_ngd_fo_ngd_fo,
    SN_SN = rsfmri_c_ngd_sa_ngd_sa,
    DMN_SN = rsfmri_c_ngd_dt_ngd_sa,
    DMN_FPN = rsfmri_c_ngd_dt_ngd_fo,
    FPN_SN = rsfmri_c_ngd_fo_ngd_sa
  ) %>%
  mutate(
    DV = factor(DV,
                levels = c(0,1),
                labels = c('Uninitiated', 'Initiated'))
  ) %>%
  left_join(baseline_age, by = 'src_subject_id')

#ANOVA  
raw_aov <- lapply(rsfmri_vars, function(var) {
  formula <- as.formula(paste(var, "~ DV + interview_age"))
  aov(formula, data = data.1)
})
names(raw_aov) <- rsfmri_vars

raw_anova_tables <- lapply(raw_aov, summary)
raw_anova_tables

# Convert each ANOVA table to a data frame and bind them
raw_aov_results <- bind_rows(
  lapply(names(raw_anova_tables), function(net) {
    df <- as.data.frame(raw_anova_tables[[net]][[1]])  # Extract the actual table
    df$Network <- net
    df$Term <- rownames(df)
    rownames(df) <- NULL
    df
  })
)

# Reorder columns
raw_aov_results <- raw_aov_results %>%
  select(Network, Term, everything())

##ALL CONNECTIONS RESIDUAL ANOVA

# Prepare data
rsfmri_vars <- grep("DMN|FPN|SN", names(data.2), value = TRUE)

# Define within vs between network measures
within_network <- c("DMN_DMN", "FPN_FPN", "SN_SN")
between_network <- c("DMN_FPN", "DMN_SN", "FPN_SN")

alpha <- 0.05

# Function to test association between demographics and connectivity
test_demographic_connectivity <- function(connectivity_vars, analysis_name) {
  
  cat("\n========================================\n")
  cat("COVARIATE SELECTION:", analysis_name, "\n")
  cat("========================================\n\n")
  
  # Create long-format data for this subset
  data_long <- data.2 %>%
    pivot_longer(
      cols = all_of(connectivity_vars),
      names_to = 'connectivity',
      values_to = 'rsfmri_value'
    ) %>%
    rename (initiation_group = DV)
  
  # Demographics to test
  demographic_vars <- c("race_4l", "eth_hisp", "sex_2l", "income", "p_edu", "interview_age")
  
  # Test each demographic's association with connectivity
  covariate_tests <- bind_rows(
    lapply(demographic_vars, function(var) {
      x <- data.2[[var]]
      
      # For categorical predictors: ANOVA testing if connectivity differs by category
      if (is.factor(x) || is.character(x)) {
        # Test each connectivity measure
        results <- lapply(connectivity_vars, function(conn_var) {
          formula <- as.formula(paste(conn_var, "~", var))
          fit <- aov(formula, data = data.2)
          tidy_fit <- broom::tidy(fit)
          
          # Get the F-test for the demographic variable
          test_row <- tidy_fit %>% filter(term == var)
          
          if(nrow(test_row) > 0) {
            tibble(
              Demographic = var,
              connectivity = conn_var,
              F_statistic = test_row$statistic,
              df1 = test_row$df,
              p.value = test_row$p.value
            )
          } else {
            NULL
          }
        }) %>% bind_rows()
        
        # Overall significance: use minimum p-value across connectivity measures
        tibble(
          Demographic = var,
          test_type = "ANOVA",
          min_p = min(results$p.value, na.rm = TRUE),
          significant = min(results$p.value, na.rm = TRUE) < alpha,
          details = list(results)
        )
        
        # For continuous predictors: correlation with connectivity
      } else if (is.numeric(x)) {
        results <- lapply(connectivity_vars, function(conn_var) {
          y <- data.2[[conn_var]]
          
          # Remove NAs
          valid_idx <- !is.na(x) & !is.na(y)
          
          if(sum(valid_idx) > 2) {
            cor_test <- cor.test(x[valid_idx], y[valid_idx])
            
            tibble(
              Demographic = var,
              connectivity = conn_var,
              r = cor_test$estimate,
              t_statistic = cor_test$statistic,
              df = cor_test$parameter,
              p.value = cor_test$p.value
            )
          } else {
            NULL
          }
        }) %>% bind_rows()
        
        # Overall significance: use minimum p-value across connectivity measures
        tibble(
          Demographic = var,
          test_type = "Correlation",
          min_p = min(results$p.value, na.rm = TRUE),
          significant = min(results$p.value, na.rm = TRUE) < alpha,
          details = list(results)
        )
      } else {
        NULL
      }
    })
  )
  
  # Identify significant covariates
  sig_covariates <- covariate_tests %>%
    filter(significant) %>%
    pull(Demographic)
  
  cat("Testing association between demographics and connectivity:\n")
  print(covariate_tests %>% select(Demographic, test_type, min_p, significant))
  cat("\nSignificant covariates (p < .05):", paste(sig_covariates, collapse = ", "), "\n")
  if(length(sig_covariates) == 0) cat("No significant covariates found.\n")
  
  list(
    covariate_tests = covariate_tests,
    sig_covariates = sig_covariates,
    analysis_name = analysis_name
  )
}

# Step 1: Identify covariates for within and between network analyses
within_covariates <- test_demographic_connectivity(within_network, "Within-Network connectivity")
between_covariates <- test_demographic_connectivity(between_network, "Between-Network connectivity")
all_covariates <- test_demographic_connectivity(rsfmri_vars, "Triple-Network connectivity")


# Save covariate selection results
write.csv(within_covariates$covariate_tests %>% select(-details), 
          'output/study3/Within-Network Covariate Selection.csv', row.names = FALSE)
write.csv(between_covariates$covariate_tests %>% select(-details), 
          'output/study3/Between-Network Covariate Selection.csv', row.names = FALSE)

# Step 2: Function to run ANOVA with selected covariates
run_connectivity_anova <- function(connectivity_vars, sig_vars, analysis_name, incl_covariate) {
  
  cat("\n\n========================================\n")
  cat("ANOVA:", analysis_name, "\n")
  cat("========================================\n\n")
  
  # Create long-format data
  data_long <- data.2 %>%
    pivot_longer(
      cols = all_of(connectivity_vars),
      names_to = "connectivity",
      values_to = "rsfmri_value"
    ) %>%
    rename(initiation_group = DV) %>%
    mutate(
      connectivity = factor(connectivity),
      initiation_group = factor(initiation_group)
    )
  
  if(!incl_covariate) {
    sig_vars <- list()
  }
  
  # Build the ANOVA formula with the appropriate covariates
  if(length(sig_vars) > 0) {
    formula_str <- paste("rsfmri_value ~ initiation_group * connectivity +", paste(sig_vars, collapse = " + "))
    cat("Including covariates:", paste(sig_vars, collapse = ", "), "\n")
  } else {
    formula_str <- "rsfmri_value ~ initiation_group * connectivity"
    cat("No covariates included (none were significant)\n")
  }
  
  anova_formula <- as.formula(formula_str)
  cat("Formula:", formula_str, "\n\n")
  
  # Fit the ANOVA model
  fit_anova <- aov(anova_formula, data = data_long)
  
  # Get ANOVA table
  anova_results <- broom::tidy(fit_anova) %>%
    mutate(
      eta_squared = sumsq / sum(sumsq),
      sig = ifelse(p.value < 0.05, "*", "")
    )
  
  cat("ANOVA Results:\n")
  print(anova_results)
  cat("Levene's Test (assumption met with p > .05):\n")
  print(data_long %>%
          group_by(connectivity) %>%
          levene_test(rsfmri_value ~ initiation_group))
  cat("Box's M (assumption met with p > .0001):\n")
  print(box_m(data_long[, "rsfmri_value", drop = FALSE], data_long$initiation_group))
  
  # Check if we need post-hoc tests
  needs_posthoc <- anova_results %>%
    filter(term %in% c("initiation_group", "connectivity", "initiation_group:connectivity"),
           p.value < alpha)
  
  posthoc_results <- list()
  
  # Post-hoc for main effect of initiation_group (if significant)
  if("initiation_group" %in% needs_posthoc$term) {
    cat("\nMain effect of initiation_group (Initiation) is significant - descriptive statistics:\n")
    dv_means <- data_long %>%
      group_by(initiation_group) %>%
      summarise(
        M = mean(rsfmri_value, na.rm = TRUE),
        SD = sd(rsfmri_value, na.rm = TRUE),
        N = n(),
        .groups = "drop"
      )
    print(dv_means)
    posthoc_results$dv_means <- dv_means
  }
  
  # Post-hoc for main effect of connectivity (if significant)
  if("connectivity" %in% needs_posthoc$term) {
    cat("\nMain effect of connectivity is significant - running pairwise comparisons:\n")
    
    # Tukey HSD for pairwise comparisons
    library(emmeans)
    emm_connectivity <- emmeans(fit_anova, ~ connectivity)
    connectivity_posthoc <- pairs(emm_connectivity, adjust = "tukey")
    connectivity_posthoc_df <- as.data.frame(connectivity_posthoc)
    
    posthoc_results$connectivity_pairwise <- connectivity_posthoc_df
    print(connectivity_posthoc_df)
    
    # Descriptive statistics by connectivity
    connectivity_means <- data_long %>%
      group_by(connectivity) %>%
      summarise(
        M = mean(rsfmri_value, na.rm = TRUE),
        SD = sd(rsfmri_value, na.rm = TRUE),
        N = n(),
        .groups = "drop"
      )
    print(connectivity_means)
    posthoc_results$connectivity_means <- connectivity_means
  }
  
  # Post-hoc for interaction (if significant)
  if("initiation_group:connectivity" %in% needs_posthoc$term) {
    cat("\nInteraction effect is significant - running simple effects:\n")
    
    # Simple effects: effect of initiation_group at each level of connectivity
    simple_effects <- lapply(levels(data_long$connectivity), function(conn) {
      df_subset <- data_long %>% filter(connectivity == conn)
      t_test <- t.test(rsfmri_value ~ initiation_group, data = df_subset)
      
      # Cohen's d calculation
      pooled_sd <- sqrt(((sum(!is.na(df_subset$rsfmri_value[df_subset$initiation_group == levels(df_subset$initiation_group)[1]])) - 1) * 
                           var(df_subset$rsfmri_value[df_subset$initiation_group == levels(df_subset$initiation_group)[1]], na.rm = TRUE) +
                           (sum(!is.na(df_subset$rsfmri_value[df_subset$initiation_group == levels(df_subset$initiation_group)[2]])) - 1) * 
                           var(df_subset$rsfmri_value[df_subset$initiation_group == levels(df_subset$initiation_group)[2]], na.rm = TRUE)) /
                          (sum(!is.na(df_subset$rsfmri_value)) - 2))
      
      tibble(
        connectivity = conn,
        statistic = t_test$statistic,
        df = t_test$parameter,
        p.value = t_test$p.value,
        mean_diff = diff(t_test$estimate),
        cohen_d = diff(t_test$estimate) / pooled_sd
      )
    }) %>% bind_rows() %>%
      mutate(
        p.adj = p.adjust(p.value, method = "none"),
        sig = ifelse(p.adj < 0.05, "*", "")
      )
    
    print(simple_effects)
    posthoc_results$simple_effects <- simple_effects
    
    # Descriptive statistics for interaction
    interaction_means <- data_long %>%
      group_by(initiation_group, connectivity) %>%
      summarise(
        M = mean(rsfmri_value, na.rm = TRUE),
        SD = sd(rsfmri_value, na.rm = TRUE),
        N = n(),
        .groups = "drop"
      )
    print(interaction_means)
    posthoc_results$interaction_means <- interaction_means
  }
  
  # Return results
  list(
    analysis_name = analysis_name,
    anova_table = anova_results,
    posthoc = posthoc_results,
    model = fit_anova,
    formula = anova_formula,
    data = data_long,
    covariates = sig_vars
  )
}

# Step 3: Run ANOVAs with their respective significant covariates
within_anova <- run_connectivity_anova(
  within_network, 
  within_covariates$sig_covariates, 
  "Within-Network connectivity",
  FALSE
)

between_anova <- run_connectivity_anova(
  between_network, 
  between_covariates$sig_covariates, 
  "Between-Network connectivity",
  FALSE
)

full_anova <- run_connectivity_anova(
  rsfmri_vars, 
  all_covariates$sig_covariates, 
  "Triple-Network connectivity",
  FALSE
)

# Step 4: Save ANOVA results
write.csv(within_anova$anova_table, 'output/study3/ANOVA_2x3 Within-Network Results.csv', row.names = FALSE)
write.csv(between_anova$anova_table, 'output/study3/ANOVA_2x3 Between-Network Results.csv', row.names = FALSE)
write.csv(full_anova$anova_table, 'output/study3/ANOVA_2x6 Triple-Network Results.csv', row.names = FALSE)


# Save post-hoc results if they exist
if(length(within_anova$posthoc) > 0) {
  if(!is.null(within_anova$posthoc$simple_effects)) {
    write.csv(within_anova$posthoc$simple_effects, 
              'output/study3/ANOVA Within-Network PostHoc Simple Effects.csv', row.names = FALSE)
  }
  if(!is.null(within_anova$posthoc$interaction_means)) {
    write.csv(within_anova$posthoc$interaction_means, 
              'output/study3/ANOVA Within-Network Interaction Means.csv', row.names = FALSE)
  }
  if(!is.null(within_anova$posthoc$connectivity_means)) {
    write.csv(within_anova$posthoc$connectivity_means, 
              'output/study3/ANOVA Within-Network connectivity Means.csv', row.names = FALSE)
  }
}

if(length(between_anova$posthoc) > 0) {
  if(!is.null(between_anova$posthoc$simple_effects)) {
    write.csv(between_anova$posthoc$simple_effects, 
              'output/study3/ANOVA Between-Network PostHoc Simple Effects.csv', row.names = FALSE)
  }
  if(!is.null(between_anova$posthoc$interaction_means)) {
    write.csv(between_anova$posthoc$interaction_means, 
              'output/study3/ANOVA Between-Network Interaction Means.csv', row.names = FALSE)
  }
  if(!is.null(between_anova$posthoc$connectivity_means)) {
    write.csv(between_anova$posthoc$connectivity_means, 
              'output/study3/ANOVA Between-Network connectivity Means.csv', row.names = FALSE)
  }
}

# Step 5: Print comprehensive summary
cat("\n\n=== COMPREHENSIVE SUMMARY ===\n")

cat("\n--- WITHIN-NETWORK CONNECTIVITY ---\n")
cat("Covariates included:", paste(within_anova$covariates, collapse = ", "), "\n")
if(length(within_anova$covariates) == 0) cat("(No significant covariates)\n")
cat("\nSignificant ANOVA effects (p < .05):\n")
sig_within <- within_anova$anova_table %>% 
  filter(p.value < 0.05) %>% 
  select(term, statistic, df, p.value, eta_squared)
if(nrow(sig_within) > 0) {
  print(sig_within)
} else {
  cat("No significant effects found.\n")
}

cat("\n--- BETWEEN-NETWORK CONNECTIVITY ---\n")
cat("Covariates included:", paste(between_anova$covariates, collapse = ", "), "\n")
if(length(between_anova$covariates) == 0) cat("(No significant covariates)\n")
cat("\nSignificant ANOVA effects (p < .05):\n")
sig_between <- between_anova$anova_table %>% 
  filter(p.value < 0.05) %>% 
  select(term, statistic, df, p.value, eta_squared)
if(nrow(sig_between) > 0) {
  print(sig_between)
} else {
  cat("No significant effects found.\n")
}

