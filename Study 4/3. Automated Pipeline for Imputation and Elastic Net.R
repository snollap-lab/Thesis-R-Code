# Fully automated imputation and elastic net pipeline --------------------------
system("powercfg -change standby-timeout-ac 0")  # Disables sleep


tryCatch({
  # ==============================================================================
  # MICE PREPROCESSING FOR 4-LEVEL OUTCOME (PAIRWISE COMPARISONS)
  # ==============================================================================
  # This script automates MICE imputation for all pairwise comparisons of a
  # 4-level outcome variable (6 total comparisons: 4 choose 2)
  # ==============================================================================
  
  # Load Required Packages -------------------------------------------------------
  library(dplyr)
  library(tidyr)
  library(mice)
  library(rsample)
  library(recipes)
  library(doParallel)
  library(doRNG)
  library(rstatix)
  library(fastDummies)
  library(blastula)
  library(keyring)
  
  
  
  # Configure Parallelization ----------------------------------------------------
  nCores <- min(parallel::detectCores()) - 1
  options(mc.cores = nCores, cores = nCores)
  doParallel::registerDoParallel(cores = nCores)
  cat("### Using", foreach::getDoParWorkers(), "cores\n")
  
  # Helper Functions -------------------------------------------------------------
  factor_sum <- function(x) {
    n <- table(x)
    proportion <- round(prop.table(n), 4)
    percentage <- (round(prop.table(n), 4)) * 100
    OUT <- cbind(n, proportion, percentage)
  }
  
  # Generate All Pairwise Comparisons --------------------------------------------
  ng4_levels <- 1:4
  pairs <- combn(ng4_levels, 2, simplify = FALSE)
  
  cat("\n### Pairwise Comparisons to Process:\n")
  for (i in seq_along(pairs)) {
    cat(sprintf("Comparison %d: Level %d vs Level %d\n", 
                i, pairs[[i]][1], pairs[[i]][2]))
  }
  
  # Master Function: Process One Split for One Comparison -----------------------
  process_split_comparison <- function(data, split_seed, split_num, 
                                       comparison_num, level1, level2,
                                       min_class_size = 30) {
    
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Processing Split %d, Comparison %d (%d vs %d)\n", 
                split_num, comparison_num, level1, level2))
    cat(sprintf("========================================\n"))
    
    # Subset data for this comparison
    data_subset <- data %>%
      filter(ng4_class %in% c(level1, level2)) %>%
      mutate(
        # Create binary outcome: level1 = 1, level2 = 0
        # Keep as numeric throughout to avoid factor conversion issues
        ng4_class_binary = ifelse(ng4_class == level2, 1, 0),
        ng4_class = ng4_class_binary
      ) %>%
      select(-ng4_class_binary)
    
    # Check total sample size
    n_level2 <- sum(data_subset$ng4_class == 1)
    n_level1 <- sum(data_subset$ng4_class == 0)
    
    cat(sprintf("Total sample sizes - Level %d: %d, Level %d: %d\n",
                level1, n_level1, level2, n_level2))
    cat(sprintf("Binary coding: Level %d = 0, Level %d = 1\n", level1, level2))
    
    # Check if either class is too small
    if (n_level1 < min_class_size || n_level2 < min_class_size) {
      warning(sprintf("⚠ Insufficient sample size for comparison %d vs %d (n1=%d, n2=%d). Skipping.\n",
                      level1, level2, n_level1, n_level2))
      return(NULL)
    }
    
    # Calculate minimum expected training size
    min_train_expected <- min(n_level1, n_level2) * 0.75
    
    if (min_train_expected < 20) {
      warning(sprintf("⚠ Expected training class size too small (%d). Adjusting split proportion.\n",
                      floor(min_train_expected)))
      prop_train <- max(0.85, 25 / min(n_level1, n_level2))
      prop_train <- min(prop_train, 0.95)
    } else {
      prop_train <- 0.75
    }
    
    cat(sprintf("Using train proportion: %.2f\n", prop_train))
    
    # Create train/test split with adjusted proportion
    # Create a temporary factor for stratification only
    data_subset_temp <- data_subset %>%
      mutate(ng4_class_factor = factor(ng4_class))
    
    set.seed(split_seed)
    split_obj <- initial_split(data_subset_temp, prop = prop_train, strata = ng4_class_factor)
    
    trn <- rsample::training(split_obj) %>% select(-ng4_class_factor)
    tst <- rsample::testing(split_obj) %>% select(-ng4_class_factor)
    
    # Verify minimum class sizes in training set
    trn_level1 <- sum(trn$ng4_class == 0)
    trn_level2 <- sum(trn$ng4_class == 1)
    
    if (trn_level1 < 10 || trn_level2 < 10) {
      warning(sprintf("⚠ Training set has insufficient observations (n1=%d, n2=%d). Skipping.\n",
                      trn_level1, trn_level2))
      return(NULL)
    }
    
    cat(sprintf("Training sizes - Level %d (coded as 1): %d, Level %d (coded as 0): %d\n",
                level1, trn_level1, level2, trn_level2))
    
    # Validation checks
    validate_split(trn, tst, level1, level2)
    
    # Run MICE on training and test
    trn_imputed <- run_mice_pipeline(trn, split_num, comparison_num, "trn")
    tst_imputed <- run_mice_pipeline(tst, split_num, comparison_num, "tst")
    
    # Return results
    list(
      trn_obs = trn,
      tst_obs = tst,
      trn_imp = trn_imputed$imp_df,
      tst_imp = tst_imputed$imp_df,
      trn_mice = trn_imputed$mice_obj,
      tst_mice = tst_imputed$mice_obj,
      sample_sizes = list(
        total = c(n_level1, n_level2),
        train = c(trn_level1, trn_level2),
        test = c(sum(tst$ng4_class == 0), sum(tst$ng4_class == 1))
      )
    )
  }
  
  # Validation Function ----------------------------------------------------------
  validate_split <- function(trn, tst, level1, level2) {
    
    # ID check
    id_overlap <- any(trn$participant_id %in% tst$participant_id)
    if (!id_overlap) {
      cat("✓ No ID overlap between train and test\n")
    } else {
      stop("✗ ID overlap detected between train and test")
    }
    
    # Outcome distribution check
    trn_prop <- mean(trn$ng4_class == level1)
    tst_prop <- mean(tst$ng4_class == level1)
    
    if (trn_prop >= 0.4 && trn_prop <= 0.6 && tst_prop >= 0.4 && tst_prop <= 0.6) {
      cat(sprintf("✓ Outcome proportions reasonable - Train: %.2f, Test: %.2f\n",
                  trn_prop, tst_prop))
    } else {
      warning(sprintf("⚠ Outcome proportions imbalanced - Train: %.2f, Test: %.2f\n",
                      trn_prop, tst_prop))
    }
    
    # Factor variable check
    factor_vars <- c('race_4l', 'eth_hisp', 'income', 'religion', 'p_edu', 
                     'sex_2l', 'rec_bin', 'exp_sub', 'tbi_injury', 
                     'mh_p_kbi__school_006', 'det_susp', 'se_services', 'cct')
    
    all_represented <- TRUE
    for (var in factor_vars) {
      if (var %in% names(trn)) {
        trn_levels <- length(unique(trn[[var]]))
        tst_levels <- length(unique(tst[[var]]))
        if (trn_levels != tst_levels) {
          warning(sprintf("⚠ %s has different levels: Train=%d, Test=%d\n",
                          var, trn_levels, tst_levels))
          all_represented <- FALSE
        }
      }
    }
    
    if (all_represented) {
      cat("✓ All factor levels represented in both sets\n")
    }
  }
  
  # MICE Pipeline Function -------------------------------------------------------
  run_mice_pipeline <- function(data, split_num, comparison_num, dataset_type) {
    
    cat(sprintf("\n--- MICE for %s dataset ---\n", toupper(dataset_type)))
    
    # ============================================================================
    # FIX: Drop unused factor levels to prevent dimension mismatches
    # ============================================================================
    factor_cols <- sapply(data, is.factor)
    if (any(factor_cols)) {
      data[factor_cols] <- lapply(data[factor_cols], droplevels)
      cat("✓ Dropped unused factor levels\n")
    }
    
    # ============================================================================
    # FIX: Check for and remove constant/near-constant variables
    # ============================================================================
    n_unique <- sapply(data, function(x) {
      if (is.factor(x)) {
        length(levels(x))
      } else {
        length(unique(na.omit(x)))
      }
    })
    
    constant_vars <- names(n_unique)[n_unique <= 1]
    
    if (length(constant_vars) > 0) {
      cat(sprintf("⚠ Removing %d constant/single-level variables: %s\n", 
                  length(constant_vars), 
                  paste(constant_vars, collapse = ", ")))
      data <- data[, !names(data) %in% constant_vars, drop = FALSE]
    }
    
    # Extract variable names
    all_vars <- names(data)
    miss_vars <- names(data)[colSums(is.na(data)) > 0]
    
    cat(sprintf("Variables with missingness: %d\n", length(miss_vars)))
    
    # Build predictor matrix
    predictor_matrix <- matrix(0, ncol = length(all_vars), nrow = length(all_vars))
    rownames(predictor_matrix) <- all_vars
    colnames(predictor_matrix) <- all_vars
    
    # Specify imputer variables
    imputer_vars <- c('age_baseline', 'sex_2l', 'race_4l', 'eth_hisp', 
                      'income', 'religion', 'p_edu')
    imputer_vars <- intersect(unique(imputer_vars), all_vars)
    
    # Specify imputed variables
    imputed_only_vars <- miss_vars
    imputed_vars <- intersect(unique(c(imputed_only_vars, imputer_vars)), miss_vars)
    
    # ============================================================================
    # FIX: Check for collinearity and reduce predictor matrix if needed
    # ============================================================================
    # Compute correlation matrix for numeric variables
    numeric_vars <- names(data)[sapply(data, is.numeric)]
    numeric_vars <- setdiff(numeric_vars, "ng4_class")  # Exclude outcome
    
    if (length(numeric_vars) > 1) {
      cor_matrix <- cor(data[, numeric_vars], use = "pairwise.complete.obs")
      
      # Find highly correlated pairs (|r| > 0.95)
      high_cor <- which(abs(cor_matrix) > 0.95 & abs(cor_matrix) < 1, arr.ind = TRUE)
      
      if (nrow(high_cor) > 0) {
        high_cor_pairs <- unique(t(apply(high_cor, 1, sort)))
        high_cor_pairs <- high_cor_pairs[!duplicated(high_cor_pairs), , drop = FALSE]
        
        cat(sprintf("⚠ Found %d highly correlated variable pairs (|r| > 0.95)\n", 
                    nrow(high_cor_pairs)))
        
        # Remove one variable from each highly correlated pair
        vars_to_remove <- c()
        for (i in 1:nrow(high_cor_pairs)) {
          var1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
          var2 <- rownames(cor_matrix)[high_cor_pairs[i, 2]]
          
          # Keep the variable with less missingness
          miss1 <- sum(is.na(data[[var1]]))
          miss2 <- sum(is.na(data[[var2]]))
          
          if (miss1 > miss2) {
            vars_to_remove <- c(vars_to_remove, var1)
          } else {
            vars_to_remove <- c(vars_to_remove, var2)
          }
        }
        
        vars_to_remove <- unique(vars_to_remove)
        
        if (length(vars_to_remove) > 0) {
          cat(sprintf("  Removing %d collinear variables: %s\n", 
                      length(vars_to_remove),
                      paste(vars_to_remove, collapse = ", ")))
          data <- data[, !names(data) %in% vars_to_remove, drop = FALSE]
          
          # Update variable lists
          all_vars <- names(data)
          miss_vars <- names(data)[colSums(is.na(data)) > 0]
          imputer_vars <- intersect(imputer_vars, all_vars)
          imputed_vars <- intersect(imputed_vars, all_vars)
          
          # Rebuild predictor matrix with new dimensions
          predictor_matrix <- matrix(0, ncol = length(all_vars), nrow = length(all_vars))
          rownames(predictor_matrix) <- all_vars
          colnames(predictor_matrix) <- all_vars
        }
      }
    }
    
    # Construct predictor matrix
    imputer_matrix <- predictor_matrix
    imputer_matrix[, imputer_vars] <- 1
    
    imputed_matrix <- predictor_matrix
    imputed_matrix[imputed_vars, ] <- 1
    
    predictor_matrix <- imputer_matrix * imputed_matrix
    diag(predictor_matrix) <- 0
    
    # ============================================================================
    # FIX: Limit number of predictors per variable to avoid over-specification
    # ============================================================================
    max_predictors <- floor(nrow(data) / 3)  # Rule of thumb: n/3 predictors max
    
    for (var in rownames(predictor_matrix)) {
      n_predictors <- sum(predictor_matrix[var, ])
      
      if (n_predictors > max_predictors) {
        cat(sprintf("⚠ Variable '%s' has %d predictors (max: %d), reducing...\n", 
                    var, n_predictors, max_predictors))
        
        # Keep only the most important predictors (imputer_vars + highest correlation)
        predictors_for_var <- names(which(predictor_matrix[var, ] == 1))
        
        # Always keep imputer_vars
        keep_vars <- intersect(predictors_for_var, imputer_vars)
        
        # Add other predictors based on correlation
        other_predictors <- setdiff(predictors_for_var, imputer_vars)
        
        if (length(other_predictors) > 0 && length(keep_vars) < max_predictors) {
          if (is.numeric(data[[var]])) {
            # Calculate correlations with numeric predictors
            numeric_predictors <- intersect(other_predictors, numeric_vars)
            
            if (length(numeric_predictors) > 0) {
              cors <- abs(cor(data[[var]], data[, numeric_predictors], 
                              use = "pairwise.complete.obs"))
              top_predictors <- names(sort(cors[1,], decreasing = TRUE)[1:(max_predictors - length(keep_vars))])
              keep_vars <- c(keep_vars, top_predictors)
            }
          }
          
          # Fill remaining slots with other predictors if needed
          if (length(keep_vars) < max_predictors) {
            remaining_slots <- max_predictors - length(keep_vars)
            keep_vars <- c(keep_vars, other_predictors[1:min(remaining_slots, length(other_predictors))])
          }
        }
        
        # Update predictor matrix
        predictor_matrix[var, ] <- 0
        predictor_matrix[var, keep_vars] <- 1
      }
    }
    
    # ============================================================================
    # FIX: Wrap dry run in tryCatch to catch dimension errors early
    # ============================================================================
    dry_mice <- tryCatch({
      mice(data = data, m = 1, predictorMatrix = predictor_matrix, maxit = 0)
    }, error = function(e) {
      cat("✗ Error in dry MICE run:\n")
      cat(sprintf("  %s\n", e$message))
      
      # Try with default predictor matrix as fallback
      cat("  Attempting with quickpred predictor matrix...\n")
      pred_quick <- quickpred(data, mincor = 0.1, minpuc = 0.1)
      mice(data = data, m = 1, predictorMatrix = pred_quick, maxit = 0)
    })
    
    predictor_matrix <- dry_mice$predictorMatrix
    dry_mice$method[setdiff(all_vars, imputed_vars)] <- ""
    
    # ============================================================================
    # FIX: Add verbose output for first imputation to catch issues
    # ============================================================================
    M <- 5
    set.seed(808915)
    
    # Run first imputation with verbose to diagnose issues
    cat("Running first imputation (verbose)...\n")
    mice_first <- tryCatch({
      mice::mice(data = data, m = 1, print = FALSE,
                 predictorMatrix = predictor_matrix,
                 method = dry_mice$method, maxit = 20)
    }, error = function(e) {
      cat(sprintf("✗ Error in first MICE imputation: %s\n", e$message))
      
      # Print diagnostic information
      cat("\nDiagnostic Information:\n")
      cat(sprintf("  Dataset dimensions: %d rows x %d cols\n", nrow(data), ncol(data)))
      cat("\n  Factor variables and their levels:\n")
      for (var in names(data)[sapply(data, is.factor)]) {
        cat(sprintf("    %s: %d levels - %s\n", 
                    var, 
                    nlevels(data[[var]]),
                    paste(levels(data[[var]]), collapse = ", ")))
      }
      stop(e)
    })
    
    # Continue with remaining imputations in parallel
    if (M > 1) {
      cat(sprintf("Running remaining %d imputations in parallel...\n", M - 1))
      mice_rest <- foreach(i = seq_len(M - 1), .combine = ibind) %dorng% {
        mice::mice(data = data, m = 1, print = FALSE,
                   predictorMatrix = predictor_matrix,
                   method = dry_mice$method, maxit = 20)
      }
      
      mice_out <- ibind(mice_first, mice_rest)
    } else {
      mice_out <- mice_first
    }
    
    # Check for logged events
    if (is.null(mice_out$loggedEvents)) {
      cat("✓ No logged events during imputation\n")
    } else {
      warning("⚠ Logged events detected:\n")
      print(mice_out$loggedEvents)
    }
    
    # Create long format dataframe
    imp_df <- mice::complete(mice_out, "long", include = TRUE)
    
    cat(sprintf("✓ MICE completed: %d imputations\n", M))
    
    list(
      mice_obj = mice_out,
      imp_df = imp_df
    )
  }
  # Center and Scale Function ----------------------------------------------------
  center_scale_data <- function(imp_df) {
    
    cat("\n--- Center and Scale ---\n")
    
    # IMPORTANT: ng4_class should already be numeric (0/1) at this point
    # Verify this
    cat(sprintf("ng4_class unique values before C&S: %s\n", 
                paste(unique(imp_df$ng4_class[imp_df$.imp != 0]), collapse=", ")))
    
    # Store ng4_class and remove observed data
    imp_df_cs <- imp_df %>% subset(.imp != 0)
    
    # Store observed data separately
    imp_df_imp0 <- imp_df %>% subset(.imp == 0) %>% select(-participant_id)
    
    # Store .imp, .id, pds, and ng4_class (don't normalize these)
    imp_df_imp15 <- imp_df %>% 
      subset(.imp != 0) %>% 
      select(.imp, .id, pds, ng4_class)
    
    # Build recipe - exclude ng4_class from normalization
    recipe_prep <- imp_df_cs %>%
      recipe(ng4_class ~ .) %>%
      step_rm(ends_with('.id') | ends_with('.imp') | participant_id | pds) %>%
      update_role(ng4_class, new_role = "outcome") %>%
      step_normalize(all_numeric_predictors()) %>%
      prep()
    
    # Apply recipe
    imp_df_cs_2 <- recipe_prep %>% bake(new_data = imp_df_cs)
    
    # Remove ng4_class from centered/scaled data (we'll add it back from imp_df_imp15)
    imp_df_cs_2 <- imp_df_cs_2 %>% select(-ng4_class)
    
    # Combine: add back .imp, .id, pds, and ng4_class
    imp_df_cs_3 <- cbind(imp_df_imp15, imp_df_cs_2)
    
    # Add back observed data
    imp_df_cs_final <- rbind(imp_df_imp0, imp_df_cs_3)
    
    # Validate ng4_class is still 0/1
    unique_outcome <- unique(imp_df_cs_final$ng4_class[imp_df_cs_final$.imp != 0])
    cat(sprintf("ng4_class unique values after C&S: %s\n", 
                paste(sort(unique_outcome), collapse=", ")))
    
    if (!all(unique_outcome %in% c(0, 1))) {
      stop("ERROR: ng4_class is not binary (0/1) after center and scale!")
    }
    
    # Validate centering and scaling of predictors
    cs_check <- imp_df_cs_final %>%
      subset(.imp != 0) %>%
      select_if(is.numeric) %>%
      select(-.imp, -.id, -ng4_class) %>%
      get_summary_stats(type = "common") %>%
      select(variable, mean, sd)
    
    mean_ok <- all(cs_check$mean >= -0.10 & cs_check$mean <= 0.10)
    sd_ok <- all(cs_check$sd >= 0.95 & cs_check$sd <= 1.05)
    
    if (mean_ok && sd_ok) {
      cat("✓ Center and scale validated\n")
    } else {
      warning("⚠ Center and scale validation failed\n")
    }
    
    imp_df_cs_final
  }
  
  # Create Dummy Variables -------------------------------------------------------
  create_dummy_variables <- function(imp_df_cs_final) {
    
    cat("\n--- Creating Dummy Variables ---\n")
    
    # Extract factor variables
    factor_vars <- imp_df_cs_final %>% select_if(is.factor)
    
    # Create dummy codes for ALL categories
    factor_dummy <- fastDummies::dummy_cols(factor_vars, remove_first_dummy = FALSE)
    
    # Define expected dummy variable names (excluding reference groups)
    expected_dummies <- c(
      #Baseline Use
      'alc0_Yes', 'mj0_Yes', 'nic0_Yes', 'oth0_Yes',
      # Demographics
      "race_4l_Asian", "race_4l_Black", "race_4l_Other_MultiRacial",
      "eth_hisp_Hispanic",
      "income_inc_1", "income_inc_2", "income_inc_3",
      "religion_rp_1", "religion_rp_2", "religion_rp_3", "religion_rp_4",
      "religion_rp_5", "religion_rp_6", "religion_rp_7", "religion_rp_8",
      "religion_rp_9", "religion_rp_10", "religion_rp_11", "religion_rp_12",
      "religion_rp_13", "religion_rp_14", "religion_rp_15", "religion_rp_16",
      "p_edu_Less_than_HS_Degree_GED_Equivalent",
      "p_edu_HS_Graduate_GED_Equivalent",
      "p_edu_Some_College_or_Associates_Degree",
      "p_edu_Masters_Degree",
      "p_edu_Professional_School_or_Doctoral_Degree",
      "sex_2l_Female",
      # Physical Health
      "rec_bin_Yes", "exp_sub_Yes", "tbi_injury_Yes",
      # Culture & Environment
      "mh_p_kbi__school_006_Grade_B", "mh_p_kbi__school_006_Grade_C",
      "mh_p_kbi__school_006_Grade_Fail",
      "det_susp_Yes",
      "se_services_Emotion_or_Learning_Support", "se_services_Gifted",
      "se_services_Other", "se_services_Combined_Services",
      # Neurocognitive
      "cct_Immediate"
    )
    
    # Select only the expected dummies that exist, create zeros for missing ones
    dummy_selected <- data.frame(matrix(0, nrow = nrow(factor_dummy), ncol = length(expected_dummies)))
    colnames(dummy_selected) <- expected_dummies
    
    # Fill in values for dummies that exist
    for (col in expected_dummies) {
      if (col %in% colnames(factor_dummy)) {
        dummy_selected[[col]] <- factor_dummy[[col]]
      } else {
        cat(sprintf("  ⚠ Creating zero column for missing dummy: %s\n", col))
      }
    }
    
    # Remove original factor variables
    imp_df_final <- imp_df_cs_final %>%
      select(-race_4l, -eth_hisp, -income, -religion, -p_edu, -sex_2l,
             -rec_bin, -exp_sub, -tbi_injury,
             -mh_p_kbi__school_006, -det_susp, -se_services, -cct)
    
    # Add dummy codes
    imp_df_final <- cbind(imp_df_final, dummy_selected)
    
    cat(sprintf("✓ Total variables after dummy coding: %d\n", ncol(imp_df_final)))
    cat(sprintf("✓ Total dummy variables: %d\n", length(expected_dummies)))
    
    imp_df_final
  }
  
  
  # Create MICE Object -----------------------------------------------------------
  create_mice_object <- function(imp_df_final) {
    
    obs_df <- imp_df_final %>% subset(.imp == 0) %>% select(-.imp, -.id)
    
    where_matrix <- matrix(TRUE, nrow = nrow(obs_df), ncol = ncol(obs_df))
    colnames(where_matrix) <- colnames(obs_df)
    
    mice_obj <- as.mids(imp_df_final, where = where_matrix, .imp = ".imp", .id = ".id")
    
    mice_obj
  }
  
  # Master Processing Loop -------------------------------------------------------
  # This processes all splits for all comparisons
  
  # Define split seeds
  split_seeds <- c(11897, 95039, 48378, 21773, 48339, 
                   88593, 87388, 12993, 52701, 34885)
  
  # Initialize storage
  results <- list()
  data.3 <- readRDS('output/study4/data_management/wBL_fulldata.rds')
  table_names_all <- readRDS('output/study4/data_management/table_names_wBL.rds')
  
  for (i in seq_along(pairs)) {
    
  }
  
  # Loop through each comparison
  for (comp_idx in seq_along(pairs)) {
    job::job({
      p <- pairs[[comp_idx]]
      level_order <- p[order(-table(data.3[data.3$ng4_class %in% p,]$ng4_class)[p])]
      level1 <- level_order[1]
      level2 <- level_order[2]
      
      cat(sprintf("\n\n"))
      cat(sprintf("================================================================================\n"))
      cat(sprintf("COMPARISON %d: LEVEL %d vs LEVEL %d\n", comp_idx, level1, level2))
      cat(sprintf("================================================================================\n"))
      
      # Loop through each split
      for (split_idx in 1:10) {
        
        tryCatch({
          
          # Process split
          split_results <- process_split_comparison(
            data = data.3 ,  # original dataset
            split_seed = split_seeds[split_idx],
            split_num = split_idx,
            comparison_num = comp_idx,
            level1 = level1,
            level2 = level2
          )
          
          # Center and scale
          split_results$trn_imp_cs <- center_scale_data(split_results$trn_imp)
          split_results$tst_imp_cs <- center_scale_data(split_results$tst_imp)
          
          # Create dummy variables
          split_results$trn_final <- create_dummy_variables(split_results$trn_imp_cs)
          split_results$tst_final <- create_dummy_variables(split_results$tst_imp_cs)
          
          # Create final MICE objects
          split_results$trn_mice_final <- create_mice_object(split_results$trn_final)
          split_results$tst_mice_final <- create_mice_object(split_results$tst_final)
          
          # Store results with naming convention
          comp_name <- sprintf("comp%d_split%d", comp_idx, split_idx)
          results[[comp_name]] <- split_results
          
          # Assign to global environment with meaningful names
          assign(sprintf("trn_%d_%dvs%d_obs", split_idx, level1, level2),
                 split_results$trn_obs, envir = .GlobalEnv)
          assign(sprintf("tst_%d_%dvs%d_obs", split_idx, level1, level2),
                 split_results$tst_obs, envir = .GlobalEnv)
          assign(sprintf("trn_%d_%dvs%d_final", split_idx, level1, level2),
                 split_results$trn_final, envir = .GlobalEnv)
          assign(sprintf("tst_%d_%dvs%d_final", split_idx, level1, level2),
                 split_results$tst_final, envir = .GlobalEnv)
          assign(sprintf("trn_%d_%dvs%d_mice", split_idx, level1, level2),
                 split_results$trn_mice_final, envir = .GlobalEnv)
          assign(sprintf("tst_%d_%dvs%d_mice", split_idx, level1, level2),
                 split_results$tst_mice_final, envir = .GlobalEnv)
          
          cat(sprintf("\n✓✓✓ Successfully completed Split %d, Comparison %d ✓✓✓\n",
                      split_idx, comp_idx))
          
        }, error = function(e) {
          cat(sprintf("\n✗✗✗ ERROR in Split %d, Comparison %d: %s ✗✗✗\n",
                      split_idx, comp_idx, e$message))
        })}
    })
  }
  
  # Save Results -----------------------------------------------------------------
  result_list <- c(
    paste0("trn_", 1:10, "_1vs2_mice"),
    paste0("trn_", 1:10, "_1vs3_mice"),
    paste0("trn_", 1:10, "_2vs3_mice"),
    paste0("trn_", 1:10, "_1vs4_mice"),
    paste0("trn_", 1:10, "_2vs4_mice"),
    paste0("trn_", 1:10, "_3vs4_mice")
  )
  while (any(!sapply(result_list, exists)) ||
         any(sapply(result_list, function(x) exists(x) && length(get(x)) == 0))) {
    message("Imputation incomplete. Waiting 30 minutes...")
    Sys.sleep(30 * 60)  # sleep for 30 minutes
    # refresh list here if needed, e.g.
    # my_list <- get_new_data()
  }
  
  cat("\n\n=== Saving Results ===\n")
  save(list = ls(pattern = "^(trn_|tst_)"), file = 'output/study4/MICE_prep/imputed_data_wBL.RData')
  cat("✓ Saved to mice_pairs.RData\n")
  
  
  
  
  email <- compose_email(
    body = "Imputation complete."
  )
  
  smtp_send(
    email,
    from = " email here",
    to   = "destination email here",
    subject = "Imputation step done",
    credentials = creds_key("gmail_smtp")
  )
  
  # Summary ----------------------------------------------------------------------
  cat("\n\n=== PROCESSING COMPLETE ===\n")
  cat(sprintf("Total comparisons processed: %d\n", length(pairs)))
  cat(sprintf("Total splits per comparison: 10\n"))
  cat(sprintf("Total datasets created: %d\n", length(results)))
  cat("\nObject naming convention:\n")
  cat("  - trn_<split>_<level1>vs<level2>_obs    : Training observed data\n")
  cat("  - tst_<split>_<level1>vs<level2>_obs    : Test observed data\n")
  cat("  - trn_<split>_<level1>vs<level2>_final  : Training imputed/processed\n")
  cat("  - tst_<split>_<level1>vs<level2>_final  : Test imputed/processed\n")
  cat("  - trn_<split>_<level1>vs<level2>_mice   : Training MICE object\n")
  cat("  - tst_<split>_<level1>vs<level2>_mice   : Test MICE object\n")
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(glmnet)
  library(rsample)
  library(purrr)
  library(pROC)
  library(tidyverse)
  library(job)
  library(ggsci)
  
  # #if necessary:
  #load("~/output/study4/MICE_prep/imputed_data_wMRI.RData")
  
  # Function: Run Elastic Net on MICE Data ---------------------------------------
  run_elastic_net_mice <- function(mice_obj, obs_df, alpha = 0.5, nfolds = 10) {
    
    cat("\n--- Running Elastic Net with MICE ---\n")
    
    # Extract number of imputations
    m <- mice_obj$m
    
    # Check sample size for CV folds
    n <- nrow(obs_df)
    min_class_size <- min(table(obs_df$ng4_class))
    
    # Adjust nfolds if necessary
    if (min_class_size < nfolds * 2) {
      nfolds_adjusted <- max(3, floor(min_class_size / 2))
      cat(sprintf("⚠ Adjusting nfolds from %d to %d due to small class size (%d)\n",
                  nfolds, nfolds_adjusted, min_class_size))
      nfolds <- nfolds_adjusted
    }
    
    # Further check: minimum 2 observations per fold per class
    if (floor(min_class_size / nfolds) < 2) {
      nfolds <- max(2, floor(min_class_size / 2))
      cat(sprintf("⚠ Further adjusting nfolds to %d to ensure min 2 obs/fold\n", nfolds))
    }
    
    # Storage for models and predictions
    models <- vector("list", m)
    lambdas <- numeric(m)
    max_aucs <- numeric(m)
    kept_columns <- vector("list", m)  # NEW: Track which columns are kept
    
    # Run elastic net on each imputation
    for (i in 1:m) {
      
      tryCatch({
        # Complete data for this imputation
        complete_data <- mice::complete(mice_obj, action = i)
        
        # Prepare X and Y
        y <- complete_data$ng4_class
        X <- complete_data %>% 
          select(-ng4_class, -age_baseline) %>%
          as.matrix()
        
        #Check for zero variance columns
        col_vars <- apply(X, 2, var)
        zero_var_cols <- which(col_vars == 0 | is.na(col_vars))
        
        if (length(zero_var_cols) > 0) {
          cat(sprintf("  Removing %d zero-variance columns\n", length(zero_var_cols)))
          X <- X[, -zero_var_cols, drop = FALSE]
        }
        
        # Store which columns were kept
        kept_columns[[i]] <- colnames(X)
        
        # Run cross-validated elastic net
        set.seed(123 + i)
        cv_model <- cv.glmnet(
          x = X,
          y = y,
          family = "binomial",
          alpha = alpha,
          nfolds = nfolds,
          type.measure = "auc",
          parallel = FALSE,  # Disable parallel to avoid nested parallelization issues
          standardize = FALSE  # Already standardized
        )
        
        models[[i]] <- cv_model
        lambdas[i] <- cv_model$lambda.min
        max_aucs[i] <- max(cv_model$cvm)
        
        cat(sprintf("  Imputation %d: Lambda = %.6f, Max AUC = %.4f\n",
                    i, cv_model$lambda.min, max(cv_model$cvm)))
        
      }, error = function(e) {
        cat(sprintf("  ✗ Error in imputation %d: %s\n", i, e$message))
        models[[i]] <- NULL
        lambdas[i] <- NA
        max_aucs[i] <- NA
      })
    }
    
    # Remove failed models
    valid_idx <- !sapply(models, is.null)
    
    if (sum(valid_idx) == 0) {
      stop("All imputations failed during elastic net fitting")
    }
    
    if (sum(valid_idx) < m) {
      cat(sprintf("⚠ Warning: %d out of %d imputations failed\n", m - sum(valid_idx), m))
      models <- models[valid_idx]
      lambdas <- lambdas[valid_idx]
      max_aucs <- max_aucs[valid_idx]
    }
    
    # Use median lambda across successful imputations
    lambda_final <- median(lambdas, na.rm = TRUE)
    cat(sprintf("\n✓ Final lambda (median): %.6f\n", lambda_final))
    cat(sprintf("✓ Mean training AUC: %.4f (SD: %.4f)\n", 
                mean(max_aucs, na.rm = TRUE), sd(max_aucs, na.rm = TRUE)))
    
    list(
      models = models,
      lambda_final = lambda_final,
      lambdas = lambdas,
      max_aucs = max_aucs,
      nfolds_used = nfolds,
      kept_columns = kept_columns  # NEW: Return column info
    )
  }
  
  # Function: Get Pooled Coefficients (Rubin's Rules) ---------------------------
  pool_coefficients <- function(enet_results, mice_obj, lambda = NULL) {
    
    if (is.null(lambda)) {
      lambda <- enet_results$lambda_final
    }
    
    m <- length(enet_results$models)
    
    # Extract coefficients from each imputation
    coef_list <- lapply(enet_results$models, function(model) {
      as.matrix(coef(model, s = lambda))
    })
    
    # Ensure all have same variables (elastic net may select different variables)
    all_vars <- unique(unlist(lapply(coef_list, rownames)))
    
    # Create matrix with all coefficients
    coef_matrix <- matrix(0, nrow = length(all_vars), ncol = m)
    rownames(coef_matrix) <- all_vars
    
    for (i in 1:m) {
      vars_i <- rownames(coef_list[[i]])
      coef_matrix[vars_i, i] <- coef_list[[i]][, 1]
    }
    
    # Pool using Rubin's rules (simple average for elastic net)
    pooled_coef <- rowMeans(coef_matrix)
    
    # Calculate stability (% of imputations where coefficient is non-zero)
    stability <- rowMeans(coef_matrix != 0) * 100
    
    # Create results dataframe
    results_df <- data.frame(
      variable = names(pooled_coef),
      coefficient = pooled_coef,
      stability_pct = stability,
      stringsAsFactors = FALSE
    ) %>%
      filter(coefficient != 0) %>%
      arrange(desc(abs(coefficient)))
    
    list(
      pooled_coef = pooled_coef,
      coef_matrix = coef_matrix,
      results_df = results_df
    )
  }
  
  # Function: Make Predictions on Test Data -------------------------------------
  predict_elastic_net_mice <- function(enet_results, pooled_coef, 
                                       test_mice_obj, test_obs_df) {
    
    m <- test_mice_obj$m
    lambda <- enet_results$lambda_final
    
    # Storage for predictions from each imputation
    pred_list <- vector("list", m)
    
    for (i in 1:m) {
      
      # Complete test data
      complete_test <- mice::complete(test_mice_obj, action = i)
      
      # Prepare X
      X_test <- complete_test %>%
        select(-ng4_class) %>%
        as.matrix()
      
      # Get columns that were kept during training
      kept_cols <- enet_results$kept_columns[[i]]
      
      # Find missing and extra columns
      missing_cols <- setdiff(kept_cols, colnames(X_test))
      extra_cols <- setdiff(colnames(X_test), kept_cols)
      
      # Handle missing columns
      if (length(missing_cols) > 0) {
        cat(sprintf("  ⚠ Test imputation %d missing %d columns, filling with zeros: %s\n",
                    i, length(missing_cols), 
                    paste(head(missing_cols, 3), collapse = ", ")))
        
        # Create matrix with missing columns (filled with 0)
        missing_matrix <- matrix(0, 
                                 nrow = nrow(X_test), 
                                 ncol = length(missing_cols))
        colnames(missing_matrix) <- missing_cols
        
        # Combine with existing test data
        X_test <- cbind(X_test, missing_matrix)
      }
      
      # Remove extra columns and reorder to match training
      X_test <- X_test[, kept_cols, drop = FALSE]
      
      # Verify dimensions match
      if (ncol(X_test) != length(kept_cols)) {
        stop(sprintf("Dimension mismatch after alignment: expected %d cols, got %d",
                     length(kept_cols), ncol(X_test)))
      }
      
      # Make predictions
      pred_list[[i]] <- predict(
        enet_results$models[[i]],
        newx = X_test,
        s = lambda,
        type = "response"
      )[, 1]
    }
    
    # Average predictions across imputations
    pred_matrix <- do.call(cbind, pred_list)
    pooled_predictions <- rowMeans(pred_matrix)
    
    # Get true outcomes from observed data
    y_true <- test_obs_df$ng4_class
    
    list(
      predictions = pooled_predictions,
      pred_matrix = pred_matrix,
      y_true = y_true
    )
  }
  
  # Function: Calculate Performance Metrics --------------------------------------
  calculate_metrics <- function(predictions, y_true) {
    
    # ROC curve and AUC
    roc_obj <- pROC::roc(y_true, predictions, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    # Find optimal threshold using Youden's index
    coords_obj <- pROC::coords(roc_obj, "best", best.method = "youden")
    threshold <- coords_obj$threshold
    
    # Predictions at optimal threshold
    pred_class <- ifelse(predictions >= threshold, 1, 0)
    
    # Confusion matrix
    cm <- table(Predicted = pred_class, Actual = y_true)
    
    # Calculate metrics
    accuracy <- sum(diag(cm)) / sum(cm)
    sensitivity <- cm[2, 2] / sum(cm[, 2])
    specificity <- cm[1, 1] / sum(cm[, 1])
    ppv <- cm[2, 2] / sum(cm[2, ])
    npv <- cm[1, 1] / sum(cm[1, ])
    
    list(
      auc = auc_val,
      threshold = threshold,
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      ppv = ppv,
      npv = npv,
      confusion_matrix = cm,
      roc_obj = roc_obj
    )
  }
  
  # Master Function: Complete Elastic Net Pipeline ------------------------------
  run_complete_pipeline <- function(split_num, comp_idx, level1, level2, 
                                    alpha = 0.5, nfolds = 10) {
    
    
    cat("\n")
    cat("========================================================================\n")
    cat(sprintf("ELASTIC NET: Split %d, Comparison %d (%d vs %d)\n", 
                split_num, comp_idx, level1, level2))
    cat("========================================================================\n")
    
    # Get object names
    trn_mice_name <- sprintf("trn_%d_%dvs%d_mice", split_num, level1, level2)
    tst_mice_name <- sprintf("tst_%d_%dvs%d_mice", split_num, level1, level2)
    trn_obs_name <- sprintf("trn_%d_%dvs%d_obs", split_num, level1, level2)
    tst_obs_name <- sprintf("tst_%d_%dvs%d_obs", split_num, level1, level2)
    
    # Get objects from global environment
    trn_mice <- get(trn_mice_name, envir = .GlobalEnv)
    tst_mice <- get(tst_mice_name, envir = .GlobalEnv)
    trn_obs <- get(trn_obs_name, envir = .GlobalEnv)
    tst_obs <- get(tst_obs_name, envir = .GlobalEnv)
    
    # 1. Train elastic net
    enet_results <- run_elastic_net_mice(trn_mice, trn_obs, alpha, nfolds)
    
    # 2. Pool coefficients
    pooled <- pool_coefficients(enet_results, trn_mice)
    
    cat("\n--- Selected Variables ---\n")
    print(head(pooled$results_df, 20))
    
    # 3. Predict on test set
    predictions <- predict_elastic_net_mice(enet_results, pooled$pooled_coef,
                                            tst_mice, tst_obs)
    
    # 4. Calculate metrics
    metrics <- calculate_metrics(predictions$predictions, predictions$y_true)
    
    cat("\n--- Test Set Performance ---\n")
    cat(sprintf("AUC: %.4f\n", metrics$auc))
    cat(sprintf("Accuracy: %.4f\n", metrics$accuracy))
    cat(sprintf("Sensitivity: %.4f\n", metrics$sensitivity))
    cat(sprintf("Specificity: %.4f\n", metrics$specificity))
    cat(sprintf("PPV: %.4f\n", metrics$ppv))
    cat(sprintf("NPV: %.4f\n", metrics$npv))
    
    # Return all results
    list(
      split = split_num,
      comparison = comp_idx,
      levels = c(level1, level2),
      enet_results = enet_results,
      pooled_coef = pooled,
      predictions = predictions,
      metrics = metrics
    )
  }
  
  # Run All Elastic Net Models ---------------------------------------------------
  cat("\n\n")
  cat("================================================================================\n")
  cat("RUNNING ELASTIC NET ON ALL COMPARISONS\n")
  cat("================================================================================\n")
  
  # Storage for all results
  enet_all_results <- list()
  
  # Loop through comparisons and splits
  for (comp_idx in seq_along(pairs)) {
    
    p <- pairs[[comp_idx]]
    level_order <- p[order(-table(data.3[data.3$ng4_class %in% p,]$ng4_class)[p])]
    level1 <- level_order[1]
    level2 <- level_order[2]
    
    job({
      for (split_idx in 1:10) {
        
        tryCatch({
          
          # Run complete pipeline
          result <- run_complete_pipeline(
            split_num = split_idx,
            comp_idx = comp_idx,
            level1 = level1,
            level2 = level2,
            alpha = 0.5,  # Elastic net (0 = ridge, 1 = lasso)
            nfolds = 10
          )
          
          # Store result
          result_name <- sprintf("comp%d_split%d", comp_idx, split_idx)
          enet_all_results[[result_name]] <- result
          
          # Assign to global environment
          assign(sprintf("enet_%d_%dvs%d", split_idx, level1, level2),
                 result, envir = .GlobalEnv)
          
        }, error = function(e) {
          cat(sprintf("\n✗ ERROR in Split %d, Comparison %d: %s\n",
                      split_idx, comp_idx, e$message))
        })
      }
      
      assign(sprintf('enet_results_comp%d', comp_idx), enet_all_results,
             envir = .GlobalEnv)}, title = sprintf("GLMNET_Comp%d", comp_idx),
      import = c(ls(pattern = paste0("^(trn_|tst_)\\d+_", level1, "vs", level2, "_(mice|obs)$")),
                 'comp_idx', 'level1', 'level2', 'calculate_metrics', 'pool_coefficients',
                 'predict_elastic_net_mice', 'run_complete_pipeline', 'run_elastic_net_mice',
                 'enet_all_results', 'pairs'))
  }
  
  # Save all elastic net results
  result_list <- c('enet_results_comp1', 'enet_results_comp2', 'enet_results_comp3',
                   'enet_results_comp4', 'enet_results_comp5', 'enet_results_comp6')
  while (any(!sapply(result_list, exists)) ||
         any(sapply(result_list, function(x) exists(x) && length(get(x)) == 0))) {
    message("Regressions incomplete. Waiting 30 minutes...")
    Sys.sleep(30 * 60)  # sleep for 30 minutes
    # refresh  list here if needed, e.g.
    # my_list <- get_new_data() 
  }
  enet_all_results = c(enet_results_comp1, enet_results_comp2, enet_results_comp3,
                       enet_results_comp4, enet_results_comp5, enet_results_comp6)
  save(list = c(ls(pattern = "^(enet_)"),'table_names_all', 'data.3'), file = "output/study4/en_results/claude_elastic_net_results_wBL.RData")
  cat("\n✓ All elastic net results saved to elastic_net_results.RData\n")
  
  email <- compose_email(
    body = "eNet complete."
  )
  
  smtp_send(
    email,
    from = " email here",
    to   = "destination email here",
    subject = "Elastic Net step 1 done",
    credentials = creds_key("gmail_smtp")
  )
  
  # Summary Statistics -----------------------------------------------------------
  cat("\n\n=== ELASTIC NET SUMMARY ===\n")
  
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
  
  # ==============================================================================
  # EXTRACT AND SUMMARIZE ELASTIC NET RESULTS FROM ENVIRONMENT
  # ==============================================================================
  # Extracts results from enet_X_XvsX objects in the global environment
  # and creates comprehensive summaries and visualizations
  # ==============================================================================
  
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
  
  # ============================================================================
  # 1. PERFORMANCE SUMMARY BY SPLIT
  # ============================================================================
  
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
          comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
          comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
          comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
          comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
          comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
          comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
        ))
    }
  }
  
  performance_df <- bind_rows(performance_by_split) %>%
    arrange(comparison, split)
  
  cat("\n--- Performance by Split ---\n")
  print(head(performance_df, 20))
  
  write.csv(performance_df,
            "output/study4/en_results/performance_by_split_wMRI.csv",
            row.names = FALSE)
  
  # ============================================================================
  # 2. AGGREGATED PERFORMANCE SUMMARY (MEAN ACROSS SPLITS)
  # ============================================================================
  
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
  
  write.csv(performance_summary,
            "output/study4/en_results/performance_summary_aggregated_wMRI.csv",
            row.names = FALSE)
  
  # ============================================================================
  # 3. TOP PREDICTORS BY COMPARISON (ACROSS ALL SPLITS)
  # ============================================================================
  
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
      comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
      comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
      comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
      comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
      comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
      comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
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
    arrange(comparison, desc(n_splits_selected), desc(abs_mean_coef))
  
  write.csv(variable_importance,
            "output/study4/en_results/variable_importance_by_comparison_wMRI.csv",
            row.names = FALSE)
  
  # Top predictors (selected in >=7 splits with high coefficients)
  top_predictors <- variable_importance %>%
    filter(n_splits_selected >= 7) %>%
    filter(variable != '(Intercept)') %>%
    group_by(comparison) %>%
    slice_max(order_by = abs_mean_coef, n = 20) %>%
    ungroup() %>%
    left_join(table_names_all[,c('variable', 'table_name', 'domain_name')], by = 'variable') %>%
    select(comparison, domain_name, table_name, variable, everything()) %>%
    arrange(comparison, rev(abs_mean_coef))
  
  cat("\n--- Top Predictors (selected in >=7 splits) ---\n")
  print(head(top_predictors, 30))
  
  write.csv(top_predictors,
            "output/study4/en_results/top_predictors_consistent_wMRI.csv",
            row.names = FALSE)
  
  
  # ============================================================================
  # 4. COEFFICIENT STABILITY HEATMAP DATA
  # ============================================================================
  
  # Create wide format: variables x comparisons, showing mean coefficient
  coef_wide <- variable_importance %>%
    filter(n_splits_selected >= 5) %>%  # Only variables selected in >=5 splits
    select(comparison, variable, mean_coefficient) %>%
    pivot_wider(
      names_from = comparison,
      values_from = mean_coefficient,
      values_fill = 0
    )
  
  write.csv(coef_wide,
            "output/study4/en_results/coefficients_wide_stable_wMRI.csv",
            row.names = FALSE)
  
  cat(sprintf("\nStable variables (selected in >=5 splits): %d\n", nrow(coef_wide)))
  
  # ============================================================================
  # 5. ROC CURVES VISUALIZATION
  # ============================================================================
  
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
            comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
            comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
            comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
            comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
            comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
            comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
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
    theme_minimal() +
    theme(legend.position = 'none')
  
  ggsave("output/study4/en_results/roc_curves_by_comparison_all_splits_wMRI.png",
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
  
  p_roc_mean <- ggplot(mean_roc_curves, 
                       aes(x = fpr_grid, y = mean_sensitivity, color = comparison_auc)) +
    geom_line(linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    labs(
      title = "Mean ROC Curves Across Splits",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Comparison"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 9)
    )
  
  ggsave("output/study4/en_results/roc_curves_mean_by_comparison_wMRI.png",
         p_roc_mean, width = 10, height = 8, dpi = 300)
  
  cat("ROC curves saved\n")
  
  # ============================================================================
  # 6. PERFORMANCE BOXPLOTS
  # ============================================================================
  
  cat("\n--- Creating Performance Boxplots ---\n")
  
  p_auc <- ggplot(performance_df, aes(x = comparison, y = auc, fill = comparison)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      title = "AUC Distribution Across Splits",
      x = "Comparison",
      y = "AUC"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(0.5, 1))
  
  ggsave("output/study4/en_results/auc_boxplot_by_comparison.png",
         p_auc, width = 10, height = 6, dpi = 300)
  
  p_accuracy <- ggplot(performance_df, aes(x = comparison, y = accuracy, fill = comparison)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      title = "Accuracy Distribution Across Splits",
      x = "Comparison",
      y = "Accuracy"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(0.5, 1))
  
  ggsave("output/study4/en_results/accuracy_boxplot_by_comparison_wMRI.png",
         p_accuracy, width = 10, height = 6, dpi = 300)
  
  cat("Performance boxplots saved\n")
  
  # ============================================================================
  # 7. VARIABLE SELECTION FREQUENCY
  # ============================================================================
  
  cat("\n--- Analyzing Variable Selection Frequency ---\n")
  
  var_selection_freq <- all_coefs_df %>%
    group_by(comparison, variable) %>%
    summarise(
      times_selected = n(),
      .groups = "drop"
    ) %>%
    arrange(comparison, desc(times_selected))
  
  # Variables selected in all 10 splits
  always_selected <- var_selection_freq %>%
    filter(times_selected == 10)
  
  cat(sprintf("Variables selected in all 10 splits: %d\n", nrow(always_selected)))
  
  write.csv(always_selected,
            "output/study4/en_results/variables_always_selected_wMRI.csv",
            row.names = FALSE)
  
  # Histogram of selection frequency
  p_freq <- ggplot(var_selection_freq %>%
                     mutate(comparison = case_when(
                       comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
                       comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
                       comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
                       comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
                       comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
                       comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
                     )), 
                   aes(x = times_selected, fill = comparison)) +
    geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
    facet_wrap(~ comparison, ncol = 2) +
    labs(
      title = "Variable Selection Frequency Across Splits",
      x = "Number of Splits Variable Was Selected",
      y = "Count of Variables"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave("output/study4/en_results/variable_selection_frequency_wMRI.png",
         p_freq, width = 10, height = 8, dpi = 300)
  
  cat("Variable selection frequency saved\n")
  
  # ============================================================================
  # 8. CONFUSION MATRIX SUMMARY
  # ============================================================================
  
  cat("\n--- Summarizing Confusion Matrices ---\n")
  
  confusion_summary <- list()
  
  for (comp_name in names(all_results)) {
    comp_results <- all_results[[comp_name]]
    
    for (split_name in names(comp_results)) {
      split_result <- comp_results[[split_name]]
      
      cm <- split_result$confusion_matrix
      
      confusion_summary[[paste(comp_name, split_name, sep = "_")]] <- data.frame(
        comparison = split_result$comparison,
        split = split_result$split,
        true_neg = cm[1, 1],
        false_pos = cm[2, 1],
        false_neg = cm[1, 2],
        true_pos = cm[2, 2],
        stringsAsFactors = FALSE
      ) %>%
        mutate(comparison = case_when(
          comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
          comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
          comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
          comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
          comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
          comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
        ))
    }
  }
  
  confusion_df <- bind_rows(confusion_summary)
  
  write.csv(confusion_df,
            "output/study4/en_results/confusion_matrices_all_splits_wMRI.csv",
            row.names = FALSE)
  
  # ============================================================================
  # FINAL SUMMARY
  # ============================================================================
  
  cat("\n\n")
  cat("================================================================================\n")
  cat("=== ELASTIC NET ANALYSIS COMPLETE ===\n")
  cat("================================================================================\n\n")
  
  cat("Files saved to output/study4/en_results_claude/:\n\n")
  cat("Performance Metrics:\n")
  cat("  - performance_by_split.csv\n")
  cat("  - performance_summary_aggregated.csv\n")
  cat("  - confusion_matrices_all_splits.csv\n\n")
  
  cat("Variable/Coefficient Information:\n")
  cat("  - variable_importance_by_comparison.csv\n")
  cat("  - top_predictors_consistent.csv\n")
  cat("  - coefficients_wide_stable.csv\n")
  cat("  - variables_always_selected.csv\n\n")
  
  cat("Visualizations:\n")
  cat("  - roc_curves_by_comparison_all_splits.png\n")
  cat("  - roc_curves_mean_by_comparison.png\n")
  cat("  - auc_boxplot_by_comparison.png\n")
  cat("  - accuracy_boxplot_by_comparison.png\n")
  cat("  - variable_selection_frequency.png\n\n")
  
  cat("Summary Statistics:\n")
  cat(sprintf("  Total comparisons analyzed: %d\n", length(all_results)))
  cat(sprintf("  Total splits analyzed: %d\n", nrow(performance_df)))
  cat(sprintf("  Mean AUC (across all): %.3f (SD: %.3f)\n", 
              mean(performance_df$auc, na.rm = TRUE),
              sd(performance_df$auc, na.rm = TRUE)))
  cat(sprintf("  Mean Accuracy (across all): %.3f (SD: %.3f)\n",
              mean(performance_df$accuracy, na.rm = TRUE),
              sd(performance_df$accuracy, na.rm = TRUE)))
  cat(sprintf("  Total unique variables selected: %d\n",
              length(unique(all_coefs_df$variable))))
  cat(sprintf("  Variables selected in all 10 splits: %d\n",
              nrow(always_selected)))
  
  cat("\n")
  cat("Best performing comparison (by mean AUC):\n")
  best_comp <- performance_summary %>% slice_max(mean_auc, n = 1)
  cat(sprintf("  %s: Mean AUC = %.3f (SD: %.3f)\n",
              best_comp$comparison, best_comp$mean_auc, best_comp$sd_auc))
  
  email <- compose_email(
    body = "eNet complete."
  )
  
  smtp_send(
    email,
    from = " email here",
    to   = "destination email here",
    subject = "Four Class Model Complete",
    credentials = creds_key("gmail_smtp")
  )
  
  
  
  rm(list = ls())
  
  
  
  
}, error = function(e) {
  # Send email on error
  smtp_send(
    compose_email(
      body = paste("Error occurred:\n\n", e$message, "\n\nTraceback:\n",
                   paste(capture.output(traceback()), collapse = "\n"))
    ),
    from = "your email here",
    to   = "destination email here",
    subject = "Error Occurred",
    credentials = creds_key("gmail_smtp")
  )
  
  stop(e)  # Re-throw the error if you want it to still fail
})

rm(list = ls())
system("powercfg -change standby-timeout-ac 15")  # Resets sleep timeout
system('shutdown -l')

