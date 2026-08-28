# Configure Parallelization ----------------------------------------------------

## Detect core count
nCores <- min(parallel::detectCores())

## Used by parallel::mclapply() as default
options(mc.cores = nCores)

## Used by doParallel as default
options(cores = nCores)

## Register doParallel as the parallel backend with foreach
doParallel::registerDoParallel(cores = nCores-1)

## Report multicore use
cat("### Using", foreach::getDoParWorkers(), "cores\n")
cat("### Using", foreach::getDoParName(), "as backend\n")

# Assign IDs to Folds -----------------------------------------------------

nfolds <- 10
set.seed(48109)
foldid <- sample(rep(seq(nfolds), length.out = nrow(trn_1_mice_m1[[1]]))) 
foldid # n = 5121 for training dataset
table(foldid) 
# n = 632 in fold 1-4, n = 631 in folds 5 - 10 (aligns w/final sample size)

# Elastic Net: Model 1 -----------------------------------------------------------------

#------------------------------------------------------------------------------#
# (A) run elastic net on training dataset
#------------------------------------------------------------------------------#

### split (1)
job::job(m1_split1 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_1_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_1_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (1): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_1_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_1_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_1_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})

### split (2)
job::job(m1_split2 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_2_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_2_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (2): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_2_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_2_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_2_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})

### split (3)
job::job(m1_split3 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_3_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_3_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (3): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_3_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_3_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_3_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})

### split (4)
job::job(m1_split4 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_4_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_4_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (4): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_4_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_4_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_4_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})

### split (5)
job::job(m1_split5 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_5_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_5_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (5): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_5_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_5_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_5_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})

### split (6)
job::job(m1_split6 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_6_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_6_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (6): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_6_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_6_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_6_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})

### split (7)
job::job(m1_split7 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_7_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_7_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (7): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_7_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_7_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_7_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})

### split (8)
job::job(m1_split8 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_8_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_8_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (8): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_8_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_8_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_8_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})

### split (9)
job::job(m1_split9 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_9_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_9_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (9): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_9_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_9_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_9_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})

### split (10)
job::job(m1_split10 ={
  # matrix of imputed values 
  dfs_trn_m1 <- lapply(1:5, function(i) complete(trn_10_mice_m1, action = i)) 
  
  # generate list of imputed design matrices and imputed responses
  x_trn_m1 <- list()
  y_trn_m1 <- list()
  for (i in 1:5) {
    x_trn_m1[[i]] <- as.matrix(dfs_trn_m1[[i]][,(2:length(trn_10_obs_df_m1))]) # drop column 1 DV
    y_trn_m1[[i]] <- dfs_trn_m1[[i]]$DV # save DV column separately
  }
  
  cat(dim(x_trn_m1[[1]]), 'split (10): IDs by predictors  \n')
  
  # calculate observational weights
  weights_trn_m1 <- 1 - rowMeans(is.na(trn_10_obs_df_m1))
  pf_trn_m1 <- rep(1, length(trn_10_obs_df_m1)-1) # repetitions = total number of predictors
  adWeight_trn_m1 <- rep(1, length(trn_10_obs_df_m1)-1) # repetitions = total number of predictors
  
  # parallelized execution of cross-validated stacked adaptive elastic net reg.
  
  kind <- RNGkind()
  kind
  
  ## Register doParallel as the parallel backend with foreach
  doParallel::registerDoParallel(cores = nCores-1)
  
  fit_m1 <- foreach(i = seq_len(1), .combine = ibind) %dorng% {
    set.seed(48109, kind = kind[1]) # set to default 'Mersenne-Twister'
    fit_m1 <- miselect::cv.saenet(
      x_trn_m1, y_trn_m1, pf_trn_m1, adWeight_trn_m1, weights_trn_m1, 
      family = "binomial",
      alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      nfolds = nfolds, foldid = foldid)
    fit_m1
  }
  
})


#------------------------------------------------------------------------------#
# (B) extract coefficients retained and associated output
#------------------------------------------------------------------------------#

# objects retained from each split: 
# - model fit details, AUC, sensitivity and specificity 

### split (1)
job::job(m1_split1_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split1$fit_m1$alpha.1se
  m1_alpha.min <- m1_split1$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split1$fit_m1$lambda.1se
  m1_lambda.min <- m1_split1$fit_m1$lambda.min
  
  m1_split1_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split1_fit) <- c('split 1')
  m1_split1_fit <- m1_split1_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split1$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split1_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split1_lambda_range) <- c('split 1')
  
  m1_df <- m1_split1$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_1/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split1$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split1$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_1/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split1$fit_m1,
                    lambda = m1_split1$fit_m1$lambda.1se, 
                    alpha = m1_split1$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split1$fit_m1, 
                    lambda = m1_split1$fit_m1$lambda.min, 
                    alpha = m1_split1$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (1)  \n') 
  
  write.csv(m1_cf,'output/study2/elastic_net/split_1/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_1_temp <- trn_1_imp_df_m1 %>% 
      subset(.imp == i)
    trn_1_temp[, names(trn_1_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_1_temp <- tst_1_imp_df_m1 %>% 
      subset(.imp == i)
    tst_1_temp[, names(tst_1_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:(length(m1_cf_min_nz)-1))]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:(length(m1_cf_min_nz)-1))]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (1)  \n') 
  
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #----------------------------------------------------------------------------#  
  #                       Calculate probabilities
  #----------------------------------------------------------------------------# 
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities across 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_1/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_1/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  
  # double check each probability across the imputations is different 
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC, 
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_1_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_1_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'training (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'training (across 5 imputations): median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'test (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'test (across 5 imputations): median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  # - training dataset
  m1_roc_trn <- roc(trn_1_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split1 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split1) <- c('training: split 1')
  
  # - test dataset
  m1_roc_tst <- roc(tst_1_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split1 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split1) <- c('test: split 1')
  
  # (2) confusion matrix
  # - threshold > 0.50
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <- m1_prob_avg_tst %>% 
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_1/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_1/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (3) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('training (split 1): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('test (split 1): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('training (split 1): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('test (split 1): best threshold')
  ss_tst_best
  
  ss_m1_split1 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})

### split (2)
job::job(m1_split2_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split2$fit_m1$alpha.1se
  m1_alpha.min <- m1_split2$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split2$fit_m1$lambda.1se
  m1_lambda.min <- m1_split2$fit_m1$lambda.min
  
  m1_split2_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split2_fit) <- c('split 2')
  m1_split2_fit <- m1_split2_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split2$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split2_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split2_lambda_range) <- c('split 2')
  
  m1_df <- m1_split2$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_2/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split2$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split2$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_2/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split2$fit_m1,
                    lambda = m1_split2$fit_m1$lambda.1se, 
                    alpha = m1_split2$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split2$fit_m1, 
                    lambda = m1_split2$fit_m1$lambda.min, 
                    alpha = m1_split2$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (2)  \n') 
  
  write.csv(m1_cf,'output/study2/elastic_net/split_2/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_2_temp <- trn_2_imp_df_m1 %>% 
      subset(.imp == i)
    trn_2_temp[, names(trn_2_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_2_temp <- tst_2_imp_df_m1 %>% 
      subset(.imp == i)
    tst_2_temp[, names(tst_2_temp) %in% names(m1_cf_min_nz)]    
    
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (2)  \n') 
  
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #------------------------------------------------------------------------------#  
  #                 Calculate AUC for training and test dataset
  #------------------------------------------------------------------------------#  
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities aross 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_2/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_2/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  # double check each probability across the imputations is different 
  
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC,
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_2_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_2_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'training (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'training (across 5 imputations): median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'test (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'test (across 5 imputations): median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  # - training dataset
  m1_roc_trn <- roc(trn_2_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split2 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split2) <- c('training: split 2')
  
  # - test dataset
  m1_roc_tst <- roc(tst_2_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split2 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split2) <- c('test: split 2')
  
  # (2) confusion matrix
  # - threshold > 0.50
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <-  m1_prob_avg_tst %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_tst = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_2/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_2/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (3) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('training (split 2): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('test (split 2): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('training (split 2): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('test (split 2): best threshold')
  ss_tst_best
  
  ss_m1_split2 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})

### split (3)
job::job(m1_split3_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split3$fit_m1$alpha.1se
  m1_alpha.min <- m1_split3$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split3$fit_m1$lambda.1se
  m1_lambda.min <- m1_split3$fit_m1$lambda.min
  
  m1_split3_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split3_fit) <- c('split 3')
  m1_split3_fit <- m1_split3_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split3$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split3_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split3_lambda_range) <- c('split 3')
  
  m1_df <- m1_split3$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_3/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split3$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split3$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_3/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split3$fit_m1,
                    lambda = m1_split3$fit_m1$lambda.1se, 
                    alpha = m1_split3$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split3$fit_m1, 
                    lambda = m1_split3$fit_m1$lambda.min, 
                    alpha = m1_split3$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients 
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (3)  \n') 
  
  write.csv(m1_cf,'output/study2/elastic_net/split_3/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_3_temp <- trn_3_imp_df_m1 %>% 
      subset(.imp == i)
    trn_3_temp[, names(trn_3_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_3_temp <- tst_3_imp_df_m1 %>% 
      subset(.imp == i)
    tst_3_temp[, names(tst_3_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (3)  \n') 
  
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #------------------------------------------------------------------------------#  
  #                 Calculate AUC for training and test dataset
  #------------------------------------------------------------------------------#  
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities aross 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_3/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_3/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  # double check each probability across the imputations is different 
  
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC, 
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_3_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_3_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'training (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'training (across 5 imputations): median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'test (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'test (across 5 imputations): median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  # - training dataset
  m1_roc_trn <- roc(trn_3_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split3 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split3) <- c('training: split 3')
  
  # - test dataset
  m1_roc_tst <- roc(tst_3_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split3 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split3) <- c('test: split 3')
  
  # (2) confusion matrix
  # - threshold > 0.50
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <- m1_prob_avg_tst %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_tst = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_3/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_3/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (3) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('training (split 3): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('test (split 3): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('training (split 3): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('test (split 3): best threshold')
  ss_tst_best
  
  ss_m1_split3 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})

### split (4)
job::job(m1_split4_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split4$fit_m1$alpha.1se
  m1_alpha.min <- m1_split4$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split4$fit_m1$lambda.1se
  m1_lambda.min <- m1_split4$fit_m1$lambda.min
  
  m1_split4_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split4_fit) <- c('split 4')
  m1_split4_fit <- m1_split4_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split4$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split4_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split4_lambda_range) <- c('split 4')
  
  m1_df <- m1_split4$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_4/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split4$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split4$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_4/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split4$fit_m1,
                    lambda = m1_split4$fit_m1$lambda.1se, 
                    alpha = m1_split4$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split4$fit_m1, 
                    lambda = m1_split4$fit_m1$lambda.min, 
                    alpha = m1_split4$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients 
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (4)  \n') 
  
  write.csv(m1_cf,'output/study2/elastic_net/split_4/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_4_temp <- trn_4_imp_df_m1 %>% 
      subset(.imp == i)
    trn_4_temp[, names(trn_4_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_4_temp <- tst_4_imp_df_m1 %>% 
      subset(.imp == i)
    tst_4_temp[, names(tst_4_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (4)  \n') 
  
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #------------------------------------------------------------------------------#  
  #                 Calculate AUC for training and test dataset
  #------------------------------------------------------------------------------#  
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities aross 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_4/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_4/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  # double check each probability across the imputations is different 
  
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC, 
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_4_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_4_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'training (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'training (across 5 imputations): median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'test (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'test (across 5 imputations): median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  # - training dataset
  m1_roc_trn <- roc(trn_4_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split4 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split4) <- c('training: split 4')
  
  # - test dataset
  m1_roc_tst <- roc(tst_4_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split4 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split4) <- c('test: split 4')
  
  # (2) confusion matrix
  # - threshold > 0.50
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <- m1_prob_avg_tst %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_tst = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_4/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_4/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (2) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('training (split 4): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('test (split 4): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('training (split 4): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('test (split 4): best threshold')
  ss_tst_best
  
  ss_m1_split4 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})

### split (5)
job::job(m1_split5_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split5$fit_m1$alpha.1se
  m1_alpha.min <- m1_split5$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split5$fit_m1$lambda.1se
  m1_lambda.min <- m1_split5$fit_m1$lambda.min
  
  m1_split5_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split5_fit) <- c('split 5')
  m1_split5_fit <- m1_split5_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split5$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split5_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split5_lambda_range) <- c('split 5')
  
  m1_df <- m1_split5$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_5/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split5$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split5$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_5/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split5$fit_m1,
                    lambda = m1_split5$fit_m1$lambda.1se, 
                    alpha = m1_split5$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split5$fit_m1, 
                    lambda = m1_split5$fit_m1$lambda.min, 
                    alpha = m1_split5$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients 
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (5)  \n') 
  
  write.csv(m1_cf,'output/study2/elastic_net/split_5/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_5_temp <- trn_5_imp_df_m1 %>% 
      subset(.imp == i)
    trn_5_temp[, names(trn_5_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_5_temp <- tst_5_imp_df_m1 %>% 
      subset(.imp == i)
    tst_5_temp[, names(tst_5_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (5)  \n') 
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #------------------------------------------------------------------------------#  
  #                 Calculate AUC for training and test dataset
  #------------------------------------------------------------------------------#  
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities aross 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_5/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_5/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  # double check each probability across the imputations is different 
  
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC, 
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_5_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_5_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'training (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'training (across 5 imputations): median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'test (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'test (across 5 imputations): median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  ## training dataset
  m1_roc_trn <- roc(trn_5_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split5 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split5) <- c('training: split 5')
  
  ## test dataset
  m1_roc_tst <- roc(tst_5_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split5 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split5) <- c('test: split 5')
  
  # (2) confusion matrix
  # - threshold > 0.50
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>%  
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <- m1_prob_avg_tst %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_tst = '.') %>%  
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_5/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_5/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (3) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('training (split 5): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('test (split 5): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('training (split 5): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('test (split 5): best threshold')
  ss_tst_best
  
  ss_m1_split5 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})

### split (6)
job::job(m1_split6_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split6$fit_m1$alpha.1se
  m1_alpha.min <- m1_split6$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split6$fit_m1$lambda.1se
  m1_lambda.min <- m1_split6$fit_m1$lambda.min
  
  m1_split6_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split6_fit) <- c('split 6')
  m1_split6_fit <- m1_split6_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split6$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split6_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split6_lambda_range) <- c('split 6')
  
  m1_df <- m1_split6$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_6/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split6$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split6$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_6/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split6$fit_m1,
                    lambda = m1_split6$fit_m1$lambda.1se, 
                    alpha = m1_split6$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split6$fit_m1, 
                    lambda = m1_split6$fit_m1$lambda.min, 
                    alpha = m1_split6$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients 
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (6)  \n') 
  
  write.csv(m1_cf,'output/study2/elastic_net/split_6/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  # generate new list of imputed design matrices and imputed responses for 
  # variables retained
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_6_temp <- trn_6_imp_df_m1 %>% 
      subset(.imp == i)
    trn_6_temp[, names(trn_6_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_6_temp <- tst_6_imp_df_m1 %>% 
      subset(.imp == i)
    tst_6_temp[, names(tst_6_temp) %in% names(m1_cf_min_nz)]  
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (6)  \n') 
  
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #------------------------------------------------------------------------------#  
  #                 Calculate AUC for training and test dataset
  #------------------------------------------------------------------------------#  
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities aross 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_6/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_6/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  # double check each probability across the imputations is different 
  
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC, 
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_6_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_6_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'training (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'training (across 5 imputations): median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'test (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'test (across 5 imputations): median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  # - training dataset
  m1_roc_trn <- roc(trn_6_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split6 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split6) <- c('training: split 6')
  
  # - test dataset
  m1_roc_tst <- roc(tst_6_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split6 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split6) <- c('test: split 6')
  
  # (2) confusion matrix
  # - threshold > 0.50 
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>%  
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <- m1_prob_avg_tst %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_tst = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_6/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_6/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (3) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('training (split 6): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('test (split 6): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('training (split 6): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('test (split 6): best threshold')
  ss_tst_best
  
  ss_m1_split6 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})

### split (7)
job::job(m1_split7_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split7$fit_m1$alpha.1se
  m1_alpha.min <- m1_split7$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split7$fit_m1$lambda.1se
  m1_lambda.min <- m1_split7$fit_m1$lambda.min
  
  m1_split7_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split7_fit) <- c('Split 7')
  m1_split7_fit <- m1_split7_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split7$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split7_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split7_lambda_range) <- c('Split 7')
  
  m1_df <- m1_split7$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_7/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split7$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split7$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_7/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split7$fit_m1,
                    lambda = m1_split7$fit_m1$lambda.1se, 
                    alpha = m1_split7$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split7$fit_m1, 
                    lambda = m1_split7$fit_m1$lambda.min, 
                    alpha = m1_split7$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients 
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (7)  \n') 
  
  write.csv(m1_cf,'output/study2/elastic_net/split_7/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_7_temp <- trn_7_imp_df_m1 %>% 
      subset(.imp == i)
    trn_7_temp[, names(trn_7_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_7_temp <- tst_7_imp_df_m1 %>% 
      subset(.imp == i)
    tst_7_temp[, names(tst_7_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (7)  \n') 
  
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #------------------------------------------------------------------------------#  
  #                 Calculate AUC for training and test dataset
  #------------------------------------------------------------------------------#  
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities aross 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_7/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_7/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  # double check each probability across the imputations is different 
  
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC, 
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_7_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_7_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'Training (across 5 imputations): Median of controls < median of cases')
  } else {
    print(
      'Training (across 5 imputations): Median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'Test (across 5 imputations): Median of controls < median of cases')
  } else {
    print(
      'Test (across 5 imputations): Median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  # - training dataset
  m1_roc_trn <- roc(trn_7_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split7 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split7) <- c('Training: split 7')
  
  # - test dataset
  m1_roc_tst <- roc(tst_7_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split7 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split7) <- c('Test: split 7')
  
  # (2) confusion matrix
  # - threshold > 0.50 
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <- m1_prob_avg_tst %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_tst = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_7/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_7/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (3) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('Training (split 7): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('Test (split 7): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('Training (split 7): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('Test (split 7): best threshold method')
  ss_tst_best
  
  ss_m1_split7 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})

### split (8)
job::job(m1_split8_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split8$fit_m1$alpha.1se
  m1_alpha.min <- m1_split8$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split8$fit_m1$lambda.1se
  m1_lambda.min <- m1_split8$fit_m1$lambda.min
  
  m1_split8_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split8_fit) <- c('split 8')
  m1_split8_fit <- m1_split8_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split8$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split8_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split8_lambda_range) <- c('split 8')
  
  m1_df <- m1_split8$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_8/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split8$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split8$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_8/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split8$fit_m1,
                    lambda = m1_split8$fit_m1$lambda.1se, 
                    alpha = m1_split8$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split8$fit_m1, 
                    lambda = m1_split8$fit_m1$lambda.min, 
                    alpha = m1_split8$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients 
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (8)  \n')  
  
  write.csv(m1_cf,'output/study2/elastic_net/split_8/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_8_temp <- trn_8_imp_df_m1 %>% 
      subset(.imp == i)
    trn_8_temp[, names(trn_8_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_8_temp <- tst_8_imp_df_m1 %>% 
      subset(.imp == i)
    tst_8_temp[, names(tst_8_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (8)  \n') 
  
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #------------------------------------------------------------------------------#  
  #                 Calculate AUC for training and test dataset
  #------------------------------------------------------------------------------#  
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities aross 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_8/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_8/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  # double check each probability across the imputations is different 
  
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC, 
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_8_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_8_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'training (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'training (across 5 imputations): median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'test (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'test (across 5 imputations): median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  # - training dataset
  m1_roc_trn <- roc(trn_8_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split8 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split8) <- c('training: split 8')
  
  # - test dataset
  m1_roc_tst <- roc(tst_8_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split8 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split8) <- c('test: split 8')
  
  # (2) confusion matrix
  # - threshold > 0.50 
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>%  
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <- m1_prob_avg_tst %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_tst = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_8/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_8/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (3) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('training (split 8): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('test (split 8): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('training (split 8): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('test (split 8): best threshold')
  ss_tst_best
  
  ss_m1_split8 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})

### split (9)
job::job(m1_split9_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split9$fit_m1$alpha.1se
  m1_alpha.min <- m1_split9$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split9$fit_m1$lambda.1se
  m1_lambda.min <- m1_split9$fit_m1$lambda.min
  
  m1_split9_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split9_fit) <- c('split 9')
  m1_split9_fit <- m1_split9_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split9$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split9_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split9_lambda_range) <- c('split 9')
  
  m1_df <- m1_split9$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_9/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split9$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split9$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_9/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split9$fit_m1,
                    lambda = m1_split9$fit_m1$lambda.1se, 
                    alpha = m1_split9$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split9$fit_m1, 
                    lambda = m1_split9$fit_m1$lambda.min, 
                    alpha = m1_split9$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients 
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (9)  \n') 
  
  write.csv(m1_cf,'output/study2/elastic_net/split_9/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_9_temp <- trn_9_imp_df_m1 %>% 
      subset(.imp == i)
    trn_9_temp[, names(trn_9_temp) %in% names(m1_cf_min_nz)]  
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_9_temp <- tst_9_imp_df_m1 %>% 
      subset(.imp == i)
    tst_9_temp[, names(tst_9_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (9)  \n') 
  
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #------------------------------------------------------------------------------#  
  #                 Calculate AUC for training and test dataset
  #------------------------------------------------------------------------------#  
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities aross 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_9/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_9/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  # double check each probability across the imputations is different 
  
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC, 
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_9_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_9_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'training (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'training (across 5 imputations): median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'test (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'test (across 5 imputations): median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  # - training dataset
  m1_roc_trn <- roc(trn_9_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split9 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split9) <- c('training: split 9')
  
  # - test dataset
  m1_roc_tst <- roc(tst_9_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split9 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split9) <- c('test: split 9')
  
  # (2) confusion matrix
  # - threshold > 0.50 
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>%  
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <- m1_prob_avg_tst %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_tst = '.') %>% 
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_9/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_9/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (3) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('training (split 9): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('test (split 9): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('training (split 9): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('test (split 9): best threshold')
  ss_tst_best
  
  ss_m1_split9 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})

### split (10)
job::job(m1_split10_output ={
  
  #----------------------------------------------------------------------------#
  #                   Extract indices of model fit
  #----------------------------------------------------------------------------#
  
  # - alpha.min & lambda.min: for model w/minimum cross validation error
  # - alpha.1se & lambda.1se: for sparsest model within 1SE of the 
  #   minimum cross validation error
  m1_alpha.1se <- m1_split10$fit_m1$alpha.1se
  m1_alpha.min <- m1_split10$fit_m1$alpha.min  
  m1_lambda.1se <- m1_split10$fit_m1$lambda.1se
  m1_lambda.min <- m1_split10$fit_m1$lambda.min
  
  m1_split10_fit <- cbind(
    m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  row.names(m1_split10_fit) <- c('split 10')
  m1_split10_fit <- m1_split10_fit %>% 
    as.data.frame() %>% 
    round(., 3)
  
  rm(m1_alpha.1se, m1_alpha.min, m1_lambda.1se, m1_lambda.min)
  
  # - lambda: sequence of lambdas fit
  # - df: number of nonzero coefficients for each value of lambda & alpha
  m1_lambda <- m1_split10$fit_m1$lambda %>% 
    as.data.frame() %>% 
    rename(lambda = '.')
  
  m1_split10_lambda_range <- m1_lambda %>% 
    summarise(
      lambda_min = min(lambda),
      lambda_max = max(lambda)) %>% 
    round(., 3) 
  row.names(m1_split10_lambda_range) <- c('split 10')
  
  m1_df <- m1_split10$fit_m1$df %>% 
    as.data.frame() %>% 
    rename(alpha_0 = V1,
           alpha_0.2 = V2,
           alpha_0.4 = V3,
           alpha_0.6 = V4,
           alpha_0.8 = V5,
           alpha_1.0 = V6)
  
  m1_lambda_df <- cbind(m1_lambda, m1_df)
  
  write.csv(
    m1_lambda_df,'output/study2/elastic_net/split_10/model_1/m1_lambda_df.csv')
  rm(m1_lambda, m1_df, m1_lambda_df)
  
  # - cvm: average cross validation error for each lambda and alpha
  # - cvse: standard error of ’cvm’
  m1_cvm <- m1_split10$fit_m1$cvm %>% 
    as.data.frame() %>% 
    rename(cvm_alpha_0 = V1,
           cvm_alpha_0.2 = V2,
           cvm_alpha_0.4 = V3,
           cvm_alpha_0.6 = V4,
           cvm_alpha_0.8 = V5,
           cvm_alpha_1 = V6)
  m1_cvm
  
  m1_cvse <- m1_split10$fit_m1$cvse %>% 
    as.data.frame() %>% 
    rename(cvse_alpha_0 = V1,
           cvse_alpha_0.2 = V2,
           cvse_alpha_0.4 = V3,
           cvse_alpha_0.6 = V4,
           cvse_alpha_0.8 = V5,
           cvse_alpha_1 = V6)
  m1_cvse
  
  m1_cv <- cbind(m1_cvm, m1_cvse)
  m1_cv
  names(m1_cv)
  
  write.csv(m1_cv,'output/study2/elastic_net/split_10/model_1/m1_cv.csv')
  rm(m1_cvm, m1_cvse, m1_cv)
  
  # extract coefficients for alpha.min and alpha.1se 
  m1_cf_1se <- coef(m1_split10$fit_m1,
                    lambda = m1_split10$fit_m1$lambda.1se, 
                    alpha = m1_split10$fit_m1$alpha.1se) 
  m1_cf_min <- coef(m1_split10$fit_m1, 
                    lambda = m1_split10$fit_m1$lambda.min, 
                    alpha = m1_split10$fit_m1$alpha.min) 
  
  m1_cf <- cbind(m1_cf_1se, m1_cf_min) %>% 
    as.data.frame() %>% 
    mutate(m1_cf_min_exp = exp(abs(m1_cf_min))) %>% 
    arrange(desc(abs(m1_cf_min))) %>% 
    rownames_to_column() %>% 
    rename(variable = rowname)
  
  m1_cf <- m1_cf %>% 
    left_join(table_names, by = 'variable') %>% 
    relocate(domain_name, table_name, variable)
  
  # subset non-zero coefficients 
  m1_cf_min_nz <- m1_cf_min %>% 
    subset(. != 0)
  names(m1_cf_min_nz) 
  cat(print(length(m1_cf_min_nz) - 1), 'non-zero predictors retained in model (1) split (10)  \n') 
  
  write.csv(m1_cf,'output/study2/elastic_net/split_10/model_1/m1_cf.csv', 
            row.names = TRUE)
  
  rm(m1_cf_1se, m1_cf_min)
  
  #----------------------------------------------------------------------------#  
  #    New list of imputed design matrices & responses for variables retained
  #----------------------------------------------------------------------------#  
  
  dfs_trn_m1_AUC <- lapply(1:5, function(i){
    trn_10_temp <- trn_10_imp_df_m1 %>% 
      subset(.imp == i)
    trn_10_temp[, names(trn_10_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  dfs_tst_m1_AUC <- lapply(1:5, function(i){
    tst_10_temp <- tst_10_imp_df_m1 %>% 
      subset(.imp == i)
    tst_10_temp[, names(tst_10_temp) %in% names(m1_cf_min_nz)]    
  }
  ) 
  
  x_trn_m1_AUC <- list()
  for (i in 1:5) {
    x_trn_m1_AUC[[i]] <- as.matrix(dfs_trn_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  x_tst_m1_AUC <- list()
  for (i in 1:5) {
    x_tst_m1_AUC[[i]] <- as.matrix(dfs_tst_m1_AUC[[i]][,(1:length(m1_cf_min_nz)-1)]) 
  }
  
  cat(print(dim(x_tst_m1_AUC[[1]])), 
      'sample size and number of non-zero predictors retained in test dataset for model (1) split (10)  \n') 
  
  rm(dfs_trn_m1_AUC, dfs_tst_m1_AUC)
  
  #------------------------------------------------------------------------------#  
  #                 Calculate AUC for training and test dataset
  #------------------------------------------------------------------------------#  
  
  # assign each imputation to an individual matrix
  m1_imp1_trn <-x_trn_m1_AUC[[1]] 
  m1_imp2_trn <-x_trn_m1_AUC[[2]] 
  m1_imp3_trn <-x_trn_m1_AUC[[3]] 
  m1_imp4_trn <-x_trn_m1_AUC[[4]] 
  m1_imp5_trn <-x_trn_m1_AUC[[5]] 
  
  m1_imp1_tst <-x_tst_m1_AUC[[1]] 
  m1_imp2_tst <-x_tst_m1_AUC[[2]] 
  m1_imp3_tst <-x_tst_m1_AUC[[3]] 
  m1_imp4_tst <-x_tst_m1_AUC[[4]] 
  m1_imp5_tst <-x_tst_m1_AUC[[5]] 
  
  # create a vector of 1s for intercept
  m1_imp1_trn <- cbind(rep(1,5121),m1_imp1_trn)
  m1_imp2_trn <- cbind(rep(1,5121),m1_imp2_trn)
  m1_imp3_trn <- cbind(rep(1,5121),m1_imp3_trn)
  m1_imp4_trn <- cbind(rep(1,5121),m1_imp4_trn)
  m1_imp5_trn <- cbind(rep(1,5121),m1_imp5_trn)
  
  m1_imp1_tst <- cbind(rep(1,1708),m1_imp1_tst)
  m1_imp2_tst <- cbind(rep(1,1708),m1_imp2_tst)
  m1_imp3_tst <- cbind(rep(1,1708),m1_imp3_tst)
  m1_imp4_tst <- cbind(rep(1,1708),m1_imp4_tst)
  m1_imp5_tst <- cbind(rep(1,1708),m1_imp5_tst)
  
  # compute probabilities
  
  m1_perc_imp1_trn <- m1_imp1_trn%*%m1_cf_min_nz 
  m1_perc_imp2_trn <- m1_imp2_trn%*%m1_cf_min_nz
  m1_perc_imp3_trn <- m1_imp3_trn%*%m1_cf_min_nz
  m1_perc_imp4_trn <- m1_imp4_trn%*%m1_cf_min_nz
  m1_perc_imp5_trn <- m1_imp5_trn%*%m1_cf_min_nz
  
  m1_perc_imp1_tst <- m1_imp1_tst%*%m1_cf_min_nz 
  m1_perc_imp2_tst <- m1_imp2_tst%*%m1_cf_min_nz
  m1_perc_imp3_tst <- m1_imp3_tst%*%m1_cf_min_nz
  m1_perc_imp4_tst <- m1_imp4_tst%*%m1_cf_min_nz
  m1_perc_imp5_tst <- m1_imp5_tst%*%m1_cf_min_nz
  
  m1_imp1_prob_trn <- exp(m1_perc_imp1_trn/(1+exp(m1_perc_imp1_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_trn = V1)
  m1_imp2_prob_trn <- exp(m1_perc_imp2_trn/(1+exp(m1_perc_imp2_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_trn = V1)
  m1_imp3_prob_trn <- exp(m1_perc_imp3_trn/(1+exp(m1_perc_imp3_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_trn = V1)
  m1_imp4_prob_trn <- exp(m1_perc_imp4_trn/(1+exp(m1_perc_imp4_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_trn = V1)
  m1_imp5_prob_trn <- exp(m1_perc_imp5_trn/(1+exp(m1_perc_imp5_trn))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_trn = V1)
  
  m1_imp1_prob_tst <- exp(m1_perc_imp1_tst/(1+exp(m1_perc_imp1_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp1_prob_tst = V1)
  m1_imp2_prob_tst <- exp(m1_perc_imp2_tst/(1+exp(m1_perc_imp2_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp2_prob_tst = V1)
  m1_imp3_prob_tst <- exp(m1_perc_imp3_tst/(1+exp(m1_perc_imp3_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp3_prob_tst = V1)
  m1_imp4_prob_tst <- exp(m1_perc_imp4_tst/(1+exp(m1_perc_imp4_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp4_prob_tst = V1)
  m1_imp5_prob_tst <- exp(m1_perc_imp5_tst/(1+exp(m1_perc_imp5_tst))) %>% 
    as.data.frame() %>% 
    rename(m1_imp5_prob_tst = V1)
  
  # average probabilities aross 5 imputations
  m1_prob_trn <- cbind(
    m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn, m1_imp4_prob_trn, 
    m1_imp5_prob_trn) %>% 
    mutate(m1_prob_avg_trn = rowMeans(.))
  write.csv(m1_prob_trn,'output/study2/elastic_net/split_10/model_1/m1_prob_trn.csv', 
            row.names = TRUE)  
  
  m1_prob_tst <- cbind(
    m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, m1_imp4_prob_tst, 
    m1_imp5_prob_tst) %>% 
    mutate(m1_prob_avg_tst = rowMeans(.))
  write.csv(m1_prob_tst,'output/study2/elastic_net/split_10/model_1/m1_prob_tst.csv', 
            row.names = TRUE)
  
  m1_prob_avg_trn <- m1_prob_trn %>% 
    select(m1_prob_avg_trn)
  m1_prob_avg_trn <- as.vector(m1_prob_avg_trn$m1_prob_avg_trn)
  
  m1_prob_avg_tst <- m1_prob_tst %>% 
    select(m1_prob_avg_tst)
  m1_prob_avg_tst <- as.vector(m1_prob_avg_tst$m1_prob_avg_tst)
  
  # double check each probability across the imputations is different 
  
  dc1 <- identical(m1_imp1_prob_tst, m1_imp2_prob_tst)
  dc2 <- identical(m1_imp1_prob_tst, m1_imp3_prob_tst)
  dc3 <- identical(m1_imp1_prob_tst, m1_imp4_prob_tst)
  dc4 <- identical(m1_imp1_prob_tst, m1_imp5_prob_tst)
  dc5 <- identical(m1_imp2_prob_tst, m1_imp3_prob_tst)
  dc6 <- identical(m1_imp2_prob_tst, m1_imp4_prob_tst)
  dc7 <- identical(m1_imp2_prob_tst, m1_imp5_prob_tst)
  dc8 <- identical(m1_imp3_prob_tst, m1_imp4_prob_tst)
  dc9 <- identical(m1_imp3_prob_tst, m1_imp5_prob_tst)
  dc10 <- identical(m1_imp4_prob_tst, m1_imp5_prob_tst)
  
  dc11 <- identical(m1_imp1_prob_trn, m1_imp2_prob_trn)
  dc12 <- identical(m1_imp1_prob_trn, m1_imp3_prob_trn)
  dc13 <- identical(m1_imp1_prob_trn, m1_imp4_prob_trn)
  dc14 <- identical(m1_imp1_prob_trn, m1_imp5_prob_trn)
  dc15 <- identical(m1_imp2_prob_trn, m1_imp3_prob_trn)
  dc16 <- identical(m1_imp2_prob_trn, m1_imp4_prob_trn)
  dc17 <- identical(m1_imp2_prob_trn, m1_imp5_prob_trn)
  dc18 <- identical(m1_imp3_prob_trn, m1_imp4_prob_trn)
  dc19 <- identical(m1_imp3_prob_trn, m1_imp5_prob_trn)
  dc20 <- identical(m1_imp4_prob_trn, m1_imp5_prob_trn)
  
  dc21 <- identical(m1_imp1_prob_tst, m1_imp1_prob_trn)
  dc22 <- identical(m1_imp2_prob_tst, m1_imp2_prob_trn)
  dc23 <- identical(m1_imp3_prob_tst, m1_imp3_prob_trn)
  dc24 <- identical(m1_imp4_prob_tst, m1_imp4_prob_trn)
  dc25 <- identical(m1_imp5_prob_tst, m1_imp5_prob_trn)
  
  dc26 <- identical(m1_prob_avg_tst, m1_prob_avg_trn)
  
  if (any(dc1 == FALSE && dc2 == FALSE && dc3 == FALSE && 
          dc4 == FALSE && dc5 == FALSE && dc6 == FALSE && 
          dc7 == FALSE && dc8 == FALSE && dc9 == FALSE && 
          dc10 == FALSE && dc11 == FALSE && dc12 == FALSE && 
          dc13 == FALSE && dc14 == FALSE && dc15 == FALSE && 
          dc16 == FALSE && dc17 == FALSE && dc18 == FALSE && 
          dc19 == FALSE && dc20 == FALSE && dc21 == FALSE && 
          dc22 == FALSE && dc23 == FALSE && dc24 == FALSE && 
          dc25 == FALSE && dc26 == FALSE)) {
    print('No issues w/probability double check')
  } else {
    print('Issues w/probability double check')
  }
  
  rm(
    x_trn_m1_AUC, x_tst_m1_AUC, 
    
    m1_prob_trn, m1_prob_tst, 
    
    m1_imp1_trn, m1_imp2_trn, m1_imp3_trn, m1_imp4_trn, m1_imp5_trn,
    m1_perc_imp1_trn, m1_perc_imp2_trn, m1_perc_imp3_trn, m1_perc_imp4_trn, 
    m1_perc_imp5_trn, m1_imp1_prob_trn, m1_imp2_prob_trn, m1_imp3_prob_trn,
    m1_imp4_prob_trn, m1_imp5_prob_trn,
    
    m1_imp1_tst, m1_imp2_tst, m1_imp3_tst, m1_imp4_tst, m1_imp5_tst,
    m1_perc_imp1_tst, m1_perc_imp2_tst, m1_perc_imp3_tst, m1_perc_imp4_tst, 
    m1_perc_imp5_tst, m1_imp1_prob_tst, m1_imp2_prob_tst, m1_imp3_prob_tst, 
    m1_imp4_prob_tst, m1_imp5_prob_tst,
    
    dc1, dc2, dc3, dc4, dc5, dc6, dc7, dc8, dc9, dc10,
    dc11, dc12, dc13, dc14, dc15, dc16, dc17, dc18, dc19, dc20,
    dc21, dc22, dc23, dc24, dc25, dc26)
  
  # double check direction for roc
  
  # - subset DV 
  trn_m1_DV_actual <- trn_10_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  tst_m1_DV_actual <- tst_10_obs_df_m1 %>% 
    rename(actual = DV) %>% 
    select(actual)
  
  # - bind actual and average predicted probability values across 5 imputations
  trn_m1_actual_prob <- cbind(trn_m1_DV_actual, m1_prob_avg_trn) 
  tst_m1_actual_prob <- cbind(tst_m1_DV_actual, m1_prob_avg_tst) 
  rm(trn_m1_DV_actual, tst_m1_DV_actual)
  
  names(trn_m1_actual_prob)
  names(tst_m1_actual_prob)
  
  # - calculate median probability split by DV
  trn_median <- trn_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_trn))
  
  tst_median <- tst_m1_actual_prob %>% 
    group_by(actual) %>% 
    summarise(median_value = median(m1_prob_avg_tst))
  
  trn_median_DV0 <- trn_median %>% 
    subset(actual == 0) 
  trn_median_DV1 <- trn_median %>% 
    subset(actual == 1) 
  
  tst_median_DV0 <- tst_median %>% 
    subset(actual == 0) 
  tst_median_DV1 <- tst_median %>% 
    subset(actual == 1) 
  
  if ((trn_median_DV0$median_value < trn_median_DV1$median_value)) {
    print(
      'training (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'training (across 5 imputations): median of controls > median of cases')
  }
  
  if ((tst_median_DV0$median_value < tst_median_DV1$median_value)) {
    print(
      'test (across 5 imputations): median of controls < median of cases')
  } else {
    print(
      'test (across 5 imputations): median of controls > median of cases')
  }
  
  rm(trn_median, tst_median, 
     trn_median_DV0, trn_median_DV1, 
     tst_median_DV0, tst_median_DV1)
  
  #----------------------------------------------------------------------------#  
  #                   Calculate AUC and additional metrics
  #----------------------------------------------------------------------------#  
  
  # (1) AUC
  
  # - training dataset
  m1_roc_trn <- roc(trn_10_obs_df_m1$DV, m1_prob_avg_trn) %>% 
    as.vector()
  
  AUC_CI_trn_m1_split10 <- ci.auc(m1_roc_trn) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_trn_m1_split10) <- c('training: split 10')
  
  # - test dataset
  m1_roc_tst <- roc(tst_10_obs_df_m1$DV, m1_prob_avg_tst) %>% 
    as.vector()
  
  AUC_CI_tst_m1_split10 <- ci.auc(m1_roc_tst) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(
      lower_CI = V1,
      AUC = V2,
      upper_CI = V3) %>% 
    relocate(AUC) %>% 
    mutate(
      AUC = as.numeric(AUC),
      lower_CI = as.numeric(lower_CI),
      upper_CI = as.numeric(upper_CI)) %>% 
    round(., 3)
  row.names(AUC_CI_tst_m1_split10) <- c('test: split 10')
  
  # (2) confusion matrix
  # - threshold > 0.50 
  
  trn_m1_binary_pred <- m1_prob_avg_trn %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_trn = '.') %>%  
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_trn > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  tst_m1_binary_pred <- m1_prob_avg_tst %>%
    as.data.frame() %>% 
    rename(m1_prob_avg_tst = '.') %>%  
    mutate(predicted_avg_0.50 = ifelse(m1_prob_avg_tst > 0.50, 1, 0)) %>% 
    select(predicted_avg_0.50)
  
  trn_m1_cm_data <- cbind(trn_m1_actual_prob, trn_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  tst_m1_cm_data <- cbind(tst_m1_actual_prob, tst_m1_binary_pred) %>% 
    select(actual, predicted_avg_0.50) %>% 
    mutate_at(c('actual', 'predicted_avg_0.50'), as.factor)
  
  trn_cm_0.50 <- confusionMatrix(
    trn_m1_cm_data$predicted_avg_0.50, trn_m1_cm_data$actual, positive = c("1")) 
  tst_cm_0.50 <- confusionMatrix(
    tst_m1_cm_data$predicted_avg_0.50, tst_m1_cm_data$actual, positive = c("1")) 
  
  capture.output(trn_cm_0.50, 
                 file = 'output/study2/elastic_net/split_10/model_1/trn_cm_0.50.txt') 
  capture.output(tst_cm_0.50, 
                 file = 'output/study2/elastic_net/split_10/model_1/tst_cm_0.50.txt') 
  
  rm(trn_m1_cm_data, tst_m1_cm_data, trn_cm_0.50, tst_cm_0.50)
  
  # (3) sensitivity and specificity
  
  # - 50% threshold 
  ss_trn_0.50 <- coords(m1_roc_trn, 0.50, transpose = FALSE) 
  row.names(ss_trn_0.50) <- c('training (split 10): 50% threshold')
  ss_trn_0.50 
  
  ss_tst_0.50 <- coords(m1_roc_tst, 0.50, transpose = FALSE)
  row.names(ss_tst_0.50) <- c('test (split 10): 50% threshold')
  ss_tst_0.50
  
  # - best threshold
  trn_best <- coords(m1_roc_trn, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  tst_best <- coords(m1_roc_tst, 'best', ret = 'threshold', transpose = FALSE) %>% 
    as.numeric()
  
  ss_trn_best <- coords(
    m1_roc_trn, trn_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_trn_best) <- c('training (split 10): best threshold')
  ss_trn_best
  
  ss_tst_best <- coords(
    m1_roc_tst, tst_best, transpose = FALSE, best.method = 'youden') 
  row.names(ss_tst_best) <- c('test (split 10): best threshold')
  ss_tst_best
  
  ss_m1_split10 <- round(
    rbind(ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best), 3)
  
  rm(
    m1_cf_min_nz, 
    m1_prob_avg_trn, m1_prob_avg_tst, 
    trn_m1_actual_prob, tst_m1_actual_prob,
    trn_m1_binary_pred, tst_m1_binary_pred,
    trn_best, tst_best, 
    ss_trn_0.50, ss_tst_0.50, ss_trn_best, ss_tst_best,
    m1_roc_trn, m1_roc_tst)
  
})


# Average Results Across Splits ------------------------------------------------

# (1) fit values: alpha & lambda min and 1se
fit_m1 <- rbind(
  m1_split1_output$m1_split1_fit, m1_split2_output$m1_split2_fit,
  m1_split3_output$m1_split3_fit, m1_split4_output$m1_split4_fit,
  m1_split5_output$m1_split5_fit, m1_split6_output$m1_split6_fit,
  m1_split7_output$m1_split7_fit, m1_split8_output$m1_split8_fit,
  m1_split9_output$m1_split9_fit, m1_split10_output$m1_split10_fit)



write.csv(fit_m1, 'output/study2/elastic_net/fit_m1.csv')    

# (2) range for lambda (individual lambda exported within job)
lambda_range_m1 <- rbind(
  m1_split1_output$m1_split1_lambda_range, m1_split2_output$m1_split2_lambda_range,
  m1_split3_output$m1_split3_lambda_range, m1_split4_output$m1_split4_lambda_range,
  m1_split5_output$m1_split5_lambda_range, m1_split6_output$m1_split6_lambda_range,
  m1_split7_output$m1_split7_lambda_range, m1_split8_output$m1_split8_lambda_range,
  m1_split9_output$m1_split9_lambda_range, m1_split10_output$m1_split10_lambda_range)

write.csv(lambda_range_m1, 'output/study2/elastic_net/lambda_range_m1.csv')

# (3) AUC
AUC_m1 <- rbind(
  m1_split1_output$AUC_CI_trn_m1_split1, m1_split2_output$AUC_CI_trn_m1_split2,
  m1_split3_output$AUC_CI_trn_m1_split3, m1_split4_output$AUC_CI_trn_m1_split4,
  m1_split5_output$AUC_CI_trn_m1_split5, m1_split6_output$AUC_CI_trn_m1_split6,
  m1_split7_output$AUC_CI_trn_m1_split7, m1_split8_output$AUC_CI_trn_m1_split8,
  m1_split9_output$AUC_CI_trn_m1_split9, m1_split10_output$AUC_CI_trn_m1_split10,
  
  m1_split1_output$AUC_CI_tst_m1_split1, m1_split2_output$AUC_CI_tst_m1_split2,
  m1_split3_output$AUC_CI_tst_m1_split3, m1_split4_output$AUC_CI_tst_m1_split4,
  m1_split5_output$AUC_CI_tst_m1_split5, m1_split6_output$AUC_CI_tst_m1_split6,
  m1_split7_output$AUC_CI_tst_m1_split7, m1_split8_output$AUC_CI_tst_m1_split8,
  m1_split9_output$AUC_CI_tst_m1_split9, m1_split10_output$AUC_CI_tst_m1_split10)


write.csv(AUC_m1, 'output/study2/elastic_net/AUC_m1.csv')

# (4) sensitivity and specificity 
ss_m1 <- rbind(
  m1_split1_output$ss_m1_split1, m1_split2_output$ss_m1_split2,
  m1_split3_output$ss_m1_split3, m1_split4_output$ss_m1_split4,
  m1_split5_output$ss_m1_split5, m1_split6_output$ss_m1_split6,
  m1_split7_output$ss_m1_split7, m1_split8_output$ss_m1_split8,
  m1_split9_output$ss_m1_split9, m1_split10_output$ss_m1_split10)



write.csv(ss_m1, 'output/study2/elastic_net/ss_m1.csv') 


# (5) coefficients for variables retained + n splits retained 

## model 1
m1_split1_coef <- m1_split1_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split1 = m1_cf_min,
         m1_split1_exp = m1_cf_min_exp) 

m1_split2_coef <- m1_split2_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split2 = m1_cf_min,
         m1_split2_exp = m1_cf_min_exp) 

m1_split3_coef <- m1_split3_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split3 = m1_cf_min,
         m1_split3_exp = m1_cf_min_exp) 

m1_split4_coef <- m1_split4_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split4 = m1_cf_min,
         m1_split4_exp = m1_cf_min_exp) 

m1_split5_coef <- m1_split5_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split5 = m1_cf_min,
         m1_split5_exp = m1_cf_min_exp) 

m1_split6_coef <- m1_split6_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split6 = m1_cf_min,
         m1_split6_exp = m1_cf_min_exp) 

m1_split7_coef <- m1_split7_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split7 = m1_cf_min,
         m1_split7_exp = m1_cf_min_exp) 

m1_split8_coef <- m1_split8_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split8 = m1_cf_min,
         m1_split8_exp = m1_cf_min_exp) 

m1_split9_coef <- m1_split9_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split9 = m1_cf_min,
         m1_split9_exp = m1_cf_min_exp) 

m1_split10_coef <- m1_split10_output$m1_cf %>% 
  subset(m1_cf_min != 0) %>% 
  select(domain_name, table_name, variable, m1_cf_min, m1_cf_min_exp) %>% 
  rename(m1_split10 = m1_cf_min,
         m1_split10_exp = m1_cf_min_exp) 


m1_coef <- m1_split1_coef %>% 
  full_join(m1_split2_coef, by = c('variable', 'domain_name', 'table_name')) %>% 
  full_join(m1_split3_coef, by = c('variable', 'domain_name', 'table_name')) %>% 
  full_join(m1_split4_coef, by = c('variable', 'domain_name', 'table_name')) %>% 
  full_join(m1_split5_coef, by = c('variable', 'domain_name', 'table_name')) %>% 
  full_join(m1_split6_coef, by = c('variable', 'domain_name', 'table_name')) %>% 
  full_join(m1_split7_coef, by = c('variable', 'domain_name', 'table_name')) %>% 
  full_join(m1_split8_coef, by = c('variable', 'domain_name', 'table_name')) %>% 
  full_join(m1_split9_coef, by = c('variable', 'domain_name', 'table_name')) %>% 
  full_join(m1_split10_coef, by = c('variable', 'domain_name', 'table_name')) 


rm(
  m1_split1_coef, m1_split2_coef, m1_split3_coef, m1_split4_coef, m1_split5_coef, 
  m1_split6_coef, m1_split7_coef, m1_split8_coef, m1_split9_coef, m1_split10_coef,
  
  m2_split1_coef, m2_split2_coef, m2_split3_coef, m2_split4_coef, m2_split5_coef, 
  m2_split6_coef, m2_split7_coef, m2_split8_coef, m2_split9_coef, m2_split10_coef,
  
  m3_split1_coef, m3_split2_coef, m3_split3_coef, m3_split4_coef, m3_split5_coef, 
  m3_split6_coef, m3_split7_coef, m3_split8_coef, m3_split9_coef, m3_split10_coef)

### average across splits (supplementary materials)

# - create dataframe with coefficients and exponentiated coefficients 
m1_coef_supp <- m1_coef %>% 
  select(
    domain_name, table_name, variable,
    ends_with('split1') | ends_with('split2') | ends_with('split3') | 
      ends_with('split4') | ends_with('split5') | ends_with('split6') |
      ends_with('split7') | ends_with('split8') | ends_with('split9') | 
      ends_with('split10'))

m1_coef_supp_exp <- m1_coef %>% 
  select(
    domain_name, table_name, variable, ends_with('exp'))


# - averages for coefficients and exponentiated coefficients
m1_coef_supp <- m1_coef_supp %>% 
  slice(-1) %>% # drop intercept
  mutate(n_split = rowSums(!is.na(select(., -domain_name, -table_name, -variable)))) %>% 
  mutate(coef_avg = rowMeans(select(., starts_with('m1')), na.rm = TRUE))  %>% 
  mutate(order_variable = abs(coef_avg)) %>% # create absolute value indicator for arranging variables 
  arrange(., desc(order_variable)) %>% 
  select(-order_variable) 

m1_coef_supp_exp <- m1_coef_supp_exp %>% 
  slice(-1) %>% # drop intercept
  mutate(n_split = rowSums(!is.na(select(., -domain_name, -table_name, -variable)))) %>% 
  mutate(coef_avg = rowMeans(select(., starts_with('m1')), na.rm = TRUE))  %>% 
  mutate(order_variable = abs(coef_avg)) %>% # create absolute value indicator for arranging variables 
  arrange(., desc(order_variable)) %>% 
  select(-order_variable) 


# add min and max for coefficients 
m1_coef_supp$coef_max <- pmax(
  m1_coef_supp$m1_split1, m1_coef_supp$m1_split2, m1_coef_supp$m1_split3, m1_coef_supp$m1_split4, 
  m1_coef_supp$m1_split5, m1_coef_supp$m1_split6, m1_coef_supp$m1_split7, m1_coef_supp$m1_split8,
  m1_coef_supp$m1_split9, m1_coef_supp$m1_split10, na.rm = TRUE)
m1_coef_supp$coef_min <- pmin(
  m1_coef_supp$m1_split1, m1_coef_supp$m1_split2, m1_coef_supp$m1_split3, m1_coef_supp$m1_split4, 
  m1_coef_supp$m1_split5, m1_coef_supp$m1_split6, m1_coef_supp$m1_split7, m1_coef_supp$m1_split8,
  m1_coef_supp$m1_split9, m1_coef_supp$m1_split10, na.rm = TRUE)

m1_coef_supp_exp$coef_max <- pmax(
  m1_coef_supp$m1_split1, m1_coef_supp$m1_split2, m1_coef_supp$m1_split3, m1_coef_supp$m1_split4, 
  m1_coef_supp$m1_split5, m1_coef_supp$m1_split6, m1_coef_supp$m1_split7, m1_coef_supp$m1_split8,
  m1_coef_supp$m1_split9, m1_coef_supp$m1_split10, na.rm = TRUE)
m1_coef_supp_exp$coef_min <- pmin(
  m1_coef_supp$m1_split1, m1_coef_supp$m1_split2, m1_coef_supp$m1_split3, m1_coef_supp$m1_split4, 
  m1_coef_supp$m1_split5, m1_coef_supp$m1_split6, m1_coef_supp$m1_split7, m1_coef_supp$m1_split8,
  m1_coef_supp$m1_split9, m1_coef_supp$m1_split10, na.rm = TRUE)


# rename variables
m1_coef_supp <- m1_coef_supp %>% 
  rename(
    'Domain' = domain_name,
    'Variable' = table_name,
    'Split 1' = m1_split1,
    'Split 2' = m1_split2,
    'Split 3' = m1_split3,
    'Split 4' = m1_split4,
    'Split 5' = m1_split5,
    'Split 6' = m1_split6,
    'Split 7' = m1_split7,
    'Split 8' = m1_split8,
    'Split 9' = m1_split9,
    'Split 10' = m1_split10,
    'Number of Splits' = n_split,
    'Average Across Splits' = coef_avg,
    'Minimum Across Splits' = coef_min,
    'Maximum Across Splits' = coef_max) %>% 
  mutate(across(where(is.numeric), round, 3)) 

m1_coef_supp_exp <- m1_coef_supp_exp %>% 
  rename(
    'Domain' = domain_name,
    'Variable' = table_name,
    'Split 1' = m1_split1_exp,
    'Split 2' = m1_split2_exp,
    'Split 3' = m1_split3_exp,
    'Split 4' = m1_split4_exp,
    'Split 5' = m1_split5_exp,
    'Split 6' = m1_split6_exp,
    'Split 7' = m1_split7_exp,
    'Split 8' = m1_split8_exp,
    'Split 9' = m1_split9_exp,
    'Split 10' = m1_split10_exp,
    'Number of Splits' = n_split,
    'Average Across Splits' = coef_avg,
    'Minimum Across Splits' = coef_min,
    'Maximum Across Splits' = coef_max) %>% 
  mutate(across(where(is.numeric), round, 3)) 






write.csv(m1_coef_supp, 'output/study2/elastic_net/m1_coef_supp.csv') # n = 85 variables

write.csv(m1_coef_supp_exp, 'output/study2/elastic_net/m1_coef_supp_exp.csv') # n = 85 variables

# indicator of number of variables retained within each domain
m1_coef_domain_n <- table(m1_coef_supp$Domain)


capture.output(m1_coef_domain_n, 
               file = 'output/study2/elastic_net/m1_coef_domain_n.txt')

### average across splits for variables retained in >5 splits  
m1_coef_top5 <- m1_coef_supp %>% 
  filter(`Number of Splits` > 5)



m1_coef_domain_n_top5 <- table(m1_coef_top5$Domain)

write.csv(m1_coef_top5, 'output/study2/elastic_net/top_5/m1_coef_top5.csv')

capture.output(m1_coef_domain_n_top5, 
               file = 'output/study2/elastic_net/top_5/m1_coef_domain_n_top5.txt')

save.image("output/study2/data_management/elastic_net_data.RData")

rm(
  m1_coef, m2_coef, m3_coef,
  m1_coef_supp, m2_coef_supp, m3_coef_supp,
  m1_coef_domain_n, m2_coef_domain_n, m3_coef_domain_n,
  m1_coef_top5, m2_coef_top5, m3_coef_top5,
  m1_coef_domain_n_top5, m2_coef_domain_n_top5, m3_coef_domain_n_top5)
