# stratified PU CV

pu_cv_xgboost <- function(data, features, label_col = "label",
                          n_folds = 5,
                          n_iterations = 50,
                          sample_frac = 0.5,
                          nrounds = 100,
                          seed = 123) {
  library(xgboost)
  library(dplyr)
  library(rsample)

  set.seed(seed)

  # Extract labels and feature matrix
  y <- data[[label_col]]
  X <- as.matrix(data[, features])
  n <- nrow(X)

  # Stratified K-fold CV
  folds <- vfold_cv(data.frame(y), v = n_folds, strata = y)
  all_cv_preds <- rep(0, n)

  for (fold in folds$splits) {
    train_idx <- as.integer(fold$in_id)
    val_idx <- setdiff(1:n, train_idx)

    y_train <- y[train_idx]
    X_train <- X[train_idx, , drop = FALSE]
    X_val <- X[val_idx, , drop = FALSE]

    positive_idx <- which(y_train == 1)
    unlabeled_idx <- which(y_train == 0)

    val_preds_matrix <- matrix(0, nrow = length(val_idx), ncol = n_iterations)

    for (i in 1:n_iterations) {
      sampled_unlabeled_idx <- sample(unlabeled_idx, size = floor(sample_frac * length(unlabeled_idx)))
      train_sub_idx <- c(positive_idx, sampled_unlabeled_idx)

      dtrain <- xgb.DMatrix(data = X_train[train_sub_idx, ],
                            label = y_train[train_sub_idx])
      pos_weight <- sum(y_train[train_sub_idx] == 0) / sum(y_train[train_sub_idx] == 1)

      model <- xgboost(
        data = dtrain,
        nrounds = nrounds,
        params = list(
          objective = "binary:logistic",
          eval_metric = "logloss",
          scale_pos_weight = pos_weight
        ),
        verbose = 0
      )

      val_preds_matrix[, i] <- predict(model, X_val)
    }

    # Average predictions across bagging iterations
    all_cv_preds[val_idx] <- rowMeans(val_preds_matrix)
  }

  return(all_cv_preds)
}

# utility functions


run_downsampled_pu_cv <- function(data, label_col, id_col = NULL,
                                  n_repeats = 5,
                                  pos_neg_ratio = 1.0,
                                  k = 5,
                                  mtry_neg_frac = 0.2,
                                  normalize = "zscore",
                                  seed = 42,
                                  verbose = TRUE) {
  set.seed(seed)
  all_results <- vector("list", n_repeats)

  for (r in seq_len(n_repeats)) {
    if (verbose) cat(sprintf("Repeat %d/%d:\n", r, n_repeats))

    positives <- data[data[[label_col]] == 1, ]
    negatives <- data[data[[label_col]] == 0, ]

    n_pos <- nrow(positives)
    n_neg_sample <- floor(pos_neg_ratio * n_pos)

    if (nrow(negatives) < n_neg_sample) {
      warning("Not enough negatives to sample desired ratio. Using all.")
      neg_sample <- negatives
    } else {
      neg_sample <- negatives[sample(nrow(negatives), n_neg_sample), ]
    }

    down_data <- rbind(positives, neg_sample)
    down_data <- down_data[sample(nrow(down_data)), , drop = FALSE]  # shuffle

    res <- pu_cv_xgboost(
      data = down_data,
      label_col = label_col,
      id_col = id_col,
      k = k,
      mtry_neg_frac = mtry_neg_frac,
      normalize = normalize,
      seed = seed + r,
      verbose = verbose
    )

    all_results[[r]] <- res
  }

  # Combine summaries
  aucs <- unlist(lapply(all_results, function(x) x$summary$aucs))
  loglosses <- unlist(lapply(all_results, function(x) x$summary$loglosses))

  list(
    mean_auc = mean(aucs, na.rm = TRUE),
    sd_auc = sd(aucs, na.rm = TRUE),
    mean_logloss = mean(loglosses, na.rm = TRUE),
    sd_logloss = sd(loglosses, na.rm = TRUE),
    results = all_results
  )
}



