whatif_algo = function(predictor, n_cfactuals, x_interest, pred_column, desired_y_hat_range, X_search, distance_function,
                       fixed_features, epsilon) {
  y_hat = setDT(predictor$predict(X_search))[[pred_column]]
  if (!is.null(fixed_features)) {
    cat("Immutable features are considered:", paste(fixed_features, collapse = ", "), ".\n")
    X_fixed <- X_search |> select(fixed_features)
    num_features <- names(X_fixed)[sapply(X_fixed, is.numeric)]
    cat_features <- setdiff(names(X_fixed), num_features)
    if (length(cat_features) != 0L) {
      filter_idx = 
        (Reduce(`&`, Map(function(col) X_search[[col]] == x_interest[[col]], cat_features))) &
        (y_hat %between% desired_y_hat_range)
      X_search = X_search[filter_idx,]
      y_hat = y_hat[filter_idx]
      
      if (nrow(X_search) == 0L) {
        warning("Impossible to find training instances for categ. variable(s). \n")
      }
      
    }
    if (length(num_features) != 0L) {
      val_features <- x_interest[, ..num_features]
      lower_bounds <- val_features - epsilon*abs(val_features)
      upper_bounds <- val_features + epsilon*abs(val_features)
      filter_idx = 
        (Reduce(`&`, 
                Map(function(col) X_search[[col]] %between% c(lower_bounds[[col]], upper_bounds[[col]]), num_features))
        ) &
        (y_hat %between% desired_y_hat_range)
      X_search = X_search[filter_idx,]
      y_hat = y_hat[filter_idx]
      
      if (nrow(X_search) == 0L) {
        warning("Impossible to find training instances for num. variable(s) under epsilon neighborhood. \n")
      }
      
      # We assign to X_search the numerical feature values of x_interest
      X_search[, (num_features) := val_features]
      y_hat <- setDT(predictor$predict(X_search))[[pred_column]]
      X_search = X_search |> filter(y_hat %between% desired_y_hat_range)
      
      if (nrow(X_search) == 0L) {
        warning("Impossible to find training instances for num. variable(s) after setting epsilon to 0. \n")
      }
      
    }
  } else {
    X_search = X_search[y_hat %between% desired_y_hat_range]
  }
  X_search = unique(X_search)
  if (nrow(X_search) < n_cfactuals) {
    warning(sprintf("Could only find %s counterfactual(s)", nrow(X_search)))
  }
  
  if (nrow(X_search) == 0L) {
    return(X_search)
  }
  
  dist_matrix = eval_distance(distance_function, x_interest, X_search, predictor$data$X)
  if ("topn" %in% class(distance_function)) {
    idx = c(dist_matrix)
  } else {
    idx = smallest_n_indices(as.vector(dist_matrix), n = n_cfactuals) 
  }
  X_search[idx]
}
