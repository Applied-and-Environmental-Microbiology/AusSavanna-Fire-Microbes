
# Bootstrap models -------------------------------------------------------------

# Load parallel processing libraries
library(parallel)
library(doParallel)

# Set up parallel processing
n_cores <- detectCores() - 1
message(paste0("Using ", n_cores, " cores for parallel processing"))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Function to bootstrap a single model
bootstrap_model <- function(model, data, formula, n_boot = 10000, min_n = 8) {
  
  # Get unique blocks and treatment combinations
  blocks <- unique(data$block)
  
  # Count samples per treatment combination
  treatment_counts <- data %>%
    group_by(regime, litter_type, block) %>%
    summarise(n = n(), .groups = "drop")
  
  # Get proper column names from original model
  original_emmeans <- emmeans(model, ~ regime * litter_type)
  original_effects <- contrast(original_emmeans, "eff")
  original_effects_df <- as.data.frame(original_effects)
  proper_names <- as.character(original_effects_df$contrast)
  
  # Export necessary objects to cluster
  clusterExport(cl, c("data", "formula", "blocks", "treatment_counts", "min_n"), 
                envir = environment())
  
  # Progress tracking
  cat("Bootstrapping", n_boot, "iterations...\n")
  
  # Parallel bootstrap
  boot_results <- foreach(
    i = 1:n_boot,
    .packages = c("lme4", "emmeans", "dplyr"),
    .errorhandling = "pass"
  ) %dopar% {
    
    # Stratified resampling: resample within each regime x litter_type combination
    boot_data_list <- list()
    
    for (reg in levels(data$regime)) {
      for (lit in levels(data$litter_type)) {
        # Get data for this treatment combination
        subset_data <- data[data$regime == reg & data$litter_type == lit, ]
        
        if (nrow(subset_data) > 0) {
          # Get blocks for this combination
          combo_blocks <- unique(subset_data$block)
          
          # Resample blocks with replacement
          sampled_blocks <- sample(combo_blocks, length(combo_blocks), replace = TRUE)
          
          # Get resampled data
          combo_boot_data <- do.call(rbind, lapply(sampled_blocks, function(b) {
            subset_data[subset_data$block == b, ]
          }))
          
          boot_data_list[[paste(reg, lit, sep = "_")]] <- combo_boot_data
        }
      }
    }
    
    # Combine all resampled data
    boot_data <- do.call(rbind, boot_data_list)
    
    # Check if all treatment combinations have at least min_n samples
    boot_counts <- boot_data %>%
      group_by(regime, litter_type) %>%
      summarise(n = n(), .groups = "drop")
    
    if (any(boot_counts$n < min_n)) {
      return(NULL)  # Skip this iteration
    }
    
    # Refit model
    boot_model <- tryCatch({
      lmer(formula, data = boot_data)
    }, error = function(e) NULL)
    
    if (!is.null(boot_model)) {
      # Extract intercept
      intercept <- fixef(boot_model)[1]
      
      # Extract effect sizes (deviations from grand mean)
      boot_emmeans <- emmeans(boot_model, ~ regime * litter_type)
      boot_effects <- contrast(boot_emmeans, "eff")
      boot_effects_df <- as.data.frame(boot_effects)
      
      # Store results as unnamed vector (we'll add names later)
      list(
        intercept = intercept,
        effects = boot_effects_df$estimate
      )
    } else {
      NULL
    }
  }
  
  cat("Completed!\n")
  
  # Clean results (remove failed iterations)
  boot_results_clean <- Filter(Negate(is.null), boot_results)
  n_successful <- length(boot_results_clean)
  cat("Successful iterations:", n_successful, "out of", n_boot, "\n")
  
  # Extract intercepts
  intercepts_raw <- sapply(boot_results_clean, function(x) x$intercept)
  
  # Extract effect sizes - create matrix
  effects_matrix <- matrix(NA, nrow = length(boot_results_clean), ncol = length(proper_names))
  
  for (i in seq_along(boot_results_clean)) {
    effects_matrix[i, ] <- boot_results_clean[[i]]$effects
  }
  
  # Convert to data frame with proper column names
  effects_df_raw <- as.data.frame(effects_matrix)
  colnames(effects_df_raw) <- proper_names
  
  # Trim intercepts
  intercept_limits <- quantile(intercepts_raw, c(0.005, 0.995), na.rm = TRUE)
  intercepts_keep <- intercepts_raw >= intercept_limits[1] & intercepts_raw <= intercept_limits[2]
  intercepts <- intercepts_raw[intercepts_keep]
  
  # Trim effects - keep rows where ALL effects are within 99% range
  effects_keep <- rep(TRUE, nrow(effects_df_raw))
  
  for (col in colnames(effects_df_raw)) {
    col_limits <- quantile(effects_df_raw[[col]], c(0.005, 0.995), na.rm = TRUE)
    effects_keep <- effects_keep & 
      (effects_df_raw[[col]] >= col_limits[1]) & 
      (effects_df_raw[[col]] <= col_limits[2])
  }
  
  effects_df <- effects_df_raw[effects_keep, ]
  n_retained <- nrow(effects_df)
  
  cat("Retained", n_retained, "iterations after trimming (", 
      round(100 * n_retained / n_successful, 1), "%)\n", sep = "")
  
  # Function for 95% percentile CI
  ci <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
  
  # Intercept estimate and 95% CI
  intercept_mean <- mean(intercepts)
  intercept_ci <- ci(intercepts)
  intercept_formatted <- sprintf(
    "%.3f [%.3f, %.3f]",
    intercept_mean,
    intercept_ci[1],
    intercept_ci[2]
  )
  
  # Effect size summary for each factor combination
  effect_summary <- data.frame(
    contrast = colnames(effects_df),
    effect_mean = colMeans(effects_df),
    lower_ci = apply(effects_df, 2, function(x) ci(x)[1]),
    upper_ci = apply(effects_df, 2, function(x) ci(x)[2])
  ) %>%
    mutate(
      regime = sapply(strsplit(contrast, " "), function(x) x[1]),
      litter_type = sapply(strsplit(contrast, " "), function(x) x[2])
    )
  
  # Density data for plotting
  density_list <- lapply(colnames(effects_df), function(combo) {
    dens <- density(effects_df[[combo]], adjust = 3.0)
    tibble(
      Effect = dens$x,
      Density = dens$y,
      contrast = combo
    )
  })
  
  density_data <- bind_rows(density_list) %>%
    mutate(
      regime = sapply(strsplit(contrast, " "), function(x) x[1]),
      litter_type = sapply(strsplit(contrast, " "), function(x) x[2])
    )
  
  return(list(
    intercept = intercept_formatted,
    intercept_distribution = intercepts,
    effect_summary = effect_summary,
    density_data = density_data,
    effects_distribution = effects_df,
    n_successful = n_successful
  ))
}
