
# Helper function --------------------------------------------------------------

# Function to bootstrap effect sizes
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
  cat("Bootstrapping", n_boot, "iterations in parallel with stratified sampling...\n")
  
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
  
  # Trim to 99% quantile (remove extreme 0.5% on each tail)
  cat("Trimming to 99% quantile range...\n")
  
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

# (1) Setup --------------------------------------------------------------------

# Required libraries
library(emmeans)
library(lme4)
library(DHARMa)
library(performance)
library(parameters)
library(parallel)
library(doParallel)
library(tidyverse)

# Soil bacteria data
data_bacteria_soil <- data.table::fread("data/sample_metadata_bacteria.txt") %>%
  filter(community == "Soil") %>%
  # Define factors and levels
  mutate(
    regime = factor(regime, levels = c("E2", "L2", "U")),
    litter_type = factor(litter_type, levels = c("Grass", "Leaf")),
    block = factor(block)
  )

# Litter bacteria data
data_bacteria_litter <- data.table::fread("data/sample_metadata_bacteria.txt") %>%
  filter(community == "Leaf") %>%
  # Define factors and levels
  mutate(
    regime = factor(regime, levels = c("E2", "L2", "U")),
    litter_type = factor(litter_type, levels = c("Grass", "Leaf")),
    block = factor(block)
  )

# Soil fungi data
data_fungi_soil <- data.table::fread("data/sample_metadata_fungi.txt") %>%
  filter(community == "Soil") %>%
  # Define factors and levels
  mutate(
    regime = factor(regime, levels = c("E2", "L2", "U")),
    litter_type = factor(litter_type, levels = c("Grass", "Leaf")),
    block = factor(block)
  )

# Litter fungi data
data_fungi_litter <- data.table::fread("data/sample_metadata_fungi.txt") %>%
  filter(community == "Leaf") %>%
  # Define factors and levels
  mutate(
    regime = factor(regime, levels = c("E2", "L2", "U")),
    litter_type = factor(litter_type, levels = c("Grass", "Leaf")),
    block = factor(block)
  )

# Set sum contrasts
contrasts(data_bacteria_soil$regime) <- contr.sum(levels(data_bacteria_soil$regime))
contrasts(data_bacteria_soil$litter_type) <- contr.sum(levels(data_bacteria_soil$litter_type))
contrasts(data_bacteria_litter$regime) <- contr.sum(levels(data_bacteria_litter$regime))
contrasts(data_bacteria_litter$litter_type) <- contr.sum(levels(data_bacteria_litter$litter_type))
contrasts(data_fungi_soil$regime) <- contr.sum(levels(data_fungi_soil$regime))
contrasts(data_fungi_soil$litter_type) <- contr.sum(levels(data_fungi_soil$litter_type))
contrasts(data_fungi_litter$regime) <- contr.sum(levels(data_fungi_litter$regime))
contrasts(data_fungi_litter$litter_type) <- contr.sum(levels(data_fungi_litter$litter_type))

# Count the number of samples for each treatment combination
data_bacteria_soil %>%
  group_by(regime, litter_type) %>%
  summarise(n = n())
data_bacteria_litter %>%
  group_by(regime, litter_type) %>%
  summarise(n = n())
data_fungi_soil %>%
  group_by(regime, litter_type) %>%
  summarise(n = n())
data_fungi_litter %>%
  group_by(regime, litter_type) %>%
  summarise(n = n())

# (2) Fit models ---------------------------------------------------------------

#### (2a) Fit models ####

# Soil bacteria community
model_richness_bacteria_soil <- lmer(
  richness_bacteria ~ regime * litter_type + (1 | block),
  data = data_bacteria_soil
)
model_evenness_bacteria_soil <- lmer(
  evenness_bacteria ~ regime * litter_type + (1 | block),
  data = data_bacteria_soil
)
model_shannon_bacteria_soil <- lmer(
  shannon_bacteria ~ regime * litter_type + (1 | block),
  data = data_bacteria_soil
)

# Litter bacteria community
model_richness_bacteria_litter <- lmer(
  richness_bacteria ~ regime * litter_type + (1 | block),
  data = data_bacteria_litter
)
model_evenness_bacteria_litter <- lmer(
  evenness_bacteria ~ regime * litter_type + (1 | block),
  data = data_bacteria_litter
)
model_shannon_bacteria_litter <- lmer(
  shannon_bacteria ~ regime * litter_type + (1 | block),
  data = data_bacteria_litter
)

# Soil fungi community
model_richness_fungi_soil <- lmer(
  richness_fungi ~ regime * litter_type + (1 | block),
  data = data_fungi_soil
)
model_evenness_fungi_soil <- lmer(
  evenness_fungi ~ regime * litter_type + (1 | block),
  data = data_fungi_soil
)
model_shannon_fungi_soil <- lmer(
  shannon_fungi ~ regime * litter_type + (1 | block),
  data = data_fungi_soil
)

# Litter fungi community
model_richness_fungi_litter <- lmer(
  richness_fungi ~ regime * litter_type + (1 | block),
  data = data_fungi_litter
)
model_evenness_fungi_litter <- lmer(
  evenness_fungi ~ regime * litter_type + (1 | block),
  data = data_fungi_litter
)
model_shannon_fungi_litter <- lmer(
  shannon_fungi ~ regime * litter_type + (1 | block),
  data = data_fungi_litter
)

#### (2b) Check assumptions ####

# Soil bacteria community
simulateResiduals(model_richness_bacteria_soil, plot = TRUE)
simulateResiduals(model_evenness_bacteria_soil, plot = TRUE)
simulateResiduals(model_shannon_bacteria_soil, plot = TRUE)

# Litter bacteria community
simulateResiduals(model_richness_bacteria_litter, plot = TRUE)
simulateResiduals(model_evenness_bacteria_litter, plot = TRUE)
simulateResiduals(model_shannon_bacteria_litter, plot = TRUE)

# Soil fungi community
simulateResiduals(model_richness_fungi_soil, plot = TRUE)
simulateResiduals(model_evenness_fungi_soil, plot = TRUE)
simulateResiduals(model_shannon_fungi_soil, plot = TRUE)

# Litter fungi community  
simulateResiduals(model_richness_fungi_litter, plot = TRUE)
simulateResiduals(model_evenness_fungi_litter, plot = TRUE)
simulateResiduals(model_shannon_fungi_litter, plot = TRUE)

# (3) Model summaries ----------------------------------------------------------

#### (3a) Model parameters ####

# Soil bacteria community
parameters_richness_bacteria_soil <- parameters(model_richness_bacteria_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil bacteria richness")
parameters_evenness_bacteria_soil <- parameters(model_evenness_bacteria_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil bacteria evenness")
parameters_shannon_bacteria_soil <- parameters(model_shannon_bacteria_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil bacteria shannon")

# Litter bacteria community
parameters_richness_bacteria_litter <- parameters(model_richness_bacteria_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter bacteria richness")
parameters_evenness_bacteria_litter <- parameters(model_evenness_bacteria_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter bacteria evenness")
parameters_shannon_bacteria_litter <- parameters(model_shannon_bacteria_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter bacteria shannon")

# Soil fungi community
parameters_richness_fungi_soil <- parameters(model_richness_fungi_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil fungi richness")
parameters_evenness_fungi_soil <- parameters(model_evenness_fungi_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil fungi evenness")
parameters_shannon_fungi_soil <- parameters(model_shannon_fungi_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil fungi shannon")

# Litter fungi community
parameters_richness_fungi_litter <- parameters(model_richness_fungi_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter fungi richness")
parameters_evenness_fungi_litter <- parameters(model_evenness_fungi_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter fungi evenness")
parameters_shannon_fungi_litter <- parameters(model_shannon_fungi_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter fungi shannon")

#### (3b) Marginal means ####

# Soil bacteria soil richness
emmeans_richness_bacteria_soil <- emmeans(
  model_richness_bacteria_soil,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil bacteria richness"
  )

# Soil bacteria community evenness
emmeans_evenness_bacteria_soil <- emmeans(
  model_evenness_bacteria_soil,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil bacteria evenness"
  )

# Soil bacteria community shannon
emmeans_shannon_bacteria_soil <- emmeans(
  model_shannon_bacteria_soil,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil bacteria shannon"
  )

# Litter bacteria community richness
emmeans_richness_bacteria_litter <- emmeans(
  model_richness_bacteria_litter,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter bacteria richness"
  )

# Litter bacteria community evenness
emmeans_evenness_bacteria_litter <- emmeans(
  model_evenness_bacteria_litter,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter bacteria evenness"
  )

# Litter bacteria community shannon
emmeans_shannon_bacteria_litter <- emmeans(
  model_shannon_bacteria_litter,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter bacteria shannon"
  )

# Soil fungi community richness
emmeans_richness_fungi_soil <- emmeans(
  model_richness_fungi_soil,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil fungi richness"
  )

# Soil fungi community evenness
emmeans_evenness_fungi_soil <- emmeans(
  model_evenness_fungi_soil,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil fungi evenness"
  )

# Soil fungi community shannon
emmeans_shannon_fungi_soil <- emmeans(
  model_shannon_fungi_soil,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil fungi shannon"
  )

# Litter fungi community
emmeans_richness_fungi_litter <- emmeans(
  model_richness_fungi_litter,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter fungi richness"
  )

# Litter fungi community evenness
emmeans_evenness_fungi_litter <- emmeans(
  model_evenness_fungi_litter,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter fungi evenness"
  )

# Litter fungi community shannon
emmeans_shannon_fungi_litter <- emmeans(
  model_shannon_fungi_litter,
  ~ regime * litter_type
) %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter fungi shannon"
  )

## (3c) Effect sizes ####

# Soil bacteria community richness
effectsize_richness_bacteria_soil <- emmeans(
  model_richness_bacteria_soil,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil bacteria richness"
  )

# Soil bacteria community evenness
effectsize_evenness_bacteria_soil <- emmeans(
  model_evenness_bacteria_soil,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil bacteria evenness"
  )

# Soil bacteria community shannon
effectsize_shannon_bacteria_soil <- emmeans(
  model_shannon_bacteria_soil,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil bacteria shannon"
  )

# Litter bacteria community richness
effectsize_richness_bacteria_litter <- emmeans(
  model_richness_bacteria_litter,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter bacteria richness"
  )

#  Litter bacteria community evenness
effectsize_evenness_bacteria_litter <- emmeans(
  model_evenness_bacteria_litter,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter bacteria evenness"
  )

# Litter bacteria community shannon
effectsize_shannon_bacteria_litter <- emmeans(
  model_shannon_bacteria_litter,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter bacteria shannon"
  )

# Soil fungi community richness
effectsize_richness_fungi_soil <- emmeans(
  model_richness_fungi_soil,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil fungi richness"
  )

# Soil fungi community evenness
effectsize_evenness_fungi_soil <- emmeans(
  model_evenness_fungi_soil,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil fungi evenness"
  )

# Soil fungi community shannon
effectsize_shannon_fungi_soil <- emmeans(
  model_shannon_fungi_soil,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Soil fungi shannon"
  )

# Litter fungi community richness
effectsize_richness_fungi_litter <- emmeans(
  model_richness_fungi_litter,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter fungi richness"
  )

# Litter fungi community evenness
effectsize_evenness_fungi_litter <- emmeans(
  model_evenness_fungi_litter,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter fungi evenness"
  )

# Litter fubngi community shannon
effectsize_shannon_fungi_litter <- emmeans(
  model_shannon_fungi_litter,
  ~ regime * litter_type
) %>%
  contrast("eff") %>%
  as.data.frame() %>%
  mutate(
    Response = "Litter fungi shannon"
  )

#### (3d) Model R2 ####

# Soil bacteria community
r2_richness_bacteria_soil <- model_performance(model_richness_bacteria_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil bacteria richness")
r2_evenness_bacteria_soil <- model_performance(model_evenness_bacteria_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil bacteria evenness")
r2_shannon_bacteria_soil <- model_performance(model_shannon_bacteria_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil bacteria shannon")

# Litter bacteria community
r2_richness_bacteria_litter <- model_performance(model_richness_bacteria_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter bacteria richness")
r2_evenness_bacteria_litter <- model_performance(model_evenness_bacteria_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter bacteria evenness")
r2_shannon_bacteria_litter <- model_performance(model_shannon_bacteria_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter bacteria shannon")
# Soil fungi community
r2_richness_fungi_soil <- model_performance(model_richness_fungi_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil fungi richness")
r2_evenness_fungi_soil <- model_performance(model_evenness_fungi_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil fungi evenness")
r2_shannon_fungi_soil <- model_performance(model_shannon_fungi_soil) %>%
  as_tibble() %>%
  mutate(Response = "Soil fungi shannon")
# Litter fungi community
r2_richness_fungi_litter <- model_performance(model_richness_fungi_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter fungi richness")
r2_evenness_fungi_litter <- model_performance(model_evenness_fungi_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter fungi evenness")
r2_shannon_fungi_litter <- model_performance(model_shannon_fungi_litter) %>%
  as_tibble() %>%
  mutate(Response = "Litter fungi shannon")

#### (3e) Join summaries ####

# Model parameters
my_parameters <- bind_rows(
  parameters_richness_bacteria_soil,
  parameters_evenness_bacteria_soil,
  parameters_shannon_bacteria_soil,
  parameters_richness_bacteria_litter,
  parameters_evenness_bacteria_litter,
  parameters_shannon_bacteria_litter,
  parameters_richness_fungi_soil,
  parameters_evenness_fungi_soil,
  parameters_shannon_fungi_soil,
  parameters_richness_fungi_litter,
  parameters_evenness_fungi_litter,
  parameters_shannon_fungi_litter
  ) %>%
  filter(Effects == "fixed") %>%
  mutate(
    Parameter = case_when(
      Parameter == "(Intercept)" ~ "Intercept (Grand Mean)",
      Parameter == "regime1" ~ "Regime: E2 effect",
      Parameter == "regime2" ~ "Regime: L2 effect",
      Parameter == "litter_type1" ~ "Litter type: Grass effect",
      Parameter == "regime1:litter_type1" ~ "Interaction: E2 x Grass",
      Parameter == "regime2:litter_type1" ~ "Interaction: L2 x Grass",
      TRUE ~ Parameter
    ),
    `95% CI` = paste0('[', round(CI_low, 2), ', ', round(CI_high, 2), ']')
    ) %>%
  select(Response, Parameter, Coefficient, `95% CI`, df = df_error, `t-value` = t, `p-value` = p)

# R2 values
my_r2 <- bind_rows(
  r2_richness_bacteria_soil,
  r2_evenness_bacteria_soil,
  r2_shannon_bacteria_soil,
  r2_richness_bacteria_litter,
  r2_evenness_bacteria_litter,
  r2_shannon_bacteria_litter,
  r2_richness_fungi_soil,
  r2_evenness_fungi_soil,
  r2_shannon_fungi_soil,
  r2_richness_fungi_litter,
  r2_evenness_fungi_litter,
  r2_shannon_fungi_litter
) %>%
  select(Response, R2_marginal = R2_marginal)

# Marginal means
my_emmeans <- bind_rows(
  emmeans_richness_bacteria_soil,
  emmeans_evenness_bacteria_soil,
  emmeans_shannon_bacteria_soil,
  emmeans_richness_bacteria_litter,
  emmeans_evenness_bacteria_litter,
  emmeans_shannon_bacteria_litter,
  emmeans_richness_fungi_soil,
  emmeans_evenness_fungi_soil,
  emmeans_shannon_fungi_soil,
  emmeans_richness_fungi_litter,
  emmeans_evenness_fungi_litter,
  emmeans_shannon_fungi_litter
) %>%
  select(Response, regime, litter_type, emmean, lower_CI = `lower.CL`, upper_CI = `upper.CL`)

# Effect sizes
my_effectsize <- bind_rows(
  effectsize_richness_bacteria_soil,
  effectsize_evenness_bacteria_soil,
  effectsize_shannon_bacteria_soil,
  effectsize_richness_bacteria_litter,
  effectsize_evenness_bacteria_litter,
  effectsize_shannon_bacteria_litter,
  effectsize_richness_fungi_soil,
  effectsize_evenness_fungi_soil,
  effectsize_shannon_fungi_soil,
  effectsize_richness_fungi_litter,
  effectsize_evenness_fungi_litter,
  effectsize_shannon_fungi_litter
) %>%
  mutate(
    regime = str_extract(contrast, "E2|L2|U"),
    litter_type = str_extract(contrast, "Grass|Leaf")
  ) %>%
  select(Response, regime, litter_type, estimate)

# (4) Bootstrap models --------------------------------------------------------

# Set up parallel processing
n_cores <- detectCores() - 1
message(paste0("Using ", n_cores, " cores for parallel processing"))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Bootstrap all models --------------------------------------------------------

cat("\n=== Soil Bacteria Richness ===\n")
boot_richness_bacteria_soil <- bootstrap_model(
  model = model_richness_bacteria_soil,
  data = data_bacteria_soil,
  formula = richness_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Soil Bacteria Evenness ===\n")
boot_evenness_bacteria_soil <- bootstrap_model(
  model = model_evenness_bacteria_soil,
  data = data_bacteria_soil,
  formula = evenness_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Soil Bacteria Shannon ===\n")
boot_shannon_bacteria_soil <- bootstrap_model(
  model = model_shannon_bacteria_soil,
  data = data_bacteria_soil,
  formula = shannon_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Litter Bacteria Richness ===\n")
boot_richness_bacteria_litter <- bootstrap_model(
  model = model_richness_bacteria_litter,
  data = data_bacteria_litter,
  formula = richness_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Litter Bacteria Evenness ===\n")
boot_evenness_bacteria_litter <- bootstrap_model(
  model = model_evenness_bacteria_litter,
  data = data_bacteria_litter,
  formula = evenness_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Litter Bacteria Shannon ===\n")
boot_shannon_bacteria_litter <- bootstrap_model(
  model = model_shannon_bacteria_litter,
  data = data_bacteria_litter,
  formula = shannon_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Soil Fungi Richness ===\n")
boot_richness_fungi_soil <- bootstrap_model(
  model = model_richness_fungi_soil,
  data = data_fungi_soil,
  formula = richness_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Soil Fungi Evenness ===\n")
boot_evenness_fungi_soil <- bootstrap_model(
  model = model_evenness_fungi_soil,
  data = data_fungi_soil,
  formula = evenness_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Soil Fungi Shannon ===\n")
boot_shannon_fungi_soil <- bootstrap_model(
  model = model_shannon_fungi_soil,
  data = data_fungi_soil,
  formula = shannon_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Litter Fungi Richness ===\n")
boot_richness_fungi_litter <- bootstrap_model(
  model = model_richness_fungi_litter,
  data = data_fungi_litter,
  formula = richness_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Litter Fungi Evenness ===\n")
boot_evenness_fungi_litter <- bootstrap_model(
  model = model_evenness_fungi_litter,
  data = data_fungi_litter,
  formula = evenness_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

cat("\n=== Litter Fungi Shannon ===\n")
boot_shannon_fungi_litter <- bootstrap_model(
  model = model_shannon_fungi_litter,
  data = data_fungi_litter,
  formula = shannon_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Compile all bootstrap results -----------------------------------------------

# Intercepts
boot_intercepts <- tibble(
  Response = c(
    "Soil bacteria richness", "Soil bacteria evenness", "Soil bacteria shannon",
    "Litter bacteria richness", "Litter bacteria evenness", "Litter bacteria shannon",
    "Soil fungi richness", "Soil fungi evenness", "Soil fungi shannon",
    "Litter fungi richness", "Litter fungi evenness", "Litter fungi shannon"
  ),
  Intercept = c(
    boot_richness_bacteria_soil$intercept,
    boot_evenness_bacteria_soil$intercept,
    boot_shannon_bacteria_soil$intercept,
    boot_richness_bacteria_litter$intercept,
    boot_evenness_bacteria_litter$intercept,
    boot_shannon_bacteria_litter$intercept,
    boot_richness_fungi_soil$intercept,
    boot_evenness_fungi_soil$intercept,
    boot_shannon_fungi_soil$intercept,
    boot_richness_fungi_litter$intercept,
    boot_evenness_fungi_litter$intercept,
    boot_shannon_fungi_litter$intercept
  )
)

# Effect sizes
boot_effects <- bind_rows(
  boot_richness_bacteria_soil$effect_summary %>% mutate(Response = "Soil bacteria richness"),
  boot_evenness_bacteria_soil$effect_summary %>% mutate(Response = "Soil bacteria evenness"),
  boot_shannon_bacteria_soil$effect_summary %>% mutate(Response = "Soil bacteria shannon"),
  boot_richness_bacteria_litter$effect_summary %>% mutate(Response = "Litter bacteria richness"),
  boot_evenness_bacteria_litter$effect_summary %>% mutate(Response = "Litter bacteria evenness"),
  boot_shannon_bacteria_litter$effect_summary %>% mutate(Response = "Litter bacteria shannon"),
  boot_richness_fungi_soil$effect_summary %>% mutate(Response = "Soil fungi richness"),
  boot_evenness_fungi_soil$effect_summary %>% mutate(Response = "Soil fungi evenness"),
  boot_shannon_fungi_soil$effect_summary %>% mutate(Response = "Soil fungi shannon"),
  boot_richness_fungi_litter$effect_summary %>% mutate(Response = "Litter fungi richness"),
  boot_evenness_fungi_litter$effect_summary %>% mutate(Response = "Litter fungi evenness"),
  boot_shannon_fungi_litter$effect_summary %>% mutate(Response = "Litter fungi shannon")
) %>%
  select(Response, regime, litter_type, contrast, effect_mean, lower_ci, upper_ci)%>%
  mutate(
    regime = case_when(
      regime == "E2" ~ "Early",
      regime == "L2" ~ "Late",
      regime == "U" ~ "Unburnt",
      TRUE ~ regime
    ),
    litter_type = ifelse(
      litter_type == "Leaf", "Broadleaf", "Grass"
    )
  )

# Density data for all models
boot_densities <- bind_rows(
  boot_richness_bacteria_soil$density_data %>% mutate(Response = "Soil bacteria richness"),
  boot_evenness_bacteria_soil$density_data %>% mutate(Response = "Soil bacteria evenness"),
  boot_shannon_bacteria_soil$density_data %>% mutate(Response = "Soil bacteria shannon"),
  boot_richness_bacteria_litter$density_data %>% mutate(Response = "Litter bacteria richness"),
  boot_evenness_bacteria_litter$density_data %>% mutate(Response = "Litter bacteria evenness"),
  boot_shannon_bacteria_litter$density_data %>% mutate(Response = "Litter bacteria shannon"),
  boot_richness_fungi_soil$density_data %>% mutate(Response = "Soil fungi richness"),
  boot_evenness_fungi_soil$density_data %>% mutate(Response = "Soil fungi evenness"),
  boot_shannon_fungi_soil$density_data %>% mutate(Response = "Soil fungi shannon"),
  boot_richness_fungi_litter$density_data %>% mutate(Response = "Litter fungi richness"),
  boot_evenness_fungi_litter$density_data %>% mutate(Response = "Litter fungi evenness"),
  boot_shannon_fungi_litter$density_data %>% mutate(Response = "Litter fungi shannon")
) %>%
  mutate(
    regime = case_when(
      regime == "E2" ~ "Early",
      regime == "L2" ~ "Late",
      regime == "U" ~ "Unburnt",
      TRUE ~ regime
    ),
    litter_type = ifelse(
      litter_type == "Leaf", "Broadleaf", "Grass"
    )
  )

# Print summaries
cat("\n=== Bootstrap Intercepts ===\n")
print(boot_intercepts)

cat("\n=== Bootstrap Effect Sizes ===\n")
print(boot_effects)

# Stop cluster when done
stopCluster(cl)
message("Parallel processing complete!")

# (5) Save the data ------------------------------------------------------------

#### (5a) Save the data for figures ####

# Organise the raw data for plotting
raw_data <- bind_rows(
  data_bacteria_soil %>%
    mutate(Response = "Soil bacteria richness") %>% 
    select(Response, observed_value = richness_bacteria, regime, litter_type),
  data_bacteria_soil %>%
    mutate(Response = "Soil bacteria evenness") %>%
    select(Response, observed_value = evenness_bacteria, regime, litter_type),
  data_bacteria_soil %>%
    mutate(Response = "Soil bacteria shannon") %>%
    select(Response, observed_value = shannon_bacteria, regime, litter_type),
  data_bacteria_litter %>%
    mutate(Response = "Litter bacteria richness") %>%
    select(Response, observed_value = richness_bacteria, regime, litter_type),
  data_bacteria_litter %>%
    mutate(Response = "Litter bacteria evenness") %>%
    select(Response, observed_value = evenness_bacteria, regime, litter_type),
  data_bacteria_litter %>%
    mutate(Response = "Litter bacteria shannon") %>%
    select(Response, observed_value = shannon_bacteria, regime, litter_type),
  data_fungi_soil %>%
    mutate(Response = "Soil fungi richness") %>%
    select(Response, observed_value = richness_fungi, regime, litter_type),
  data_fungi_soil %>%
    mutate(Response = "Soil fungi evenness") %>%
    select(Response, observed_value = evenness_fungi, regime, litter_type),
  data_fungi_soil %>%
    mutate(Response = "Soil fungi shannon") %>%
    select(Response, observed_value = shannon_fungi, regime, litter_type),
  data_fungi_litter %>%
    mutate(Response = "Litter fungi richness") %>%
    select(Response, observed_value = richness_fungi, regime, litter_type),
  data_fungi_litter %>%
    mutate(Response = "Litter fungi evenness") %>%
    select(Response, observed_value = evenness_fungi, regime, litter_type),
  data_fungi_litter %>%
    mutate(Response = "Litter fungi shannon") %>%
    select(Response, observed_value = shannon_fungi, regime, litter_type)
  ) %>%
  mutate(
    regime = case_when(
      regime == "U" ~ "Unburnt",
      regime == "E2" ~ "Early",
      regime == "L2" ~ "Late",
      TRUE ~ regime
    ),
    litter_type = case_when(
      litter_type == "Leaf" ~ "Broadleaf",
      litter_type == "1" ~ "Grass",
      TRUE ~ litter_type
    )
  ) %>%
  # Level the factors
  mutate(
    regime = factor(regime, levels = c("Early", "Late", "Unburnt")),
    litter_type = factor(litter_type, levels = c("Broadleaf", "Grass"))
  )

# Organise the emmean data for plotting
emmean <- my_emmeans %>%
  mutate(
    regime = case_when(
      regime == "U" ~ "Unburnt",
      regime == "E2" ~ "Early",
      regime == "L2" ~ "Late",
      TRUE ~ regime
    ),
    litter_type = case_when(
      litter_type == "Leaf" ~ "Broadleaf",
      litter_type == "1" ~ "Grass",
      TRUE ~ litter_type
    )
  ) %>%
  # Level the factors
  mutate(
    regime = factor(regime, levels = c("Early", "Late", "Unburnt")),
    litter_type = factor(litter_type, levels = c("Broadleaf", "Grass"))
  )

# Model summaries
model_summaries <- inner_join(
  my_parameters,
  my_r2,
  by = "Response"
)

# Create data and output folders if they do not exist
if (!dir.exists("data")) {
  dir.create("data")
}
if (!dir.exists("output")) {
  dir.create("output")
}

# Save the data
save(
  raw_data,
  emmean,
  boot_densities,
  boot_effects,
  boot_intercepts,
  file = "data/alpha_diversity_plot_data.RData"
)
data.table::fwrite(
  model_summaries,
  file = "output/alpha_diversity_model_summaries.txt",
  sep = "\t"
)
