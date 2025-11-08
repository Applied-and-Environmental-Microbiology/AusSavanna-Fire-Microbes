
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
source("code/bootstrap_effects.R")

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

# (4) Bootstrap effect sizes and intercepts ------------------------------------

# Set up parallel processing
n_cores <- detectCores() - 1
message(paste0("Using ", n_cores, " cores for parallel processing"))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

#### (4a) Bootstrap models ####

# Soil bacterial richness
boot_richness_bacteria_soil <- bootstrap_effects(
  model = model_richness_bacteria_soil,
  data = data_bacteria_soil,
  formula = richness_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Soil bacterial evenness
boot_evenness_bacteria_soil <- bootstrap_effects(
  model = model_evenness_bacteria_soil,
  data = data_bacteria_soil,
  formula = evenness_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Soil bacterial Shannon
boot_shannon_bacteria_soil <- bootstrap_effects(
  model = model_shannon_bacteria_soil,
  data = data_bacteria_soil,
  formula = shannon_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Litter bacterial richness
boot_richness_bacteria_litter <- bootstrap_effects(
  model = model_richness_bacteria_litter,
  data = data_bacteria_litter,
  formula = richness_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Litter bacterial evenness
boot_evenness_bacteria_litter <- bootstrap_effects(
  model = model_evenness_bacteria_litter,
  data = data_bacteria_litter,
  formula = evenness_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Litter bacterial Shannon
boot_shannon_bacteria_litter <- bootstrap_effects(
  model = model_shannon_bacteria_litter,
  data = data_bacteria_litter,
  formula = shannon_bacteria ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Soil fungi richness
boot_richness_fungi_soil <- bootstrap_effects(
  model = model_richness_fungi_soil,
  data = data_fungi_soil,
  formula = richness_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Soil fungi evenness
boot_evenness_fungi_soil <- bootstrap_effects(
  model = model_evenness_fungi_soil,
  data = data_fungi_soil,
  formula = evenness_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Soil fungi Shannon
boot_shannon_fungi_soil <- bootstrap_effects(
  model = model_shannon_fungi_soil,
  data = data_fungi_soil,
  formula = shannon_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Litter fungi richness
boot_richness_fungi_litter <- bootstrap_effects(
  model = model_richness_fungi_litter,
  data = data_fungi_litter,
  formula = richness_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Litter fungi evenness
boot_evenness_fungi_litter <- bootstrap_effects(
  model = model_evenness_fungi_litter,
  data = data_fungi_litter,
  formula = evenness_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

# Litter fungi Shannon
boot_shannon_fungi_litter <- bootstrap_effects(
  model = model_shannon_fungi_litter,
  data = data_fungi_litter,
  formula = shannon_fungi ~ regime * litter_type + (1 | block),
  n_boot = 10000
)

#### (4b) Compile results ####

# Intercepts
boot_intercepts <- tibble(
  Response = c(
    "Soil bacteria richness", "Soil bacteria evenness", "Soil bacteria shannon",
    "Litter bacteria richness", "Litter bacteria evenness", "Litter bacteria shannon",
    "Soil fungi richness", "Soil fungi evenness", "Soil fungi shannon",
    "Litter fungi richness", "Litter fungi evenness", "Litter fungi shannon"
  ),
  Intercept = c(
    sprintf(
      "&beta;<sub>0</sub> = %.0f [%.0f, %.0f]",
      mean(boot_richness_bacteria_soil$intercept_distribution),
      ci(boot_richness_bacteria_soil$intercept_distribution)[1],
      ci(boot_richness_bacteria_soil$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.2f [%.2f, %.2f]",
      mean(boot_evenness_bacteria_soil$intercept_distribution),
      ci(boot_evenness_bacteria_soil$intercept_distribution)[1],
      ci(boot_evenness_bacteria_soil$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.2f [%.2f, %.2f]",
      mean(boot_shannon_bacteria_soil$intercept_distribution),
      ci(boot_shannon_bacteria_soil$intercept_distribution)[1],
      ci(boot_shannon_bacteria_soil$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.0f [%.0f, %.0f]",
      mean(boot_richness_bacteria_litter$intercept_distribution),
      ci(boot_richness_bacteria_litter$intercept_distribution)[1],
      ci(boot_richness_bacteria_litter$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.2f [%.2f, %.2f]",
      mean(boot_evenness_bacteria_litter$intercept_distribution),
      ci(boot_evenness_bacteria_litter$intercept_distribution)[1],
      ci(boot_evenness_bacteria_litter$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.2f [%.2f, %.2f]",
      mean(boot_shannon_bacteria_litter$intercept_distribution),
      ci(boot_shannon_bacteria_litter$intercept_distribution)[1],
      ci(boot_shannon_bacteria_litter$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.0f [%.0f, %.0f]",
      mean(boot_richness_fungi_soil$intercept_distribution),
      ci(boot_richness_fungi_soil$intercept_distribution)[1],
      ci(boot_richness_fungi_soil$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.2f [%.2f, %.2f]",
      mean(boot_evenness_fungi_soil$intercept_distribution),
      ci(boot_evenness_fungi_soil$intercept_distribution)[1],
      ci(boot_evenness_fungi_soil$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.2f [%.2f, %.2f]",
      mean(boot_shannon_fungi_soil$intercept_distribution),
      ci(boot_shannon_fungi_soil$intercept_distribution)[1],
      ci(boot_shannon_fungi_soil$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.0f [%.0f, %.0f]",
      mean(boot_richness_fungi_litter$intercept_distribution),
      ci(boot_richness_fungi_litter$intercept_distribution)[1],
      ci(boot_richness_fungi_litter$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.2f [%.2f, %.2f]",
      mean(boot_evenness_fungi_litter$intercept_distribution),
      ci(boot_evenness_fungi_litter$intercept_distribution)[1],
      ci(boot_evenness_fungi_litter$intercept_distribution)[2]
    ),
    sprintf(
      "&beta;<sub>0</sub> = %.2f [%.2f, %.2f]",
      mean(boot_shannon_fungi_litter$intercept_distribution),
      ci(boot_shannon_fungi_litter$intercept_distribution)[1],
      ci(boot_shannon_fungi_litter$intercept_distribution)[2]
    )
  )
)

# Effect sizes
effect_mean <- bind_rows(
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
  select(Response, regime, litter_type, contrast, effect_mean, lower_ci, upper_ci) %>%
  mutate(
    regime = case_when(
      regime == "E2" ~ "Early",
      regime == "L2" ~ "Late",
      regime == "U" ~ "Unburnt",
      TRUE ~ regime
    ),
    litter_type = case_when(
      litter_type == "Leaf" ~ "Broadleaf",
      TRUE ~ litter_type
    )
  ) %>%
  # Level the factors
  mutate(
    regime = factor(regime, levels = c("Early", "Late", "Unburnt")),
    litter_type = factor(litter_type, levels = c("Broadleaf", "Grass"))
  )

# Density data for all models
effect_density <- bind_rows(
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
    litter_type = case_when(
      litter_type == "Leaf" ~ "Broadleaf",
      TRUE ~ litter_type
    )
  ) %>%
  # Level the factors
  mutate(
    regime = factor(regime, levels = c("Early", "Late", "Unburnt")),
    litter_type = factor(litter_type, levels = c("Broadleaf", "Grass"))
  )

# Print summaries
print(boot_intercepts)
print(boot_effects)

# Stop cluster when done
stopCluster(cl)


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

# Save the data
save(
  raw_data,
  emmean,
  file = "data/alpha_diversity_plot_data.RData"
)
fwrite(
  model_summaries,
  file = "output/alpha_diversity_model_summaries.txt",
  sep = "\t"
)
