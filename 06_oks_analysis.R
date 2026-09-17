################################################################################
##
## 06_oks_analysis.R
##
## Empirical analysis: longitudinal measurement invariance in the Oxford Knee
## Score (OKS), pre- and 1-year post-total knee arthroplasty (TKA).
##
## Analyses:
##   1. Descriptive summary of the analytic cohort
##   2. Baseline GRM fit and fit indices (RMSEA, CFI, SRMR)
##   3. GRMTree with time as the sole partitioning covariate
##   4. Two-factor longitudinal GRM + LRT (global and item-level with Bonferroni
##      correction)
##
## Input:  cleaned_oks_data.xlsx  (not shared publicly; contact registry)
## Output: printed tables + saved fitted models as .rds
##
## Data note. The cleaned OKS dataset is derived from the Winnipeg Regional
## Health Authority (WRHA) Joint Replacement Registry and cannot be shared with
## this repository due to data-use restrictions. Investigators seeking access
## should contact the registry directly.
##
################################################################################

suppressPackageStartupMessages({
  library(grmtree)
  library(mirt)
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(EnvStats)
})

################################################################################
## 1. Load cleaned data
##
## Expected columns:
##   - Q1_Pre  ... Q12_Pre   (OKS items pre-surgery, 0-4)
##   - Q1_Post ... Q12_Post  (OKS items post-surgery, 0-4)
##   - Patient covariates (age, Sex, BMI, ComorbidityCount, ...)
##
## Scoring direction: OKS is scored 0 (worst function) to 48 (best function).
## Item-level responses in this file follow the same direction: 0 = worst,
## 4 = best.
################################################################################

oks_data <- read_xlsx("cleaned_oks_data.xlsx")

cat("Loaded", nrow(oks_data), "patients.\n")
cat("Column names:\n")
print(names(oks_data))

################################################################################
## 2. Descriptive summary of the analytic cohort
################################################################################

cat("\n============================================================\n")
cat("DESCRIPTIVE SUMMARY (Table 3 and supplementary file in manuscript)\n")
cat("============================================================\n")

# Pre-/post-OKS totals
oks_data <- oks_data %>%
  mutate(
    OKS_Pre  = rowSums(across(Q1_Pre:Q12_Pre)),
    OKS_Post = rowSums(across(Q1_Post:Q12_Post))
  )

cat("\nBaseline OKS total:\n")
print(summary(oks_data$OKS_Pre))
cat("Median (IQR) =", median(oks_data$OKS_Pre, na.rm = TRUE),
    "(", IQR(oks_data$OKS_Pre, na.rm = TRUE), ")\n")

cat("\nPost-surgery OKS total:\n")
print(summary(oks_data$OKS_Post))
cat("Median (IQR) =", median(oks_data$OKS_Post, na.rm = TRUE),
    "(", IQR(oks_data$OKS_Post, na.rm = TRUE), ")\n")

cat("\nAge (years):\n")
cat("Median (IQR) =", median(oks_data$age, na.rm = TRUE),
    "(", IQR(oks_data$age, na.rm = TRUE), ")\n")

cat("\nBMI:\n")
cat("Median (IQR) =", median(oks_data$BMI, na.rm = TRUE),
    "(", IQR(oks_data$BMI, na.rm = TRUE), ")\n")

cat("\nSex:\n")
print(table(oks_data$Sex, useNA = "ifany"))


## Item summary (Table A5)
oks_data %>% 
  group_by(Q1_Pre) %>% 
  summarise(n = n()) %>% 
  mutate(percent = round(100 * n/sum(n), 1))


################################################################################
## 3. Baseline GRM fit and fit indices
##
## We fit a unidimensional GRM to the pre-surgery OKS data and evaluate model
## fit using RMSEA, CFI, and SRMR (M2 procedure). Benchmarks: RMSEA < 0.08,
## CFI > 0.90, SRMR < 0.08.
################################################################################

cat("\n============================================================\n")
cat("BASELINE GRM FIT (pre-surgery)\n")
cat("============================================================\n")

resp_pre <- oks_data %>%
  select(Q1_Pre:Q12_Pre) %>%
  drop_na() %>%
  as.matrix()

oks_grm_pre <- mirt(
  data     = resp_pre,
  model    = 1,               # unidimensional
  itemtype = "graded",
  SE       = TRUE,
  verbose  = FALSE
)

cat("\nGRM item parameters (baseline):\n")
print(coef(oks_grm_pre, IRTpars = TRUE, simplify = TRUE))

cat("\nM2-based fit statistics (RMSEA, CFI, SRMR):\n")
print(M2(oks_grm_pre, type = "C2"))

################################################################################
## 4. GRMTree with time as sole partitioning covariate (main analysis)
##
## The wide-format data are reshaped into person-period long format. Each
## individual contributes two rows (one for pre-surgery, one for post-surgery).
## The tree is fit with time as the only candidate splitting variable; a
## split on time indicates longitudinal measurement non-invariance.
##
## The minimum terminal node size is set to 600 -- 10 times the number of item
## parameters (12 items x 5 params = 60; 10 x 60 = 600) -- following the rule
## of thumb for stable GRM estimation at each terminal node.
################################################################################

cat("\n============================================================\n")
cat("GRMTREE ANALYSIS (time as sole covariate)\n")
cat("============================================================\n")

# Reshape to person-period format
n_patients <- nrow(oks_data)
item_names <- paste0("Q", 1:12, "_oks")

df_pre <- oks_data %>%
  select(Q1_Pre:Q12_Pre) %>%
  setNames(item_names) %>%
  mutate(id = 1:n_patients, time = "Pre-surgery")

df_post <- oks_data %>%
  select(Q1_Post:Q12_Post) %>%
  setNames(item_names) %>%
  mutate(id = 1:n_patients, time = "Post-surgery")

df_long <- bind_rows(df_pre, df_post) %>%
  mutate(time = factor(time, levels = c("Pre-surgery", "Post-surgery"))) %>%
  drop_na()

##---------------------------------------------------------------------------##
## Create a plot comparing the total scores for two times
##---------------------------------------------------------------------------##

## Create smmary scores at baseline and 1 year
oks_data_sum <- df_long %>% 
  rowwise() %>% 
  dplyr::mutate(total_oks = sum(c_across(Q1_oks:Q12_oks))) 

## Summary scores
oks_data_sum %>% 
  group_by(time) %>% 
  summarise(median = median(total_oks),
            IQ_oks = IQR(total_oks))

# Calculating the mean of each group to annotate
mean_data <- oks_data_sum %>% 
  group_by(time) %>% 
  summarise(mean_cond = round(mean(total_oks),2)) %>% 
  .$mean_cond

# Visualisation
ggplot(oks_data_yr_sum,
       aes(x = time, y = total_oks)) + # Setup x and y axis
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(width = 0.25) + # Draw default boxplot first
  
  geom_point(aes(color = time), alpha = 0.8, size = 1) + # Adding points and coloring them
  
  #geom_line(aes(group = Patient, col = slope), alpha = 0.8, linetype = 2) + # Connection lines
  
  # Drawing two dots which are group means and a connection line  
  stat_summary(fun=mean, color = "#880000", 
               geom="line",linetype = 1, aes(group=1), 
               size = 1, alpha = 0.5) + 
  stat_summary(fun=mean, color = "#880000", 
               geom="point", size = 2, alpha = 1, aes(group=1)) + 
  
  # Control dot colors and slope colors  
  # Group 1, Group 2 - dots and lines in that order
  #scale_colour_manual(values = c("#007020", "#bc7a00","#bb6688","black")) + 
  
  # Add n number above x axis
  stat_n_text(size= 3.5,
              color = "steelblue",
              text.box = TRUE) + 
  
  # Labels
  labs(x = "Time", y = "Oxford knee score") + 
  
  # Theming
  theme_bw() + 
  
  # Annotate means. Control the location using x and y
  annotate("text", label = paste("mu == ",mean_data[1]),
           size = 3.25, x = 0.75, y = mean_data[1], parse =TRUE) + 
  annotate("text", label = paste("mu == ",mean_data[2]),
           size = 3.25, x = 2.25, y = mean_data[2], parse = TRUE) + 
  
  # Remove all legends for clean look
  theme(legend.position = "none")


# Prepare data as a matrix
df_long$resp <- as.matrix(df_long[, item_names])

cat("Person-period dataset:", nrow(df_long), "rows (",
    nrow(df_long) / 2, "patients).\n")

# Fit the GRMTree
oks_tree <- grmtree(
  resp ~ time,
  data    = df_long,
  control = grmtree.control(p_adjust = "bonferroni", minbucket = 600)
)

cat("\nGRMTree summary:\n")
print(oks_tree)
plot(oks_tree)

# Extract item parameters at each terminal node (Table 4 in manuscript)
cat("\nItem parameters by terminal node (Table 4):\n")
cat("\nThresholds:\n")
print(threshpar_grmtree(oks_tree))

cat("\nAverage thresholds (b_bar):\n")
print(itempar_grmtree(oks_tree))

cat("\nDiscriminations (a):\n")
print(discrpar_grmtree(oks_tree))

# Save the fitted tree
saveRDS(oks_tree, file = "oks_grmtree_fit.rds")

################################################################################
## 5. Two-factor longitudinal GRM + LRT
##
## Global test: constrained model (item parameters equal across time) vs.
## unconstrained model (freely estimated at each time). Both models estimate
## the T2 latent mean and the T1-T2 covariance freely, so impact is absorbed
## by the structural parameters rather than by the item parameters.
##
## Item-level test: release the constraints for one item at a time and
## compare to the fully constrained model. Bonferroni correction over 12 items.
################################################################################

cat("\n============================================================\n")
cat("LONGITUDINAL GRM + LRT (global and item-level)\n")
cat("============================================================\n")

# Wide-format response matrix: Q1_Pre..Q12_Pre, Q1_Post..Q12_Post
oks_wide <- oks_data %>%
  select(Q1_Pre:Q12_Pre, Q1_Post:Q12_Post) %>%
  drop_na() %>%
  as.matrix()

n_items <- 12

# --- Fully constrained model (invariance)

# Build the constraint string programmatically
# For discrimination: constrain a1 of item i to a2 of item i+12
discr_constraints <- paste0("(", 1:n_items, ", a1, ",
                            (n_items + 1):(2 * n_items), ", a2)",
                            collapse = ", ")

# For thresholds: constrain d1-d4 across time (same parameter name, different items)
thresh_parts <- c()
for (d in 1:4) {
  thresh_parts <- c(thresh_parts,
                    paste0("(", 1:n_items, ",",
                           (n_items + 1):(2 * n_items), ", d", d, ")",
                           collapse = ", "))
}

thresh_constraints <- paste(thresh_parts, collapse = ", ")

# Combine all constraints
all_constraints    <- paste(discr_constraints, thresh_constraints, sep = ", ")

# Build the full model specification
model_constrained_spec <- paste0("
  Theta_T1 = 1-", n_items, "
  Theta_T2 = ", n_items + 1, "-", 2 * n_items, "
  COV      = Theta_T1*Theta_T2
  MEAN     = Theta_T2
  CONSTRAIN = ", all_constraints)

model_unconstrained_spec <- paste0("
  Theta_T1 = 1-", n_items, "
  Theta_T2 = ", n_items + 1, "-", 2 * n_items, "
  COV      = Theta_T1*Theta_T2
  MEAN     = Theta_T2")

# Print to verify
cat(model_constrained_spec)
cat(model_unconstrained_spec)

cat("\nFitting constrained model...\n")
mod_constrained <- mirt(
  data     = oks_wide,
  model    = mirt.model(model_constrained_spec),
  itemtype = "graded",
  SE       = FALSE,
  verbose  = TRUE
)

cat("\nFitting unconstrained model...\n")
mod_unconstrained <- mirt(
  data     = oks_wide,
  model    = mirt.model(model_unconstrained_spec),
  itemtype = "graded",
  SE       = FALSE,
  verbose  = TRUE
)

cat("\nGlobal LRT (constrained vs. unconstrained):\n")
print(anova(mod_constrained, mod_unconstrained))

# --- Item-level LRT (release one item's constraints at a time)

# Function to build constraint string EXCLUDING a specific item
build_constraints_exclude_item <- function(exclude_item, n_items = 12) {

  items_to_constrain <- setdiff(1:n_items, exclude_item)
  discr <- paste0("(", items_to_constrain, ", a1, ",
                  items_to_constrain + n_items, ", a2)", collapse = ", ")
  parts <- c()
  for (d in 1:4) {
    parts <- c(parts,
               paste0("(", items_to_constrain, ",",
                      items_to_constrain + n_items, ", d", d, ")",
                      collapse = ", "))
  }
  thresh <- paste(parts, collapse = ", ")
  paste(discr, thresh, sep = ", ")
}

# Function to test RS for a specific item
test_item_rs <- function(item_num, baseline_model, data, n_items = 12) {
  
  # Build constraints excluding this item
  partial_constraints <- build_constraints_exclude_item(item_num, n_items)

  # Build model spec with partial constraints
  spec_partial <- paste0("
    Theta_T1 = 1-", n_items, "
    Theta_T2 = ", n_items + 1, "-", 2 * n_items, "
    COV      = Theta_T1*Theta_T2
    MEAN     = Theta_T2
    CONSTRAIN = ", partial_constraints)

  # Fit model with this item's parameters free
  mod_free_item <- mirt(
    data     = data,
    model    = mirt.model(spec_partial),
    itemtype = "graded",
    SE       = FALSE,
    verbose  = FALSE
  )

  # Compare to fully constrained model
  comparison <- anova(baseline_model, mod_free_item)

  data.frame(
    item      = item_num,
    item_name = paste0("Q", item_num),
    X2        = comparison$X2[2],
    df        = comparison$df[2],
    p         = comparison$p[2]
  )
}

# Test each item (this will take a few minutes)
cat("\nItem-level LRT (releasing one item at a time)...\n")
item_results <- list()
for (i in 1:n_items) {
  cat("  Testing item", i, "of", n_items, "...\n")
  item_results[[i]] <- test_item_rs(i, mod_constrained, oks_wide, n_items)
}

item_lrt_table <- do.call(rbind, item_results)
item_lrt_table$p_adj       <- p.adjust(item_lrt_table$p, method = "bonferroni")
item_lrt_table$significant <- ifelse(item_lrt_table$p_adj < 0.05, "Yes", "No")

cat("\nItem-level LRT results (Bonferroni-corrected) -- Supplementary Table A5:\n")
print(item_lrt_table, row.names = FALSE)

# Save
saveRDS(list(
  global_constrained   = mod_constrained,
  global_unconstrained = mod_unconstrained,
  item_lrt_table       = item_lrt_table
), file = "oks_lrt_results.rds")

#------------------------------------------------------------------
# TEST RESPONSE SHIFT BY TYPE
# Uniform RS = thresholds differ (d parameters)
# Non-uniform RS = discrimination differs (a parameters)
#------------------------------------------------------------------

# 1. TEST UNIFORM RS (thresholds free, discriminations constrained)
cat("\n=== Testing UNIFORM Response Shift ===\n")

# Only constrain discriminations
discr_only_constraints <- paste0("(", 1:12, ", a1, ", 13:24, ", a2)", collapse = ", ")

model_uniform_rs_spec <- paste0('
  Theta_T1 = 1-12
  Theta_T2 = 13-24
  COV = Theta_T1*Theta_T2
  MEAN = Theta_T2
  CONSTRAIN = ', discr_only_constraints
)

mod_uniform_rs <- mirt(
  data = oks_wide,
  model = mirt.model(model_uniform_rs_spec),
  itemtype = 'graded',
  SE = TRUE,
  verbose = TRUE
)

# Compare: Is there uniform RS?
cat("\nUniform RS Test (constrained vs. thresholds-free):\n")
uniform_test <- anova(mod_constrained, mod_uniform_rs)
print(uniform_test)

#------------------------------------------------------------------
# 2. TEST NON-UNIFORM RS (discriminations free, thresholds constrained)
cat("\n=== Testing NON-UNIFORM Response Shift ===\n")

# Only constrain thresholds
thresh_parts <- c()
for (d in 1:4) {
  thresh_parts <- c(thresh_parts,
                    paste0("(", 1:12, ",", 13:24, ", d", d, ")", collapse = ", "))
}
thresh_only_constraints <- paste(thresh_parts, collapse = ", ")

model_nonuniform_rs_spec <- paste0('
  Theta_T1 = 1-12
  Theta_T2 = 13-24
  COV = Theta_T1*Theta_T2
  MEAN = Theta_T2
  CONSTRAIN = ', thresh_only_constraints
)

mod_nonuniform_rs <- mirt(
  data = oks_wide,
  model = mirt.model(model_nonuniform_rs_spec),
  itemtype = 'graded',
  SE = TRUE,
  verbose = TRUE
)

# Compare: Is there non-uniform RS?
cat("\nNon-uniform RS Test (constrained vs. discriminations-free):\n")
nonuniform_test <- anova(mod_constrained, mod_nonuniform_rs)
print(nonuniform_test)

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Saved fitted models to: oks_grmtree_fit.rds, oks_lrt_results.rds\n")
