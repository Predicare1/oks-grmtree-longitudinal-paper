################################################################################
##
## 01_helpers.R
##
## Shared helper functions for the GRMTree longitudinal simulation study.
## Sourced by 02_grmtree_type1.R, 03_grmtree_power.R,
## 04_lrt_type1.R, and 05_lrt_power.R.
##
## Authors: O. I. Arimoro, T. T. Sajobi, E. Bohm, L. M. Lix
## Contact: olayinka.arimoro@ucalgary.ca
##
################################################################################

suppressPackageStartupMessages({
  library(grmtree)
  library(mirt)
  library(tidyverse)
})

################################################################################
## Baseline item parameters
##
## For each scale length (6 or 12 items), returns the T1 discrimination and
## threshold parameters, and the indices of items to receive non-invariance.
##
## Discriminations are in the range 1.1-1.9 (typical of PROM applications).
## Thresholds span approximately -2 to +2 so each response category has
## non-trivial probability across plausible latent trait values.
################################################################################

get_item_params <- function(n_items) {

  if (n_items == 6) {
    # Discrimination parameters (a): between 0.8 and 2.0
    # These stay fixed across all replications
    a_params <- c(1.2, 1.5, 1.1, 1.8, 1.3, 1.6)
    
    # Threshold/difficulty parameters (d) for mirt::simdata()
    # For graded response: need n_categories - 1 = 4 thresholds per item
    # Thresholds must be DESCENDING within each item (mirt convention)
    # These are the "d" parameters in mirt's parameterization
    d_params <- matrix(c(
      2.0,  0.7, -0.7, -2.0,
      1.8,  0.5, -0.5, -1.8,
      2.2,  0.8, -0.6, -2.1,
      1.9,  0.6, -0.8, -1.9,
      2.1,  0.7, -0.5, -2.0,
      1.7,  0.4, -0.6, -1.7
    ), nrow = n_items, ncol = 4, byrow = TRUE)

    dif_items <- c(1, 3, 6)

  } else if (n_items == 12) {

    a_params <- c(1.2, 1.5, 1.1, 1.8, 1.3, 1.6,
                  1.4, 1.7, 1.2, 1.5, 1.3, 1.9)

    d_params <- matrix(c(
      2.0,  0.7, -0.7, -2.0,
      1.8,  0.5, -0.5, -1.8,
      2.2,  0.8, -0.6, -2.1,
      1.9,  0.6, -0.8, -1.9,
      2.1,  0.7, -0.5, -2.0,
      1.7,  0.4, -0.6, -1.7,
      2.0,  0.6, -0.7, -1.8,
      1.8,  0.5, -0.5, -1.9,
      2.2,  0.8, -0.6, -2.2,
      1.9,  0.6, -0.7, -1.8,
      2.0,  0.7, -0.8, -2.0,
      1.8,  0.5, -0.6, -1.7
    ), nrow = n_items, ncol = 4, byrow = TRUE)

    dif_items <- c(1, 6, 12)

  } else {
    stop("n_items must be 6 or 12.")
  }

  list(a_params = a_params, d_params = d_params, dif_items = dif_items)
}

################################################################################
## Create T2 item parameters with non-invariance applied
##
## For the specified items, applies a constant shift Delta to the discrimination
## (non-uniform), the thresholds (uniform), or both.
################################################################################

create_t2_params <- function(a_params_t1, d_params_t1, dif_items,
                             dif_type = c("nonuniform", "uniform", "both"),
                             dif_magnitude) {
  
  #' Create T2 item parameters with DIF applied to specified items
  #' 
  
  #' @param a_params_t1 Discrimination parameters at T1
  #' @param d_params_t1 Threshold parameters at T1 (matrix)
  #' @param dif_items Vector of item indices with DIF
  #' @param dif_type "nonuniform", "uniform", or "both"
  #' @param dif_magnitude Size of DIF effect
  #' @return List with a_params_t2 and d_params_t2
  
  dif_type <- match.arg(dif_type)

  # Start with T1 parameter
  a_params_t2 <- a_params_t1
  d_params_t2 <- d_params_t1

  # Apply DIF to specified items
  for (item in dif_items) {

    if (dif_type %in% c("nonuniform", "both")) {
      # Non-uniform DIF: increase discrimination at T2
      a_params_t2[item] <- a_params_t1[item] + dif_magnitude
    }
    if (dif_type %in% c("uniform", "both")) {
      # Uniform DIF: shift all thresholds down at T2
      # (makes item "easier" at T2 - higher probability of higher responses)
      d_params_t2[item, ] <- d_params_t1[item, ] + dif_magnitude
    }
  }

  list(a_params_t2 = a_params_t2, d_params_t2 = d_params_t2)
}

################################################################################
## Generate correlated latent trait values across two occasions
##
## theta_T2 = rho * theta_T1 + sqrt(1 - rho^2) * eps + mean_t2
##
## Ensures Var(theta_T2) = 1 and Cor(theta_T1, theta_T2) = rho.
################################################################################

generate_correlated_theta <- function(n, rho, mean_t1 = 0, mean_t2 = 0) {
  #' Generate correlated theta values for two time points
  
  #' Uses the formula: theta_t2 = rho * theta_t1 + sqrt(1-rho^2) * error
  #' This ensures theta_t2 has variance = 1 and cor(theta_t1, theta_t2) = rho
  #' 
  #' @param n Sample size (number of individuals)
  #' @param rho Correlation between theta at T1 and T2
  #' @param mean_t1 Mean of theta at time 1 (default = 0)
  #' @param mean_t2 Mean of theta at time 2 (default = 0)
  #' @return Data frame with theta_t1 and theta_t2
  
  # Generate theta at T1
  theta_t1 <- rnorm(n, mean = mean_t1, sd = 1)
  
  # Generate theta at T2 with specified correlation
  # This formula ensures: Var(theta_t2) = rho^2 * 1 + (1-rho^2) * 1 = 1
  # and Cor(theta_t1, theta_t2) = rho
  theta_t2 <- rho * theta_t1 + sqrt(1 - rho^2) * rnorm(n, mean = 0, sd = 1) + mean_t2

  data.frame(theta_t1 = theta_t1, theta_t2 = theta_t2)
}

################################################################################
## Reshape T1 and T2 response matrices into person-period long format
##
## Each individual contributes two rows: one for T1, one for T2. A binary
## factor `time` distinguishes the two occasions and serves as the sole
## partitioning covariate in the GRMTree fit.
################################################################################

create_person_period_data <- function(resp_t1, resp_t2, n_items) {
 
  #' Stack T1 and T2 responses into person-period format
  #' 
  #' @param resp_t1 Response MATRIX at time 1 (n x n_items) from simdata()
  #' @param resp_t2 Response MATRIX at time 2 (n x n_items) from simdata()
  #' @param n_items Number of items
  #' @return Data frame in person-period format with time indicator and resp matrix
  
  n <- nrow(resp_t1)
  
  # Create item names
  item_names <- paste0("Item_", 1:n_items)

  # Convert matrices to data frames with proper names
  df_t1 <- as.data.frame(resp_t1)
  names(df_t1) <- item_names
  df_t1$id   <- 1:n
  df_t1$time <- "T1"

  df_t2 <- as.data.frame(resp_t2)
  names(df_t2) <- item_names
  df_t2$id   <- 1:n
  df_t2$time <- "T2"

  # Stack into person-period format
  df_long <- rbind(df_t1, df_t2)
  
  # Convert time to factor (T1 as reference)
  df_long$time <- factor(df_long$time, levels = c("T1", "T2"))
  
  # Create response matrix for grmtree
  df_long$resp <- as.matrix(sapply(df_long[, item_names], as.numeric))

  df_long
}

################################################################################
## Build mirt model specifications for the longitudinal GRM-based LRT
##
## The constrained model holds all item parameters (discriminations and
## thresholds) equal across the two occasions. The unconstrained model
## estimates them freely at each occasion. In both models, the T2 latent mean
## and the T1-T2 covariance are freely estimated (i.e., impact is absorbed by
## the structural parameters rather than by the item parameters).
################################################################################

build_lrt_model_specs <- function(n_items) {

  t1_items <- 1:n_items
  t2_items <- (n_items + 1):(2 * n_items)

  # Cross-time equality constraints on discriminations
  discr_constraints <- paste0("(", t1_items, ", a1, ", t2_items, ", a2)",
                              collapse = ", ")

  # Cross-time equality constraints on thresholds (d1-d4)
  thresh_parts <- c()
  for (d in 1:4) {
    thresh_parts <- c(thresh_parts,
                      paste0("(", t1_items, ",", t2_items, ", d", d, ")",
                             collapse = ", "))
  }
  thresh_constraints <- paste(thresh_parts, collapse = ", ")

  all_constraints <- paste(discr_constraints, thresh_constraints, sep = ", ")

  constrained_spec <- paste0(
    "Theta_T1 = 1-", n_items, "\n",
    "Theta_T2 = ", n_items + 1, "-", 2 * n_items, "\n",
    "COV  = Theta_T1*Theta_T2\n",
    "MEAN = Theta_T2\n",
    "CONSTRAIN = ", all_constraints
  )

  unconstrained_spec <- paste0(
    "Theta_T1 = 1-", n_items, "\n",
    "Theta_T2 = ", n_items + 1, "-", 2 * n_items, "\n",
    "COV  = Theta_T1*Theta_T2\n",
    "MEAN = Theta_T2"
  )

  list(constrained_spec = constrained_spec,
       unconstrained_spec = unconstrained_spec)
}

################################################################################
## Single-replication fit functions
################################################################################

# GRMTree: returns 1 if the tree splits on time, 0 if it does not, NA on error.
fit_grmtree_replication <- function(n, rho, n_items,
                                    a_params_t1, d_params_t1,
                                    a_params_t2, d_params_t2,
                                    theta_mean_change,
                                    minbucket = 50) {

  theta_df <- generate_correlated_theta(n, rho, mean_t1 = 0,
                                        mean_t2 = theta_mean_change)

  resp_t1 <- simdata(a = a_params_t1, d = d_params_t1, N = n,
                     itemtype = "graded",
                     Theta = matrix(theta_df$theta_t1, ncol = 1))

  resp_t2 <- simdata(a = a_params_t2, d = d_params_t2, N = n,
                     itemtype = "graded",
                     Theta = matrix(theta_df$theta_t2, ncol = 1))

  df_long <- create_person_period_data(resp_t1, resp_t2, n_items)

  tryCatch({
    tree_model <- grmtree(
      resp ~ time,
      data = df_long,
      control = grmtree.control(
        minbucket = minbucket,
        p_adjust  = "bonferroni",
        alpha     = 0.05
      )
    )
    n_nodes <- length(nodeids(tree_model, terminal = TRUE))
    ifelse(n_nodes > 1, 1, 0)
  }, error = function(e) NA_real_)
}

# LRT via a two-factor longitudinal GRM: returns 1 if the LRT p < 0.05, 0
# otherwise, NA on error.
fit_lrt_replication <- function(n, rho, n_items,
                                a_params_t1, d_params_t1,
                                a_params_t2, d_params_t2,
                                theta_mean_change,
                                model_specs) {

  theta_df <- generate_correlated_theta(n, rho, mean_t1 = 0,
                                        mean_t2 = theta_mean_change)

  resp_t1 <- simdata(a = a_params_t1, d = d_params_t1, N = n,
                     itemtype = "graded",
                     Theta = matrix(theta_df$theta_t1, ncol = 1))

  resp_t2 <- simdata(a = a_params_t2, d = d_params_t2, N = n,
                     itemtype = "graded",
                     Theta = matrix(theta_df$theta_t2, ncol = 1))

  # Wide format for the two-factor longitudinal specification
  colnames(resp_t1) <- paste0("T1_Item", 1:n_items)
  colnames(resp_t2) <- paste0("T2_Item", 1:n_items)
  data_wide <- cbind(resp_t1, resp_t2)

  tryCatch({

    mod_constrained <- mirt(
      data     = data_wide,
      model    = mirt.model(model_specs$constrained_spec),
      itemtype = "graded",
      SE       = FALSE,
      verbose  = FALSE
    )

    mod_unconstrained <- mirt(
      data     = data_wide,
      model    = mirt.model(model_specs$unconstrained_spec),
      itemtype = "graded",
      SE       = FALSE,
      verbose  = FALSE
    )

    lrt_result <- anova(mod_constrained, mod_unconstrained)
    p_value <- lrt_result$p[2]
    ifelse(p_value < 0.05, 1, 0)

  }, error = function(e) NA_real_)
}
