################################################################################
##
## 03_grmtree_power.R
##
## Power simulation for GRMTree in the longitudinal setting.
##
## Design: non-invariance is induced on three items at T2 by adding a constant
## Delta to their discrimination (non-uniform), their thresholds (uniform), or
## both. Statistical power is the proportion of replications in which GRMTree
## correctly splits on time.
##
## Factors varied (500 replications per cell):
##
##   - scale length (n_items):        6, 12
##   - sample size (n):               200, 500, 750, 1500
##   - latent trait correlation (rho): 0.3, 0.6, 0.9
##   - non-invariance type:           uniform, non-uniform, both
##   - non-invariance magnitude:      0.5 (medium), 1.0 (large)
##   - latent mean shift at T2:       0, 0.5
##
## Output: power_grmtree_results.rds
##
################################################################################

source("01_helpers.R")

set.seed(2026)

################################################################################
## Simulation design
################################################################################

sim_conditions <- expand.grid(
  n_items           = c(6, 12),
  n                 = c(200, 500, 750, 1500),
  rho               = c(0.3, 0.6, 0.9),
  dif_type          = c("nonuniform", "uniform", "both"),
  dif_magnitude     = c(0.5, 1.0),
  theta_mean_change = c(0, 0.5),
  stringsAsFactors  = FALSE
)

n_replications <- 500

cat("\n============================================================\n")
cat("POWER SIMULATION FOR GRMTree (LONGITUDINAL)\n")
cat("============================================================\n")
cat("Total conditions:  ", nrow(sim_conditions), "\n")
cat("Replications/cond: ", n_replications, "\n")
cat("Total fits:        ", nrow(sim_conditions) * n_replications, "\n")
cat("============================================================\n\n")

################################################################################
## Main loop
################################################################################

overall_start <- Sys.time()
results_list  <- vector("list", nrow(sim_conditions))

for (cond_idx in seq_len(nrow(sim_conditions))) {

  cond <- sim_conditions[cond_idx, ]

  cat("------------------------------------------------------------\n")
  cat("Condition ", cond_idx, " of ", nrow(sim_conditions),
      " | items=", cond$n_items, " | n=", cond$n, " | rho=", cond$rho,
      " | type=", cond$dif_type, " | mag=", cond$dif_magnitude,
      " | mu_T2=", cond$theta_mean_change, "\n", sep = "")

  # Baseline (T1) item parameters
  params      <- get_item_params(cond$n_items)
  a_params_t1 <- params$a_params
  d_params_t1 <- params$d_params
  dif_items   <- params$dif_items

  # Induce non-invariance at T2
  t2_params <- create_t2_params(
    a_params_t1   = a_params_t1,
    d_params_t1   = d_params_t1,
    dif_items     = dif_items,
    dif_type      = cond$dif_type,
    dif_magnitude = cond$dif_magnitude
  )

  cond_start  <- Sys.time()
  rep_results <- numeric(n_replications)

  for (rep in seq_len(n_replications)) {

    if (rep %% 100 == 0) cat("    rep ", rep, " / ", n_replications, "\n", sep = "")

    rep_results[rep] <- fit_grmtree_replication(
      n                 = cond$n,
      rho               = cond$rho,
      n_items           = cond$n_items,
      a_params_t1       = a_params_t1,
      d_params_t1       = d_params_t1,
      a_params_t2       = t2_params$a_params_t2,
      d_params_t2       = t2_params$d_params_t2,
      theta_mean_change = cond$theta_mean_change,
      minbucket         = 50
    )
  }

  valid_idx <- !is.na(rep_results)
  n_valid   <- sum(valid_idx)
  n_splits  <- sum(rep_results[valid_idx] == 1)
  power_pct <- (n_splits / n_valid) * 100

  results_list[[cond_idx]] <- data.frame(
    n_items           = cond$n_items,
    n                 = cond$n,
    rho               = cond$rho,
    dif_type          = cond$dif_type,
    dif_magnitude     = cond$dif_magnitude,
    theta_mean_change = cond$theta_mean_change,
    n_valid           = n_valid,
    n_splits          = n_splits,
    power_pct         = round(power_pct, 2),
    time_mins         = round(as.numeric(difftime(Sys.time(), cond_start, units = "mins")), 2)
  )

  cat("  Power: ", round(power_pct, 2), "% (", n_splits, "/", n_valid, ")\n",
      sep = "")
}

overall_end <- Sys.time()
total_mins  <- as.numeric(difftime(overall_end, overall_start, units = "mins"))

################################################################################
## Summarize and save
################################################################################

summary_table <- do.call(rbind, results_list)

cat("\n============================================================\n")
cat("POWER RESULTS\n")
cat("============================================================\n\n")
print(summary_table, row.names = FALSE)

output <- list(
  summary    = summary_table,
  parameters = list(
    n_replications = n_replications,
    seed           = 2026,
    minbucket      = 50
  ),
  timing = list(
    start      = overall_start,
    end        = overall_end,
    total_mins = total_mins
  )
)

saveRDS(output, file = "power_grmtree_results.rds")

cat("\nResults saved to: power_grmtree_results.rds\n")
cat("Total runtime:  ", round(total_mins, 1), " minutes\n", sep = "")
