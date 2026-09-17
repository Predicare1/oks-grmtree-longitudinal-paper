################################################################################
##
## 02_grmtree_type1.R
##
## Type I error simulation for GRMTree in the longitudinal setting.
##
## Design: item parameters are identical at both occasions (null of invariance
## holds). We record whether GRMTree spuriously splits on time (a false
## positive) across 500 replications for each combination of:
##
##   - scale length (n_items): 6, 12
##   - sample size (n):        200, 500, 750, 1500
##   - latent trait correlation across occasions (rho): 0.3, 0.6, 0.9
##
## Output: type1_grmtree_results.rds
##
################################################################################

# Load helpers and packages
source("01_helpers.R")

set.seed(2026)

################################################################################
## Simulation design
################################################################################

sim_conditions <- expand.grid(
  n_items = c(6, 12),
  n       = c(200, 500, 750, 1500),
  rho     = c(0.3, 0.6, 0.9),
  stringsAsFactors = FALSE
)

n_replications <- 500

cat("\n============================================================\n")
cat("TYPE I ERROR SIMULATION FOR GRMTree (LONGITUDINAL)\n")
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

  n_items <- sim_conditions$n_items[cond_idx]
  n       <- sim_conditions$n[cond_idx]
  rho     <- sim_conditions$rho[cond_idx]

  cat("------------------------------------------------------------\n")
  cat("Condition ", cond_idx, " of ", nrow(sim_conditions),
      " | items = ", n_items, " | n = ", n, " | rho = ", rho, "\n", sep = "")

  # Get baseline item parameters. For the null condition, T2 parameters equal
  # T1 parameters (no non-invariance is induced).
  params      <- get_item_params(n_items)
  a_params_t1 <- params$a_params
  d_params_t1 <- params$d_params

  a_params_t2 <- a_params_t1   # identical -> null of invariance
  d_params_t2 <- d_params_t1

  cond_start  <- Sys.time()
  rep_results <- numeric(n_replications)

  for (rep in seq_len(n_replications)) {

    if (rep %% 100 == 0) cat("    rep ", rep, " / ", n_replications, "\n", sep = "")

    rep_results[rep] <- fit_grmtree_replication(
      n                 = n,
      rho               = rho,
      n_items           = n_items,
      a_params_t1       = a_params_t1,
      d_params_t1       = d_params_t1,
      a_params_t2       = a_params_t2,
      d_params_t2       = d_params_t2,
      theta_mean_change = 0,           # no mean shift in the null
      minbucket         = 50
    )
  }

  valid_idx <- !is.na(rep_results)
  n_valid   <- sum(valid_idx)
  n_splits  <- sum(rep_results[valid_idx] == 1)
  type1     <- (n_splits / n_valid) * 100

  results_list[[cond_idx]] <- data.frame(
    n_items    = n_items,
    n          = n,
    rho        = rho,
    n_valid    = n_valid,
    n_splits   = n_splits,
    type1_pct  = round(type1, 2),
    time_mins  = round(as.numeric(difftime(Sys.time(), cond_start, units = "mins")), 2)
  )

  cat("  Type I error: ", round(type1, 2), "% (", n_splits, "/", n_valid, ")\n",
      sep = "")
}

overall_end <- Sys.time()
total_mins  <- as.numeric(difftime(overall_end, overall_start, units = "mins"))

################################################################################
## Summarize and save
################################################################################

summary_table <- do.call(rbind, results_list)

cat("\n============================================================\n")
cat("TYPE I ERROR RESULTS\n")
cat("============================================================\n\n")
print(summary_table, row.names = FALSE)

# Save results
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

saveRDS(output, file = "type1_grmtree_results.rds")

cat("\nResults saved to: type1_grmtree_results.rds\n")
cat("Total runtime:  ", round(total_mins, 1), " minutes\n", sep = "")
