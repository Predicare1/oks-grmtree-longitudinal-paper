################################################################################
##
## 06_oks_analysis.R
##
## Empirical analysis: longitudinal measurement invariance in the Oxford Knee
## Score (OKS), before and 1 year after total knee arthroplasty (TKA).
##
## Analyses:
##   1. Load data
##   2. Descriptive summary (Table 3, Supplementary Table A5, Supplementary
##      Figure A1)
##   3. Baseline GRM fit and fit indices (RMSEA, CFI, SRMSR)
##   4. GRMTree with time as the only partitioning covariate (Figure 1, Table 4)
##   5. Two-factor longitudinal GRM and likelihood ratio tests (LRT): global
##      test and item-level tests (Supplementary Table A6). The post-surgery
##      latent mean and variance are freely estimated.
##
## Input:  cleaned_oks_data.xlsx (not shared publicly; see data note)
## Output: printed tables, Supplementary Figure A1, and fitted models (.rds)
##
## Data note: the cleaned OKS dataset is derived from the Winnipeg Regional
## Health Authority (WRHA) Joint Replacement Registry and cannot be shared due
## to data-use restrictions. Investigators seeking access should contact the
## registry directly.
##
## Authors: O. I. Arimoro, T. T. Sajobi, E. Bohm, L. M. Lix
## Contact: olayinka.arimoro@ucalgary.ca
##
################################################################################

suppressPackageStartupMessages({
  library(grmtree)
  library(mirt)
  library(tidyverse)
  library(readxl)
})

# EnvStats is not attached because it masks base print.default, which breaks
# printing of the GRMTree object; stat_n_text() is called as EnvStats::stat_n_text()

################################################################################
## 1. Load data
##
## Expected columns:
##   - Q1_Pre  ... Q12_Pre   (OKS items before surgery, scored 0-4)
##   - Q1_Post ... Q12_Post  (OKS items 1 year after surgery, scored 0-4)
##   - Patient characteristics (age, Sex, BMI, ...)
##
## Scoring direction: items are scored 0 (worst) to 4 (best); the total OKS
## ranges from 0 (worst function) to 48 (best function).
################################################################################

oks_data <- read_xlsx("cleaned_oks_data.xlsx")

cat("Loaded", nrow(oks_data), "patients.\n")

n_items <- 12

################################################################################
## 2. Descriptive summary
################################################################################

cat("\n============================================================\n")
cat("DESCRIPTIVE SUMMARY\n")
cat("============================================================\n")

oks_data <- oks_data %>%
  mutate(
    OKS_Pre  = rowSums(across(Q1_Pre:Q12_Pre)),
    OKS_Post = rowSums(across(Q1_Post:Q12_Post))
  )

# Table 3: patient characteristics and OKS totals
median_iqr <- function(x) {
  sprintf("%.1f (%.1f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
}

cat("\nTable 3: median (IQR) or n (%)\n")
cat("  Age:                 ", median_iqr(oks_data$age), "\n")
cat("  BMI:                 ", median_iqr(oks_data$BMI), "\n")
cat("  OKS before surgery:  ", median_iqr(oks_data$OKS_Pre), "\n")
cat("  OKS after surgery:   ", median_iqr(oks_data$OKS_Post), "\n")
cat("  Sex:\n")
print(oks_data %>% count(Sex) %>% mutate(percent = round(100 * n / sum(n), 1)))

# Supplementary Table A5: distribution of item responses before and after surgery
item_distribution <- oks_data %>%
  dplyr::select(Q1_Pre:Q12_Pre, Q1_Post:Q12_Post) %>%
  pivot_longer(everything(), names_to = c("item", "time"),
               names_sep = "_", values_to = "response") %>%
  count(time, item, response) %>%
  group_by(time, item) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  mutate(item = factor(item, levels = paste0("Q", 1:n_items)),
         time = factor(time, levels = c("Pre", "Post"))) %>%
  arrange(item, time, response)

cat("\nSupplementary Table A5: item response distribution, n (%)\n")
print(as.data.frame(item_distribution), row.names = FALSE)

################################################################################
## 3. Baseline GRM fit and fit indices
##
## Unidimensional GRM fitted to the pre-surgery data by marginal maximum
## likelihood. Benchmarks: RMSEA and SRMSR < 0.08, CFI > 0.90.
################################################################################

cat("\n============================================================\n")
cat("BASELINE GRM FIT (before surgery)\n")
cat("============================================================\n")

resp_pre <- oks_data %>%
  dplyr::select(Q1_Pre:Q12_Pre) %>%
  drop_na() %>%
  as.matrix()

oks_grm_pre <- mirt(resp_pre, model = 1, itemtype = "graded",
                    SE = TRUE, verbose = FALSE)

cat("\nGRM item parameters (baseline):\n")
print(coef(oks_grm_pre, IRTpars = TRUE, simplify = TRUE)$items)

cat("\nM2-based fit statistics (RMSEA, CFI, SRMSR):\n")
print(M2(oks_grm_pre, type = "C2"))

################################################################################
## 4. GRMTree with time as the only partitioning covariate
##
## Data are restructured into person-period format: each patient contributes
## one row before and one row after surgery, and time is the partitioning
## covariate. A split on time indicates longitudinal measurement
## non-invariance.
##
## Minimum terminal node size = 600, i.e., 10 times the number of item
## parameters (12 items x 5 parameters = 60).
################################################################################

cat("\n============================================================\n")
cat("GRMTREE (time as the only covariate)\n")
cat("============================================================\n")

item_names <- paste0("Q", 1:n_items, "_oks")
n_patients <- nrow(oks_data)

df_pre <- oks_data %>%
  dplyr::select(Q1_Pre:Q12_Pre) %>%
  setNames(item_names) %>%
  mutate(id = 1:n_patients, time = "Pre-surgery")

df_post <- oks_data %>%
  dplyr::select(Q1_Post:Q12_Post) %>%
  setNames(item_names) %>%
  mutate(id = 1:n_patients, time = "Post-surgery")

df_long <- bind_rows(df_pre, df_post) %>%
  mutate(time = factor(time, levels = c("Pre-surgery", "Post-surgery"))) %>%
  drop_na()

cat("Person-period dataset:", nrow(df_long), "rows (",
    nrow(df_long) / 2, "patients).\n")

# Supplementary Figure A1: total OKS before and after surgery
df_long <- df_long %>%
  mutate(total_oks = rowSums(across(all_of(item_names))))

group_means <- df_long %>%
  group_by(time) %>%
  summarise(mean_oks = round(mean(total_oks), 2), .groups = "drop")

fig_a1 <- ggplot(df_long, aes(x = time, y = total_oks)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(width = 0.25, outlier.shape = NA) +
  geom_point(aes(color = time), alpha = 0.8, size = 1,
             position = position_jitter(width = 0.05, height = 0)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1),
               color = "#880000", linewidth = 1, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", aes(group = 1),
               color = "#880000", size = 2) +
  EnvStats::stat_n_text(size = 3.5, color = "steelblue", text.box = TRUE) +
  annotate("text", label = paste("mu ==", group_means$mean_oks[1]),
           size = 3.25, x = 0.75, y = group_means$mean_oks[1], parse = TRUE) +
  annotate("text", label = paste("mu ==", group_means$mean_oks[2]),
           size = 3.25, x = 2.25, y = group_means$mean_oks[2], parse = TRUE) +
  labs(x = "Time", y = "Oxford Knee Score") +
  theme_bw() +
  theme(legend.position = "none")

print(fig_a1)
ggsave("supplementary_figure_A1.png", fig_a1, width = 6, height = 5, dpi = 300)

# Fit the GRMTree
df_long$resp <- as.matrix(df_long[, item_names])

oks_tree <- grmtree(
  resp ~ time,
  data    = df_long,
  control = grmtree.control(p_adjust = "bonferroni", minbucket = 600)
)

cat("\nGRMTree (Figure 1):\n")
print(oks_tree)
plot(oks_tree)

# Item parameters by terminal node (Table 4)
cat("\nTable 4: thresholds by terminal node\n")
print(threshpar_grmtree(oks_tree))

cat("\nTable 4: discriminations by terminal node\n")
print(discrpar_grmtree(oks_tree))

saveRDS(oks_tree, file = "oks_grmtree_fit.rds")

################################################################################
## 5. Two-factor longitudinal GRM and likelihood ratio tests
##
## Items before surgery load on Theta_T1 and items after surgery load on
## Theta_T2. Theta_T1 is fixed at N(0, 1) to set the scale. In the fully and
## partially constrained models, the post-surgery latent mean and variance and
## the covariance between occasions are freely estimated, so a change in the
## latent trait distribution is not attributed to the item parameters.
##
## The unconstrained model frees all post-surgery item parameters. The
## post-surgery variance is not freed there because the free post-surgery
## discriminations already allow for it.
##
## Global test: fully constrained vs. unconstrained model.
## Item-level tests (Supplementary Table A6): release one item at a time and
## compare with the fully constrained model.
##   - overall:     discrimination and 4 thresholds released (df = 5)
##   - uniform:     4 thresholds released (df = 4)
##   - non-uniform: discrimination released (df = 1)
## P-values are Bonferroni-adjusted for 12 comparisons within each test type.
################################################################################

cat("\n============================================================\n")
cat("LONGITUDINAL GRM AND LIKELIHOOD RATIO TESTS\n")
cat("============================================================\n")

oks_wide <- oks_data %>%
  dplyr::select(Q1_Pre:Q12_Pre, Q1_Post:Q12_Post) %>%
  drop_na() %>%
  as.matrix()

# CONSTRAIN string: discriminations (a_items) and thresholds (d_items) equal
# across time for the listed items
build_constraints <- function(a_items = 1:n_items, d_items = 1:n_items) {
  parts <- c()
  if (length(a_items) > 0) {
    parts <- c(parts, paste0("(", a_items, ", a1, ", a_items + n_items, ", a2)",
                             collapse = ", "))
  }
  if (length(d_items) > 0) {
    for (k in 1:4) {
      parts <- c(parts, paste0("(", d_items, ",", d_items + n_items, ", d", k, ")",
                               collapse = ", "))
    }
  }
  paste(parts, collapse = ", ")
}

factor_lines <- paste0("
  Theta_T1 = 1-", n_items, "
  Theta_T2 = ", n_items + 1, "-", 2 * n_items, "
  MEAN = Theta_T2")

# Constrained or partially constrained model: post-surgery variance free
spec_constrained <- function(constraints) {
  paste0(factor_lines, "
  COV = Theta_T1*Theta_T2, Theta_T2*Theta_T2
  CONSTRAIN = ", constraints)
}

# Unconstrained model: all post-surgery item parameters free
spec_unconstrained <- paste0(factor_lines, "
  COV = Theta_T1*Theta_T2")

fit_grm <- function(spec) {
  mirt(oks_wide, mirt.model(spec), itemtype = "graded",
       SE = FALSE, verbose = FALSE)
}

# --- Global test
cat("\nFitting fully constrained and unconstrained models...\n")
mod_constrained   <- fit_grm(spec_constrained(build_constraints()))
mod_unconstrained <- fit_grm(spec_unconstrained)

cat("\nGlobal LRT (fully constrained vs. unconstrained):\n")
print(anova(mod_constrained, mod_unconstrained))

cat("\nEstimated latent means (post-surgery relative to baseline):\n")
print(coef(mod_constrained, simplify = TRUE)$means)
cat("\nEstimated latent covariance matrix:\n")
print(coef(mod_constrained, simplify = TRUE)$cov)

# --- Item-level tests
others <- function(i) setdiff(1:n_items, i)

item_test <- function(i, type) {
  constraints <- switch(type,
    "Overall"     = build_constraints(a_items = others(i), d_items = others(i)),
    "Uniform"     = build_constraints(a_items = 1:n_items, d_items = others(i)),
    "Non-uniform" = build_constraints(a_items = others(i), d_items = 1:n_items))
  cmp <- anova(mod_constrained, fit_grm(spec_constrained(constraints)))
  data.frame(item = paste0("Q", i), type = type,
             X2 = cmp$X2[2], df = cmp$df[2], p = cmp$p[2])
}

cat("\nItem-level LRTs (36 model fits; this takes a while)...\n")
item_lrt_table <- do.call(rbind, lapply(c("Overall", "Uniform", "Non-uniform"),
  function(type) {
    cat("  ", type, "\n")
    res <- do.call(rbind, lapply(1:n_items, item_test, type = type))
    res$p_adj <- p.adjust(res$p, method = "bonferroni")
    res
  }))

# Supplementary Table A6 layout: one row per item, chi-square and adjusted p by type
table_a6 <- item_lrt_table %>%
  mutate(X2 = round(X2, 1), p_adj = round(p_adj, 3)) %>%
  dplyr::select(item, type, X2, p_adj) %>%
  pivot_wider(names_from = type, values_from = c(X2, p_adj),
              names_glue = "{type}_{.value}") %>%
  dplyr::select(item,
                Overall_X2, Overall_p_adj,
                Uniform_X2, Uniform_p_adj,
                `Non-uniform_X2`, `Non-uniform_p_adj`)

cat("\nSupplementary Table A6: item-level LRTs (Bonferroni-adjusted p)\n")
print(as.data.frame(table_a6), row.names = FALSE)

saveRDS(list(global_constrained   = mod_constrained,
             global_unconstrained = mod_unconstrained,
             item_lrt_table       = item_lrt_table),
        file = "oks_lrt_results.rds")

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("Saved: oks_grmtree_fit.rds, oks_lrt_results.rds, supplementary_figure_A1.png\n")
cat("============================================================\n")
