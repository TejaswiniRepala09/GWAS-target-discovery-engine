#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/run_susie_diagnostics_pilot.R <locus_id> [summary_path] [ld_path] [fit_path] [z_path] [output_tag]")
}

locus_id <- args[[1]]

summary_path <- if (length(args) >= 2) args[[2]] else file.path("results", "fine_mapping", "loci_inputs", locus_id, "summary_stats.tsv")
ld_path <- if (length(args) >= 3) args[[3]] else file.path("results", "fine_mapping", "loci_inputs", locus_id, "ld_matrix.tsv")
fit_path <- if (length(args) >= 4) args[[4]] else file.path("results", "fine_mapping", "susie_raw", locus_id, "susie_fit.rds")
z_path <- if (length(args) >= 5) args[[5]] else file.path("results", "fine_mapping", "susie_raw", locus_id, "zscores.tsv")
output_tag <- if (length(args) >= 6) args[[6]] else locus_id
diag_dir <- file.path("results", "fine_mapping", "diagnostics")
plot_dir <- file.path("results", "plots", "phase2")
diag_txt <- file.path(diag_dir, paste0(output_tag, "_diagnostics.txt"))
eig_plot <- file.path(plot_dir, paste0(output_tag, "_ld_eigenvalues.png"))
resid_plot <- file.path(plot_dir, paste0(output_tag, "_residuals_plot.png"))

dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(summary_path, ld_path, fit_path, z_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(paste("Missing required input(s):", paste(missing_paths, collapse = ", ")))
}

ss <- fread(summary_path, sep = "\t")
ld <- fread(ld_path, sep = "\t")
zdt <- fread(z_path, sep = "\t")
fit <- readRDS(fit_path)

if (!("SNP" %in% names(ld))) stop("ld_matrix.tsv must include SNP as first column")
if (!("SNP" %in% names(zdt) && "Z" %in% names(zdt))) stop("zscores.tsv must include SNP and Z columns")

snps <- ld$SNP
R <- as.matrix(ld[, -1, with = FALSE])
storage.mode(R) <- "double"
rownames(R) <- snps
colnames(R) <- snps

setkey(zdt, SNP)
z <- zdt[snps]$Z
if (any(is.na(z))) stop("Z-score alignment contains NA values after matching SNP order")

n_eff <- suppressWarnings(as.numeric(median(ss$N, na.rm = TRUE)))
if (!is.finite(n_eff) || is.na(n_eff) || n_eff <= 0) {
  stop("Could not infer valid effective sample size from summary_stats.tsv column N")
}

# LD matrix quality checks
max_abs_asym <- max(abs(R - t(R)))
R_sym <- (R + t(R)) / 2
diag_before <- diag(R_sym)
diag_dev_max <- max(abs(diag_before - 1))
diag(R_sym) <- 1
diag_after <- diag(R_sym)

eig <- eigen(R_sym, symmetric = TRUE, only.values = TRUE)$values
min_eig <- min(eig)
max_eig <- max(eig)
neg_eig_count <- sum(eig < -1e-8)
near_zero_count <- sum(abs(eig) < 1e-8)
pos_eig <- eig[eig > 1e-8]
cond_num <- if (length(pos_eig) > 0) max(pos_eig) / min(pos_eig) else NA_real_

eig_df <- data.frame(index = seq_along(eig), eigenvalue = sort(eig, decreasing = TRUE))
p_eig <- ggplot(eig_df, aes(x = index, y = eigenvalue)) +
  geom_line(color = "#264653", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#d62828") +
  labs(
    title = paste0(locus_id, ": LD Eigenvalue Spectrum"),
    x = "Eigenvalue rank",
    y = "Eigenvalue"
  ) +
  theme_minimal(base_size = 12)
ggsave(eig_plot, p_eig, width = 8, height = 5, dpi = 220)

# RSS consistency diagnostics
s_est <- estimate_s_rss(z = z, R = R_sym, n = n_eff)
kr <- kriging_rss(z = z, R = R_sym, n = n_eff, s = s_est)
resid_df <- as.data.frame(kr$conditional_dist)

if (!all(c("z", "condmean", "z_std_diff", "logLR") %in% colnames(resid_df))) {
  stop("kriging_rss did not return expected conditional diagnostics columns")
}

resid_mean_abs <- mean(abs(resid_df$z_std_diff))
resid_q95_abs <- quantile(abs(resid_df$z_std_diff), probs = 0.95, names = FALSE)
resid_gt2 <- sum(abs(resid_df$z_std_diff) > 2)
resid_gt3 <- sum(abs(resid_df$z_std_diff) > 3)
resid_loglr_max <- max(resid_df$logLR, na.rm = TRUE)

p_resid <- ggplot(resid_df, aes(x = condmean, y = z)) +
  geom_point(alpha = 0.45, size = 0.8, color = "#3a86ff") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#d62828") +
  labs(
    title = paste0(locus_id, ": RSS Observed vs Expected z"),
    x = "Expected z (conditional mean from LD)",
    y = "Observed z"
  ) +
  theme_minimal(base_size = 12)
ggsave(resid_plot, p_resid, width = 8, height = 6, dpi = 220)

# Convergence diagnostics from stored fit
niter <- if (!is.null(fit$niter)) as.integer(fit$niter) else NA_integer_
converged <- if (!is.null(fit$converged)) as.logical(fit$converged) else NA
elbo <- if (!is.null(fit$elbo)) fit$elbo else numeric(0)
elbo_len <- length(elbo)
elbo_delta <- if (elbo_len >= 2) elbo[elbo_len] - elbo[1] else NA_real_
elbo_tail_delta <- if (elbo_len >= 2) elbo[elbo_len] - elbo[max(1, elbo_len - 1)] else NA_real_
sigma2 <- if (!is.null(fit$sigma2)) as.numeric(fit$sigma2) else NA_real_
hit_max_iter <- isTRUE(!is.na(niter) && niter >= 200)

well_conditioned <- (max_abs_asym < 1e-6) && (neg_eig_count == 0) && is.finite(cond_num) && (cond_num < 1e8)
rss_consistent <- is.finite(resid_q95_abs) && (resid_q95_abs < 2.5) && (resid_gt3 < 0.05 * length(z))

lines_out <- c(
  paste0("locus_id: ", locus_id),
  paste0("timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  "[input_files]",
  paste0("summary_stats: ", summary_path),
  paste0("ld_matrix: ", ld_path),
  paste0("susie_fit: ", fit_path),
  paste0("zscores: ", z_path),
  "",
  "[alignment]",
  paste0("variant_count: ", length(snps)),
  paste0("effective_n_median: ", format(n_eff, scientific = FALSE)),
  "",
  "[ld_quality]",
  paste0("max_abs_asymmetry: ", signif(max_abs_asym, 6)),
  paste0("diag_max_abs_deviation_before_reset: ", signif(diag_dev_max, 6)),
  paste0("diag_min_after_reset: ", signif(min(diag_after), 6)),
  paste0("diag_max_after_reset: ", signif(max(diag_after), 6)),
  paste0("min_eigenvalue: ", signif(min_eig, 6)),
  paste0("max_eigenvalue: ", signif(max_eig, 6)),
  paste0("negative_eigenvalue_count(<-1e-8): ", neg_eig_count),
  paste0("near_zero_eigenvalue_count(|eig|<1e-8): ", near_zero_count),
  paste0("condition_number_pos_eigs: ", signif(cond_num, 6)),
  paste0("ld_well_conditioned_flag: ", well_conditioned),
  "",
  "[rss_diagnostics]",
  paste0("estimate_s_rss: ", signif(s_est, 6)),
  paste0("mean_abs_z_std_diff: ", signif(resid_mean_abs, 6)),
  paste0("q95_abs_z_std_diff: ", signif(resid_q95_abs, 6)),
  paste0("count_abs_z_std_diff_gt_2: ", resid_gt2),
  paste0("count_abs_z_std_diff_gt_3: ", resid_gt3),
  paste0("max_logLR: ", signif(resid_loglr_max, 6)),
  paste0("rss_consistency_flag: ", rss_consistent),
  "",
  "[susie_fit]",
  paste0("niter: ", niter),
  paste0("converged: ", converged),
  paste0("hit_max_iter_assumed_200: ", hit_max_iter),
  paste0("sigma2(residual_variance_proxy): ", signif(sigma2, 6)),
  paste0("elbo_length: ", elbo_len),
  paste0("elbo_delta_first_to_last: ", signif(elbo_delta, 6)),
  paste0("elbo_delta_last_step: ", signif(elbo_tail_delta, 6)),
  "",
  "[interpretation]",
  "If ld_well_conditioned_flag=FALSE, LD matrix regularization/shrinkage is recommended.",
  "If rss_consistency_flag=FALSE, summary-statistics / LD mismatch is likely and SuSiE outputs are exploratory.",
  "If hit_max_iter_assumed_200=TRUE and converged=FALSE, increase max_iter only after LD/summary harmonization checks."
)

writeLines(lines_out, con = diag_txt)
cat("Wrote diagnostics report:", diag_txt, "\n")
cat("Wrote plot:", eig_plot, "\n")
cat("Wrote plot:", resid_plot, "\n")
