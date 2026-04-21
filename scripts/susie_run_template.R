#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(susieR)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('Usage: Rscript susie_run_template.R <summary_stats.tsv> <ld_matrix.tsv> <out_dir>')
summary_path <- args[[1]]
ld_path <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(summary_path, sep='	', header=TRUE)
if (!all(c('SNP','BETA','SE') %in% names(dat))) stop('summary_stats.tsv must include SNP,BETA,SE')
dat <- dat[!is.na(SNP) & !is.na(BETA) & !is.na(SE) & SE > 0]
dat[, Z := BETA/SE]
fwrite(dat[, .(SNP, Z)], file.path(out_dir, 'zscores.tsv'), sep='	')

ld <- fread(ld_path, sep='	', header=TRUE)
if (!('SNP' %in% names(ld))) stop('ld_matrix.tsv must include SNP first column')
ld_snps <- ld$SNP
R <- as.matrix(ld[, setdiff(names(ld),'SNP'), with=FALSE])
if (nrow(R) != ncol(R)) stop('LD matrix must be square')
common <- intersect(dat$SNP, ld_snps)
if (length(common) < 10) stop('Too few overlapping SNPs for susie_rss')
dat2 <- dat[match(common, SNP), ]
idx <- match(common, ld_snps)
R2 <- R[idx, idx, drop=FALSE]
diag(R2) <- 1
n_eff <- if ('N' %in% names(dat2)) median(dat2$N, na.rm=TRUE) else 100000
if (!is.finite(n_eff)) n_eff <- 100000
fit <- susie_rss(z=dat2$Z, R=R2, n=as.numeric(n_eff), max_iter=200, estimate_residual_variance=TRUE)

pip <- data.table(SNP=dat2$SNP, PIP=fit$pip)
fwrite(pip[order(-PIP)], file.path(out_dir, 'pip.tsv'), sep='	')

cs <- susie_get_cs(fit)
if (length(cs$cs) > 0) {
  cs_rows <- rbindlist(lapply(seq_along(cs$cs), function(i) {
    ids <- cs$cs[[i]]
    data.table(credible_set_id=paste0('CS', i), SNP=dat2$SNP[ids], cs_coverage=cs$coverage[i], cs_log10bf=cs$cs_log10bf[i])
  }))
} else {
  cs_rows <- data.table(credible_set_id=character(), SNP=character(), cs_coverage=numeric(), cs_log10bf=numeric())
}
fwrite(cs_rows, file.path(out_dir, 'credible_sets.tsv'), sep='	')
saveRDS(fit, file.path(out_dir, 'susie_fit.rds'))
