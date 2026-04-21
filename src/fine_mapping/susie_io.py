"""SuSiE script scaffolding and raw output parsing."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import polars as pl

from src.utils.path_utils import ensure_dir

LOGGER = logging.getLogger(__name__)


def write_susie_r_template(template_path: str | Path) -> Path:
    """Write exact susieR command template script for locus-level execution."""
    path = Path(template_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    script = """#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(susieR)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('Usage: Rscript susie_run_template.R <summary_stats.tsv> <ld_matrix.tsv> <out_dir>')
summary_path <- args[[1]]
ld_path <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(summary_path, sep='\t', header=TRUE)
if (!all(c('SNP','BETA','SE') %in% names(dat))) stop('summary_stats.tsv must include SNP,BETA,SE')
dat <- dat[!is.na(SNP) & !is.na(BETA) & !is.na(SE) & SE > 0]
dat[, Z := BETA/SE]
fwrite(dat[, .(SNP, Z)], file.path(out_dir, 'zscores.tsv'), sep='\t')

ld <- fread(ld_path, sep='\t', header=TRUE)
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
fwrite(pip[order(-PIP)], file.path(out_dir, 'pip.tsv'), sep='\t')

cs <- susie_get_cs(fit)
if (length(cs$cs) > 0) {
  cs_rows <- rbindlist(lapply(seq_along(cs$cs), function(i) {
    ids <- cs$cs[[i]]
    data.table(credible_set_id=paste0('CS', i), SNP=dat2$SNP[ids], cs_coverage=cs$coverage[i], cs_log10bf=cs$cs_log10bf[i])
  }))
} else {
  cs_rows <- data.table(credible_set_id=character(), SNP=character(), cs_coverage=numeric(), cs_log10bf=numeric())
}
fwrite(cs_rows, file.path(out_dir, 'credible_sets.tsv'), sep='\t')
saveRDS(fit, file.path(out_dir, 'susie_fit.rds'))
"""
    path.write_text(script, encoding="utf-8")
    return path


def generate_susie_commands(bundle_manifest: pl.DataFrame, fine_cfg: dict[str, Any]) -> pl.DataFrame:
    """Create runnable command table for each locus bundle."""
    template = fine_cfg["susie"].get("rscript_template", "scripts/susie_run_template.R")
    write_susie_r_template(template)
    out_root = ensure_dir(fine_cfg["paths"]["susie_raw_dir"])

    rows: list[dict[str, Any]] = []
    for row in bundle_manifest.iter_rows(named=True):
        locus_id = str(row["locus_id"])
        bundle_dir = Path(fine_cfg["paths"]["bundles_dir"]) / locus_id
        summary = bundle_dir / fine_cfg["output_contract"]["summary_stats_filename"]
        ld = bundle_dir / fine_cfg["output_contract"]["ld_matrix_filename"]
        out_dir = out_root / locus_id
        cmd = f"Rscript {template} {summary} {ld} {out_dir}"
        rows.append({"locus_id": locus_id, "susie_command": cmd, "susie_output_dir": str(out_dir)})
    return pl.from_dicts(rows) if rows else pl.DataFrame()


def parse_susie_outputs(fine_cfg: dict[str, Any]) -> dict[str, pl.DataFrame]:
    """Parse all locus outputs under susie_raw into aggregate tables."""
    raw_root = Path(fine_cfg["paths"]["susie_raw_dir"])
    raw_root.mkdir(parents=True, exist_ok=True)

    pip_list: list[pl.DataFrame] = []
    cs_list: list[pl.DataFrame] = []
    missing_rows: list[dict[str, Any]] = []

    for locus_dir in sorted([p for p in raw_root.iterdir() if p.is_dir()]):
        locus_id = locus_dir.name
        pip = locus_dir / "pip.tsv"
        cs = locus_dir / "credible_sets.tsv"
        if not pip.exists() or not cs.exists():
            missing_rows.append({"locus_id": locus_id, "missing": ",".join([str(x) for x in [pip, cs] if not x.exists()])})
            continue
        pip_df = pl.read_csv(pip, separator="\t").with_columns(pl.lit(locus_id).alias("locus_id"))
        cs_df = pl.read_csv(cs, separator="\t").with_columns(pl.lit(locus_id).alias("locus_id"))
        pip_list.append(pip_df)
        cs_list.append(cs_df)

    pip_all = pl.concat(pip_list, how="diagonal_relaxed") if pip_list else pl.DataFrame()
    cs_all = pl.concat(cs_list, how="diagonal_relaxed") if cs_list else pl.DataFrame()

    return {
        "pip_summary": pip_all,
        "credible_sets": cs_all,
        "susie_missing": pl.from_dicts(missing_rows) if missing_rows else pl.DataFrame(),
    }
