#!/usr/bin/env bash
set -euo pipefail

LOG_DIR="results/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/pilot_chr7_locus10.log"

log() { printf "%s\n" "$1" | tee -a "$LOG_FILE"; }
run_cmd() {
  log "\n>>> $*"
  "$@" 2>&1 | tee -a "$LOG_FILE"
}

LOCUS_ID="chr7_locus10"
BUNDLE_DIR="results/fine_mapping/loci_inputs/${LOCUS_ID}"
SUSIE_RAW_DIR="results/fine_mapping/susie_raw/${LOCUS_ID}"
VEP_RAW_DIR="results/annotations/vep_raw"
VEP_LOCUS_OUT="${VEP_RAW_DIR}/vep_output_${LOCUS_ID}.tsv"

log "Starting pilot locus Phase 2 run for ${LOCUS_ID}"

# 1) activate env and verify python deps (no reinstall here)
run_cmd source .venv/bin/activate
run_cmd python3 -c "import importlib; [importlib.import_module(m) for m in ['polars','yaml','duckdb','numpy','pandas','pyarrow','httpx']]; print('Python deps OK')"

# 2) strict external tool checks (manual install expected)
for tool in docker bcftools tabix plink Rscript; do
  log "\n>>> command -v ${tool}"
  if ! command -v "$tool" >/dev/null 2>&1; then
    log "ERROR: Missing required binary on PATH: ${tool}"
    exit 1
  fi
  command -v "$tool" | tee -a "$LOG_FILE"
done

# 3) ensure cache and pilot directories
run_cmd mkdir -p "$VEP_RAW_DIR" "$SUSIE_RAW_DIR"
run_cmd test -d data/reference/vep
run_cmd test -f "${BUNDLE_DIR}/summary_stats.tsv"
run_cmd test -f "${BUNDLE_DIR}/metadata.yaml"
run_cmd sed -n 1,120p "${BUNDLE_DIR}/metadata.yaml"

# 4) extract pilot region from metadata
log "\n>>> python3 (extract region from metadata.yaml)"
python3 - <<'PY' 2>&1 | tee -a "$LOG_FILE"
import yaml
from pathlib import Path
locus='chr7_locus10'
meta=yaml.safe_load(Path(f'results/fine_mapping/loci_inputs/{locus}/metadata.yaml').read_text())
Path(f'results/fine_mapping/susie_raw/{locus}/region.txt').write_text(f"{meta['chr']}:{meta['window_start']}-{meta['window_end']}\\n")
print('Region:', Path(f'results/fine_mapping/susie_raw/{locus}/region.txt').read_text().strip())
PY

# 5) VEP per-locus run (Docker-based Ensembl VEP)
run_cmd awk -F '\t' 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="CHR")c=i; if($i=="BP")b=i; if($i=="A1")a1=i; if($i=="A2")a2=i; if($i=="SNP")s=i;} next} {id=($s!=""?$s:$c":"$b); print $c,$b,$b,$a1"/"$a2,"+",id}' "${BUNDLE_DIR}/summary_stats.tsv"
run_cmd bash -lc "awk -F '\t' 'BEGIN{OFS=\"\\t\"} NR==1{for(i=1;i<=NF;i++){if(\$i==\"CHR\")c=i; if(\$i==\"BP\")b=i; if(\$i==\"A1\")a1=i; if(\$i==\"A2\")a2=i; if(\$i==\"SNP\")s=i;} next} {id=(\$s!=""?\$s:\$c\":\"\$b); print \$c,\$b,\$b,\$a1\"/\"\$a2,\"+\",id}' ${BUNDLE_DIR}/summary_stats.tsv > ${VEP_RAW_DIR}/${LOCUS_ID}_vep_input_ensembl.tsv"
run_cmd docker run --rm \
  -v "$(pwd)":/work \
  -v "$(pwd)/data/reference/vep":/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_115.2 \
  vep \
  --input_file "/work/${VEP_RAW_DIR}/${LOCUS_ID}_vep_input_ensembl.tsv" \
  --format ensembl \
  --output_file "/work/${VEP_LOCUS_OUT}" \
  --tab \
  --everything \
  --cache \
  --offline \
  --dir_cache /opt/vep/.vep \
  --assembly GRCh37 \
  --force_overwrite \
  --fork 4 \
  --stats_file "/work/${VEP_RAW_DIR}/vep_summary_${LOCUS_ID}.html"
run_cmd test -s "$VEP_LOCUS_OUT"

# 6) parse per-locus VEP output
run_cmd python3 scripts/parse_vep_results.py --vep-raw-output "$VEP_LOCUS_OUT"
run_cmd test -s results/annotations/variant_annotations.tsv
run_cmd test -s results/annotations/variant_consequence_scores.tsv

# 7) region-based VCF extraction
run_cmd bash -lc "bcftools view -r \$(cat ${SUSIE_RAW_DIR}/region.txt) data/reference/Ld/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz -o ${SUSIE_RAW_DIR}/${LOCUS_ID}.region.vcf.gz"
run_cmd tabix -f -p vcf "${SUSIE_RAW_DIR}/${LOCUS_ID}.region.vcf.gz"

# 8) LD matrix generation with PLINK
run_cmd plink --vcf "${SUSIE_RAW_DIR}/${LOCUS_ID}.region.vcf.gz" --double-id --allow-extra-chr --snps-only just-acgt --r square --out "${SUSIE_RAW_DIR}/${LOCUS_ID}"
run_cmd test -s "${SUSIE_RAW_DIR}/${LOCUS_ID}.ld"
run_cmd test -s "${SUSIE_RAW_DIR}/${LOCUS_ID}.bim"

# 9) convert + validate LD matrix
log "\n>>> python3 (LD validation/alignment/symmetry/diag/mismatch log)"
python3 - <<'PY' 2>&1 | tee -a "$LOG_FILE"
from pathlib import Path
import numpy as np
import pandas as pd

locus='chr7_locus10'
bundle=Path(f'results/fine_mapping/loci_inputs/{locus}')
raw=Path(f'results/fine_mapping/susie_raw/{locus}')
summary=pd.read_csv(bundle/'summary_stats.tsv', sep='\t')
req={'SNP','CHR','BP','A1','A2'}
missing=req-set(summary.columns)
if missing:
    raise SystemExit(f'summary_stats.tsv missing required columns: {sorted(missing)}')

bim=pd.read_csv(raw/f'{locus}.bim', sep='\s+', header=None, names=['CHR','ID','CM','BP','A1','A2'])
M=np.loadtxt(raw/f'{locus}.ld')
if M.shape[0]!=M.shape[1]:
    raise SystemExit(f'LD not square: {M.shape}')
if M.shape[0]!=len(bim):
    raise SystemExit(f'LD dim and BIM mismatch: {M.shape[0]} vs {len(bim)}')

summary=summary.copy()
bim=bim.copy()
summary['CHR']=summary['CHR'].astype(str).str.replace('chr','', regex=False)
bim['CHR']=bim['CHR'].astype(str).str.replace('chr','', regex=False)
summary['BP']=summary['BP'].astype(int)
bim['BP']=bim['BP'].astype(int)
summary['A1u']=summary['A1'].astype(str).str.upper()
summary['A2u']=summary['A2'].astype(str).str.upper()
bim['A1u']=bim['A1'].astype(str).str.upper()
bim['A2u']=bim['A2'].astype(str).str.upper()
summary['PAIR']=summary.apply(lambda r: '|'.join(sorted([r['A1u'], r['A2u']])), axis=1)
bim['PAIR']=bim.apply(lambda r: '|'.join(sorted([r['A1u'], r['A2u']])), axis=1)
bim['LD_INDEX']=np.arange(len(bim), dtype=int)

# Keep one summary row per coordinate+allele pair, preferring strongest p-value.
summary['P_NUM']=pd.to_numeric(summary['P'], errors='coerce')
summary_keyed=summary.sort_values(['P_NUM','SNP'], na_position='last').drop_duplicates(subset=['CHR','BP','PAIR'], keep='first')
joined=summary_keyed.merge(bim[['CHR','BP','PAIR','A1u','A2u','LD_INDEX']], on=['CHR','BP','PAIR'], how='inner', suffixes=('_SUM','_LD'))
if joined.empty:
    raise SystemExit('No overlapping variants after CHR/BP/allele harmonization')

joined['ALLELE_ORIENTATION']=np.where(
    (joined['A1u_SUM']==joined['A1u_LD']) & (joined['A2u_SUM']==joined['A2u_LD']),
    'same',
    np.where(
        (joined['A1u_SUM']==joined['A2u_LD']) & (joined['A2u_SUM']==joined['A1u_LD']),
        'swapped',
        'other'
    )
)

mm=raw/'ld_snp_mismatch.tsv'
rows=[]
summary_keys=set(summary_keyed[['CHR','BP','PAIR']].astype(str).agg(':'.join, axis=1).tolist())
bim_keys=set(bim[['CHR','BP','PAIR']].astype(str).agg(':'.join, axis=1).tolist())
rows += [{'type':'summary_not_in_ld','key':x} for x in sorted(summary_keys-bim_keys)]
rows += [{'type':'ld_not_in_summary','key':x} for x in sorted(bim_keys-summary_keys)]
rows += [{'type':'allele_orientation', 'key':k, 'status':s} for k,s in joined[['SNP','ALLELE_ORIENTATION']].astype(str).itertuples(index=False, name=None) if s!='same']
pd.DataFrame(rows).to_csv(mm, sep='\t', index=False)
print('Mismatch log:', mm, 'rows=', len(rows))
print('Overlap variants after harmonization:', len(joined))

joined=joined.sort_values(['BP','SNP','LD_INDEX']).drop_duplicates(subset=['LD_INDEX'], keep='first')
common=[f"{c}:{bp}:{a1}:{a2}" for c,bp,a1,a2 in joined[['CHR','BP','A1u_SUM','A2u_SUM']].itertuples(index=False, name=None)]
if len(common)<10:
    raise SystemExit(f'Too few overlapping SNPs after harmonization: {len(common)}')
ord_idx=joined['LD_INDEX'].astype(int).tolist()
A=M[np.ix_(ord_idx, ord_idx)]
A=(A + A.T)/2.0
np.fill_diagonal(A, 1.0)
max_asym=np.max(np.abs(A-A.T))
if max_asym>1e-8:
    raise SystemExit(f'LD matrix still asymmetric after symmetrization: {max_asym}')
out=pd.DataFrame(A, columns=common, index=common)
out.insert(0,'SNP',common)
out_path=bundle/'ld_matrix.tsv'
out.to_csv(out_path, sep='\t', index=False)
check=pd.read_csv(out_path, sep='\t')['SNP'].tolist()
if check!=common:
    raise SystemExit('LD SNP order mismatch after write')
summary_cols=[c for c in summary_keyed.columns if c not in {'A1u','A2u','PAIR','P_NUM'}]
summary_aligned=joined[summary_cols].copy()
summary_aligned['SNP']=common
summary_aligned_out=bundle/'summary_stats.tsv'
summary_aligned.to_csv(summary_aligned_out, sep='\t', index=False)
print('Wrote aligned summary stats:', summary_aligned_out, 'rows=', len(summary_aligned))
print('Wrote aligned LD matrix:', out_path, 'shape=', out.shape)
print('LD SNP order validation OK')
PY
run_cmd test -s "${BUNDLE_DIR}/ld_matrix.tsv"

# 10) SuSiE run + parse
run_cmd Rscript scripts/susie_run_template.R "${BUNDLE_DIR}/summary_stats.tsv" "${BUNDLE_DIR}/ld_matrix.tsv" "${SUSIE_RAW_DIR}"
run_cmd test -s "${SUSIE_RAW_DIR}/pip.tsv"
run_cmd test -s "${SUSIE_RAW_DIR}/credible_sets.tsv"
run_cmd python3 scripts/parse_susie_results.py

# 11) full phase2 rerun + validation
run_cmd ln -sf "$(basename "$VEP_LOCUS_OUT")" "${VEP_RAW_DIR}/vep_output.tsv"
run_cmd python3 scripts/run_phase2.py
run_cmd test -s results/target_prioritization/variant_priority_scores.tsv
run_cmd test -s results/target_prioritization/gene_prioritization.tsv

log "\n>>> python3 (DuckDB validation)"
python3 - <<'PY' 2>&1 | tee -a "$LOG_FILE"
import duckdb
con=duckdb.connect('data/processed/ckd_target_discovery.duckdb')
print('tables=', [x[0] for x in con.execute('SHOW TABLES').fetchall()])
print('pilot variant rows=', con.execute("SELECT count(*) FROM variant_priority_scores WHERE locus_id='chr7_locus10'").fetchone()[0])
print('pilot pip rows=', con.execute("SELECT count(*) FROM pip_summary WHERE locus_id='chr7_locus10'").fetchone()[0])
con.close()
PY

log "Pilot run completed successfully for ${LOCUS_ID}"
