.PHONY: install pipeline phase2 manhattan qq loci susie vep \
	fineprep susie_scaffold susie_parse vep_parse ot rank structure \
	benchmark setup_finemap_docker docs_index

install:
	python -m pip install -r requirements.txt

pipeline:
	python scripts/run_pipeline.py

phase2:
	python scripts/run_phase2.py

manhattan:
	python scripts/plot_manhattan.py

qq:
	python scripts/plot_qq.py

loci:
	python scripts/extract_lead_loci.py

susie:
	python scripts/export_for_susie.py

vep:
	python scripts/export_for_vep.py

fineprep:
	python scripts/prepare_fine_mapping_inputs.py

susie_scaffold:
	python scripts/run_susie_scaffold.py

susie_parse:
	python scripts/parse_susie_results.py

vep_parse:
	python scripts/parse_vep_results.py

ot:
	python scripts/fetch_open_targets.py

rank:
	python scripts/rank_targets.py

structure:
	python scripts/prepare_structure_candidates.py

setup_finemap_docker:
	bash scripts/setup_finemap_docker.sh

benchmark:
	bash scripts/run_benchmark.sh

docs_index:
	@echo "Docs index: docs/README.md"
