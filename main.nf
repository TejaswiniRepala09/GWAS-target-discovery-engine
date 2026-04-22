nextflow.enable.dsl = 2

include { PHASE2_ORCHESTRATION } from './workflows/phase2_orchestration'

workflow {
  /*
   * Stage mapping to existing project scripts:
   * 1) locus selection / batch setup           -> scripts/run_phase2_multi_locus.py
   * 2) VEP annotation                           -> internal to run_phase2_multi_locus.py
   * 3) LD extraction + harmonization            -> internal to run_phase2_multi_locus.py
   * 4) SuSiE fine-mapping                       -> internal to run_phase2_multi_locus.py
   * 5) prioritization                           -> internal to run_phase2_multi_locus.py
   * 6) integration / reporting                  -> scripts/integrate_system_layer.py
   * 7) benchmarking (optional)                  -> scripts/run_benchmark.sh
   */
  start = Channel.value('start')
  PHASE2_ORCHESTRATION(start)
}

