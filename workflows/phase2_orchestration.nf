include { PHASE2_BATCH; INTEGRATE_SYSTEM; RUN_BENCHMARK } from '../modules/pipeline_processes'

workflow PHASE2_ORCHESTRATION {
  take:
  trigger

  main:
  batch_done = PHASE2_BATCH(trigger)
  integrated = INTEGRATE_SYSTEM(batch_done)

  if (params.run_benchmark) {
    bench = RUN_BENCHMARK(integrated)
    done = bench
  } else {
    done = integrated
  }

  emit:
  completed = done
}

