process PHASE2_BATCH {
  tag "${params.mode}:${params.locus_ids ?: params.chromosomes}"

  input:
  val trigger

  output:
  val('phase2_batch_complete'), emit: done

  script:
  def continueFlag = params.continue_on_error ? '--continue-on-error' : ''
  def locusArg = params.locus_ids ? "--locus-ids '${params.locus_ids}'" : ''
  def chrArg = params.chromosomes ? "--chromosomes '${params.chromosomes}'" : ''
  """
  cd '${params.repo_root}'
  ${params.python_bin} scripts/run_phase2_multi_locus.py \
    --mode '${params.mode}' \
    --max-loci ${params.max_loci} \
    ${chrArg} \
    ${locusArg} \
    --docker-image '${params.docker_image}' \
    ${continueFlag}
  """
}

process INTEGRATE_SYSTEM {
  tag 'duckdb+postgres'

  input:
  val trigger

  output:
  val('integration_complete'), emit: done

  script:
  def pgArg = params.postgres_dsn ? "--postgres-dsn '${params.postgres_dsn}'" : ''
  """
  cd '${params.repo_root}'
  ${params.python_bin} scripts/integrate_system_layer.py ${pgArg}
  """
}

process RUN_BENCHMARK {
  tag 'susie-vs-finemap'

  input:
  val trigger

  output:
  val('benchmark_complete'), emit: done

  script:
  """
  cd '${params.repo_root}'
  bash scripts/run_benchmark.sh
  """
}

