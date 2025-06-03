nextflow.enable.dsl = 2

// Workflow
workflow {

  // Step 1: Import all long-read FASTQ files as samples
  Channel.fromPath(params.longReadsGlob)
    .map { file ->
      def sampleId = file.baseName.replaceFirst(/\.fastq(\.gz)?$/, '')
      tuple(sampleId, file)
    }
    .set { longReadSamples }

  // Step 2: Perform quality control on long-reads
  qcLongReadsChannel = QC_LONGREADS(longReadSamples)

  // Step 3: For each sample, estimate genome size via Raven helper script
  qcLongReadsChannel
    .map { sampleId, fastq -> tuple(sampleId, fastq, file('bin/genome_size_raven.sh')) }
    .set { genomeSizeInputs }

  genomeSizeChannel = GENOME_SIZE_ESTIMATE(genomeSizeInputs)

  // Step 4: Subsample reads for each sample
  subsampleInputs = qcLongReadsChannel
    .join(genomeSizeChannel.map { sampleId, _originalReads, genomeSize -> tuple(sampleId, genomeSize) })
    .map { sampleId, qcReads, genomeSize -> tuple(sampleId, qcReads, genomeSize) }

  subsampledReadsChannel = AUTOCYCLER_SUBSAMPLE(subsampleInputs)

  // Step 5: Assemble each subsample with each assembler script in /bin

  // Note: canu takes a lot longer than the other assemblers, so only run when everything else is debugged
  // assemblers = ['canu', 'flye', 'miniasm', 'raven']
  assemblers = ['flye', 'miniasm', 'raven']

  assemblyInputs = subsampledReadsChannel.flatMap { sampleId, subsampledReads, genomeSize ->
    subsampledReads.collectMany { subsampledRead ->
      def subsampleReadId = subsampledRead.baseName.replaceFirst(/\.fastq$/, '')
      assemblers.collect { assembler ->
        tuple(
          sampleId,
          subsampledRead,
          subsampleReadId,
          genomeSize,
          assembler,
          file("bin/${assembler}.sh"),
          file('bin/canu_trim.py'),
        )
      }
    }
  }

  assembledChannel = AUTOCYCLER_ASSEMBLY(assemblyInputs)

  // Step 6: Compress assembled contigs into unitig graphs
  compressInputs = assembledChannel
    .groupTuple(by: 0)
    .map { sampleId, _qcReads -> tuple(sampleId, file("${params.outdir}/${sampleId}/assemblies")) }

  compressedChannel = AUTOCYCLER_COMPRESS(compressInputs)

  // Step 7: Cluster unitigs per sample
  clusterInputs = compressedChannel.map { sampleId, _autocycler_dir ->
    tuple(sampleId, file("${params.outdir}/${sampleId}/autocycler_outputs"))
  }

  clusteredChannel = AUTOCYCLER_CLUSTER(clusterInputs)

  // Step 8: Trim clusters in series
  trimInputs = clusteredChannel.map { sampleId, _autocycler_dir ->
    tuple(sampleId, file("${params.outdir}/${sampleId}/autocycler_outputs"))
  }

  trimmedChannel = AUTOCYCLER_TRIM(trimInputs)

  // Step 9: Resolve clusters in series
  resolveInputs = trimmedChannel.map { sampleId, _autocycler_dir ->
    tuple(sampleId, file("${params.outdir}/${sampleId}/autocycler_outputs"))
  }

  resolvedChannel = AUTOCYCLER_RESOLVE(resolveInputs)

  // Step 10: Combine final GFA graphs per sample
  combineInputs = resolvedChannel.map { sampleId, _autocycler_dir ->
    tuple(sampleId, file("${params.outdir}/${sampleId}/autocycler_outputs"))
  }

  // NB underscore
  _combinedChannel = AUTOCYCLER_COMBINE(combineInputs)

  if (params.runShortReads) {

    Channel.fromFilePairs(params.shortReadsGlob, flat: true)
      .map { sampleId, pair -> tuple(sampleId, pair[0], pair[1]) }
      .set { shortReadSamples }

    shortReadSamples.view()
  }
  else {

    log.info("Short-read processing is disabled (params.runShortReads = false)")
  }
}

// Processes

process QC_LONGREADS {
  tag "QC_LONGREADS:${sampleId}"

  container 'quay.io/biocontainers/filtlong:0.2.1--hdcf5f25_4'

  memory '32GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/qced_longreads", mode: 'copy', pattern: "${sampleId}_qc.fastq"

  input:
  tuple val(sampleId), path(rawLongReads)

  output:
  tuple val(sampleId), path("${sampleId}_qc.fastq")

  script:
  // Perform quality control on long-reads using Filtlong
  """
  filtlong --min_length 1000 --keep_percent 95 ${rawLongReads} > ${sampleId}_qc.fastq
  """
}

process GENOME_SIZE_ESTIMATE {
  tag "GENOME_SIZE_ESTIMATE:${sampleId}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/genome_size_estimate", mode: 'copy', pattern: "${sampleId}_genome_size.txt"

  input:
  tuple val(sampleId), path(qcLongReads), path(script)

  output:
  tuple val(sampleId), path(qcLongReads), path("${sampleId}_genome_size.txt")

  script:
  // Compute genome size using Raven helper script (skipping polishing step)
  """
  chmod +x "${script}"
  "${script}" "${qcLongReads}" "${params.threads}" > ${sampleId}_genome_size.txt
  """
}

process AUTOCYCLER_SUBSAMPLE {
  tag "AUTOCYCLER_SUBSAMPLE:${sampleId}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: 'copy', pattern: 'subsampled_longreads/*.fastq'

  input:
  tuple val(sampleId), path(qcLongReads), path(genomeSize)

  output:
  tuple val(sampleId), path('subsampled_longreads/*.fastq'), path(genomeSize)

  script:
  // Subsample reads using the genome size
  """
  mkdir -p subsampled_longreads
  autocycler subsample \
    --reads "${qcLongReads}" \
    --out_dir subsampled_longreads \
    --genome_size \$(head -n1 "${genomeSize}")
  echo -e "${sampleId}\tsubsampled_longreads/${sampleId}\t${genomeSize}"
  """
}

process AUTOCYCLER_ASSEMBLY {
  tag { "AUTOCYCLER_ASSEMBLY:${sampleId}:${assembler}:${subsampleId}" }

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/assemblies", mode: 'copy', pattern: '*.fasta'

  input:
  tuple val(sampleId), path(subsampledLongReads), val(subsampleId), path(genomeSize), val(assembler), path(script), path(canuTrim)

  output:
  tuple val(sampleId), path('*.fasta')

  script:
  // Assemble reads using the specified assembler script
  // Canu requires an additional script
  """
  chmod +x "${script}"
  chmod +x "${canuTrim}"
  "${script}" \\
    ${subsampledLongReads} \\
    "${sampleId}_${assembler}_${subsampleId}" \\
    "${params.threads}" \\
    \$(head -n1 "${genomeSize}")
  """
}

process AUTOCYCLER_COMPRESS {
  tag "AUTOCYCLER_COMPRESS:${sampleId}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/autocycler_outputs", mode: 'copy', pattern: '**'

  input:
  tuple val(sampleId), path(assembliesDir)

  output:
  tuple val(sampleId), path('**')

  script:
  // Merge assembled contigs into a unitig graph
  """
  autocycler compress \
    -i ${assembliesDir} \
    -a . \
    -t ${params.threads}
  """
}

process AUTOCYCLER_CLUSTER {
  tag "AUTOCYCLER_CLUSTER:${sampleId}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: 'copy', pattern: '**', overwrite: true

  input:
  tuple val(sampleId), path(autocyclerDir)

  output:
  tuple val(sampleId), path('**')

  script:
  // Group unitigs into putative chromosomes/plasmids
  """
  autocycler cluster -a ${autocyclerDir}
  """
}

process AUTOCYCLER_TRIM {
  tag "AUTOCYCLER_TRIM:${sampleId}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: 'copy', pattern: '**', overwrite: true

  input:
  tuple val(sampleId), path(autocyclerDir)

  output:
  tuple val(sampleId), path('**')

  script:
  // Remove low-quality tips from each cluster graph
  """
  for c in ${autocyclerDir}/clustering/qc_pass/cluster_*; do
    autocycler trim -c "\$c"
  done
  """
}

process AUTOCYCLER_RESOLVE {
  tag "AUTOCYCLER_RESOLVE:${sampleId}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: 'copy', pattern: '**', overwrite: true

  input:
  tuple val(sampleId), path(autocyclerDir)

  output:
  tuple val(sampleId), path('**')

  script:
  // Resolve repeats and ambiguities for each cluster graph
  """
  for c in ${autocyclerDir}/clustering/qc_pass/cluster_*; do
    autocycler resolve -c "\$c"
  done
  """
}

process AUTOCYCLER_COMBINE {
  tag "AUTOCYCLER_COMBINE:${sampleId}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: 'copy', pattern: '**', overwrite: true

  input:
  tuple val(sampleId), path(autocyclerDir)

  output:
  tuple val(sampleId), path('**')

  script:
  // Merge all resolved graphs into one final assembly per sample
  """
  autocycler combine -a ${autocyclerDir} -i ${autocyclerDir}/clustering/qc_pass/cluster_*/5_final.gfa
  """
}

process REORIENT_ASSEMBLY {
  tag "REORIENT_ASSEMBLY:${sampleId}"

  container 'quay.io/biocontainers/dnaapler:1.2.0--pyhdfd78af_0'

  memory '32 GB'

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/reoriented_assembly", mode: 'copy', pattern: '*.fasta'

  input:
  tuple val(sampleId), path(autocyclerDir)

  output:
  tuple val(sampleId), path('*.fasta')

  script:
  // Reorient the final assembly using dnaapler
  """
  dnaapler all -i ${autocyclerDir}/consensus_assembly.fasta \\
   -o reoriented_assembly \\
   -p "${sampleId}"
  """
}
