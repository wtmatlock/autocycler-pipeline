nextflow.enable.dsl = 2

// Workflow
workflow {

  // Step 1: Import all long-read FASTQ files as samples
  Channel.fromPath(params.longReadsGlob)
    .map { file ->
      def sampleId = file.baseName.replaceFirst(/\.fastq(\.gz)?$/, "")
      tuple(sampleId, file)
    }
    .set { longReadSamples }

  // Step 2: Perform quality control on long-reads
  qcLongReadsChannel = QC_LONGREADS(longReadSamples)

  // Step 3: For each sample, estimate genome size via Raven helper script
  qcLongReadsChannel
    .map { sampleId, fastq -> tuple(sampleId, fastq, file("bin/genome_size_raven.sh")) }
    .set { genomeSizeInputs }

  genomeSizeChannel = GENOME_SIZE_ESTIMATE(genomeSizeInputs)

  // Step 4: Subsample reads for each sample
  subsampleInputs = qcLongReadsChannel
    .join(genomeSizeChannel.map { sampleId, _originalReads, genomeSize -> tuple(sampleId, genomeSize) })
    .map { sampleId, qcReads, genomeSize -> tuple(sampleId, qcReads, genomeSize) }

  subsampledReadsChannel = AUTOCYCLER_SUBSAMPLE(subsampleInputs)

  // Step 5: Assemble each subsample with each assembler script in /bin
  plassemblerDatabaseChannel = DOWNLOAD_PLASSEMBLER_DB()

  assemblers = params.assemblers

  assemblyInputs = subsampledReadsChannel.flatMap { sampleId, subsampledReads, genomeSize ->
    subsampledReads.collectMany { subsampledRead ->
      def subsampleReadId = subsampledRead.baseName.replaceFirst(/\.fastq$/, "")
      assemblers.collect { assembler ->
        tuple(
          [
            sampleId,
            subsampledRead,
            subsampleReadId,
            genomeSize,
            assembler,
            file("bin/${assembler}.sh"),
            file("bin/canu_trim.py"),
          ]
        )
      }
    }
  }

  assemblyInputsWithDB = assemblyInputs.combine(plassemblerDatabaseChannel)

  assembledChannel = AUTOCYCLER_ASSEMBLY(assemblyInputsWithDB)

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

  combinedChannel = AUTOCYCLER_COMBINE(combineInputs)

  // Optional step: QC short-reads and polish final assembly
  if (params.runShortReads) {

    Channel.fromFilePairs(params.shortReadsGlob, flat: true)
      .set { shortReadSamples }

    // Perform quality control on short-reads
    qcShortReadsChannel = QC_SHORTREADS(shortReadSamples)

    // Index the combined assembly
    indexInputs = qcShortReadsChannel
      .join(combinedChannel, by: 0)
      .map { sampleId, qc1, qc2, assembly ->
        tuple(sampleId, assembly, qc1, qc2)
      }

    indexedAssemblyChannel = INDEX_ASSEMBLY(indexInputs)

    // Polish the combined assembly using short-reads
    polishInputs = indexedAssemblyChannel.map { sampleId, assembly, aln1, aln2 ->
      tuple(sampleId, assembly, aln1, aln2)
    }

    polishedAssemblyChannel = POLISH_ASSEMBLY(polishInputs)

    // Step 11: Reorient the final assembly using dnaapler
    _reorientedAssemblyChannel = REORIENT_ASSEMBLY(polishedAssemblyChannel)
  }
  else {

    // Step 11: Reorient the final assembly using dnaapler
    _reorientedAssemblyChannel = REORIENT_ASSEMBLY(combinedChannel)
  }
}

// Processes

process QC_LONGREADS {
  tag "QC_LONGREADS:${sampleId}"

  container "quay.io/biocontainers/filtlong:0.2.1--hdcf5f25_4"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/qced_longreads", mode: "copy", pattern: "${sampleId}_qc.fastq"

  input:
  tuple val(sampleId), path(rawLongReads)

  output:
  tuple val(sampleId), path("${sampleId}_qc.fastq")

  script:
  // Perform quality control on long-reads using Filtlong
  """
  filtlong \\
    --min_length 1000 \\
    --keep_percent 95 ${rawLongReads} \\
    > ${sampleId}_qc.fastq
  """
}

process GENOME_SIZE_ESTIMATE {
  tag "GENOME_SIZE_ESTIMATE:${sampleId}"

  container "wtmatlock/autocycler-suite:0.4.0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/genome_size_estimate", mode: "copy", pattern: "${sampleId}_genome_size.txt"

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

  container "wtmatlock/autocycler-suite:0.4.0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: "copy", pattern: "subsampled_longreads/*.fastq"

  input:
  tuple val(sampleId), path(qcLongReads), path(genomeSize)

  output:
  tuple val(sampleId), path("subsampled_longreads/*.fastq"), path(genomeSize)

  script:
  // Subsample reads using the genome size
  """
  mkdir -p subsampled_longreads
  autocycler subsample \\
    --reads "${qcLongReads}" \\
    --out_dir subsampled_longreads \\
    --genome_size \$(head -n1 "${genomeSize}")
  """
}

process DOWNLOAD_PLASSEMBLER_DB {
  tag "DOWNLOAD_PLASSEMBLER_DB"

  container "wtmatlock/autocycler-suite:0.4.0"

  memory params.memory

  cpus params.threads

  publishDir "./", mode: "copy"

  output:
  path "plassembler_db"

  script:
  """
    plassembler download -d plassembler_db
  """
}

process AUTOCYCLER_ASSEMBLY {
  tag { "AUTOCYCLER_ASSEMBLY:${sampleId}:${assembler}:${subsampleId}" }

  container "wtmatlock/autocycler-suite:0.4.0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/assemblies", mode: "copy", pattern: "*.fasta"

  input:
  tuple val(sampleId), path(subsampledLongReads), val(subsampleId), path(genomeSize), val(assembler), path(script), path(canuTrim), path(plassemblerDB)

  output:
  tuple val(sampleId), path("*.fasta")

  script:
  // Assemble reads using the specified assembler script
  // Note: Canu requires an additional script and Plassembler requires an additional database
  // I set 'TERM=xterm-256color' so Plassembler can find Unicycler (otherwise it fails to run)
  """
  export TERM=xterm-256color
  export PLASSEMBLER_DB="${plassemblerDB}"
  chmod +x "${script}"
  chmod +x "${canuTrim}"
  ./${script} \\
    ${subsampledLongReads} \\
    "${sampleId}_${assembler}_${subsampleId}" \\
    "${params.threads}" \\
    \$(head -n1 "${genomeSize}")
  """
}

process AUTOCYCLER_COMPRESS {
  tag "AUTOCYCLER_COMPRESS:${sampleId}"

  container "wtmatlock/autocycler-suite:0.4.0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/autocycler_outputs", mode: "copy", pattern: "**"

  input:
  tuple val(sampleId), path(assembliesDir)

  output:
  tuple val(sampleId), path("**")

  script:
  // Merge assembled contigs into a unitig graph
  // Remove empty FASTA files before compression 
  // i.e. from Plassambler when no plasmids are recovered
  """
  find ${assembliesDir} -name '*.fasta' -size 0 -delete
  autocycler compress \\
    -i ${assembliesDir} \\
    -a . \\
    -t ${params.threads}
  """
}

process AUTOCYCLER_CLUSTER {
  tag "AUTOCYCLER_CLUSTER:${sampleId}"

  container "wtmatlock/autocycler-suite:0.4.0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: "copy", pattern: "**", overwrite: true

  input:
  tuple val(sampleId), path(autocyclerDir)

  output:
  tuple val(sampleId), path("**")

  script:
  // Group unitigs into putative chromosomes/plasmids
  """
  autocycler cluster -a ${autocyclerDir}
  """
}

process AUTOCYCLER_TRIM {
  tag "AUTOCYCLER_TRIM:${sampleId}"

  container "wtmatlock/autocycler-suite:0.4.0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: "copy", pattern: "**", overwrite: true

  input:
  tuple val(sampleId), path(autocyclerDir)

  output:
  tuple val(sampleId), path("**")

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

  container "wtmatlock/autocycler-suite:0.4.0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: "copy", pattern: "**", overwrite: true

  input:
  tuple val(sampleId), path(autocyclerDir)

  output:
  tuple val(sampleId), path("**")

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

  container "wtmatlock/autocycler-suite:0.4.0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: "copy", pattern: "**", overwrite: true

  input:
  tuple val(sampleId), path(autocyclerDir)

  output:
  tuple val(sampleId), path("${autocyclerDir}/consensus_assembly.fasta")

  script:
  // Merge all resolved graphs into one final assembly per sample
  """
  autocycler combine -a ${autocyclerDir} \\
                     -i ${autocyclerDir}/clustering/qc_pass/cluster_*/5_final.gfa
  """
}

process QC_SHORTREADS {
  tag "QC_SHORTREADS:${sampleId}"

  container "quay.io/biocontainers/fastp:0.24.2--heae3180_0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/qced_shortreads", mode: "copy", pattern: "${sampleId}_qc_*.fastq.gz"

  input:
  tuple val(sampleId), path(rawShortReads1), path(rawShortReads2)

  output:
  tuple val(sampleId), path("${sampleId}_qc_1.fastq.gz"), path("${sampleId}_qc_2.fastq.gz")

  script:
  // Perform quality control on short-reads using Filtlong
  """
  fastp --in1 ${rawShortReads1} \\
        --in2 ${rawShortReads2}  \\
        --out1 ${sampleId}_qc_1.fastq.gz \\
        --out2 ${sampleId}_qc_2.fastq.gz
  """
}

process INDEX_ASSEMBLY {
  tag "INDEX_ASSEMBLY:${sampleId}"

  container "quay.io/biocontainers/bwa:0.7.19--h577a1d6_1"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/indexed_assembly", mode: "copy", pattern: "${sampleId}_alignments_*.sam"

  input:
  tuple val(sampleId), path(assembly), path(qcShortReads1), path(qcShortReads2)

  output:
  tuple val(sampleId), path(assembly), path("${sampleId}_alignments_1.sam"), path("${sampleId}_alignments_2.sam")

  script:
  // Index assembly using BWA
  """
  bwa index ${assembly}
  bwa mem -t ${params.threads} \\
          -a ${assembly} ${qcShortReads1} \\
          > ${sampleId}_alignments_1.sam
  bwa mem -t ${params.threads} \\
          -a ${assembly} ${qcShortReads2} \\
          > ${sampleId}_alignments_2.sam
  """
}

process POLISH_ASSEMBLY {
  tag "POLISH_ASSEMBLY:${sampleId}"

  container "quay.io/biocontainers/polypolish:0.6.0--h3ab6199_3"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}/polished_assembly", mode: "copy", pattern: "${sampleId}_polished_assembly.fasta"

  input:
  tuple val(sampleId), path(assembly), path(alignments1), path(alignments2)

  output:
  tuple val(sampleId), path("${sampleId}_polished_assembly.fasta")

  script:
  // Perform quality control on short-reads using Filtlong
  """
  polypolish filter --in1 ${alignments1} \\
                    --in2 ${alignments2} \\
                    --out1 "${sampleId}_filtered_1.sam" \\
                    --out2 "${sampleId}_filtered_2.sam"
  polypolish polish ${assembly} \\
                    "${sampleId}_filtered_1.sam" \\
                    "${sampleId}_filtered_2.sam" \\
                    > "${sampleId}_polished_assembly.fasta"
  """
}

process REORIENT_ASSEMBLY {
  tag "REORIENT_ASSEMBLY:${sampleId}"

  container "wtmatlock/dnaapler:1.2.0"

  memory params.memory

  cpus params.threads

  publishDir "${params.outdir}/${sampleId}", mode: "copy", pattern: "reoriented_assembly/*.fasta"

  input:
  tuple val(sampleId), path(assembly)

  output:
  tuple val(sampleId), path("reoriented_assembly/*.fasta")

  script:
  // Reorient the final assembly using dnaapler
  """
  dnaapler all -i ${assembly} \\
   -o reoriented_assembly \\
   -t ${params.threads} \\
   -p "${sampleId}"
  """
}
