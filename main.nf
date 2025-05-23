nextflow.enable.dsl=2

// Workflow

workflow {

  // Step 1: Import all FASTQ files as samples
Channel
  .fromPath(params.reads_glob)
  .map { file ->
    def sample_id = file.getBaseName().replaceFirst(/\.fastq(\.gz)?$/, '')
    tuple(sample_id, file)
  }
  .set { samples }

  // Step 2: For each sample, compute genome size via Raven helper script
  samples
    .map { sample_id, fastq -> tuple(sample_id, fastq, params.threads, file('bin/genome_size_raven.sh')) }
    .set { genome_inputs }

  genome_size_ch = genome_size(genome_inputs)

  // Step 3: Perform quality control on long-reads
  qc_longreads_ch = samples
    .map { sample_id, fastq -> tuple(sample_id, fastq, params.threads) }
    .set { qc_inputs }

  qc_longreads_ch = qc_longreads(qc_inputs)

  // Step 4: Subsample reads for each sample
  subsample_inputs = qc_longreads_ch
    .join(genome_size_ch.map{ sample_id, orig_reads, gsize -> tuple(sample_id, gsize) })
    .map { sample_id, qc_reads, gsize -> tuple(sample_id, qc_reads, gsize) }

  subsample_inputs.view()
  
  subsampled_reads_ch = autocycler_subsample(subsample_inputs)

  // Step 5: Assemble each subsample with each assembler script in /bin

  // NB canu takes a lot longer than the other assemblers, so only run when everything else is debugged
  // assemblers = ['canu', 'flye', 'miniasm', 'raven']

  assemblers = ['flye', 'miniasm', 'raven']

  assembly_tasks = subsampled_reads_ch
    .flatMap { sample_id, fastq_files, gsize ->
      fastq_files.collectMany { fastq_file ->
        assemblers.collect { asm ->
          tuple(sample_id, fastq_file, gsize, params.threads, asm, file("bin/${asm}.sh"), file("bin/canu_trim.py"))
        }
      }
    }

  assembly_tasks.view()
  assembled_ch = autocycler_assemble(assembly_tasks)
  assembled_ch.view()

  // Step 6: Compress assembled contigs into unitig graphs
  // Group assemblies by sample_id
  //assembled_grouped = assembled_ch.groupTuple(by: 0)
  //compressed_ch = autocycler_compress(assembled_grouped)

  // Step 7: Cluster unitigs per sample
  //clustered_ch = autocycler_cluster(compressed_ch)

  // Step 8: Trim and resolve clusters in series
  //trimmed_ch = autocycler_trim(
  //  clustered_ch.flatMap { sample_id, clusters ->
  //  clusters.collect { cluster_dir -> tuple(sample_id, cluster_dir) }
  //  }
  //  )


  //resolved_ch = autocycler_resolve(trimmed_ch)


  // Step 9: Combine final GFA graphs per sample
  // Group resolved graphs by sample_id
  // combined = resolved_ch.groupBy { it[0] }
  // combined.subscribe { sample_id, records ->
  // def gfas = records.collect { it[2] }
  // autocycler_combine(sample_id, gfas)
  // }

}

// Processes

process genome_size {

  tag "genome_size:${sample_id}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}/genome_size", mode: 'copy', pattern: "${sample_id}.gsize.txt"

  input:
    tuple val(sample_id), path(long_reads), val(threads), path(script)

  output:
    tuple val(sample_id), path(long_reads), path("${sample_id}.gsize.txt")

script:
// Compute genome size using Raven helper script (skipping polishing step)
"""
chmod +x "$script"
"$script" "$long_reads" "$threads" > ${sample_id}.gsize.txt
"""
}

process qc_longreads {

    tag "qc_longreads:${sample_id}"

    container 'quay.io/biocontainers/filtlong:0.2.1--hdcf5f25_4'
    
    memory '32GB'

    publishDir "${params.outdir}/${sample_id}/qced_reads", mode: 'copy', pattern: "${sample_id}_qc.fastq"

    input:
       tuple val(sample_id), path(long_reads), val(threads)

    output:
    tuple val(sample_id), path("${sample_id}_qc.fastq")

    script:
    // Perform quality control on long-reads using Filtlong
    """
    filtlong --min_length 1000 --keep_percent 95 $long_reads > ${sample_id}_qc.fastq      
    """
}

process autocycler_subsample {

  tag "subsample:${sample_id}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "subsampled_reads/*.fastq"

  input:
    tuple val(sample_id), path(fastq), path(gsize)

  output:
    tuple val(sample_id), path("subsampled_reads/*.fastq"), path(gsize)

  script:
  // Subsample reads using the genome size
  """
  mkdir -p subsampled_reads
  autocycler subsample \
    --reads "$fastq" \
    --out_dir subsampled_reads \
    --genome_size \$(head -n1 "$gsize")
  echo -e "${sample_id}\tsubsampled_reads/${sample_id}\t$gsize"
  """
}

process autocycler_assemble {

  tag {"assemble:${sample_id}:${asm}"}

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}/assemblies", mode: 'copy', pattern: "*.fasta"

  input:
    tuple val(sample_id), path(reads), path(gsize), val(threads), val(asm), path(script), path(canu_trim)

  output:
    tuple val(sample_id), path("*.fasta")
  
  script:
  // Assemble reads using the specified assembler script
  // Canu requires an additional script
  """
  read_id=\$(basename ${reads} | sed 's/\\.fastq\\(\\.gz\\)\\?\$//')
  chmod +x "$script"
  "$script" \\
    ${reads} \\
    "${sample_id}_${asm}_\${read_id}" \\
    "$threads" \\
    \$(head -n1 "$gsize")
  """
}

process autocycler_compress {
  
  container 'wtmatlock/autocycler-suite'

  tag "compress:${sample_id}"

  input:
    tuple val(sample_id), path(asm_dirs)

  output:
    tuple val(sample_id), path("${params.outdir}/${sample_id}")

  script:
  // Merge assembled contigs into a unitig graph
  """
  mkdir -p ${params.outdir}/${sample_id}
  autocycler compress \
    -i ${asm_dirs.join(' ')} \
    -a ${params.outdir}/${sample_id}
  echo -e "${sample_id}\t${params.outdir}/${sample_id}"
  """
}

process autocycler_cluster {

  container 'wtmatlock/autocycler-suite'

  tag "cluster:${sample_id}"

  input:
    tuple val(sample_id), path(params.outdir)

  output:
    tuple val(sample_id), path(params.outdir + '/clustering/qc_pass/cluster_*')

  script:
  // Group unitigs into putative chromosomes/plasmids
  """
  autocycler cluster -a ${params.outdir}
  echo -n "${sample_id}"
  """
}

process autocycler_trim {

  container 'wtmatlock/autocycler-suite'

  tag "trim:${sample_id}"

  input:
    tuple val(sample_id), path(params.outdir)

  output:
    tuple val(sample_id), path(params.outdir)

  script:
  // Remove low-quality tips from each cluster graph
  """
  autocycler trim -c ${params.outdir}
  echo -e "${sample_id}\t${params.outdir}"
  """
}

process autocycler_resolve {

  container 'wtmatlock/autocycler-suite'

  tag "resolve:${sample_id}"

  input:
    tuple val(sample_id), path(params.outdir)

  output:
    tuple val(sample_id), path(params.outdir + '/5_final.gfa')

  script:
  // Resolve repeats and bubbles to finalize each graph
  """
  autocycler resolve -c ${params.outdir}
  echo -e "${sample_id}\t${params.outdir}/5_final.gfa"
  """
}

process autocycler_combine {

  container 'wtmatlock/autocycler-suite'

  tag "combine:${sample_id}"

  input:
    tuple val(sample_id), path(final_gfas)

  output:
    path("${params.outdir}/${sample_id}/final_assembly.gfa")

  script:
  // Merge all resolved graphs into one final assembly per sample
  """
  autocycler combine \
    -a ${params.outdir}/${sample_id} \
    -i ${final_gfas.join(' ')}
  """
}


