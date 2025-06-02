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

  assembled_ch = autocycler_assemble(assembly_tasks)

  // Step 6: Compress assembled contigs into unitig graphs

  // Group all assembly FASTAs by sample_id
  assemblies_grouped = assembled_ch
    .groupTuple(by: 0)
    .map { sample_id, fasta_files -> tuple(sample_id, file("${params.outdir}/${sample_id}/assemblies")) }

  compressed_ch = autocycler_compress(assemblies_grouped)

  // Step 7: Cluster unitigs per sample

  cluster_inputs = compressed_ch.map { sample_id, autocycler_dir  ->
    tuple(sample_id, file("${params.outdir}/${sample_id}/autocycler_outputs"))
    }

  clustered_ch = autocycler_cluster(cluster_inputs)

  // Step 8: Trim clusters in series

  trimmed_inputs = clustered_ch.map { sample_id, autocycler_dir  ->
    tuple(sample_id, file("${params.outdir}/${sample_id}/autocycler_outputs"))
  }
  
  trimmed_ch = autocycler_trim(trimmed_inputs)

  // Step 8: Resolve clusters in series

  resolved_inputs = trimmed_ch.map { sample_id, autocycler_dir  ->
    tuple(sample_id, file("${params.outdir}/${sample_id}/autocycler_outputs"))
  }

  resolved_ch = autocycler_resolve(resolved_inputs)

  // Step 9: Combine final GFA graphs per sample
  combined_inputs = resolved_ch.map { sample_id, autocycler_dir  ->
    tuple(sample_id, file("${params.outdir}/${sample_id}/autocycler_outputs"))
  }

  combined_ch = autocycler_combine(combined_inputs)

  // Step 10: Reorient the final assembly using dnaapler
  reorient_inputs = combined_inputs.map { sample_id, autocycler_dir  ->
    tuple(sample_id, file("${params.outdir}/${sample_id}/autocycler_outputs"))
  }

  //reorient_ch = reorient_assembly(reorient_inputs)

}

// Processes

process genome_size {

  tag "genome_size:${sample_id}"

  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}/genome_size_estimate", mode: 'copy', pattern: "${sample_id}.gsize.txt"

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

    publishDir "${params.outdir}/${sample_id}/qced_longreads", mode: 'copy', pattern: "${sample_id}_qc.fastq"

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

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "subsampled_longreads/*.fastq"

  input:
    tuple val(sample_id), path(fastq), path(gsize)

  output:
    tuple val(sample_id), path("subsampled_longreads/*.fastq"), path(gsize)

  script:
  // Subsample reads using the genome size
  """
  mkdir -p subsampled_longreads
  autocycler subsample \
    --reads "$fastq" \
    --out_dir subsampled_longreads \
    --genome_size \$(head -n1 "$gsize")
  echo -e "${sample_id}\tsubsampled_longreads/${sample_id}\t$gsize"
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

  tag "compress:${sample_id}"
  
  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}/autocycler_outputs", mode: 'copy', pattern: "**"

  input:
    tuple val(sample_id), path(assemblies_dir)

  output:
    tuple val(sample_id), path("**")

  script:
  // Merge assembled contigs into a unitig graph
  """
  autocycler compress \
    -i ${assemblies_dir} \
    -a . \
    -t ${params.threads}
  """
}

process autocycler_cluster {

  tag "cluster:${sample_id}"
  
  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "**", overwrite: true

  input:
    tuple val(sample_id), path(autocycler_dir)

  output:
    tuple val(sample_id), path("**")

  script:
  // Group unitigs into putative chromosomes/plasmids
  """
  autocycler cluster -a ${autocycler_dir}
  """
}

process autocycler_trim {

  tag "trim:${sample_id}"
  
  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "**", overwrite: true

  input:
    tuple val(sample_id), path(autocycler_dir)

  output:
    tuple val(sample_id), path("**")

  script:
  // Remove low-quality tips from each cluster graph
  """
  for c in ${autocycler_dir}/clustering/qc_pass/cluster_*; do
    autocycler trim -c "\$c"
  done
  """
}

process autocycler_resolve {

  tag "resolve:${sample_id}"
  
  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "**", overwrite: true

  input:
    tuple val(sample_id), path(autocycler_dir)

  output:
    tuple val(sample_id), path("**")

  script:
  // Resolve repeats and ambiguities for each cluster graph
  """
  for c in ${autocycler_dir}/clustering/qc_pass/cluster_*; do
    autocycler resolve -c "\$c"
  done
  """
}

process autocycler_combine {

  tag "resolve:${sample_id}"
  
  container 'wtmatlock/autocycler-suite:linux_amd64'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "**", overwrite: true

  input:
    tuple val(sample_id), path(autocycler_dir)

  output:
    tuple val(sample_id), path("**")

  script:
  // Merge all resolved graphs into one final assembly per sample
  """
  autocycler combine -a ${autocycler_dir} -i ${autocycler_dir}/clustering/qc_pass/cluster_*/5_final.gfa
  """
}

process reorient_assembly {

  tag "reorient:${sample_id}"
  
  container 'quay.io/biocontainers/dnaapler:1.2.0--pyhdfd78af_0'

  memory '32 GB'

  publishDir "${params.outdir}/${sample_id}/reoriented_assembly", mode: 'copy', pattern: "*.fasta"

  input:
    tuple val(sample_id), path(autocycler_dir)

  output:
    tuple val(sample_id), path("*.fasta")

  script:
  // Reorient the final assembly using dnaapler
  """
  dnaapler all -i ${autocycler_dir}/consensus_assembly.fasta \\
   -o reoriented_assembly \\
   -p "${sample_id}"
  """
}


