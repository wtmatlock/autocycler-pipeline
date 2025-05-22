**Directory layout**

```bash
autocycler_pipeline/
├── nextflow.config      # Pipeline config
├── main.nf              # Nextflow script
├── bin/                 # Assembler helper scripts from Autocycler
│   ├── canu_trim.py
│   ├── canu.sh
│   ├── flye.sh
│   ├── genome_size_raven.sh
│   ├── miniasm.sh
│   └── raven.sh
├── long_reads/               # Input long-reads
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   └── ...
└── results/             # (Optional) empty; Nextflow writes outputs here or under params.outdir
```

**Running the pipeline**

1. Make the wrapper scripts executable:
   ```bash
   chmod +x bin/*.sh
   ```

2. Launch Nextflow from the top-level directory:
   ```bash
   nextflow run main.nf \
     -c nextflow.config \
     --reads_glob "long_reads/*.fastq.gz" \
     --threads 16 \
     --outdir results
   ```

3. You can monitor progress in the built-in Nextflow UI:
   ```bash
   nextflow logs
   nextflow trace
   ```

4. Final assemblies will be in:
   ```bash
   results/<sample_id>/final_assembly.gfa
   ```
