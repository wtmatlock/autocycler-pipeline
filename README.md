**Directory layout**

```bash
autocycler_pipeline/
├── nextflow.config           # Pipeline config
├── main.nf                   # Nextflow script
├── bin/                      # Assembler helper scripts from Autocycler
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
└── autocycler_out/           # Output directories are created
    ├── sample1/
    ├── sample2/
    └── ...
```
