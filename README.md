**Directory layout**

```bash
./
├── nextflow.config           # Nextflow config
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
├── short_reads/              # Input short-reads
│   ├── sample1_1.fastq.gz
│   ├── sample1_2.fastq.gz
│   ├── sample2_1.fastq.gz
│   ├── sample2_2.fastq.gz
│   └── ...
└── outputs/                  # Output directories are created
    ├── sample1/
    ├── sample2/
    └── ...
```