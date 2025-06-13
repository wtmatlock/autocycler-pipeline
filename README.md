# An automated Nextflow pipeline for Autocycler using Docker
## Directory layout

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

## Usage

1. You will need to install [Nextflow](https://nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/engine/install/).

2. Clone this repository to your local machine:

```
git clone https://github.com/wtmatlock/autocycler_pipeline.git
cd autocycler_pipeline
```

3. Add your input data to the appropriate directories:

- Place your long-read FASTQ files in a `long_reads` directory, naming each file with its sample label (e.g. `sample1.fastq.gz`, `sample2.fastq.gz`).
- If you have paired-end short reads, add them to a `short_reads` directory using the convention `sample1_1.fastq.gz` and `sample1_2.fastq.gz` for each sample. In the `nextflow.config`, ensure that `runShortReads = true`. 
- Review and adjust the `nextflow.config` file to set resource requirements e.g. `threads = 8` and `memory    = '16 GB'`.

4. You can then launch the pipeline with:
```
nextflow run main.nf 
```
> Note: if you are running on linux/arm64 (Apple Silicon), you will need to use Docker's emulation for linux/amd64.
