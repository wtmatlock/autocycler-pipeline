# A Dockerised Nextflow pipeline for Autocycler

## Citation

If you use this pipeline, please cite:

Matlock, W. (2025). A Dockerised Nextflow pipeline for Autocycler (v. 0.4.0) (Version 1.0.2) [Computer software]. https://doi.org/10.5281/zenodo.15706447

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15706447.svg)](https://doi.org/10.5281/zenodo.15706447)

Additionally, please cite the following:

- [Filtlong (v. 0.2.1)](https://github.com/rrwick/Filtlong)
- [Autocycler (v. 0.4.0)](https://github.com/rrwick/Autocycler)
- [Dnaapler (v. 1.2.0)](https://github.com/gbouras13/dnaapler)

Then, depending on which assemblers you used:

- [Canu (v. 2.3)](https://github.com/marbl/canu)
- [Flye (v. 2.9.5)](https://github.com/mikolmogorov/Flye)
- [Minisam (v. 0.3)](https://github.com/lh3/miniasm)
- [plassembler (v. 1.7.0)](https://github.com/gbouras13/plassembler)
- [Raven (v. 1.8.3)](https://github.com/lbcb-sci/raven)

And, if you used the short-read functionality:

- [fastp (v. 0.24.2)](https://github.com/OpenGene/fastp)
- [BWA (v. 0.7.19)](https://github.com/lh3/bwa)
- [Polypolish (v. 0.6.0)](https://github.com/rrwick/Polypolish)

## Usage

1. You will need to install [Nextflow](https://nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/engine/install/).

2. Clone this repository to your machine:

   ```bash
   git clone https://github.com/wtmatlock/autocycler-pipeline.git
   cd autocycler-pipeline
   ```

3. Add your input data to the appropriate directories:

- Place your long-read FASTQ files in a `long_reads` directory, naming each file with its sample label (e.g. `sample1.fastq.gz`, `sample2.fastq.gz`).
- If you have paired-end short reads, add them to a `short_reads` directory using the convention `sample1_1.fastq.gz` and `sample1_2.fastq.gz` for each sample. In the `nextflow.config`, ensure that `runShortReads = true`. 
- Review and adjust the `nextflow.config` file to set resource requirements (e.g. `threads = 8` and `memory    = '16 GB'`), and assemblers (e.g. `assemblers = ["canu", "flye", "miniasm", "plassembler", "raven"]`). Currently, the pipeline supports Canu, Flye, Miniasm, Plassembler, and Raven.

4. You can then launch the pipeline with:
   ```bash
   nextflow run main.nf -with-report -with-trace
   ```
   I recommend including the `-with-report` and `-with-trace` flags to make it easier to locate relevant logs in work directories and to report computational resource usage.

5. Outputs will be organised by sample label in the `outputs` directory. The final assembly FASTA can be found in `${params.outdir}/${sampleId}/reoriented_assembly`. The output from `autocycler table` can be found in `${params.outdir}/autocycler_metrics.tsv`. See the below for the full directory organisation.

> If you are running on linux/arm64 (Apple Silicon), you will need to use Docker's emulation for linux/amd64. This is currently the default in the `nextflow.config` with `runOptions = '--platform=linux/amd64'`.

> The Dockerfiles for the Autocycler suite (Autocycler plus assemblers) and Dnaapler are included in this repository but are not required to run the pipeline. All other Docker containers are pulled from [BioContainers](https://biocontainers.pro).

## Workflow

```mermaid
flowchart TB
    subgraph Inputs
        v0["Long-reads"]
        v34["Short-reads"]
        v12["Assemblers"]
    end

    v2([QC_LONGREADS])
    v5([GENOME_SIZE_ESTIMATE])
    v8([AUTOCYCLER_SUBSAMPLE])
    v10([DOWNLOAD_PLASSEMBLER_DB])
    v16([AUTOCYCLER_ASSEMBLY])
    v19([AUTOCYCLER_COMPRESS])
    v22([AUTOCYCLER_CLUSTER])
    v25([AUTOCYCLER_TRIM])
    v28([AUTOCYCLER_RESOLVE])
    v31([AUTOCYCLER_COMBINE])
    v45([AUTOCYCLER_TABLE])
    v36([QC_SHORTREADS])
    v39([INDEX_ASSEMBLY])
    v42([POLISH_ASSEMBLY])
    v44([REORIENT_ASSEMBLY])

    decision1{"Run short-reads?"}

    %% Long read flow
    v0 --> v2
    v2 --> v5
    v2 --> v8
    v5 --> v8
    v8 --> v16
    v10 --> v16
    v12 --> v16
    v16 --> v19
    v19 --> v22
    v22 --> v25
    v25 --> v28
    v28 --> v31
    v31 --> v45

    %% Conditional branching after combine
    v45 --> decision1
    decision1 -- Yes --> v36
    decision1 -- No --> v44
    v34 --> v36
    v36 --> v39

    %% Continue pipeline

    v39 --> v42
    v42 --> v44
```

## Directory organisation

```bash
./
├── nextflow.config              # Nextflow config
├── main.nf                      # Nextflow script
├── bin/                         # Assembler helper scripts from Autocycler
│   ├── canu_trim.py
│   ├── canu.sh
│   ├── flye.sh
│   ├── genome_size_raven.sh
│   ├── miniasm.sh
│   ├── plassembler.sh
│   └── raven.sh
├── long_reads/                  # Input long-reads
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   └── ...
├── short_reads/                 # Input short-reads
│   ├── sample1_1.fastq.gz
│   ├── sample1_2.fastq.gz
│   ├── sample2_1.fastq.gz
│   ├── sample2_2.fastq.gz
│   └── ...
└── outputs/                     # Output directories are created
    ├── autocycler_metrics.tsv   # Autocycler table output
    ├── sample1/
    ├── sample2/
    └── ...
```
