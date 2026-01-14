# pacbio-kinnex-nf

## Background

`pacbio-kinnex-nf` is an implementation of the PacBio Single-Cell Iso-Seq pipeline in Nextflow. It takes Circular Consensus Sequencing (CCS) reads, the output of the SMRTLINK software, and processes them using the CLI Workflow from [isoseq.how](https://isoseq.how).

## Usage

Clone the repository

```bash
git clone https://github.com/kevinpryan/nf-hlamajority.git
```

Run the pipeline on test data with:

```bash
nextflow run main.nf \
       -profile test,<singularity/docker/...>
```

The `PBBM2_ALIGN` step is computationally intensive, even when running on the test dataset. The test profile requests 8 CPUs and 64 GB of RAM for this step.

To run on your own data, place all your CCS BAMs in the same directory and run:

```bash
nextflow run main.nf \
       -profile <singularity/docker/...> \
       -outdir <outdir_name> \
       --bams "/path/to/bam/dir/*.bam"
```

Expected output from running test profile:

```bash
├── multiqc
│   ├── multiqc_data
│   │   ├── isoseq_refine_boxplot_fivelen.txt
│   │   ├── isoseq_refine_boxplot_insertlen.txt
│   │   ├── isoseq_refine_boxplot_polyAlen.txt
│   │   ├── isoseq_refine_boxplot_threelen.txt
│   │   ├── llms-full.txt
│   │   ├── multiqc_citations.txt
│   │   ├── multiqc_data.json
│   │   ├── multiqc_isoseq_refine_csv.txt
│   │   ├── multiqc_isoseq_refine_json.txt
│   │   ├── multiqc.log
│   │   ├── multiqc.parquet
│   │   └── multiqc_sources.txt
│   └── multiqc_report.html
├── pigeon_seurat
│   └── ccs_pigeon_seurat
│       ├── ccs.annotated.info.csv
│       ├── ccs.info.csv
│       ├── genes_seurat
│       │   ├── barcodes.tsv
│       │   ├── genes.tsv
│       │   └── matrix.mtx
│       └── isoforms_seurat
│           ├── barcodes.tsv
│           ├── genes.tsv
│           └── matrix.mtx
├── pipeline_info
│   ├── execution_report.html
│   ├── execution_timeline.html
│   └── execution_trace.txt
└── reference_cache
    └── reference.gtf
```

Intermediate files can be found in the relevant process `work` directory.

## Pipeline description

The steps are as follows:

### Primer removal and identification of barcodes (LIMA)

Inputs: 

- CCS BAM
- Primers fasta (default is primers from 10x 3' kit)

`primers.fasta`

```bash
>5p
AAGCAGTGGTATCAACGCAGAGTACATGGG
>3p
AGATCGGAAGAGCGTCGTGTAG
```

Outputs: 

- BAM with primers removed, correctly oriented sequences. Name of primer sequences will be inserted into BAM name (e.g. for test data, BAM output name is `ccs.lima.output.5p--3p.bam`)
- BAM index (`*.bam.pbi`)
- Consensus read set XML (`*.consensusreadset.xml`)
- Counts (`*.output.lima.counts`)
- Summary (`*.output.lima.summary`)


### Tag reads with UMIs and barcode information (ISOSEQ_TAG)

Clip UMIs and cell barcodes from reads and associate with the reads for deduplication.

Inputs:

- `LIMA.out.bam`
- UMI design (default `T-12U-16B`, customise by passing `--design my-design` when running pipeline)

Outputs:
- Tagged BAM ("*.flt.bam")
- BAM index ("*.flt.bam.pbi")

### Refine reads (ISOSEQ_REFINE)

Trims poly(A) tails and removes unintended concatemer

Inputs:

- `ISOSEQ_TAG.out.bam`
- Primers fasta

Outputs:

- Full-length non-concatemer reads ("*.fltnc.bam")
- Consensus read set XML ("*.fltnc.consensusreadset.xml")
- Summary report json, used in `MULTIQC` ("*.fltnc.filter_summary.report.json")
- Summary csv, used in `MULTIQC` ("*.fltnc.report.csv")

### Correct cell barcode errors and identify real cells (ISOSEQ_CORRECT)

Identify and correct errors in cell barcodes. Requires a cell barcode whitelist, which can likely be found [here](https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/REF-10x_barcodes/). 

Inputs:

- `ISOSEQ_REFINE.out.bam`
- Barcodes, default: [3M-february-2018-REVERSE-COMPLEMENTED.txt.gz](https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/REF-10x_barcodes/3M-february-2018-REVERSE-COMPLEMENTED.txt.gz)

Outputs:

- Corrected BAM (`*.corrected.bam`)
- Corrected BAM index (`*.corrected.bam.pbi`)
- Intermediate BAM index (`*.corrected_intermediate.bam.pbi`)
- Corrected report json (`*.corrected.report.json`)


