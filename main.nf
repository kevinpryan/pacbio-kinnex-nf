nextflow.enable.dsl = 2

/*
 * -------------------------------------------------
 *  PACBIO KINNEX PREPROCESSING PIPELINE 
 * -------------------------------------------------
 *  Inputs: CCS reads
 *  
 * 
 * -------------------------------------------------
*/

// Default parameters (can be overwritten via command line)
params.bams   = "data/*.bam"            // Path to input BAMs
params.pbi    = "data/*.bam.pbi"
params.primers = ""
params.barcodes = ""
params.outdir = "results"               // Output directory
params.design = "T-12U-16B"

log.info """
P A C B I O   K I N N E X    S C   P R E P R O C E S S I N G   P I P E L I N E
===============================================================================
BAMs        : ${params.bams}
BAM Indices   : ${params.pbi}
cDNA Primers  : ${params.primers}
Barcodes : ${params.barcodes}
Output      : ${params.outdir}
"""

include { LIMA } from "./modules/local/lima"
include { ISOSEQ_TAG } from "./modules/local/isoseq/tag"
include { ISOSEQ_REFINE } from "./modules/local/isoseq/refine"
include { ISOSEQ_BCSTATS } from "./modules/local/isoseq/bcstats"
include { GENERATE_KNEE_PLOT } from "./modules/local/knee"
//include { PREPROCESS } from "./subworkflows/preprocess"

process SAMTOOLS_SORT {
    label "SAMTOOLS_CONTAINER"
    tag "$bam_id"

    input:
    tuple val(bam_id), path(bam)

    output:
    tuple val(bam_id), path("${bam_id}.corrected.sorted.bam"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -t CB $bam -o ${bam_id}.corrected.sorted.bam
    """
}

/*
process LIMA {
    label "LIMA_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    path primers
    output:
    tuple val(bam_id), path("${bam_id}.lima.output.*.bam"), emit: bam
    path "${bam_id}.lima.output.*.bam.pbi", emit: index
    tuple val(bam_id), path("*.output.lima.counts"), path("*.output.lima.summary"), emit: stats
    script:
    """
    lima --per-read --isoseq $bam $primers ${bam_id}.lima.output.bam
    """
}
*/

/*
// Clip UMI and cell barcode

process ISOSEQ_TAG {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    val(design)

    output:
    tuple val(bam_id), path("${bam_id}.flt.bam"), emit: bam

    script:
    """
    isoseq tag $bam ${bam_id}.flt.bam --design $design --num-threads ${task.cpus}
    """
}
*/

/*
// Remove poly(A) tails and concatemer

process ISOSEQ_REFINE {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    path primers
    output:
    tuple val(bam_id), path("${bam_id}.fltnc.bam"), emit: bam
    //tuple val(bam_id), path("*.report.json"), path("*.report.csv"), emit: reports
    path("*.report.json"), emit: reports_json
    path("*.report.csv"), emit: reports_csv
    tuple val(bam_id), path("*.report.json"), path("*.report.csv"), emit: reports
    """
    isoseq refine $bam \
                  $primers \
                  ${bam_id}.fltnc.bam \
                  --require-polya \
                  --num-threads ${task.cpus}
    """
}
*/
// # Correct single cell barcodes based on an include list

process ISOSEQ_CORRECT {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    path barcodes
    val(opts)
    output:
    tuple val(bam_id), path("${bam_id}.corrected.bam"), emit: bam
 
    script:
    """
    isoseq correct \
          --num-threads ${task.cpus} \
          -B $barcodes \
          ${opts} \
          $bam ${bam_id}.corrected.bam
    """
}

/*
process ISOSEQ_CORRECT_KNEE {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    path barcodes
    output:
    tuple val(bam_id), path("${bam_id}.corrected.bam"), emit: bam

    script:
    """
    isoseq correct \
          --num-threads ${task.cpus} \
          -B $barcodes \
          $bam ${bam_id}.corrected.bam
    """
}

process ISOSEQ_CORRECT_PERCENTILE {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    path barcodes
    val(percentile)
    output:
    tuple val(bam_id), path("${bam_id}.corrected.bam"), emit: bam

    script:
    """
    isoseq correct \
          --num-threads ${task.cpus} \
          -B $barcodes \
          --method percentile \
          --percentile $percentile \
          $bam ${bam_id}.corrected.bam
    """
}
*/

process ISOSEQ_GROUPDEDUP {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    val parstr

    output:
    tuple val(bam_id), path("${bam_id}.dedup.bam"), emit: bam
    tuple val(bam_id), path("${bam_id}.dedup.fasta"), emit: fasta
    script:
    """
    isoseq groupdedup \
           --num-threads ${task.cpus} \
           ${parstr} \
           $bam \
           ${bam_id}.dedup.bam 
    """
}

process PBMM2_ALIGN {
    label "PBMM2_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    path ref

    output:
    tuple val(bam_id), path("${bam_id}.aligned.sorted.bam"), emit: bam

    script:
    """
    pbmm2 align --preset ISOSEQ \
                --num-threads ${task.cpus} \
                --sort \
                $bam \
                $ref \
                ${bam_id}.aligned.sorted.bam
    """
}

process ISOSEQ_COLLAPSE {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    output:
    tuple val(bam_id), path("${bam_id}.collapsed.gff"), emit: gff
    tuple val(bam_id), path("${bam_id}.collapsed.abundance.txt"), emit: abundance
    tuple val(bam_id), path("${bam_id}.collapsed.group.txt"), emit: group

    script:
    """
    isoseq collapse $bam ${bam_id}.collapsed.gff   
    """
}

// CAGE peak bed and intropolis tsv are optional inputs to this, not included here
process PIGEON_PREPARE {
    label "PIGEON_CONTAINER"
    tag "GTF"
    input:
    path gtf
    path reference
    output:
    path "*.{fasta,fa}.fai", emit: fai
    path "*.sorted.gtf", emit: gtf_sorted    
    path "*.sorted.gtf.pgi", emit: gtf_sorted_pgi
    script:
    """
    pigeon prepare *.gtf $reference 
    """
}

process PIGEON_PREPARE_GFF {
    label "PIGEON_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(gff)
    output:
    tuple val(bam_id), path("${bam_id}.collapsed.sorted.gff"), emit: gff 

    script:
    """
    pigeon prepare $gff
    """
}

/*
process PIGEON_CLASSIFY {
    label "PIGEON_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(gff), path(abundance)
    path reference
    path reference_index
    path annotation
    path annotation_index
    output:
    tuple val(bam_id), path("${bam_id}_junctions.txt"), emit: junctions
    tuple val(bam_id), path("${bam_id}_classification.txt"), emit: classification
    tuple val(bam_id), path("${bam_id}.summary.txt"), emit: summary
    tuple val(bam_id), path("${bam_id}.report.json"), emit: report
    script:
    """
    pigeon classify $gff \
                    $annotation \
                    $reference \
                    --fl $abundance \
                    --num-threads ${task.cpus}
    """
}
*/
process PIGEON_CLASSIFY {
    label "PIGEON_CONTAINER"
    tag "$bam_id"

    input:
    tuple val(bam_id), path(gff), path(abundance)
    path reference
    path reference_index
    path annotation
    path annotation_index
    path cage_bed
    path cage_index
    path poly_a
    path intropolis
    path intropolis_index

    output:
    tuple val(bam_id), path("${bam_id}_junctions.txt"), emit: junctions
    tuple val(bam_id), path("${bam_id}_classification.txt"), emit: classification
    tuple val(bam_id), path("${bam_id}.summary.txt"), emit: summary
    tuple val(bam_id), path("${bam_id}.report.json"), emit: report

    script:
    """
    pigeon classify $gff \
                    --ref $reference
                    --fl $abundance \
                    --cage-peak $cage_bed \
                    --poly-a $poly_a \
                    --intropolis $intropolis \
                    --num-threads ${task.cpus}
    """
}

// uses collapsed sorted gff
// classification is the isoforms coming from pigeon classify

process PIGEON_FILTER {
    label "PIGEON_CONTAINER"   
    tag "$bam_id"

    input:
    tuple val(bam_id), path(classification), path(junctions), path(gff)
    output:
    tuple val(bam_id), path("${bam_id}_classification.filtered_lite_classification.txt"), path("${bam_id}_classification.filtered_lite_junctions.txt"), path("${bam_id}_classification.filtered_lite_reasons.txt"), path("${bam_id}.collapsed.filtered_lite.gff"), emit: all_outputs
    tuple val(bam_id), path("${bam_id}_classification.filtered_lite_classification.txt"), emit: classification
    script:
    """
    pigeon filter $classification \
                  --isoforms $gff \
                  --num-threads ${task.cpus}
    """ 
}

process PIGEON_REPORT {
    label "PIGEON_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(classification)

    output:
    tuple val(bam_id), path("${bam_id}_saturation.txt"), emit: saturation_report

    script:
    """
    pigeon report --num-threads ${task.cpus} \
                  $classification \
                  ${bam_id}_saturation.txt
    """
}

process PREPARE_REFERENCE_GTF {
    tag "Unzipping Reference"
    // Use a small container with basic tools (gzip/bash)
    container 'spvensko/cloud-sdk-ps:495.0.0' 
    storeDir "${params.outdir}/reference_cache" 
 
    input:
    path raw_file

    output:
    path "reference.gtf", emit: gtf

    script:
    """
    if [[ "$raw_file" == *.gz ]]; then
        gunzip -c $raw_file > reference.gtf
    else
        ln -s $raw_file reference.gtf
    fi
    """
}

process PIGEON_SEURAT {
    tag "$bam_id"
    label "PIGEON_CONTAINER"
    publishDir "${params.outdir}/pigeon_seurat", mode: 'copy'

    input:
    tuple val(bam_id), path(dedup_fasta), path(collapse_group), path(classification_filtered_lite)
    val parstr

    output:
    tuple val(bam_id), path("${bam_id}_pigeon_seurat"), emit: seurat

    script:
    """
    pigeon make-seurat \
           --dedup $dedup_fasta \
           --group $collapse_group \
           -d ${bam_id}_pigeon_seurat \
           ${parstr} \
           $classification_filtered_lite
    """
}

process MULTIQC {
    tag "ALL_SAMPLES"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    //path isoseq_refine_reports_json
    //path isoseq_refine_reports_csv
    path multiqc_files

    output:
    path "*.html"
    path "*data"

    script:
    """
    multiqc .
    """
}

process SAMTOOLS_MERGE {
    label "SAMTOOLS_CONTAINER"
    tag "$bam_id"
    
    input:
    tuple val(bam_id), path(bams)

    output:
    tuple val(bam_id), path("${bam_id}.merged.aligned.bam"), emit: bam

    script:
    """
    samtools merge -@ ${task.cpus} -o ${bam_id}.merged.aligned.bam $bams
    samtools index ${bam_id}.merged.aligned.bam
    """
}

/*
workflow qc {
    // 1. Setup Channels
    bams_ch = Channel.fromPath(params.bams).map{ file -> tuple(file.simpleName, file) }
    primers_ch = Channel.value(file(params.primers))
    design_ch = Channel.value(params.design)

    // 2. Run Common Steps
    PREPROCESS_READS(bams_ch, primers_ch, design_ch)

    // 3. Run Plotting
    GENERATE_KNEE_PLOT(PREPROCESS_READS.out.fltnc_bam)
}
*/

process SKERA_SPLIT {
    label "SKERA_CONTAINER" 
    tag "$bam_id"
    
    input:
    tuple val(bam_id), path(bam)
    path mas_adapters 

    output:
    tuple val(bam_id), path("${bam_id}.skera.bam"), emit: bam
    path "${bam_id}.skera.summary.json", emit: summary_json
    path "${bam_id}.skera.summary.csv", emit: summary_csv
    path "${bam_id}.skera.read_lengths.csv", emit: read_lengths_csv

    script:
    """
    skera split \
          --num-threads ${task.cpus} \
          $bam \
          $mas_adapters \
          ${bam_id}.skera.bam
    """
}

process SAMTOOLS_SUBSAMPLE {
    label "SAMTOOLS_CONTAINER"
    tag "$bam_id"

    input: 
    tuple val(bam_id), path(bam)
    val frac_reads_subsample

    output:
    tuple val(bam_id), path("${bam_id}_subsampled.bam"), emit: bam
    script:
    """
    samtools view -b \
             --subsample ${frac_reads_subsample} \
             --subsample-seed 42 \
             -@ ${task.cpus} \
             $bam > ${bam_id}_subsampled.bam
    """
}


// WORKFLOW
workflow {
    // Channel setup
    bams_ch = Channel.fromPath(params.bams)
        .map { file -> tuple(file.simpleName, file) }
        .view { "Found sample: ${it[0]}" }

    primers_ch = Channel.value(file(params.primers))
    design_ch = Channel.value(params.design)
    barcodes_ch = Channel.value(file(params.barcodes))
    reference_ch = Channel.value(file(params.reference))
    reference_fai_ch = Channel.value(file(params.reference_fai))
    gtf_sorted_ch = Channel.value(file(params.gtf_sorted))
    gtf_sorted_pgi_ch = Channel.value(file(params.gtf_sorted_pgi))
    //annotation_ch_gzipped = Channel.value(file(params.gtf))

    //PREPARE_REFERENCE_GTF(
    //                      annotation_ch_gzipped
     //                     )

    // annotation_ch = PREPARE_REFERENCE_GTF.out.gtf
    cage_bed_ch = Channel.value(file(params.cage_bed))
    cage_index_ch = Channel.value(file(params.cage_bed_index)) // The .pgi file
    poly_a_ch = Channel.value(file(params.poly_a))
    intropolis_ch = Channel.value(file(params.intropolis))
    intropolis_index_ch = Channel.value(file(params.intropolis_index)) // The .pgi file

    //reference_set_ch = Channel.value(file(params.reference_set))
    if ( params.run_skera ) {
        mas_adapters_ch = Channel.value(file(params.adapters))
        if ( params.subsample_hifi_bam ) {
            frac_subsample_ch = Channel.value(params.frac_subsample)
            SAMTOOLS_SUBSAMPLE(bams_ch, frac_subsample_ch)
            SKERA_INPUT_BAM = SAMTOOLS_SUBSAMPLE.out.bam
        } else {
            SKERA_INPUT_BAM = bams_ch
        }
        SKERA_SPLIT(
           SKERA_INPUT_BAM,
           mas_adapters_ch    
        )
        bams_lima_ch = SKERA_SPLIT.out.bam
    } else {
        bams_lima_ch = bams_ch
    }

    LIMA(
        bams_lima_ch,
        primers_ch
    )

    ISOSEQ_TAG(
        LIMA.out.bam,
        design_ch
    )
    // can specify --require-polya if you data contains polyA tails - can remove too many reads if set and not all reads have this

    isoseq_refine_params_ch = Channel.value(params.isoseq_refine_params)

    ISOSEQ_REFINE(
        ISOSEQ_TAG.out.bam,
        primers_ch,
        isoseq_refine_params_ch        
    )

    def method = params.cell_calling_method ?: "knee"
    def bc_opts = ""
    def correct_opts = ""
    def knee_opts = ""
    def group_dedup_parstr = params.group_dedup_params ?: ""
    group_dedup_parstr_ch = Channel.value(group_dedup_parstr)

    if ( method == "percentile" ) {
        if ( !params.percentile ) { 
            error "Cell calling method is 'percentile' but --percentile was not provided!" 
        } else {
            bc_opts = "--method percentile --percentile ${params.percentile}"
            correct_opts = "--method percentile --percentile ${params.percentile}"
            knee_opts = "--estimate_percentile ${params.percentile}"
        }
    } else {
        bc_opts = "--method knee"
        correct_opts = "--method knee"
        knee_opts = ""
    }

    // 3. Create a Value Channel (sticky, reusable)
    bc_opts_ch = Channel.value(bc_opts)
    correct_opts_ch = Channel.value(correct_opts)
    knee_opts_ch = Channel.value(knee_opts)

    ISOSEQ_CORRECT(
        ISOSEQ_REFINE.out.bam,
        barcodes_ch,
        correct_opts_ch
    )

    ISOSEQ_BCSTATS(
        ISOSEQ_CORRECT.out.bam,
        bc_opts_ch
    )

    GENERATE_KNEE_PLOT(
        //ISOSEQ_CORRECT.out.bam,
        ISOSEQ_BCSTATS.out.tsv,
        knee_opts_ch  
    )

    // # Deduplicate reads based on UMIs
    // sort by barcode first, then groupdedup
    SAMTOOLS_SORT(
        ISOSEQ_CORRECT.out.bam
    )

// add --group_dedup_params "--keep-non-real-cells" to keep non-real cells
    ISOSEQ_GROUPDEDUP(
        SAMTOOLS_SORT.out.bam,
        group_dedup_parstr_ch
    )

    PBMM2_ALIGN(
        ISOSEQ_GROUPDEDUP.out.bam,
        reference_ch
    )

    ISOSEQ_COLLAPSE(
        PBMM2_ALIGN.out.bam
    )
    /*
    PIGEON_PREPARE(
        annotation_ch,
        reference_ch
    )
    */

    PIGEON_PREPARE_GFF(
        ISOSEQ_COLLAPSE.out.gff
    )

    classify_input_ch = PIGEON_PREPARE_GFF.out.gff
                        .join(ISOSEQ_COLLAPSE.out.abundance)
    /*
    PIGEON_CLASSIFY(
        classify_input_ch,
        reference_ch,
        PIGEON_PREPARE.out.fai,
        PIGEON_PREPARE.out.gtf_sorted,
        PIGEON_PREPARE.out.gtf_sorted_pgi
    )
    */

    PIGEON_CLASSIFY(
        classify_input_ch,
        reference_ch,
        reference_fai_ch,
        gtf_sorted_ch,
        gtf_sorted_pgi_ch,
        // Pass the new channels here
        cage_bed_ch,
        cage_index_ch,
        poly_a_ch,
        intropolis_ch,
        intropolis_index_ch
    )

    filter_input_prep_ch = PIGEON_CLASSIFY.out.classification
                           .join(PIGEON_CLASSIFY.out.junctions)

    filter_input_ch = filter_input_prep_ch
                      .join(ISOSEQ_COLLAPSE.out.gff)

    PIGEON_FILTER(
        filter_input_ch
    )

    PIGEON_REPORT(
        PIGEON_FILTER.out.classification
    )

    pigeon_seurat_in_prep_ch = ISOSEQ_GROUPDEDUP.out.fasta
                            .join(ISOSEQ_COLLAPSE.out.group)

    pigeon_seurat_in_ch = pigeon_seurat_in_prep_ch
                       .join(PIGEON_FILTER.out.classification)
    
    if ( params.pigeon_remove_ribosomal_mitochondrial ) {
        pigeon_seurat_opts = ""
    } else {
        pigeon_seurat_opts = "--keep-ribo-mito-genes"
    }

    pigeon_seurat_opts_ch = Channel.value(pigeon_seurat_opts)
    
    PIGEON_SEURAT(
                  pigeon_seurat_in_ch,
                  pigeon_seurat_opts_ch
                  )

    multiqc_input_ch = ISOSEQ_REFINE.out.reports.map{ [ it[1], it[2] ] }
                       .mix(
                           LIMA.out.stats.map{ [ it[1], it[2] ] } 
                       )
                       .collect()

    multiqc_input_ch.view()

    MULTIQC(
            multiqc_input_ch
            //ISOSEQ_REFINE.out.reports_json.collect(),
            //ISOSEQ_REFINE.out.reports_csv.collect()
            )
}
