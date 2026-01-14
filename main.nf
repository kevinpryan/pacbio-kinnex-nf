// lima
// cDNA primer removal and read orientation

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
    isoseq refine $bam $primers ${bam_id}.fltnc.bam --require-polya --num-threads ${task.cpus}
    """
}

// # Correct single cell barcodes based on an include list

process ISOSEQ_CORRECT {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    path barcodes

    output:
    tuple val(bam_id), path("${bam_id}.corrected.bam"), emit: bam
 
    script:
    """
    isoseq correct --num-threads ${task.cpus} -B $barcodes $bam ${bam_id}.corrected.bam
    """
}


process ISOSEQ_GROUPDEDUP {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)

    output:
    tuple val(bam_id), path("${bam_id}.dedup.bam"), emit: bam
    tuple val(bam_id), path("${bam_id}.dedup.fasta"), emit: fasta
    script:
    """
    isoseq groupdedup --keep-non-real-cells --num-threads ${task.cpus} $bam ${bam_id}.dedup.bam 
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
    path "*.fasta.fai", emit: fai
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

    output:
    tuple val(bam_id), path("${bam_id}_pigeon_seurat"), emit: seurat

    script:
    """
    pigeon make-seurat --dedup $dedup_fasta --group $collapse_group -d ${bam_id}_pigeon_seurat $classification_filtered_lite
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


// WORKFLOW
workflow {
    // Channel setup
    bams_ch = Channel.fromPath(params.bams)
        .map { file -> tuple(file.simpleName, file) }
        .view { "Found sample: ${it[0]}" }  // <--- ADD THIS LINE

    primers_ch = Channel.value(file(params.primers))
    design_ch = Channel.of(params.design)
    barcodes_ch = Channel.value(file(params.barcodes))
    reference_ch = Channel.value(file(params.reference))
    annotation_ch_gzipped = Channel.value(file(params.gtf))
    PREPARE_REFERENCE_GTF(annotation_ch_gzipped)
    annotation_ch = PREPARE_REFERENCE_GTF.out.gtf
    LIMA(
        bams_ch,
        primers_ch
    )
    ISOSEQ_TAG(
        LIMA.out.bam,
        design_ch
    )
    ISOSEQ_REFINE(
        ISOSEQ_TAG.out.bam,
        primers_ch
    )
    ISOSEQ_CORRECT(
        ISOSEQ_REFINE.out.bam,
        barcodes_ch
    )
    // # Deduplicate reads based on UMIs
    // sort by barcode first, then groupdedup
    SAMTOOLS_SORT(
        ISOSEQ_CORRECT.out.bam
    )
    ISOSEQ_GROUPDEDUP(
        SAMTOOLS_SORT.out.bam
    )
    PBMM2_ALIGN(
        ISOSEQ_GROUPDEDUP.out.bam,
        reference_ch
    )
    ISOSEQ_COLLAPSE(
        PBMM2_ALIGN.out.bam
    )
    PIGEON_PREPARE(
        annotation_ch,
        reference_ch
    )
    PIGEON_PREPARE_GFF(
        ISOSEQ_COLLAPSE.out.gff
    )
    classify_input_ch = PIGEON_PREPARE_GFF.out.gff
                        .join(ISOSEQ_COLLAPSE.out.abundance)
    PIGEON_CLASSIFY(
        classify_input_ch,
        reference_ch,
        PIGEON_PREPARE.out.fai,
        PIGEON_PREPARE.out.gtf_sorted,
        PIGEON_PREPARE.out.gtf_sorted_pgi
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
    PIGEON_SEURAT(pigeon_seurat_in_ch)
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

