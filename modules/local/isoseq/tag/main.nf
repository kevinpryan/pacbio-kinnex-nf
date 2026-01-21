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

