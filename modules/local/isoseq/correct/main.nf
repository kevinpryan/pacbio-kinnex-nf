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
