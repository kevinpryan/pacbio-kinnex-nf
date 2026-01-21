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

