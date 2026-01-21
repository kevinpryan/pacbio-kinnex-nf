process ISOSEQ_GROUPDEDUP {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)

    output:
    tuple val(bam_id), path("${bam_id}.dedup.bam"), emit: bam
    tuple val(bam_id), path("${bam_id}.dedup.fasta"), emit: fasta
    val(parstr)

    script:
    """
    isoseq groupdedup \
           --num-threads ${task.cpus} \
           ${parstr} \
           $bam \
           ${bam_id}.dedup.bam
    """
}

// --keep-non-real-cells
