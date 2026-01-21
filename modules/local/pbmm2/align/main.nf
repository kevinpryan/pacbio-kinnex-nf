process PBMM2_ALIGN {
    label "PBMM2_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    path ref

    output:
    tuple val(bam_id), path("${bam_id}.aligned.sorted.bam"), emit: bam

    scri:pt
    """
    pbmm2 align --preset ISOSEQ \
                --num-threads ${task.cpus} \
                --sort \
                $bam \
                $ref \
                ${bam_id}.aligned.sorted.bam
    """
}

