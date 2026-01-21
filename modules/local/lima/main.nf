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
