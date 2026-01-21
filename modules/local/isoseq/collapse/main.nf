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
