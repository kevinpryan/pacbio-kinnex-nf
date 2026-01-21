process ISOSEQ_BCSTATS {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    publishDir "${params.outdir}/isoseq/bcstats", mode: 'copy'

    input:
    tuple val(bam_id), path(bam)
    val parstr
    output:
    tuple val(bam_id), path("${bam_id}.bcstats.json"), emit: json
    tuple val(bam_id), path("${bam_id}.bcstats.tsv"), emit: tsv
    """
    isoseq bcstats \
           --json ${bam_id}.bcstats.json \
           -o ${bam_id}.bcstats.tsv \
           $parstr \
           --num-threads ${task.cpus} \
           $bam
    """
}
