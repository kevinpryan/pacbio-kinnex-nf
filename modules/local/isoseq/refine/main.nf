// Remove poly(A) tails and concatemer

process ISOSEQ_REFINE {
    label "ISOSEQ_CONTAINER"
    tag "$bam_id"
    input:
    tuple val(bam_id), path(bam)
    path primers
    val parstr
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
                  ${parstr} \
                  --num-threads ${task.cpus}
    """
}

