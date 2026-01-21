process GENERATE_KNEE_PLOT {
    tag "$bam_id"
    publishDir "${params.outdir}/knee_plots", mode: 'copy'
    label "PYTHON_KNEE_CONTAINER" 

    input:
    //tuple val(bam_id), path(bam)
    tuple val(bam_id), path(bc_stats_tsv)
    val knee_opts
    output:
    path "${bam_id}.knee.png"
    // Only include this if you correspond it to a flag in the script below
    //path "${bam_id}_reads_per_barcode.txt" 

    script:
    """
    plot_knees.py \
        --tsv ${bc_stats_tsv} \
        ${knee_opts} \
        --output ${bam_id}
    """
}
