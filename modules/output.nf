process GENERATE_SAMPLE_REPORT {
    label 'python_container'
    label 'farm_low'

    tag "$sample_id"

    publishDir "${params.output}/qc_reports", mode: 'copy'

    input:
    tuple val(sample_id), path(reports)

    output:
    tuple val(sample_id), path("${sample_id}_qc.csv"), emit: report

    script:
    sample_report="${sample_id}_qc.csv"
    """
    python3 ${projectDir}/bin/generate_sample_report.py \
        "$sample_id" \
        "$sample_report" \
        ${reports.join(' ')}
    """
}


process GENERATE_OVERALL_REPORT {
    label 'python_container'
    label 'farm_low'

    publishDir "${params.output}", mode: 'copy'

    input:
    path sample_reports

    output:
    path "summary.csv", emit: report

    script:
    """
    python3 ${projectDir}/bin/generate_overall_report.py \
        summary.csv \
        ${sample_reports.join(' ')}
    """
}