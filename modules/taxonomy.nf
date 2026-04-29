// Return Kraken2 database path, download if necessary
process GET_KRAKEN2_DB {
    label 'bash_container'
    label 'farm_low'
    label 'farm_scratchless'
    label 'farm_slow'

    input:
    val remote
    path db

    output:
    path kraken2_db, emit: path

    script:
    kraken2_db="${db}/kraken2"
    json='done_kraken.json'
    checksum='checksum.md5'
    """
    DB_REMOTE="$remote"
    DB_LOCAL="$kraken2_db"
    JSON_FILE="$json"
    CHECKSUM_FILE="$checksum"

    source check-download_kraken2_db.sh
    """
}


// Run Kraken2 to classify reads
process TAXONOMY {
    label 'kraken2_container'
    label 'farm_high'

    tag "$sample_id"

    input:
    path kraken2_db
    val kraken2_memory_mapping
    tuple val(sample_id), path(read1), path(read2), path(unpaired)

    output:
    tuple val(sample_id), path(report), emit: report

    script:
    report='kraken2_report.txt'

    if (kraken2_memory_mapping == true)
        """
        kraken2 --threads "`nproc`" --use-names --memory-mapping --db "$kraken2_db" --paired "$read1" "$read2" --report "$report" --output -
        """
    else if (kraken2_memory_mapping == false)
        """
        kraken2 --threads "`nproc`" --use-names --db "$kraken2_db" --paired "$read1" "$read2" --report "$report" --output -
        """
    else
        error "The value for --kraken2_memory_mapping is not valid."
}


// Run Bracken to estimate species abundance
process BRACKEN {
    label 'bracken_container'
    label 'farm_high'

    tag "$sample_id"

    input:
    path kraken2_db
    tuple val(sample_id), path(kraken2_report)
    val read_len
    val classification_level
    val threshold

    output:
    tuple val(sample_id), path("${sample_id}.bracken_report.txt"), emit: report

    script:
    bracken_report="${sample_id}.bracken_report.txt"
    """
    bracken -d "$kraken2_db" \
            -i "$kraken2_report" \
            -o "$bracken_report" \
            -r "$read_len" \
            -l "$classification_level" \
            -t "$threshold"
    """
}


// Extract dominant SAG species from Bracken and determine taxonomy QC
process TAXONOMY_QC {
    label 'bash_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bracken_report)
    val(qc_sag_species_percentage)
    val(qc_top_non_sag_percentage)

    output:
    tuple val(sample_id), env('TAXONOMY_QC'), emit: result
    tuple val(sample_id), env('SAG_SPECIES'), emit: species
    tuple val(sample_id), path(taxonomy_qc_report), emit: report

    script:
    taxonomy_qc_report='taxonomy_qc_report.csv'
    """
    BRACKEN_REPORT="$bracken_report"
    QC_SAG_SPECIES_PERCENTAGE="$qc_sag_species_percentage"
    QC_TOP_NON_SAG_PERCENTAGE="$qc_top_non_sag_percentage"
    TAXONOMY_QC_REPORT="$taxonomy_qc_report"

    source get_taxonomy_qc.sh
    """
}