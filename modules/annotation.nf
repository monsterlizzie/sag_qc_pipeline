// Return Bakta database path, download and create database if necessary
process GET_BAKTA_DB {
    label 'bash_container'
    label 'farm_low'
    label 'farm_scratchless'
    label 'farm_slow'

    input:
    val db_remote
    path db

    output:
    path bakta_db, emit: path

    script:
    bakta_db="${db}/bakta"
    json='done_bakta.json'
    checksum='checksum.md5'
    """
    DB_REMOTE="$db_remote"
    DB_LOCAL="$bakta_db"
    JSON_FILE="$json"
    CHECKSUM_FILE="$checksum"

    source check-download_bakta_db.sh
    """
}


// Run Bakta annotation on QC-pass SAG assemblies
process ANNOTATE {
    label 'bakta_container'
    label 'farm_high'

    tag "$sample_id"

    input:
    path bakta_db
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.gff3"), emit: gff
    tuple val(sample_id), path("${sample_id}.gbff"), emit: gbff
    tuple val(sample_id), path("${sample_id}.faa"), emit: faa
    tuple val(sample_id), path("${sample_id}.ffn"), emit: ffn
    tuple val(sample_id), path("${sample_id}.json"), emit: json

    script:
    """
    bakta --db "$bakta_db" \
          --prefix "$sample_id" \
          --threads "`nproc`" \
          --skip-plot \
          "$assembly"
    """
}