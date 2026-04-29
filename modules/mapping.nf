// Return species-specific BWA database path and prefix, construct if necessary
process GET_REF_GENOME_BWA_DB {
    label 'bwa_container'
    label 'farm_mid'
    label 'farm_scratchless'

    input:
    path reference
    path db
    val prefix

    output:
    tuple val(prefix), path(bwa_db), emit: db

    script:
    bwa_db="${db}/bwa_${prefix}"
    json="done_bwa_db_${prefix}.json"
    checksum="checksum_${prefix}.md5"
    """
    REFERENCE="$reference"
    DB_LOCAL="$bwa_db"
    PREFIX="$prefix"
    JSON_FILE="$json"
    CHECKSUM_FILE="$checksum"

    source check-create_ref_genome_bwa_db.sh
    """
}


// Map reads to species-specific reference using BWA-MEM
process MAPPING {
    label 'bwa_container'
    label 'farm_mid'

    tag "$sample_id"

    input:
    tuple val(sanginosus_prefix), path(sanginosus_bwa_db)
    tuple val(sintermedius_prefix), path(sintermedius_bwa_db)
    tuple val(sconstellatus_prefix), path(sconstellatus_bwa_db)
    tuple val(sample_id), path(read1), path(read2), path(unpaired), val(species)

    output:
    tuple val(sample_id), val(species), path(sam), emit: sam

    script:
    sam="${sample_id}_mapped.sam"
    """
    if [[ "${species}" == "Streptococcus_anginosus" ]]; then
        BWA_REF="${sanginosus_bwa_db}/${sanginosus_prefix}"
    elif [[ "${species}" == "Streptococcus_intermedius" ]]; then
        BWA_REF="${sintermedius_bwa_db}/${sintermedius_prefix}"
    elif [[ "${species}" == "Streptococcus_constellatus" ]]; then
        BWA_REF="${sconstellatus_bwa_db}/${sconstellatus_prefix}"
    else
        echo "Unknown SAG species: $species"
        exit 1
    fi

    bwa mem -t "`nproc`" "\$BWA_REF" <(zcat -f -- < "$read1") <(zcat -f -- < "$read2") > "$sam"
    """
}


// Convert mapped SAM into BAM and sort it
// Return mapped and sorted BAM, and reference coverage percentage by reads
process SAM_TO_SORTED_BAM {
    label 'samtools_container'
    label 'farm_mid'

    tag "$sample_id"

    input:
    tuple val(sample_id), val(species), path(sam)
    val lite

    output:
    tuple val(sample_id), val(species), path(sorted_bam), emit: sorted_bam
    tuple val(sample_id), val(species), env('COVERAGE'), emit: ref_coverage

    script:
    sorted_bam="${sample_id}_mapped_sorted.bam"
    max_thread="8"
    """
    SAM="$sam"
    BAM="mapped.bam"
    SORTED_BAM="$sorted_bam"
    LITE="$lite"
    MAX_THREAD="$max_thread"

    source convert_sam_to_sorted_bam.sh
    source get_ref_coverage.sh
    """
}


// Call SNPs using species-specific reference
process SNP_CALL {
    label 'bcftools_container'
    label 'farm_mid'

    tag "$sample_id"

    input:
    path sanginosus_ref
    path sintermedius_ref
    path sconstellatus_ref
    tuple val(sample_id), val(species), path(sorted_bam)
    val lite

    output:
    tuple val(sample_id), val(species), path(vcf), emit: vcf

    script:
    vcf="${sample_id}.vcf"
    """
    if [[ "${species}" == "Streptococcus_anginosus" ]]; then
        REFERENCE="$sanginosus_ref"
    elif [[ "${species}" == "Streptococcus_intermedius" ]]; then
        REFERENCE="$sintermedius_ref"
    elif [[ "${species}" == "Streptococcus_constellatus" ]]; then
        REFERENCE="$sconstellatus_ref"
    else
        echo "Unknown SAG species: $species"
        exit 1
    fi

    SORTED_BAM="$sorted_bam"
    VCF="$vcf"
    LITE="$lite"

    source call_snp.sh
    """
}


// Return non-cluster heterozygous SNP site count
process HET_SNP_COUNT {
    label 'python_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), val(species), path(vcf)

    output:
    tuple val(sample_id), val(species), env('OUTPUT'), emit: result

    script:
    het_snp_count_output='output.txt'
    """
    het_snp_count.py "$vcf" 50 "$het_snp_count_output"
    OUTPUT=`cat $het_snp_count_output`
    """
}


// Extract mapping QC information and determine QC result
process MAPPING_QC {
    label 'bash_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), val(species), val(ref_coverage), val(het_snp_count)
    val(qc_ref_coverage)
    val(qc_het_snp_site)

    output:
    tuple val(sample_id), env('MAPPING_QC'), emit: result
    tuple val(sample_id), path(mapping_qc_report), emit: report

    script:
    mapping_qc_report='mapping_qc_report.csv'
    """
    SPECIES="$species"
    COVERAGE="$ref_coverage"
    HET_SNP="$het_snp_count"
    QC_REF_COVERAGE="$qc_ref_coverage"
    QC_HET_SNP_SITE="$qc_het_snp_site"
    MAPPING_QC_REPORT="$mapping_qc_report"

    source get_mapping_qc.sh
    """
}