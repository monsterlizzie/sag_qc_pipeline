nextflow.enable.dsl=2

params.reads = null
params.output = "mapping_only_output"

params.sanginosus_ref = "$projectDir/data/S_anginosus_GCF_900636475.1_42197_F01.fna"
params.sintermedius_ref = "$projectDir/data/S_intermedius_GCF_900475975.1_46931_F01.fna"
params.sconstellatus_ref = "$projectDir/data/S_constellatus_GCF_016725005.1_ASM1672500v1.fna"

params.db = "$projectDir/db"
params.lite = false
params.ref_coverage = 70
params.het_snp_site = 40

if (!params.reads) {
    error "Please provide reads with --reads, e.g. --reads 'input/*_{R1,R2}.fastq.gz'"
}


process GET_REF_GENOME_BWA_DB {
    label 'bwa_container'
    label 'farm_mid'
    label 'farm_scratchless'

    tag "$species"

    input:
    tuple val(species), path(reference), path(db), val(prefix)

    output:
    tuple val(species), val(prefix), path(bwa_db), emit: db

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


process MAPPING {
    label 'bwa_container'
    label 'farm_mid'

    tag "${sample_id}_${species}"

    input:
    tuple val(sample_id), path(read1), path(read2), path(unpaired), val(species), val(prefix), path(bwa_db)

    output:
    tuple val(sample_id), val(species), path(sam), emit: sam

    script:
    sam="${sample_id}_${species}_mapped.sam"

    """
    BWA_REF="${bwa_db}/${prefix}"

    bwa mem -t ${task.cpus} "\$BWA_REF" \
        <(zcat -f -- "$read1") \
        <(zcat -f -- "$read2") \
        > "$sam"
    """
}


process SAM_TO_SORTED_BAM {
    label 'samtools_container'
    label 'farm_mid'

    tag "${sample_id}_${species}"

    publishDir "${params.output}/bam", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(sample_id), val(species), path(sam)
    val lite

    output:
    tuple val(sample_id), val(species), path(sorted_bam), emit: sorted_bam
    tuple val(sample_id), val(species), env('COVERAGE'), emit: ref_coverage

    script:
    sorted_bam="${sample_id}_${species}_mapped_sorted.bam"
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


process SNP_CALL {
    label 'bcftools_container'
    label 'farm_mid'

    tag "${sample_id}_${species}"

    input:
    path sanginosus_ref
    path sintermedius_ref
    path sconstellatus_ref
    tuple val(sample_id), val(species), path(sorted_bam)
    val lite

    output:
    tuple val(sample_id), val(species), path(vcf), emit: vcf

    script:
    vcf="${sample_id}_${species}.vcf"

    """
    if [[ "${species}" == "Streptococcus_anginosus" ]]; then
        REFERENCE="$sanginosus_ref"
    elif [[ "${species}" == "Streptococcus_intermedius" ]]; then
        REFERENCE="$sintermedius_ref"
    elif [[ "${species}" == "Streptococcus_constellatus" ]]; then
        REFERENCE="$sconstellatus_ref"
    else
        echo "Unknown SAG species: ${species}"
        exit 1
    fi

    SORTED_BAM="$sorted_bam"
    VCF="$vcf"
    LITE="$lite"

    source call_snp.sh
    """
}


process HET_SNP_COUNT {
    label 'python_container'
    label 'farm_low'

    tag "${sample_id}_${species}"

    input:
    tuple val(sample_id), val(species), path(vcf)

    output:
    tuple val(sample_id), val(species), env('OUTPUT'), emit: result

    script:
    het_snp_count_output='output.txt'

    """
    het_snp_count.py "$vcf" 50 "$het_snp_count_output"
    OUTPUT=\$(cat "$het_snp_count_output")
    """
}


process MAPPING_QC {
    label 'bash_container'
    label 'farm_low'

    tag "${sample_id}_${species}"

    publishDir "${params.output}/mapping_qc", mode: 'copy'

    input:
    tuple val(sample_id), val(species), val(ref_coverage), val(het_snp_count)
    val qc_ref_coverage
    val qc_het_snp_site

    output:
    tuple val(sample_id), val(species), env('MAPPING_QC'), emit: result
    tuple val(sample_id), val(species), path(mapping_qc_report), emit: report

    script:
    mapping_qc_report="${sample_id}_${species}_mapping_qc_report.csv"

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


process COMBINE_MAPPING_QC {
    label 'bash_container'
    label 'farm_low'

    publishDir "${params.output}", mode: 'copy'

    input:
    path reports

    output:
    path "mapping_only_summary.csv"

    script:
    """
    set -euo pipefail

    first=\$(ls *.csv | head -n 1)

    old_header=\$(head -n 1 "\$first")
    echo "\\"Sample_ID\\",\${old_header}" > mapping_only_summary.csv

    for f in *_mapping_qc_report.csv; do
        base=\$(basename "\$f" _mapping_qc_report.csv)
        sample_id=\${base%%_Streptococcus_*}

        tail -n +2 "\$f" | awk "NF" | while IFS= read -r line; do
            echo "\\""\${sample_id}"\\",\${line}" >> mapping_only_summary.csv
        done
    done
    """
}


workflow {

    reads_ch = Channel
        .fromFilePairs(params.reads, flat: true)
        .map { sample_id, read1, read2 ->
            tuple(sample_id, read1, read2, [], "NA")
        }

    ref_ch = Channel.of(
        tuple("Streptococcus_anginosus", file(params.sanginosus_ref), file(params.db), "sanginosus"),
        tuple("Streptococcus_intermedius", file(params.sintermedius_ref), file(params.db), "sintermedius"),
        tuple("Streptococcus_constellatus", file(params.sconstellatus_ref), file(params.db), "sconstellatus")
    )

    bwa_db_ch = GET_REF_GENOME_BWA_DB(ref_ch).db

    reads_x_db_ch = reads_ch
        .combine(bwa_db_ch)
        .map { sample_id, read1, read2, unpaired, old_species, species, prefix, bwa_db ->
            tuple(sample_id, read1, read2, unpaired, species, prefix, bwa_db)
        }

    MAPPING(reads_x_db_ch)

    SAM_TO_SORTED_BAM(
        MAPPING.out.sam,
        Channel.value(params.lite)
    )

    sanginosus_ref_ch = Channel.value(file(params.sanginosus_ref))
    sintermedius_ref_ch = Channel.value(file(params.sintermedius_ref))
    sconstellatus_ref_ch = Channel.value(file(params.sconstellatus_ref))

    SNP_CALL(
        sanginosus_ref_ch,
        sintermedius_ref_ch,
        sconstellatus_ref_ch,
        SAM_TO_SORTED_BAM.out.sorted_bam,
        Channel.value(params.lite)
    )

    HET_SNP_COUNT(SNP_CALL.out.vcf)

    mapping_qc_input_ch = SAM_TO_SORTED_BAM.out.ref_coverage
        .join(HET_SNP_COUNT.out.result)
        .map { sample_id, species, ref_coverage, species2, het_snp_count ->
            tuple(sample_id, species, ref_coverage, het_snp_count)
        }

    MAPPING_QC(
        mapping_qc_input_ch,
        Channel.value(params.ref_coverage),
        Channel.value(params.het_snp_site)
    )

    COMBINE_MAPPING_QC(
        MAPPING_QC.out.report
            .map { sample_id, species, report -> report }
            .collect()
    )
}