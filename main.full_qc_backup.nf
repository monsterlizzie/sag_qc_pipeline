
nextflow.enable.dsl=2

include { FILE_VALIDATION; PREPROCESS; READ_QC } from './modules/preprocess'
include { ASSEMBLY_SHOVILL; ASSEMBLY_ASSESS; ASSEMBLY_QC } from './modules/assembly.nf'
include { GET_KRAKEN2_DB; TAXONOMY; BRACKEN;TAXONOMY_QC } from './modules/taxonomy'
include { GET_REF_GENOME_BWA_DB as GET_SANGINOSUS_BWA_DB } from './modules/mapping'
include { GET_REF_GENOME_BWA_DB as GET_SINTERMEDIUS_BWA_DB } from './modules/mapping'
include { GET_REF_GENOME_BWA_DB as GET_SCONSTELLATUS_BWA_DB } from './modules/mapping'
include { MAPPING; SAM_TO_SORTED_BAM; SNP_CALL; HET_SNP_COUNT; MAPPING_QC } from './modules/mapping'
include { OVERALL_QC } from './modules/overall_qc.nf'
include { GENERATE_SAMPLE_REPORT; GENERATE_OVERALL_REPORT } from './modules/output.nf'

workflow {

    // Input reads

    raw_read_pairs_ch = Channel.fromFilePairs(

        "${params.reads}/*{_,.}R{1,2}{,_001}.{fq,fastq}{,.gz}",

        checkIfExists: true

    )

    // Validate input files

    FILE_VALIDATION(raw_read_pairs_ch)

    VALID_READS_ch = FILE_VALIDATION.out.result

        .join(raw_read_pairs_ch, failOnDuplicate: true)

        .filter { sample_id, file_validity, reads -> file_validity == 'PASS' }

        .map { sample_id, file_validity, reads ->

            tuple(sample_id, reads)

        }

    // Preprocess reads

    PREPROCESS(VALID_READS_ch)

    // Read QC

    READ_QC(

        PREPROCESS.out.json,

        params.length_low,

        params.depth

    )

    READ_QC_PASSED_READS_ch = READ_QC.out.result

        .join(PREPROCESS.out.processed_reads, failOnDuplicate: true)

        .filter { sample_id, read_qc, read1, read2, unpaired -> read_qc == 'PASS' }

        .map { sample_id, read_qc, read1, read2, unpaired ->

            tuple(sample_id, read1, read2, unpaired)

        }

// Assembly

ASSEMBLY_ch = ASSEMBLY_SHOVILL(READ_QC_PASSED_READS_ch, params.min_contig_length, params.assembler_thread)

ASSEMBLY_ASSESS(ASSEMBLY_ch)

ASSEMBLY_QC(
    ASSEMBLY_ASSESS.out.report
        .join(READ_QC.out.bases, failOnDuplicate: true),
    params.contigs,
    params.length_low,
    params.length_high,
    params.depth
)
// Taxonomy
// Get Kraken2 database
GET_KRAKEN2_DB(params.kraken2_db_remote, params.db)

// Run Kraken2 on READ_QC-passed reads
TAXONOMY(
    GET_KRAKEN2_DB.out.path,
    params.kraken2_memory_mapping,
    READ_QC_PASSED_READS_ch
    )

// Extract dominant SAG species and taxonomy QC
// Run Bracken on Kraken2 report
BRACKEN(
    GET_KRAKEN2_DB.out.path,
    TAXONOMY.out.report,
    params.bracken_read_len,
    params.bracken_level,
    params.bracken_threshold
)

// Extract dominant SAG species and taxonomy QC from Bracken report
TAXONOMY_QC(
    BRACKEN.out.report,
    params.sag_species_percentage,
    params.top_non_sag_percentage
)

// Keep only taxonomy-passed samples and attach species label
TAXONOMY_PASSED_SPECIES_ch = TAXONOMY_QC.out.result
    .join(TAXONOMY_QC.out.species, failOnDuplicate: true)
    .filter { sample_id, taxonomy_qc, species -> taxonomy_qc == 'PASS' }
    .map { sample_id, taxonomy_qc, species ->
        tuple(sample_id, species)
    }

SPECIES_READS_ch = READ_QC_PASSED_READS_ch
    .join(TAXONOMY_PASSED_SPECIES_ch, failOnDuplicate: true)
    .map { sample_id, read1, read2, unpaired, species ->
        tuple(sample_id, read1, read2, unpaired, species)
    }

//Mapping 
// Build species-specific BWA databases
GET_SANGINOSUS_BWA_DB(file(params.sanginosus_ref), file(params.db), "sanginosus")
GET_SINTERMEDIUS_BWA_DB(file(params.sintermedius_ref), file(params.db), "sintermedius")
GET_SCONSTELLATUS_BWA_DB(file(params.sconstellatus_ref), file(params.db), "sconstellatus")


// Map reads to species-specific reference
MAPPING(
    GET_SANGINOSUS_BWA_DB.out.db,
    GET_SINTERMEDIUS_BWA_DB.out.db,
    GET_SCONSTELLATUS_BWA_DB.out.db,
    SPECIES_READS_ch
    )

// Convert SAM to sorted BAM and calculate reference coverage
SAM_TO_SORTED_BAM(MAPPING.out.sam, params.lite)

// Call SNPs with matching species-specific reference
SNP_CALL(
    file(params.sanginosus_ref),
    file(params.sintermedius_ref),
    file(params.sconstellatus_ref),
    SAM_TO_SORTED_BAM.out.sorted_bam,
    params.lite
    )

// Count HetSNPs
HET_SNP_COUNT(SNP_CALL.out.vcf)

// Mapping QC
MAPPING_QC(
    SAM_TO_SORTED_BAM.out.ref_coverage
        .join(HET_SNP_COUNT.out.result, failOnDuplicate: true)
        .map { sample_id, species_from_cov, ref_coverage, species_from_het, het_snp_count ->
            tuple(sample_id, species_from_cov, ref_coverage, het_snp_count)
        },
    params.ref_coverage,
    params.het_snp_site
)

// Determine overall QC result based on Assembly QC, Mapping QC and Taxonomy QC
OVERALL_QC(
    raw_read_pairs_ch.map { sample_id, reads -> tuple(sample_id) }
        .join(FILE_VALIDATION.out.result, failOnDuplicate: true, remainder: true)
        .join(READ_QC.out.result, failOnDuplicate: true, remainder: true)
        .join(ASSEMBLY_QC.out.result, failOnDuplicate: true, remainder: true)
        .join(MAPPING_QC.out.result, failOnDuplicate: true, remainder: true)
        .join(TAXONOMY_QC.out.result, failOnDuplicate: true, remainder: true)
)

// Get assemblies that passed overall QC
OVERALL_QC_PASSED_ASSEMBLIES_ch = OVERALL_QC.out.result
    .join(ASSEMBLY_ch, failOnDuplicate: true)
    .filter { sample_id, overall_qc, assembly -> overall_qc == 'PASS' }
    .map { sample_id, overall_qc, assembly ->
        tuple(sample_id, assembly)
    }

GENERATE_SAMPLE_REPORT(
    raw_read_pairs_ch.map { sample_id, reads -> tuple(sample_id) }
        .join(READ_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(ASSEMBLY_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(MAPPING_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(TAXONOMY_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(OVERALL_QC.out.report, failOnDuplicate: true, remainder: true)
        .map { row ->
            tuple(row[0], row[1..-1].findAll { it != null })
        }
)

GENERATE_OVERALL_REPORT(
    GENERATE_SAMPLE_REPORT.out.report
        .map { sample_id, report -> report }
        .collect()
)
}
