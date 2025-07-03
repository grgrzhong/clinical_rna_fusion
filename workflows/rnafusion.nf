// Nextflow DSL2 configuration
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP_TRIM                } from '../modules/local/fastp/main'
include { FASTQC                    } from '../modules/local/fastqc/main'
include { STAR_ALIGN_FOR_STARFUSION } from '../modules/local/star/align/starfusion/main.nf'

workflow {

    star_index = file(params.star_index)
    // ctat_resource_lib = file(params.ctat_resource_lib)


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Prepare sample sheet
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Read input CSV
    sample_list = Channel.fromPath(params.input)
        .ifEmpty { exit(1, "Samplesheet not found: ${params.input}") }
        .splitCsv(header: true, quote: '"', strip: true)

    input_samples = sample_list.map { row ->
        // Helper function to clean csv
        def cleanValue = { value ->
            if (value == null || value == 'NA' || value == 'NULL' || value == '') {
                return null
            }
            // Remove quotes and trim whitespace
            return value.toString().replaceAll('^"(.*)"$', '$1').trim()
        }

        // Extract sample info
        def sample_id = cleanValue(row.sample)

        // Clean and validate file paths
        def fastq_1_clean = cleanValue(row.fastq_1)
        def fastq_2_clean = cleanValue(row.fastq_2)
        def fastq_trimmed_1_clean = cleanValue(row.fastq_trimmed_1)
        def fastq_trimmed_2_clean = cleanValue(row.fastq_trimmed_2)

        // Determine available input types
        def has_fastq = fastq_1_clean && fastq_2_clean
        def has_fastq_trimmed = fastq_trimmed_1_clean && fastq_trimmed_2_clean

        // Exit if no input is provided
        if (!has_fastq && !has_fastq_trimmed) {
            error("Either FASTQ or Trimmed FASTQ must be provided for sample ${sample_id}")
        }

        def meta = [
            id: sample_id,
            sample_id: sample_id,
        ]

        // Initialize all file variables with null
        def fastq_1 = null
        def fastq_2 = null
        def fastq_trimmed_1 = null
        def fastq_trimmed_2 = null


        if (has_fastq) {
            // Validate FASTQ files exist
            try {
                fastq_1 = file(fastq_1_clean, checkIfExists: true)
                fastq_2 = file(fastq_2_clean, checkIfExists: true)
            }
            catch (Exception e) {
                error("FASTQ file not found for sample ${sample_id}: ${e.message}")
            }
        }

        if (has_fastq_trimmed) {
            // Validate FASTQ files exist
            try {
                fastq_trimmed_1 = file(fastq_trimmed_1_clean, checkIfExists: true)
                fastq_trimmed_2 = file(fastq_trimmed_2_clean, checkIfExists: true)
            }
            catch (Exception e) {
                error("FASTQ file not found for sample ${sample_id}: ${e.message}")
            }
        }

        // Return unified format with all possible inputs
        return [meta, fastq_1, fastq_2, fastq_trimmed_1, fastq_trimmed_2]
    }

    fastq = input_samples
        .filter { _meta, fastq_1, fastq_2, _fastq_trimmed_1, _fastq_trimmed_2 ->
            fastq_1 != null && fastq_2 != null
        }
        .map { meta, fastq_1, fastq_2, _fastq_trimmed_1, _fastq_trimmed_2 ->
            [meta, fastq_1, fastq_2]
        }

    // fastq_trimmed = input_samples
    //     .filter {
    //         _meta, _fastq_1, _fastq_2, fastq_trimmed_1, fastq_trimmed_2 -> 
    //         fastq_trimmed_1 != null && fastq_trimmed_2 != null
    //     }
    //     .map {
    //         meta, _fastq_1, _fastq_2, fastq_trimmed_1, fastq_trimmed_2 -> 
    //         [meta, fastq_trimmed_1, fastq_trimmed_2]
    //     }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Preprocessing: FASTP trimming
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // FASTP trimming adapter sequence
    fastp_output = FASTP_TRIM(
        fastq
    )

    // Fastqc on the trimmed FASTQ files
    FASTQC(fastp_output.reads)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STAR Fusion workflow
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Run STAR alignment for starfusion
    STAR_ALIGN_FOR_STARFUSION(
        fastp_output.reads,
        star_index,
    )


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Arrib workflow
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // // Run STAR alignment for arriba
    // STAR_ALIGN_FOR_ARRIBA(
    //     fastp_output.reads,
    //     star_index,
    // )
}
