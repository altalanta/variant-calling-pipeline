//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.sample_id   = row.sample_id
    meta.single_end  = false
    meta.platform    = row.platform ?: 'ILLUMINA'
    meta.library     = row.library ?: row.sample_id
    meta.lane        = row.lane ?: '1'

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.read1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.read1}"
    }
    if (!file(row.read2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.read2}"
    }
    fastq_meta = [ meta, [ file(row.read1), file(row.read2) ] ]
    return fastq_meta
}