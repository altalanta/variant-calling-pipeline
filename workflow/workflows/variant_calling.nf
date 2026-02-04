/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowVariantCalling.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { INPUT_CHECK       } from '../subworkflows/local/input_check'

//
// MODULES: Local to the pipeline
//
include { QC_FASTQ          } from '../modules/qc_fastq'
include { ALIGN_BWA_MEM2    } from '../modules/align_bwa_mem2'
include { SORT_INDEX        } from '../modules/sort_index'
include { MARKDUP           } from '../modules/markdup'
include { BQSR              } from '../modules/bqsr'
include { CALL_GATK_HC      } from '../modules/call_gatk_hc'
include { JOINT_GENOTYPE    } from '../modules/joint_genotype'
include { FILTER_ANNOTATE   } from '../modules/filter_annotate'
include { MULTIQC           } from '../modules/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VARIANT_CALLING {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC/fastp on raw reads
    //
    QC_FASTQ (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(QC_FASTQ.out.versions)

    //
    // MODULE: Align reads with BWA-MEM2
    //
    ALIGN_BWA_MEM2 (
        QC_FASTQ.out.reads,
        params.ref_fasta
    )
    ch_versions = ch_versions.mix(ALIGN_BWA_MEM2.out.versions)

    //
    // MODULE: Sort and index BAM files
    //
    SORT_INDEX (
        ALIGN_BWA_MEM2.out.sam
    )
    ch_versions = ch_versions.mix(SORT_INDEX.out.versions)

    //
    // MODULE: Mark duplicates
    //
    MARKDUP (
        SORT_INDEX.out.bam
    )
    ch_versions = ch_versions.mix(MARKDUP.out.versions)

    //
    // MODULE: Base Quality Score Recalibration (optional)
    //
    if (params.run_bqsr && params.known_sites_vcf) {
        BQSR (
            MARKDUP.out.bam,
            params.ref_fasta,
            params.known_sites_vcf
        )
        ch_bam_final = BQSR.out.bam
        ch_versions = ch_versions.mix(BQSR.out.versions)
    } else {
        ch_bam_final = MARKDUP.out.bam
    }

    //
    // MODULE: Call variants with GATK HaplotypeCaller
    //
    CALL_GATK_HC (
        ch_bam_final,
        params.ref_fasta,
        params.targets_bed
    )
    ch_versions = ch_versions.mix(CALL_GATK_HC.out.versions)

    //
    // MODULE: Joint genotyping (optional)
    //
    if (params.run_joint) {
        JOINT_GENOTYPE (
            CALL_GATK_HC.out.gvcf.collect(),
            params.ref_fasta
        )
        ch_vcf_raw = JOINT_GENOTYPE.out.vcf
        ch_versions = ch_versions.mix(JOINT_GENOTYPE.out.versions)
    } else {
        ch_vcf_raw = CALL_GATK_HC.out.vcf
    }

    //
    // MODULE: Filter and annotate variants
    //
    FILTER_ANNOTATE (
        ch_vcf_raw,
        params.ref_fasta
    )
    ch_versions = ch_versions.mix(FILTER_ANNOTATE.out.versions)

    //
    // MODULE: MultiQC
    //
    MULTIQC (
        QC_FASTQ.out.qc.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    multiqc_report = MULTIQC.out.report.toList()
}