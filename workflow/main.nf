#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Variant Calling Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/your-org/variant-calling-pipeline
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = [:]

// Validate input parameters
if (params.input) { 
    ch_input = file(params.input) 
} else { 
    exit 1, 'Input samplesheet not specified!' 
}

if (params.ref_fasta) {
    ch_ref_fasta = file(params.ref_fasta)
} else {
    exit 1, 'Reference FASTA not specified!'
}

if (params.known_sites_vcf && params.run_bqsr) {
    ch_known_sites = file(params.known_sites_vcf)
} else {
    ch_known_sites = []
}

if (params.targets_bed) {
    ch_targets = file(params.targets_bed)
} else {
    ch_targets = []
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK       } from './subworkflows/local/input_check'
include { QC_FASTQ          } from './modules/qc_fastq'
include { ALIGN_BWA_MEM2    } from './modules/align_bwa_mem2'
include { SORT_INDEX        } from './modules/sort_index'
include { MARKDUP           } from './modules/markdup'
include { BQSR              } from './modules/bqsr'
include { CALL_GATK_HC      } from './modules/call_gatk_hc'
include { JOINT_GENOTYPE    } from './modules/joint_genotype'
include { FILTER_ANNOTATE   } from './modules/filter_annotate'
include { MULTIQC           } from './modules/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

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
    if (params.run_fastp) {
        QC_FASTQ (
            INPUT_CHECK.out.reads
        )
        ch_processed_reads = QC_FASTQ.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(QC_FASTQ.out.json.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(QC_FASTQ.out.versions.first())
    } else {
        ch_processed_reads = INPUT_CHECK.out.reads
    }

    //
    // MODULE: Align reads with BWA-MEM2
    //
    ALIGN_BWA_MEM2 (
        ch_processed_reads,
        ch_ref_fasta
    )
    ch_versions = ch_versions.mix(ALIGN_BWA_MEM2.out.versions.first())

    //
    // MODULE: Sort and index BAM files
    //
    SORT_INDEX (
        ALIGN_BWA_MEM2.out.sam
    )
    ch_versions = ch_versions.mix(SORT_INDEX.out.versions.first())

    //
    // MODULE: Mark duplicates
    //
    MARKDUP (
        SORT_INDEX.out.bam
    )
    ch_multiqc_files = ch_multiqc_files.mix(MARKDUP.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_versions = ch_versions.mix(MARKDUP.out.versions.first())

    //
    // MODULE: Base Quality Score Recalibration (optional)
    //
    if (params.run_bqsr && ch_known_sites) {
        BQSR (
            MARKDUP.out.bam,
            ch_ref_fasta,
            ch_known_sites
        )
        ch_bam_final = BQSR.out.bam
        ch_versions = ch_versions.mix(BQSR.out.versions.first())
    } else {
        ch_bam_final = MARKDUP.out.bam
    }

    //
    // MODULE: Call variants with GATK HaplotypeCaller
    //
    CALL_GATK_HC (
        ch_bam_final,
        ch_ref_fasta,
        ch_targets
    )
    ch_versions = ch_versions.mix(CALL_GATK_HC.out.versions.first())

    //
    // MODULE: Joint genotyping (optional)
    //
    if (params.run_joint) {
        JOINT_GENOTYPE (
            CALL_GATK_HC.out.gvcf.collect{it[1]},
            ch_ref_fasta
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
        ch_ref_fasta
    )
    ch_versions = ch_versions.mix(FILTER_ANNOTATE.out.versions.first())

    //
    // Collect versions and create software versions file
    //
    ch_versions
        .unique()
        .collectFile(name: 'collated_versions.yml')
        .set { ch_versions_yaml }

    //
    // MODULE: MultiQC
    //
    MULTIQC (
        ch_multiqc_files.collect(),
        [],
        [],
        ch_versions_yaml
    )

}