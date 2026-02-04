process CALL_GATK_HC {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'variantcalling/pipeline:latest' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(ref_fasta)
    path(targets_bed)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), emit: gvcf
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def targets_cmd = targets_bed ? "-L $targets_bed" : ""
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        HaplotypeCaller \\
        -R $ref_fasta \\
        -I $bam \\
        $targets_cmd \\
        -O ${prefix}.g.vcf.gz \\
        $args

    # Also create a regular VCF for single-sample analysis
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        GenotypeGVCFs \\
        -R $ref_fasta \\
        -V ${prefix}.g.vcf.gz \\
        -O ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}