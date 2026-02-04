process FILTER_ANNOTATE {
    tag "$meta.sample_id"
    label 'process_low'

    conda "bioconda::gatk4=4.4.0.0 bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b70a85927d9c2d23df14b8b3d8d5bb3b01e8b492:1.0--0' :
        'variantcalling/pipeline:latest' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(ref_fasta)

    output:
    tuple val(meta), path("*.snv_indel.vcf.gz"), path("*.snv_indel.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    
    """
    # Apply hard filters
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        VariantFiltration \\
        -R $ref_fasta \\
        -V $vcf \\
        --filter-expression "QD < 2.0" \\
        --filter-name "QD2" \\
        --filter-expression "QUAL < 30.0" \\
        --filter-name "QUAL30" \\
        --filter-expression "SOR > 3.0" \\
        --filter-name "SOR3" \\
        --filter-expression "FS > 60.0" \\
        --filter-name "FS60" \\
        --filter-expression "MQ < 40.0" \\
        --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5" \\
        --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" \\
        --filter-name "ReadPosRankSum-8" \\
        -O ${prefix}.filtered.vcf.gz

    # Normalize variants with bcftools
    bcftools norm \\
        -f $ref_fasta \\
        -m-both \\
        -O z \\
        -o ${prefix}.snv_indel.vcf.gz \\
        ${prefix}.filtered.vcf.gz

    # Index the final VCF
    tabix -p vcf ${prefix}.snv_indel.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}