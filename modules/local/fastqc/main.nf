process FASTQC {
    
    tag "${meta.id}"
    
    label 'process_medium'
    
    input:
    tuple val(meta), path(reads1), path(reads2)
    
    output:
    tuple val(meta), path("*.html"),    emit: html
    tuple val(meta), path("*.zip"),     emit: zip
    path "versions.yml",                emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    fastqc \\
        $args \\
        --threads ${task.cpus} \\
        ${reads1} \\
        ${reads2}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed -e "s/FastQC v//g")
    END_VERSIONS
    """
}