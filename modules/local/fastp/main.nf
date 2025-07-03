process FASTP_TRIM {
    
    tag "${meta.id}"
    
    label 'process_medium'
    
    input:
    tuple val(meta), path(reads1), path(reads2)
    
    output:
    tuple val(meta), path("*_trimmed_1.fastq.gz"), path("*_trimmed_2.fastq.gz"), emit: reads
    path "*.json", emit: json
    path "*.html", emit: html
    path "*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fastp \\
        --in1 ${reads1} \\
        --in2 ${reads2} \\
        --out1 ${prefix}_trimmed_1.fastq.gz \\
        --out2 ${prefix}_trimmed_2.fastq.gz \\
        --json ${prefix}.json \\
        --html ${prefix}.html \\
        --detect_adapter_for_pe \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \\
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \\
        --dont_eval_duplication \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --trim_poly_x \\
        2> ${prefix}.fastp.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}