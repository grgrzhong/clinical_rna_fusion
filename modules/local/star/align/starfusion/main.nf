process STAR_ALIGN_FOR_STARFUSION {
    
    tag "$meta.id"
    
    label 'process_high'
    
    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path star_index

    output:
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*.tab')            , optional:true, emit: tab
    tuple val(meta), path('*.SJ.out.tab')     , optional:true, emit: spl_junc_tab
    tuple val(meta), path('*.ReadsPerGene.out.tab')  , optional:true, emit: read_per_gene_tab
    tuple val(meta), path('*.out.junction')          , optional:true, emit: junction

    tuple val(meta), path("*Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("*Aligned.sortedByCoord.out.bam.bai"), emit: bai
    tuple val(meta), path("*.log")     , emit: log 
    
    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // # Clean up any existing temp directory
    // rm -rf /tmp/star_${prefix}/
    // --outTmpDir /tmp/star_${prefix}/ \\
    """
    
    export TMPDIR=/tmp

    STAR \\
        --genomeDir ${star_index} \\
        --readFilesIn ${fastq_1} ${fastq_2} \\
        --outReadsUnmapped None \\
        --runThreadN  ${task.cpus} \\
        --twopassMode Basic \\
        --readFilesCommand "gunzip -c" \\
        --outSAMstrandField intronMotif \\
        --outSAMunmapped Within \\
        --chimSegmentMin 12 \\
        --chimJunctionOverhangMin 8 \\
        --chimOutJunctionFormat 1 \\
        --alignSJDBoverhangMin 10 \\
        --alignMatesGapMax 100000 \\
        --alignIntronMax 100000 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --outSAMattrRGline ID:GRPundef SM:$prefix \\
        --chimMultimapScoreRange 3 \\
        --chimScoreJunctionNonGTAG -4 \\
        --chimMultimapNmax 20 \\
        --chimNonchimScoreDropMin 10 \\
        --peOverlapNbasesMin 12 \\
        --peOverlapMMp 0.1 \\
        --alignInsertionFlush Right \\
        --alignSplicedMateMapLminOverLmate 0 \\
        --alignSplicedMateMapLmin 30 \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${prefix}. \\
        --quantMode GeneCounts \\
        $args \\
        2> ${prefix}.staralignment.log
    
    // Index the output BAM file
    samtools index ${prefix}.Aligned.sortedByCoord.out.bam

    # Clean up temporary directory
    // rm -rf /tmp/star_${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
