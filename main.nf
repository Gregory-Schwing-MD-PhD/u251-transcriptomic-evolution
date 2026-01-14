nextflow.enable.dsl = 2

// Defaults set to null so they MUST be provided via command line (--host_fasta, etc)
params.input = null
params.outdir = "ANALYSIS"
params.host_fasta = null
params.graft_fasta = null

workflow {
    // Check if required params are provided
    if (!params.host_fasta || !params.graft_fasta) {
        error "Error: --host_fasta and --graft_fasta must be specified on the command line."
    }

    ch_input = Channel.fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> [ row.sample, [file(row.fastq_1), file(row.fastq_2)] ] }

    INDEX( file(params.host_fasta), file(params.graft_fasta) )

    SORT_READS( ch_input, INDEX.out.index_files.collect() )

    MULTIQC( SORT_READS.out.stats.collect() )
}

process INDEX {
    storeDir "${params.outdir}/xengsort_index_clean"
    cpus 16
    memory '80 GB'

    input:
        path host
        path graft
    output:
        path "index*", emit: index_files
    script:
    """
    xengsort index -k 29 -n 4500000000 --subtables 15 --fill 0.85 --index index --host $host --graft $graft --weakthreads ${task.cpus}
    """
}

process SORT_READS {
    tag "${sample_id}"
    // Pattern preserved for downstream usage
    publishDir "${params.outdir}/sorted_fastqs", mode: 'copy', pattern: "*_human_R*.fq.gz"
    publishDir "${params.outdir}/sorted_fastqs", mode: 'copy', pattern: "*_rat_R*.fq.gz"
    publishDir "${params.outdir}/xengsort_out", mode: 'copy', pattern: "*.txt"
    
    cpus 8
    memory '32 GB'

    input:
        tuple val(sample_id), path(reads)
        path index_files

    output:
        path "${sample_id}_human_R*.fq.gz", emit: human_reads
        path "${sample_id}_rat_R*.fq.gz",   emit: rat_reads
        path "${sample_id}.txt",            emit: stats

    script:
    """
    # 1. Run xengsort
    # FIXED: Removed invalid flags (--out-graft/host/both) and --classification count
    # FIXED: Redirected output to .txt to capture stats table
    
    xengsort classify \\
        --index index \\
        --fastq ${reads[0]} \\
        --pairs ${reads[1]} \\
        --mode coverage \\
        --out ${sample_id} \\
        --compression gz \\
        -T ${task.cpus} > ${sample_id}.txt 2>&1

    # 2. Merge Steps
    # Merging 'Both' (conserved) into species bins.
    # Using wildcards to safely handle .1.fq.gz vs _1.fq.gz variations

    # Human (Graft + Both)
    cat ${sample_id}*graft*1*fq.gz ${sample_id}*both*1*fq.gz > ${sample_id}_human_R1.fq.gz
    cat ${sample_id}*graft*2*fq.gz ${sample_id}*both*2*fq.gz > ${sample_id}_human_R2.fq.gz

    # Rat (Host + Both)
    cat ${sample_id}*host*1*fq.gz ${sample_id}*both*1*fq.gz > ${sample_id}_rat_R1.fq.gz
    cat ${sample_id}*host*2*fq.gz ${sample_id}*both*2*fq.gz > ${sample_id}_rat_R2.fq.gz
    """
}

process MULTIQC {
    publishDir "${params.outdir}/results_therapy", mode: 'copy'

    input:
        path xengsort_logs

    output:
        path "U251_Final_Report.html"

    script:
    """
    multiqc . \\
        --force \\
        --title "U251 Transcriptomic Evolution" \\
        --filename "U251_Final_Report.html"
    """
}
