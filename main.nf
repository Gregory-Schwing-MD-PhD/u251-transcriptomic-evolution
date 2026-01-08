nextflow.enable.dsl = 2

workflow {
    // 1. Load samples from CSV
    ch_input = Channel.fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> [ row.sample, [file(row.fastq_1), file(row.fastq_2)] ] }

    // 2. Build the Rat/Human combined index
    INDEX( file(params.host_fasta), file(params.graft_fasta) )

    // 3. Parallel classification (Human vs Rat)
    CLASSIFY( ch_input, INDEX.out.collect() )
    
    // 4. Summarize all logs
    MULTIQC( CLASSIFY.out.logs.collect() )
}

process INDEX {
    storeDir "ANALYSIS/xengsort_index"
    input:  path host; path graft
    output: path "index*"
    script: "xengsort index -n 4500000000 --index index --host $host --graft $graft"
}

process CLASSIFY {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastqs", mode: 'copy', pattern: "*.fq.gz"

    input:
    tuple val(sample_id), path(reads)
    path index_files

    output: 
    path "${sample_id}_human.fq.gz",  emit: human_fastqs
    path "${sample_id}_rat.fq.gz",    emit: rat_fastqs
    path "${sample_id}.xengsort.log", emit: logs

    script:
    """
    xengsort classify \
        --index index \
        --fastq ${reads[0]} \
        --pairs ${reads[1]} \
        --out-graft ${sample_id}_human.fq.gz \
        --out-host ${sample_id}_rat.fq.gz \
        --prefix ${sample_id}_ > ${sample_id}.xengsort.log 2>&1
    """
}

process MULTIQC {
    publishDir "${params.outdir}/reports", mode: 'copy'
    input:  path xeng_logs
    output: path "multiqc_report.html"
    script: "multiqc . --title 'Xengsort Deconvolution Report' --filename 'multiqc_report.html'"
}
