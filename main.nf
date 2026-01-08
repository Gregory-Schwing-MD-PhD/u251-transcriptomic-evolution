nextflow.enable.dsl = 2

workflow {
    // 1. Load samples from the CSV
    ch_input = Channel.fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> [ row.sample, [file(row.fastq_1), file(row.fastq_2)] ] }

    // 2. Build the combined Rat/Human index
    INDEX( file(params.host_fasta), file(params.graft_fasta) )

    // 3. Run deconvolution on all samples in parallel
    CLASSIFY( ch_input, INDEX.out.collect() )
    
    // 4. Generate the final MultiQC report
    MULTIQC( CLASSIFY.out.logs.collect() )
}

process INDEX {
    storeDir "ANALYSIS/xengsort_index"

    input:
    path host
    path graft

    output:
    path "index*"

    script:
    """
    xengsort index \
        -n 4500000000 \
        --index index \
        --host $host \
        --graft $graft
    """
}

process CLASSIFY {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastqs", mode: 'copy', pattern: "*.fq.gz"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(reads)
    path index_files

    output: 
    path "${sample_id}_human_plus_ambig.fq.gz", emit: h_ambig
    path "${sample_id}_rat_plus_ambig.fq.gz",   emit: r_ambig
    path "${sample_id}_human_pure.fq.gz",       emit: h_pure
    path "${sample_id}_rat_pure.fq.gz",         emit: r_pure
    path "${sample_id}_ambig_only.fq.gz",       emit: ambig_only
    path "${sample_id}.xengsort.log",           emit: logs

    script:
    """
    # 1. Run deconvolution to generate pure and ambiguous streams
    xengsort classify \
        --index index \
        --fastq ${reads[0]} \
        --pairs ${reads[1]} \
        --filter graft \
        --out-graft ${sample_id}_human_pure.fq.gz \
        --out-host ${sample_id}_rat_pure.fq.gz \
        --out-ambiguous ${sample_id}_ambig_only.fq.gz \
        --prefix ${sample_id}_ > ${sample_id}.xengsort.log 2>&1

    # 2. Create the "Greedy" Human set (Human + Ambiguous)
    cat ${sample_id}_human_pure.fq.gz ${sample_id}_ambig_only.fq.gz > ${sample_id}_human_plus_ambig.fq.gz

    # 3. Create the "Greedy" Rat set (Rat + Ambiguous)
    cat ${sample_id}_rat_pure.fq.gz ${sample_id}_ambig_only.fq.gz > ${sample_id}_rat_plus_ambig.fq.gz
    """
}

process MULTIQC {
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path xeng_logs

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . --title "Xengsort Deconvolution Report" --filename "multiqc_report.html"
    """
}
