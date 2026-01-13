nextflow.enable.dsl = 2

workflow {
    ch_input = Channel.fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> [ row.sample, [file(row.fastq_1), file(row.fastq_2)] ] }

    // Logic: We take whatever reference is passed via CLI params
    INDEX( file(params.host_fasta), file(params.graft_fasta) )
    SORT_READS( ch_input, INDEX.out.collect() )
}

process INDEX {
    // We save this index in a 'clean' folder to avoid mixing with the old toplevel index
    storeDir "ANALYSIS/xengsort_index_clean"
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
    publishDir "${params.outdir}/sorted_fastqs", mode: 'copy'
    cpus 8
    memory '32 GB'
    input:
        tuple val(sample_id), path(reads)
        path index_files
    output:
        // 1. Human (Tumor) Reads
        path "${sample_id}_human_R1.fq.gz", emit: human_r1
        path "${sample_id}_human_R2.fq.gz", emit: human_r2
        
        // 2. Rat (Host) Reads - CRITICAL for your new Microenvironment analysis
        path "${sample_id}_rat_R1.fq.gz",   emit: rat_r1
        path "${sample_id}_rat_R2.fq.gz",   emit: rat_r2

    script:
    """
    xengsort classify --index index --fastq ${reads[0]} --pairs ${reads[1]} --mode coverage --out ${sample_id} --prefix ${sample_id} -T ${task.cpus}

    # MERGE 1: Human (Graft + Ambiguous)
    cat ${sample_id}*graft*1*fq.gz ${sample_id}*ambig*1*fq.gz > ${sample_id}_human_R1.fq.gz
    cat ${sample_id}*graft*2*fq.gz ${sample_id}*ambig*2*fq.gz > ${sample_id}_human_R2.fq.gz

    # MERGE 2: Rat (Host + Ambiguous)
    # We include ambiguous reads so housekeeping genes are present for normalization
    cat ${sample_id}*host*1*fq.gz ${sample_id}*ambig*1*fq.gz > ${sample_id}_rat_R1.fq.gz
    cat ${sample_id}*host*2*fq.gz ${sample_id}*ambig*2*fq.gz > ${sample_id}_rat_R2.fq.gz
    """
}
