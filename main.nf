nextflow.enable.dsl = 2

workflow {
    ch_input = Channel.fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> [ row.sample, [file(row.fastq_1), file(row.fastq_2)] ] }

    // Pass the fasta files from params into the INDEX process
    INDEX( file(params.host_fasta), file(params.graft_fasta) )
    SORT_READS( ch_input, INDEX.out.collect() )
}

process INDEX {
    storeDir "ANALYSIS/xengsort_index"
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
    path "${sample_id}_human_R1.fq.gz", emit: human_r1
    path "${sample_id}_human_R2.fq.gz", emit: human_r2
    path "${sample_id}_rat_R1.fq.gz",   emit: rat_r1
    path "${sample_id}_rat_R2.fq.gz",   emit: rat_r2

    script:
    """
    xengsort classify \
        --index index \
        --fastq ${reads[0]} \
        --pairs ${reads[1]} \
        --mode coverage \
        --out ${sample_id}_ \
        --prefix ${sample_id}_ \
        -T ${task.cpus}

    # Human Paired Files (Graft + Ambiguous)
    cat ${sample_id}_graft_1.fastq.gz ${sample_id}_ambig_1.fastq.gz > ${sample_id}_human_R1.fq.gz
    cat ${sample_id}_graft_2.fastq.gz ${sample_id}_ambig_2.fastq.gz > ${sample_id}_human_R2.fq.gz

    # Rat Paired Files (Host + Ambiguous)
    cat ${sample_id}_host_1.fastq.gz ${sample_id}_ambig_1.fastq.gz > ${sample_id}_rat_R1.fq.gz
    cat ${sample_id}_host_2.fastq.gz ${sample_id}_ambig_2.fastq.gz > ${sample_id}_rat_R2.fq.gz

    rm -f ${sample_id}_*.fastq.gz
    """
}
