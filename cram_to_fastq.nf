#!/usr/bin/env nextflow
crams = Channel.fromPath(params.cramsPath)

process cram_to_fastq{
    publishDir "${params.outdir}/cram_to_fastq_results/", mode: 'copy'
    memory '7 GB'
    cpus 2

    input:
    file cram_file from crams
    path "ref.fa" from params.fa_ref

    output:
    file "${cram_file.simpleName}_1.fastq.gz"
    file "${cram_file.simpleName}_2.fastq.gz"

    script:
    """
    samtools collate $cram_file ${cram_file.simpleName}.collated
    samtools fastq -c 6 -1 ${cram_file.simpleName}_1.fastq.gz -2 ${cram_file.simpleName}_2.fastq.gz --reference ref.fa ${cram_file.simpleName}.collated.bam
    """
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops ... something went wrong" )
}