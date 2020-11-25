#!/usr/bin/env nextflow

def helpMessage(){
    log.info "help message"
}

if (params.help) {
    helpMessage()
    exit 0
}

hs2_indices = Channel.fromPath("${params.hisat2_index}*")
    .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }

Channel.fromPath(params.samples_path)
    .ifEmpty { error "Cannot find any samples_path file in: ${params.samples_path}" }
    .splitCsv(header: false, sep: '\t', strip: true)
    .map{row -> [ row[0], row[1] , row[2] ]}
    .set { align_fastq } 

process alignWithHisat2 { 
    // increased resources requirements
    // label "process_low"

    input: 
    file hs2_indices from hs2_indices.collect()
    tuple val(sample), path(read_1), path(read_2) from align_fastq

    output:
    tuple val(sample), path("${sample}.bam") into aligned_bam
    
    script:
    index_base = hs2_indices[0].toString() - ~/.\d.ht2/

    """
    hisat2 -x ${index_base} -1 ${read_1} -2 ${read_2} \\
        | samtools view -b > ${sample}.bam
    """

}

process sortBam {
    input:
    tuple val(sample), path(aligned) from aligned_bam

    output:
    tuple val(sample), path("${sample}_sorted.bam") into bin_scores

    """
    samtools sort ${aligned} -O bam -o ${sample}_sorted.bam
    """
}

process binScores {
    input:
    tuple val(sample), path(sorted) from bin_scores

    output:
    tuple val(sample), path("${sample}_qualbin.bam") into to_cram

    """
    htsbox qualbin -m 3 -b ${sorted} > ${sample}_qualbin.bam
    """
}


process bamToCram {
    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple val(sample), path(qualbin) from to_cram
    path "ref.fa" from params.fa_ref

    output:
    path "${sample}.cram"

    script:
    """
    samtools view ${qualbin} -C -T ref.fa > ${sample}.cram
    """
}