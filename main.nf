#!/usr/bin/env nextflow

def helpMessage(){
    log.info "help message"

}

if (params.help) {
    helpMessage()
    exit 0
}

params.fastq_pair_1 = "/home/peikova/thesis/fastq/data/ERR204864_1.fastq.gz"
params.fastq_pair_2 = "/home/peikova/thesis/fastq/data/ERR204864_2.fastq.gz"
params.hisat2_index = "/home/peikova/thesis/fastq/data/index/Homo_sapiens.GRCh38.dna.primary_assembly"
params.fa_ref = "/home/peikova/thesis/fastq/data/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

hs2_indices = Channel.fromPath("${params.hisat2_index}*")
    .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }

process alignWithHisat2 {
    input: 
    path "read_1.fastq.gz" from params.fastq_pair_1
    path "read_2.fastq.gz" from params.fastq_pair_2
    file hs2_indices from hs2_indices.collect()

    output:
    path "aligned.bam" into aligned_bam
    
    script:
    index_base = hs2_indices[0].toString() - ~/.\d.ht2/

    """
    hisat2 -x $index_base -1 read_1.fastq.gz -2 read_2.fastq.gz -u 100 \\
        | samtools view -bS -F 4 -F 256 - > aligned.bam
    """

}

process sortBam {
    input:
    path "aligned.bam" from aligned_bam

    output:
    path "sorted.bam" into bin_scores

    """
    samtools sort aligned.bam -O bam -o sorted.bam
    """
}

process binScores {
    input:
    path "sorted.bam" from bin_scores

    output:
    path "qualbin.bam" into to_cram

    """
    htsbox qualbin -m 3 -b sorted.bam > qualbin.bam
    """
}


process bamToCram {
    publishDir "$launchDir/result", mode: "copy"

    input:
    path "qualbin.bam" from to_cram
    path "ref.fa" from params.fa_ref

    output:
    path "file.cram"

    script:
    """
    samtools view qualbin.bam -C -T ref.fa > file.cram
    """
}