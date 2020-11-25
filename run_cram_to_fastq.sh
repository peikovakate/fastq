./nextflow run fastq/cram_to_fastq.nf\
 -profile tartu_hpc\
 --cramsPath "/gpfs/hpc/home/peikova/fastq_compress/results/*.cram"\
 -resume