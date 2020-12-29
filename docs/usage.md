# nf-core/fastq: Usage

## Parameters 

#### `--hisat2_index`
Path to hisat2_index.

#### `--fa_ref`
Path to reference fa file. 

#### `--samples_path`
TSV file with 3 columns: sample name, fastq first and fastq second reads. (Example: `data/geuvadis_eur_reads.tsv`)


## How to run

```
nextflow run main.nf -resume -profile tartu_hpc \
  --samples_path data/geuvadis_eur_reads.tsv
```

## Convert back to fastq

#### `--cramsPath` 
One or multiple cram files. 

```
nextflow run cram_to_fastq.nf -profile tartu_hpc -resume \
 --cramsPath "/gpfs/hpc/home/peikova/fastq_compress/results/*.cram"
```