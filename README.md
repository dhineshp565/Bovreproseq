# Bovreproseq
Targeted multiplex amplicon Sequencing for Bovine Reproductive Pathogens.
Requires input directory with fastq folders, reference sequence (fasta)
Outputs consensus sequences,kraken,krona and multiqc report if any reads are mapped to reference
conda or docker needs to be installed

Usage:
```
nextflow run main.nf --input path_to_input --out_dir Results --kraken_db path_to_kraken_database
```

```
Parameters:

--input      Path to input directory
--out_dir     Output directory
--kraken_db  path to kraken database 
optional
--trim_barcodes barcode and adapter trimming using porechop

```
