# Bovreproseq
Targeted multiplex amplicon Sequencing for Bovine Reproductive Pathogens.
Pipeline should work for any multiplex amplicon sequencing. Outputs consensus sequences,kraken,centrifuge, krona and multiqc report if any reads are mapped to reference. 
Requires input directory with fastq folders, reference sequence (fasta),primerbedfile (chrom name-should match header in the reference sequence, primer start and end coordinates (5'-3')), the path to Kraken2 database
conda or docker needs to be installed.

Usage:
```
nextflow run main.nf --input path_to_input --outdir Results --reference test_multi.fasta --primerbed primer.bed --kraken_db minikraken2_v1_8GB -profile docker
```

```
Parameters:

--input      Path to input directory
--outdir     Output directory
--reference  fasta or multi-fasta 
--primerbed  Bovreproseq_primer.bed
--profile    conda or docker
Taxonomic classification - one should be chosen
--kraken_db  path to kraken database 
--centri_db  path to centrifuge database
optional
--trim_barcodes barcode and adapter trimming using porechop

```
