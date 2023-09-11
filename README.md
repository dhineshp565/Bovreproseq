# Bovreproseq
Targeted multiplex amplicon Sequencing for Bovine Reproductive Pathogens.
Pipeline should work for any multiplex amplicon sequencing. Outputs consensus sequences,kraken,centrifuge, krona and multiqc report if any reads are mapped to reference. 
Requires input csv file(SampleName, Samplepath), reference sequence (fasta),primerbedfile (chrom name-should match header in the reference sequence, primer start and end coordinates (5'-3')), the path to Kraken2 and centrifuge databases.

conda or docker needs to be installed.

Usage:
```
nextflow run main.nf --input samplelist.csv --outdir Results --reference test_multi.fasta --primerbed primer.bed --kraken_db minikraken2_v1_8GB --centri_db centrifuge -profile docker
```

```
Options:

--input      csv file with two columns with headers(SampleName,SamplePath).See samplelist.csv
--outdir     Output directory
--reference  fasta or multi-fasta 
--primerbed  Bovreproseq_primer.bed 
--kraken_db  path to kraken database 
--centri_db  path to centrifuge database
--profile    conda or docker
optional
--trim_barcodes barcode and adapter trimming using porechop
```
