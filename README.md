# Bovreproseq
Targeted multiplex amplicon Sequencing for Bovine Reproductive Pathogens.
Pipeline should work for any multiplex amplicon sequencing. Outputs consensus sequences,kraken,centrifuge, krona and multiqc report if there are any reads mapped to reference. 
Requires input csv file(sample name,samplepath),reference sequence (fasta),primerbedfile (chrom name-should match header in reference sequence,primer start and end coordinates (5'-3')),path to kraken and centrifuge databases.

conda or docker needs to be installed.

Usage:
```
nextflow run main.nf --input samplelist.csv --outdir Results_reruns --reference Bovreproseq_ref_multi.fasta --primerbed Bovreproseq_primer.bed --db /data/referenceDB/kraken/minikraken2_v1_8GB/ --centri /data/referenceDB/centrifuge/ -profile docker 

Options:

--input     csv file with two columns with headers(sample,sample_path).See samplelist.csv
--outdir    Output directory
--reference fasta or multi-fasta 
--primerbed Bovreproseq_primer.bed 
--db        path to kraken database 
--centri    path to centrifuge database
--profile   conda or docker
```
