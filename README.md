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
## Dependencies
* nextflow
* docker
* wsl2
## Software Used
* nanoq (Steinig and Coin (2022). Nanoq: ultra-fast quality control for nanopore reads. Journal of Open Source Software, 7(69), 2991, https://doi.org/10.21105/joss.02991)
* porechop (https://github.com/rrwick/Porechop)
* minimap2 (Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191)
* kraken2 (https://ccb.jhu.edu/software/kraken2/)
* krona:2.7.1(Ondov, B.D., Bergman, N.H. & Phillippy, A.M. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics 12, 385 (2011). https://doi.org/10.1186/1471-2105-12-385)
* abricate (https://github.com/tseemann/abricate)
* mlst (https://github.com/tseemann/mlst,This publication made use of the PubMLST website (https://pubmlst.org/) developed by Keith Jolley (Jolley & Maiden 2010, BMC Bioinformatics, 11:595) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust)
* rmarkdown (https://rmarkdown.rstudio.com/)
