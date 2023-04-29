#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//parse input csv file with sample,path to sample fastq

params.input='./samplelist.csv'
params.outdir='Results'
params.reference='reference.fasta'

process checksamplelist {
	
	input:
	tuple val(sample),val(sample_path)
	output:
	stdout
	
	"""
	echo 'this $sample is in  $sample_path'

	"""
}
//merge fastq files for each sample
process merge_fastq {
	publishDir "${params.outdir}/merged"
	input:
	tuple val(sample),path(sample_path)
	output:
	tuple val(sample),path("${sample}.fastq")
	script:
	"""
	cat $sample_path/*.fastq > ${sample}.fastq
	"""
}

//trim barcodes and adapter using porechop

process porechop {
	
	publishDir "${params.outdir}/trimmed"
	input:
	tuple val(sample),path(sample_path)
	output:
	tuple val(sample),path ("${sample}_trimmed.fastq")
	script:
	"""
	porechop -i $sample_path -o ${sample}_trimmed.fastq
	"""
}
//getting index from reference
process samtools {
	publishDir "${params.outdir}/samtools"
	input:
	path reference
	tuple val(sample),path(sample_path)
	output:
	val(sample),emit:sample
	path ("${reference}.fai")
	path ("${sample}_sorted.bam"),emit:bam
	path ("${reference}.bed"),emit:bed
	script:
	"""
	samtools view -b ${sample_path}|samtools sort > ${sample}_sorted.bam
	samtools faidx $reference
	awk 'BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}' "${reference}.fai" > ${reference}.bed  
	"""
}
// sequence alignment
process minimap2 {
	publishDir "${params.outdir}/minimap2/"
	input:
	path (reference)	
	tuple val(sample),path(sample_path)
	output:
        tuple val(sample),path ("${sample}.sam")
	script:
	"""
	minimap2 -ax map-ont $reference $sample_path > ${sample}.sam
	"""
}
// generating bedgraph files from alignments
process bedtools {
	publishDir "${params.outdir}/bedtools/"
	input:
	val(sample)
	path (sample_path)
	
	output:
	val(sample)
	path ("${sample_path}.bedgraph"),emit:bedgraph
	script:
	"""
	bedtools genomecov -ibam ${sample_path} -bga > ${sample_path}.bedgraph
	"""
}

//igv reports from bedgraph
process igvreports {
	publishDir "${params.outdir}/igvreports/"
	input:
	path(reference)
        path(bed)
	val(sample)
	path(bedgraph)
	output:
	path{"*.html"}
	script:
	"""
	create_report $bed $reference\
	--tracks *.bedgraph\
	--output ${sample}.html 
	"""
}

workflow {
	data=Channel
	.fromPath(params.csvfile)
	.splitCsv(header:true)
        .map { row-> tuple(row.sample,row.sample_path) }
	reference=file(params.reference)	
        merge_fastq(data)
	|view()
	porechop(merge_fastq.out)
	|view ()
	minimap2(reference,porechop.out)
	samtools(reference,minimap2.out)
	bedtools(samtools.out.sample,samtools.out.bam)
	igvreports(reference,samtools.out.bed,bedtools.out)
}
