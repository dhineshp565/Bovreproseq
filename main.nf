#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//parse input csv file with sample,path to sample fastq

params.input='./samplelist.csv'
params.outdir='Results'
params.reference='reference.fasta'
params.primerbed='primer.bed'

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
	path ("${sample}.bam"),emit:bam
	path ("${reference}.bed"),emit:bed
	script:
	"""
	samtools view -b ${sample_path}|samtools sort > ${sample}.bam
	samtools faidx ${reference}
	awk 'BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}' "${reference}.fai" > ${reference}.bed  
	"""
}
//splitbam
process splitbam {
	publishDir "${params.outdir}/splitbam"
	input:
	val(sample)
	path(sample_path)
	path (primerbed)
	output:
	val(sample),emit:sample
	path("${sample}_mappedreads.txt"),emit:mapped
	path("${sample}_*.bam"),emit:bam
	path("${sample}_consensus.fasta"),emit:consensus
	shell:
	"""
	samtools ampliconclip --both-ends -b ${primerbed} ${sample_path} > ${sample}_trimmed.bam
	samtools sort "${sample}_trimmed.bam" > ${sample}_tr_sorted.bam
	samtools index "${sample}_tr_sorted.bam" > ${sample}_sorted.bai
	samtools idxstats "${sample}_tr_sorted.bam" > ${sample}_idxstats.txt
	awk '{if (\$3!=0) print \$1,\$3}' "${sample}_idxstats.txt" > ${sample}_mappedreads.txt
	while read lines
	do 
		amp=\$(echo \$lines|cut -f1 -d' ')
		samtools view -b "${sample}_tr_sorted.bam" "\${amp}" > ${sample}_\${amp}.bam
		samtools consensus -f fasta ${sample}_\${amp}.bam > ${sample}_\${amp}.fasta
		sed -i "s/>.*/>${sample}_\${amp}_consensus/" ${sample}_\${amp}.fasta
	done < "${sample}_mappedreads.txt"
	cat ${sample}_*.fasta > ${sample}_consensus.fasta
	"""
}
//split fasta for mapped reads only using seqtk

process mapped_ref {
	publishDir "${params.outdir}/mapped_ref"
	input:
	val (sample)
	path(txtfile)
	path (reference)
	output:
	val(sample),emit:sample
	path("${sample}_mapped_ref.fasta"),emit:fasta
	path("${txtfile}_nameonly.txt"),emit:nameonly
	script:
	"""
	cat ${txtfile}|cut -f 1 -d' ' > ${txtfile}_nameonly.txt
	seqtk subseq ${reference} "${txtfile}_nameonly.txt" > ${sample}_mapped_ref.fasta
	"""
}
//Generate fai and bed from mapped reference
process mapped_ref_bed {
	publishDir "${params.outdir}/mapped_ref_bed"	
	input:
	val (sample)
	path (mapped_fasta)
	output:
	path ("${sample}_mapped.bed")
	script:
	"""
	samtools faidx ${mapped_fasta}
        awk 'BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}' "${mapped_fasta}.fai" > ${sample}_mapped.bed
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
	path(sample_path)
	path(txtfile)
	output:
	val(sample)
	path ("${sample}*.bedgraph")
	script:
	"""
	while read lines;do bedtools genomecov -ibam ${sample}_\$lines.bam -bga > ${sample}_\${lines}.bedgraph;done < ${txtfile}

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
	.fromPath(params.input)
	.splitCsv(header:true)
        .map { row-> tuple(row.sample,row.sample_path) }
	reference=file(params.reference)
	primerbed=file(params.primerbed)	
        merge_fastq(data)
	|view()
	porechop(merge_fastq.out)
	|view ()
	minimap2(reference,porechop.out)
	samtools(reference,minimap2.out)
	splitbam(samtools.out.sample,samtools.out.bam,primerbed)
	mapped_ref(splitbam.out.sample,splitbam.out.mapped,reference)
	mapped_ref_bed(mapped_ref.out.sample,mapped_ref.out.fasta)
	bedtools(splitbam.out.sample,splitbam.out.bam,mapped_ref.out.nameonly)
	igvreports(mapped_ref.out.fasta,mapped_ref_bed.out,bedtools.out)
}
