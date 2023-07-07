#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//parse input csv file with sample,path to sample fastq

params.input='./samplelist.csv'
params.outdir='Results'
params.reference='reference.fasta'
params.primerbed='primer.bed'
params.trim_barcodes=null
params.centri='centrifuge_index'

//merge fastq files for each sample and create a merged file for each samples
process merge_fastq {
	publishDir "${params.outdir}/merged"
	label "low"
	input:
	tuple val(sample),path(sample_path)
	output:
	tuple val(sample),path("${sample}.{fastq,fastq.gz}")
	shell:
	"""
	count=\$(ls -1 $sample_path/*.gz 2>/dev/null | wc -l)
	
		if [[ "\${count}" != "0" ]];
		then
			cat $sample_path/*.fastq.gz > ${sample}.fastq.gz
		
		else
			cat $sample_path/*.fastq > ${sample}.fastq
		fi
	"""
}

//trim barcodes and adapter using porechop

process porechop {
	label "medium"
	publishDir "${params.outdir}/trimmed",mode:"copy",overwrite: false
	input:
	tuple val(sample),path(sample_path)
	output:
	tuple val(sample),path ("${sample}_trimmed.fastq")
	script:
	"""
	porechop -i $sample_path -o ${sample}_trimmed.fastq
	"""
}
// sequence alignment using minimap2
process minimap2 {
        publishDir "${params.outdir}/minimap2/",mode:"copy",overwrite: false
		label "low"
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

//convert minimap2 output sam to sorted bam, create reference index and bed file from the given reference input
process samtools {
	publishDir "${params.outdir}/samtools",mode:"copy",overwrite: false
	label "medium"
	input:
	path reference
	tuple val(sample),path(sample_path)
	output:
	val(sample),emit:sample
	path ("${reference}.fai")
	path ("${sample}.bam"),emit:bam
	path ("${reference}.bed"),emit:bed
	path ("${sample}_stats.txt"),emit:stats
	script:
	"""
#generate a bam file with primary alignments
	samtools view -b -F 256 ${sample_path}|samtools sort > ${sample}.bam
	samtools stats "${sample}.bam" > ${sample}_stats.txt
	samtools faidx ${reference}
	awk 'BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}' "${reference}.fai" > ${reference}.bed  
	"""
}
//split bam files and create consensus
process splitbam {
	publishDir "${params.outdir}/splitbam",mode:"copy",overwrite: false
	label "medium"
	input:
	val(sample)
	path(sample_path)
	path (primerbed)
	output:
	val(sample),emit:sample
	path("${sample}_mappedreads.txt"),emit:mapped
	path("${sample}_idxstats.txt"),emit:idxstats
	path("${sample}_*_*_idxstats.txt"),emit:full_idxstats
	path("${sample}_*.bam"),emit:bam
	path("${sample}_consensus.fasta"),emit:consensus
	script:
	"""
#trims primers from both ends of the ampliocn using primer bed file
	samtools ampliconclip --both-ends -b ${primerbed} ${sample_path} > ${sample}_trimmed.bam
# sorts the bam file for spitting	
	samtools sort "${sample}_trimmed.bam" > ${sample}_tr_sorted.bam
#index sorted bam file and generate read counts for each amplicon
	samtools index "${sample}_tr_sorted.bam" > ${sample}_sorted.bai
	samtools idxstats "${sample}_tr_sorted.bam" > ${sample}_idxstats.txt
	awk '{if (\$3!=0) print \$1,\$2,\$3}' "${sample}_idxstats.txt" > ${sample}_mappedreads.txt
#using the list of mapped amplicons from text file, bam file is split based on amplicons and consensus is generated for each mapped amplicons
	while read lines
	do 
		amp=\$(echo \$lines|cut -f1 -d' ')
		len=\$(echo \$lines|cut -f2 -d' ')
		samtools view -b "${sample}_tr_sorted.bam" "\${amp}" > ${sample}_\${amp}.bam
	    samtools view -h "${sample}_\${amp}.bam"|awk -v l=\${len} '/^@/|| length(\$10)>=l-50 && length(\$10)<=l+50'|samtools sort > ${sample}_\${len}_\${amp}.bam
		samtools index "${sample}_\${len}_\${amp}.bam" > ${sample}_\${len}_\${amp}.bai
		samtools idxstats "${sample}_\${len}_\${amp}.bam" > ${sample}_\${len}_\${amp}_idxstats.txt
		samtools consensus -f fasta "${sample}_\${len}_\${amp}.bam" > ${sample}_\${amp}.fasta
		sed -i "s/>.*/>${sample}_\${amp}_consensus/" ${sample}_\${amp}.fasta
	done < "${sample}_mappedreads.txt"
	cat ${sample}_*.fasta > ${sample}_consensus.fasta
# insert headers to mappedreads.txt
	sed -i '1i SampleID ${sample}' "${sample}_mappedreads.txt"
	"""
}
//split fasta for mapped reads only using seqtk

process mapped_ref {
	publishDir "${params.outdir}/mapped_ref",mode:"copy",overwrite: false
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
	tail -n +2  ${txtfile}|cut -f 1 -d' ' > ${txtfile}_nameonly.txt
	seqtk subseq ${reference} "${txtfile}_nameonly.txt" > ${sample}_mapped_ref.fasta
	"""
}
//Generate fasta index (fai) and bed from fasta file created in mapped ref
process mapped_ref_bed {
	publishDir "${params.outdir}/mapped_ref_bed",mode:"copy",overwrite: false
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

// generating bedgraph files from alignments
process bedtools {
	publishDir "${params.outdir}/bedtools/",mode:"copy",overwrite: false
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
	publishDir "${params.outdir}/igvreports/",mode:"copy",overwrite: false
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

//multiqc

process multiqc {
	publishDir "${params.outdir}/multiqc/",mode:"copy",overwrite: false
	label "low"
	input:
	path '*'
	output:
	file ("multiqc_report.html")
	file ("multiqc_data")
	script:
	"""
	multiqc .
	"""
}

//kraken2 for classification
process kraken2 {
	publishDir "${params.outdir}/kraken2/",mode:"copy",overwrite: false
	label "high"
	input:
	tuple val(sample),path (trimmed_fastq)
	path(db_path)
	path (consensus)
	
	output:
	path ("${sample}_kraken.csv")
	path ("${sample}_cons_kraken.csv")
	path ("${sample}_kraken_report.csv"),emit:(kraken2_raw)
	path ("${sample}_cons_kraken_report.csv"),emit:(kraken2_consensus)

	script:
	"""
	kraken2 --db $db_path --output ${sample}_kraken.csv --report ${sample}_kraken_report.csv --threads 1 $trimmed_fastq
	kraken2 --db $db_path --output ${sample}_cons_kraken.csv --report ${sample}_cons_kraken_report.csv --threads 1 $consensus
	"""
}
//centrifuge for taxonomy classification
process centrifuge {
	publishDir "${params.outdir}/centrifuge/",mode:"copy",overwrite: false
	label "medium"
	input:
	tuple val(sample),path (trimmed_fastq)
	path(db_path)
	path (consensus)
	
	output:
	path ("${sample}_cent_report.csv")
	path ("${sample}_centrifuge_kstyle.csv"),emit:(centrifuge_raw)
	path ("${sample}_consensus_kstyle.csv"),emit:(centrifuge_consensus)
	path("${sample}_centrifuge.csv")

	script:
	"""
	index=\$(find -L ${db_path} -name "*.1.cf" -not -name "._*"  | sed 's/.1.cf//')
	centrifuge -x \${index} -U $trimmed_fastq -q -S ${sample}_centrifuge.csv --report ${sample}_cent_report.csv
	centrifuge-kreport -x \${index} "${sample}_centrifuge.csv" > ${sample}_centrifuge_kstyle.csv
	centrifuge -x \${index} -U $consensus -f -S ${sample}_consensus_centrifuge.csv 
	centrifuge-kreport -x \${index} "${sample}_consensus_centrifuge.csv" > ${sample}_consensus_kstyle.csv
	"""
}
//krona plots
process krona {
	publishDir "${params.outdir}/krona/",mode:"copy",overwrite: false
	label "low"
	errorStrategy 'ignore'

	input:
	path(kraken_raw)
	path(centrifuge_raw)
	path(kraken_consensus)
	path(centrifuge_consensus)
	output:
	path ("kraken_raw.html")
	path ("centrifuge_raw.html")
	path("kraken_consensus.html")
	path("centrifuge_consensus.html")
	script:
	"""
	ktImportTaxonomy -t 5 -m 3 -o kraken_raw.html $kraken_raw
	ktImportTaxonomy -t 5 -m 3 -o centrifuge_raw.html $centrifuge_raw
	ktImportTaxonomy -t 5 -m 3 -o kraken_consensus.html $kraken_consensus
	ktImportTaxonomy -t 5 -m 3 -o centrifuge_consensus.html $centrifuge_consensus
	"""
}


workflow {
	data=Channel
	.fromPath(params.input)
	.splitCsv(header:true)
        .map { row-> tuple(row.sample,row.sample_path) }
	reference=file(params.reference)
	primerbed=file(params.primerbed)
	db=file(params.db)	
    merge_fastq(data)
//trim barcodes and adapter sequences
	if (params.trim_barcodes){
		porechop(merge_fastq.out)
		minimap2(reference,porechop.out)
		 
	} else {
                minimap2(reference,merge_fastq.out)
		
        }
	samtools(reference,minimap2.out)
	splitbam(samtools.out.sample,samtools.out.bam,primerbed)
	stats=samtools.out.stats
	idxstats=splitbam.out.idxstats
		if (params.trim_barcodes){
              kraken2(porechop.out,db,splitbam.out.consensus)          
			  centrifuge(porechop.out,params.centri,splitbam.out.consensus)
	 } else {
		kraken2(merge_fastq.out,db,splitbam.out.consensus)
		centrifuge(merge_fastq.out,params.centri,splitbam.out.consensus)
	}
	
	multiqc(stats.mix(idxstats).collect())
	kraken_raw=kraken2.out.kraken2_raw
	kraken_cons=kraken2.out.kraken2_consensus
	centri_raw=centrifuge.out.centrifuge_raw
	centri_cons=centrifuge.out.centrifuge_consensus
	krona(kraken_raw.collect(),centri_raw.collect(),kraken_cons.collect(),centri_cons.collect())


/*	mapped_ref(splitbam.out.sample,splitbam.out.mapped,reference)
	mapped_ref_bed(mapped_ref.out.sample,mapped_ref.out.fasta)
	bedtools(splitbam.out.sample,splitbam.out.bam,mapped_ref.out.nameonly)
	igvreports(mapped_ref.out.fasta,mapped_ref_bed.out,bedtools.out)*/
}
