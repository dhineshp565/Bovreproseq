#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//parse input csv file with SampleName,path to SampleName fastq

params.input='./SampleList.csv'
params.outdir='Results'
params.reference='reference.fasta'
params.primerbed='primer.bed'
params.trim_barcodes=null
params.centri='centrifuge_index'
params.db='kraken_db'

//merge fastq files for each SampleName and create a merged file for each SampleNames
process merge_fastq {
	publishDir "${params.outdir}/merged"
	label "low"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path("${SampleName}.{fastq,fastq.gz}")
	
	shell:
	"""
	count=\$(ls -1 $SamplePath/*.gz 2>/dev/null | wc -l)
	
	
		if [[ "\${count}" != "0" ]];
		then
			cat $SamplePath/*.fastq.gz > ${SampleName}.fastq.gz
		
		else
			cat $SamplePath/*.fastq > ${SampleName}.fastq
		fi
	"""
}

//trim barcodes and adapter using porechop

process porechop {
	label "medium"
	publishDir "${params.outdir}/trimmed",mode:"copy",overwrite: false
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path ("${SampleName}_trimmed.fastq")
	script:
	"""
	porechop -i $SamplePath -o ${SampleName}_trimmed.fastq
	"""
}
// sequence alignment using minimap2
process minimap2 {
        publishDir "${params.outdir}/minimap2/",mode:"copy",overwrite: false
		label "low"
        input:
        path (reference)
        tuple val(SampleName),path(SamplePath)
        output:
        tuple val(SampleName),path ("${SampleName}.sam")
        script:
        """
        minimap2 -ax map-ont $reference $SamplePath > ${SampleName}.sam
        """
}

//convert minimap2 output sam to sorted bam, create reference index and bed file from the given reference input
process samtools {
	publishDir "${params.outdir}/samtools",mode:"copy",overwrite: false
	label "medium"
	input:
	path reference
	tuple val(SampleName),path(SamplePath)
	output:
	val(SampleName),emit:SampleName
	path ("${reference}.fai")
	path ("${SampleName}.bam"),emit:bam
	path ("${reference}.bed"),emit:bed
	path ("${SampleName}_stats.txt"),emit:stats
	script:
	"""
#generate a bam file with primary alignments
	samtools view -b -F 256 ${SamplePath}|samtools sort > ${SampleName}.bam
	samtools stats "${SampleName}.bam" > ${SampleName}_stats.txt
	samtools faidx ${reference}
	awk 'BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}' "${reference}.fai" > ${reference}.bed  
	"""
}
//split bam files and create consensus
process splitbam {
	publishDir "${params.outdir}/splitbam",mode:"copy",overwrite: false
	label "medium"
	input:
	val(SampleName)
	path(SamplePath)
	path (primerbed)
	output:
	val(SampleName),emit:SampleName
	path("${SampleName}_mappedreads.txt"),emit:mapped
	path("${SampleName}_idxstats.txt"),emit:idxstats
	path("${SampleName}_*_*_idxstats.txt"),emit:full_idxstats
	path("${SampleName}_*.bam"),emit:bam
	path("${SampleName}_consensus.fasta"),emit:consensus
	script:
	"""
#trims primers from both ends of the ampliocn using primer bed file
	samtools ampliconclip --both-ends -b ${primerbed} ${SamplePath} > ${SampleName}_trimmed.bam
# sorts the bam file for spitting	
	samtools sort "${SampleName}_trimmed.bam" > ${SampleName}_tr_sorted.bam
#index sorted bam file and generate read counts for each amplicon
	samtools index "${SampleName}_tr_sorted.bam" > ${SampleName}_sorted.bai
	samtools idxstats "${SampleName}_tr_sorted.bam" > ${SampleName}_idxstats.txt
	awk '{if (\$3!=0) print \$1,\$2,\$3}' "${SampleName}_idxstats.txt" > ${SampleName}_mappedreads.txt
#using the list of mapped amplicons from text file, bam file is split based on amplicons and consensus is generated for each mapped amplicons
	while read lines
	do 
		amp=\$(echo \$lines|cut -f1 -d' ')
		len=\$(echo \$lines|cut -f2 -d' ')
		samtools view -b "${SampleName}_tr_sorted.bam" "\${amp}" > ${SampleName}_\${amp}.bam
	    samtools view -h "${SampleName}_\${amp}.bam"|awk -v l=\${len} '/^@/|| length(\$10)>=l-50 && length(\$10)<=l+50'|samtools sort > ${SampleName}_\${len}_\${amp}.bam
		samtools index "${SampleName}_\${len}_\${amp}.bam" > ${SampleName}_\${len}_\${amp}.bai
		samtools idxstats "${SampleName}_\${len}_\${amp}.bam" > ${SampleName}_\${len}_\${amp}_idxstats.txt
		samtools consensus -f fasta "${SampleName}_\${len}_\${amp}.bam" > ${SampleName}_\${amp}.fasta
		sed -i "s/>.*/>${SampleName}_\${amp}_consensus/" ${SampleName}_\${amp}.fasta
	done < "${SampleName}_mappedreads.txt"
	cat ${SampleName}_*.fasta > ${SampleName}_consensus.fasta
# insert headers to mappedreads.txt
	sed -i '1i SampleID ${SampleName}' "${SampleName}_mappedreads.txt"
	"""
}
//split fasta for mapped reads only using seqtk

process mapped_ref {
	publishDir "${params.outdir}/mapped_ref",mode:"copy",overwrite: false
	input:
	val (SampleName)
	path(txtfile)
	path (reference)
	output:
	val(SampleName),emit:SampleName
	path("${SampleName}_mapped_ref.fasta"),emit:fasta
	path("${txtfile}_nameonly.txt"),emit:nameonly
	script:
	"""
	tail -n +2  ${txtfile}|cut -f 1 -d' ' > ${txtfile}_nameonly.txt
	seqtk subseq ${reference} "${txtfile}_nameonly.txt" > ${SampleName}_mapped_ref.fasta
	"""
}
//Generate fasta index (fai) and bed from fasta file created in mapped ref
process mapped_ref_bed {
	publishDir "${params.outdir}/mapped_ref_bed",mode:"copy",overwrite: false
	input:
	val (SampleName)
	path (mapped_fasta)
	output:
	path ("${SampleName}_mapped.bed")
	script:
	"""
	samtools faidx ${mapped_fasta}
        awk 'BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}' "${mapped_fasta}.fai" > ${SampleName}_mapped.bed
	"""
}

// generating bedgraph files from alignments
process bedtools {
	publishDir "${params.outdir}/bedtools/",mode:"copy",overwrite: false
	input:
	val(SampleName)
	path(SamplePath)
	path(txtfile)
	output:
	val(SampleName)
	path ("${SampleName}*.bedgraph")
	script:
	"""
	while read lines;do bedtools genomecov -ibam ${SampleName}_\$lines.bam -bga > ${SampleName}_\${lines}.bedgraph;done < ${txtfile}

	"""
}

//igv reports from bedgraph
process igvreports {
	publishDir "${params.outdir}/igvreports/",mode:"copy",overwrite: false
	input:
	path(reference)
        path(bed)
	val(SampleName)
	path(bedgraph)
	output:
	path{"*.html"}
	script:
	"""
	create_report $bed $reference\
	--tracks *.bedgraph\
	--output ${SampleName}.html 
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
	tuple val(SampleName),path (trimmed_fastq)
	path(db_path)
	path (consensus)
	
	output:
	path ("${SampleName}_kraken.csv")
	path ("${SampleName}_cons_kraken.csv")
	path ("${SampleName}_kraken_report.csv"),emit:(kraken2_raw)
	path ("${SampleName}_cons_kraken_report.csv"),emit:(kraken2_consensus)

	script:
	"""
	kraken2 --db $db_path --output ${SampleName}_kraken.csv --report ${SampleName}_kraken_report.csv --threads 1 $trimmed_fastq
	kraken2 --db $db_path --output ${SampleName}_cons_kraken.csv --report ${SampleName}_cons_kraken_report.csv --threads 1 $consensus
	"""
}
//centrifuge for taxonomy classification
process centrifuge {
	publishDir "${params.outdir}/centrifuge/",mode:"copy",overwrite: false
	label "medium"
	errorStrategy 'ignore'
	input:
	tuple val(SampleName),path (trimmed_fastq)
	path(db_path)
	path (consensus)
	
	output:
	path ("${SampleName}_cent_report.csv")
	path ("${SampleName}_centrifuge_kstyle.csv"),emit:(centrifuge_raw)
	path ("${SampleName}_consensus_kstyle.csv"),emit:(centrifuge_consensus)
	path("${SampleName}_centrifuge.csv")

	script:
	"""
	index=\$(find -L ${db_path} -name "*.1.cf" -not -name "._*"  | sed 's/.1.cf//')
	centrifuge -x \${index} -U $trimmed_fastq -q -S ${SampleName}_centrifuge.csv --report ${SampleName}_cent_report.csv
	centrifuge-kreport -x \${index} "${SampleName}_centrifuge.csv" > ${SampleName}_centrifuge_kstyle.csv
	centrifuge -x \${index} -U $consensus -f -S ${SampleName}_consensus_centrifuge.csv 
	centrifuge-kreport -x \${index} "${SampleName}_consensus_centrifuge.csv" > ${SampleName}_consensus_kstyle.csv
	"""
}
//krona plots
process krona {
	publishDir "${params.outdir}/krona/",mode:"copy",overwrite: false
	label "low"
	

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
    .map { row-> tuple(row.SampleName,row.SamplePath) }
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
	splitbam(samtools.out.SampleName,samtools.out.bam,primerbed)
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


/*	mapped_ref(splitbam.out.SampleName,splitbam.out.mapped,reference)
	mapped_ref_bed(mapped_ref.out.SampleName,mapped_ref.out.fasta)
	bedtools(splitbam.out.SampleName,splitbam.out.bam,mapped_ref.out.nameonly)
	igvreports(mapped_ref.out.fasta,mapped_ref_bed.out,bedtools.out)*/
}
