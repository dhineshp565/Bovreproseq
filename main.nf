#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
	count=\$(ls -1 ${SamplePath}/*.gz 2>/dev/null | wc -l)
	
	
		if [[ "\${count}" != "0" ]];
		then
			cat ${SamplePath}/*.fastq.gz > ${SampleName}.fastq.gz
		
		else
			cat ${SamplePath}/*.fastq > ${SampleName}.fastq
		fi
	"""
}

//trim barcodes and adapter using porechop

process porechop {
	label "medium"
	publishDir "${params.outdir}/trimmed"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path ("${SampleName}_trimmed.fastq")
	script:
	"""
	porechop -i ${SamplePath} -o ${SampleName}_trimmed.fastq
	"""
}
// sequence alignment using minimap2
process minimap2 {
        publishDir "${params.outdir}/minimap2/",mode:"copy"
		label "low"
        input:
        path (reference)
        tuple val(SampleName),path(SamplePath)
        output:
        tuple val(SampleName),path ("${SampleName}.sam")
        script:
        """
        minimap2 -ax map-ont ${reference} ${SamplePath} > ${SampleName}.sam
        """
}

//convert minimap2 output sam to sorted bam, create reference index and bed file from the given reference input
process samtools {
	publishDir "${params.outdir}/samtools",mode:"copy"
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
	#generate a  sorted bam file with primary alignments
	samtools view -b -F 256 ${SamplePath}|samtools sort > ${SampleName}.bam
	samtools stats "${SampleName}.bam" > ${SampleName}_stats.txt  
	"""
}
//split bam files and create consensus
process splitbam {
	publishDir "${params.outdir}/splitbam",mode:"copy"
	label "medium"
	input:
	tuple val(SampleName),path(SamplePath)
	path (primerbed)
	output:
	val(SampleName),emit:SampleName
	path("${SampleName}_stats.txt"),emit:stats
	path("${SampleName}_mappedreads.txt"),emit:mapped
	path("${SampleName}_idxstats.txt"),emit:idxstats
	path("${SampleName}_*_*_idxstats.txt"),emit:full_idxstats
	path("${SampleName}_*.bam"),emit:bam
	tuple val(SampleName),path("${SampleName}_consensus.fasta"),emit:consensus
	script:
	"""
#generate a  sorted bam file with primary alignments
	samtools view -b -F 256 ${SamplePath}|samtools sort > ${SampleName}.bam	
	samtools stats "${SampleName}.bam" > ${SampleName}_stats.txt
#trims primers from both ends of the ampliocn using primer bed file
	samtools ampliconclip --both-ends -b ${primerbed} "${SampleName}.bam"	> ${SampleName}_trimmed.bam
# sorts the bam file for spitting	
	samtools sort "${SampleName}_trimmed.bam" > ${SampleName}_tr_sorted.bam
#index sorted bam file and generate read counts for each amplicon
	samtools index "${SampleName}_tr_sorted.bam" > ${SampleName}_sorted.bai
	samtools idxstats "${SampleName}_tr_sorted.bam" > ${SampleName}_idxstats.txt
	awk '{if (\$3!=0) print \$1,\$2,\$3}' "${SampleName}_idxstats.txt" > ${SampleName}_mappedreads.txt
#using the list of mapped amplicons from text file, bam file is split based on amplicons and consensus is gen
	while read lines
	do 
		amp=\$(echo \$lines|cut -f1 -d' ')
		len=\$(echo \$lines|cut -f2 -d' ')
		# split bam 
		samtools view -b "${SampleName}_tr_sorted.bam" "\${amp}" > ${SampleName}_\${amp}.bam
		# Only reads length with + or - 50 bases is used for consenus
	    samtools view -h "${SampleName}_\${amp}.bam"|awk -v l=\${len} '/^@/|| length(\$10)>=l-50 && length(\$10)<=l+50'|samtools sort > ${SampleName}_\${len}_\${amp}.bam
		# generate stats for near full length reads
		samtools index "${SampleName}_\${len}_\${amp}.bam" > ${SampleName}_\${len}_\${amp}.bai
		samtools idxstats "${SampleName}_\${len}_\${amp}.bam" > ${SampleName}_\${len}_\${amp}_idxstats.txt
		# generate consensus for full length reads
		samtools consensus -f fasta "${SampleName}_\${len}_\${amp}.bam" > ${SampleName}_\${amp}.fasta
		# change fasta header with sample and amplicon names
		sed -i "s/>.*/>${SampleName}_\${amp}_consensus/" ${SampleName}_\${amp}.fasta
	done < "${SampleName}_mappedreads.txt"
	# merge consensus from all amplicons
	cat ${SampleName}_*.fasta > ${SampleName}_consensus.fasta
	# insert headers to mappedreads.txt
	sed -i '1i Amplicon_Name Size ${SampleName}' "${SampleName}_mappedreads.txt"
	"""
}

//multiqc

process multiqc {
	publishDir "${params.outdir}/multiqc/",mode:"copy"
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
	publishDir "${params.outdir}/kraken2/",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path (SamplePath)
	path(db_path)
	
	output:
	path ("${SampleName}_kraken.csv")
	path ("${SampleName}_kraken_report.csv"),emit:(kraken2_raw)
	
	script:
	"""
	kraken2 --db $db_path --output ${SampleName}_kraken.csv --report ${SampleName}_kraken_report.csv --threads 1 ${SamplePath}
	"""
}
process kraken2_consensus {
	publishDir "${params.outdir}/kraken2_cons/",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path (SamplePath)
	path(db_path)
	
	output:
	path ("${SampleName}_cons_kraken.csv")
	path ("${SampleName}_cons_kraken_report.csv"),emit:(kraken2_cons)

	script:
	"""
	kraken2 --db $db_path --output ${SampleName}_cons_kraken.csv --report ${SampleName}_cons_kraken_report.csv --threads 1 ${SamplePath}
	"""
}

//centrifuge for taxonomy classification
process centrifuge {
	publishDir "${params.outdir}/centrifuge/",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path (SamplePath)
	path(db_path)
	
	output:
	path ("${SampleName}_cent_report.csv")
	path ("${SampleName}_centrifuge_kstyle.csv"),emit:(centrifuge_raw)

	script:
	"""
	index=\$(find -L ${db_path} -name "*.1.cf" -not -name "._*"  | sed 's/.1.cf//')
	centrifuge -x \${index} -U ${SamplePath} -q -S ${SampleName}_centrifuge.csv --report ${SampleName}_cent_report.csv
	centrifuge-kreport -x \${index} "${SampleName}_centrifuge.csv" > ${SampleName}_centrifuge_kstyle.csv
	
	"""
}
process centrifuge_consensus {
	publishDir "${params.outdir}/centrifuge_cons/",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path (SamplePath)
	path(db_path)
	
	output:
	path ("${SampleName}_consensus_centrifuge.csv")
	path ("${SampleName}_consensus_kstyle.csv"),emit:(centrifuge_cons)
	
	script:
	"""
	index=\$(find -L ${db_path} -name "*.1.cf" -not -name "._*"  | sed 's/.1.cf//')
	centrifuge -x \${index} -U ${SamplePath} -f -S ${SampleName}_consensus_centrifuge.csv 
	centrifuge-kreport -x \${index} "${SampleName}_consensus_centrifuge.csv" > ${SampleName}_consensus_kstyle.csv
	"""
}

//krona plots
process krona_kraken {
	publishDir "${params.outdir}/krona_kraken/",mode:"copy"
	label "low"
	input:
	path(raw)
	path(consensus)
	
	output:
	path ("rawreads_classified.html")
	path("consensus_classified.html")
	script:
	"""
	ktImportTaxonomy -t 5 -m 3 -o rawreads_classified.html ${raw}
	ktImportTaxonomy -t 5 -m 3 -o consensus_classified.html ${consensus}
	"""
}
process krona_centrifuge {
	publishDir "${params.outdir}/krona_centrifuge/",mode:"copy"
	label "low"
	input:
	path(raw)
	path(consensus)
	
	output:
	path ("rawreads_classified.html")
	path("consensus_classified.html")
	script:
	"""
	ktImportTaxonomy -t 5 -m 3 -o rawreads_classified.html ${raw}
	ktImportTaxonomy -t 5 -m 3 -o consensus_classified.html ${consensus}
	"""
}


workflow {
	data=Channel
	.fromPath(params.input)
	.splitCsv(header:true)
    .map { row-> tuple(row.SampleName,row.SamplePath) }
	reference=file(params.reference)
	primerbed=file(params.primerbed)
    merge_fastq(data)
//trim barcodes and adapter sequences
	if (params.trim_barcodes){
		porechop(merge_fastq.out)
		minimap2(reference,porechop.out)
		 
	} else {
            minimap2(reference,merge_fastq.out)
		
        }
	// conditional for trim barcodes option
	if (params.trim_barcodes){
		if (params.kraken_db) {
			kraken=params.kraken_db
			kraken2(porechop.out,kraken)
		}
		if (params.centri_db){
			centri=params.centri_db
			centrifuge(porechop.out,centri)
		}           
			  
	 } else {
		if (params.kraken_db){
			kraken=params.kraken_db
			kraken2(merge_fastq.out,kraken)
		}
		if (params.centri_db){
			centri=params.centri_db
			centrifuge(merge_fastq.out,centri)
		}
	}
	splitbam(minimap2.out,primerbed)
	
	//condition for kraken2 classification
	if (params.kraken_db){
		kraken=params.kraken_db
		kraken2_consensus(splitbam.out.consensus,kraken)
		kraken_raw=kraken2.out.kraken2_raw
		kraken_cons=kraken2_consensus.out.kraken2_cons
		krona_kraken(kraken_raw.collect(),kraken_cons.collect())
		
	}
	//condition for centrifuge classification
	if (params.centri_db){
		centri=params.centri_db
		centrifuge_consensus(splitbam.out.consensus,centri)
		centri_raw=centrifuge.out.centrifuge_raw
		centri_cons=centrifuge_consensus.out.centrifuge_cons
		krona_centrifuge(centri_raw.collect(),centri_cons.collect())
	}
	stats=splitbam.out.stats
	idxstats=splitbam.out.idxstats
	multiqc(stats.mix(idxstats).collect())

}
