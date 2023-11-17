#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// make csv file with headers from the given input

process make_csv {
	publishDir "${params.outdir}"
	input:
	path(fastq_input)
	output:
	path("samplelist.csv")
	
	script:
	"""
	# get list of directories and replaced slash

	ls -d */ ${fastq_input} > samples.csv

	sed -i 's#/##g' samples.csv

	# get path of each folder
	realpath ${fastq_input}/* > paths.csv

	# concatenate samplenames and path with comma as delimiter
	paste samples.csv paths.csv > samplelist.csv
	sed -i 's/	/,/g' samplelist.csv
	
	# add headers to the csv file
	sed -i '1i SampleName,SamplePath' samplelist.csv
	
	"""

}

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
	
	
		if [[ "\${count}" != "0" ]]
		then
			cat ${SamplePath}/*.fastq.gz > ${SampleName}.fastq.gz
		
		else
			count=\$(ls -1 ${SamplePath}/*.fastq 2>/dev/null | wc -l)
			if [[ "\${count}" != "0" ]]
			then
				cat ${SamplePath}/*.fastq > ${SampleName}.fastq
			fi
		fi
	"""
}

//trim barcodes and adapter using porechop

process porechop {
	label "high"
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

//convert minimap2 output sam to sorted bam and split bam files and create consensus
process splitbam {
	publishDir "${params.outdir}/splitbam",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path(SamplePath)
	path (primerbed)
	output:
	val(SampleName),emit:SampleName
	path("${SampleName}_stats.txt"),emit:stats
	path("${SampleName}_mappedreads.txt"),emit:mapped
	path("${SampleName}_idxstats.txt"),emit:idxstats
	tuple val(SampleName),path("${SampleName}_consensus.fasta"),emit:consensus
	path("${SampleName}_unfilt_stats.txt"),emit:unfilt_stats
	path("${SampleName}_flagstat.txt"),emit:flagstat
	path ("${SampleName}_unfilt_idxstats.csv"),emit:unfilt_idxstats
	script:
	"""
# generate stats prior to read filtering
	samtools view -b -h ${SamplePath}|samtools sort > ${SampleName}_unfilt.bam
	samtools stats "${SampleName}_unfilt.bam" > ${SampleName}_unfilt_stats.txt
	samtools flagstat "${SampleName}_unfilt.bam" > ${SampleName}_flagstat.txt
	samtools index "${SampleName}_unfilt.bam" > ${SampleName}_unfilt.bai
	samtools idxstats "${SampleName}_unfilt.bam" > ${SampleName}_unfilt_idxstats.csv
#generate a  sorted bam file with primary alignments
	samtools view -b -h -F 0x900 -q 30 ${SamplePath}|samtools sort > ${SampleName}.bam	
	samtools stats "${SampleName}.bam" > ${SampleName}_stats.txt
#trims primers from both ends of the ampliocn using primer bed file
	samtools ampliconclip --both-ends -b ${primerbed} "${SampleName}.bam"	> ${SampleName}_trimmed.bam
# sorts the bam file for spitting	
	samtools sort "${SampleName}_trimmed.bam" > ${SampleName}_tr_sorted.bam
#index sorted bam file and generate read counts for each amplicon
	samtools index "${SampleName}_tr_sorted.bam" > ${SampleName}_sorted.bai
	samtools idxstats "${SampleName}_tr_sorted.bam" > ${SampleName}_idxstats.txt
	awk '{if (\$3 >= 10) print \$1,\$2,\$3}' "${SampleName}_idxstats.txt" > ${SampleName}_mappedreads.txt
#conditional for empty mapped reads.txt file
	if [ \$(wc -l < "${SampleName}_mappedreads.txt") -ne 0 ]
	then 
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
	else
		echo -e ">${SampleName} No consensus\n" > ${SampleName}_consensus.fasta

	fi
	# insert headers to mappedreads.txt
	sed -i '1i Amplicon_Name Size ${SampleName}' "${SampleName}_mappedreads.txt"
	"""
}


//multiqc generate mapped read statistics from samtools output

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
	path ("${SampleName}_cons_kraken.csv"),emit:(kraken2_cons)
	path ("${SampleName}_cons_kraken_report.csv")

	script:
	"""
	kraken2 --db $db_path --output ${SampleName}_cons_kraken.csv --report ${SampleName}_cons_kraken_report.csv --threads 3 ${SamplePath} --use-names --use-mpa-style
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
	index=\$(find -L ${db_path} -name "*.2.cf" -not -name "._*"  | sed 's/.2.cf//')
	centrifuge -x \${index} -U ${SamplePath} -q -S ${SampleName}_centrifuge.csv --report ${SampleName}_cent_report.csv -p 2
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
	index=\$(find -L ${db_path} -name "*.2.cf" -not -name "._*"  | sed 's/.2.cf//')
	centrifuge -x \${index} -U ${SamplePath} -f -S ${SampleName}_consensus_centrifuge.csv -p 2
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
	path ("rawreads_classified.html"),emit:raw
	path("consensus_classified.html"),emit:cons
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
	path ("rawreads_classified.html"),emit:raw
	path("consensus_classified.html"),emit:cons
	script:
	"""
	ktImportTaxonomy -t 5 -m 3 -o rawreads_classified.html ${raw}
	ktImportTaxonomy -t 5 -m 3 -o consensus_classified.html ${consensus}
	"""
}
//make html reprot with rmarkdown
process make_report {
	publishDir "${params.outdir}/results_report/",mode:"copy"
	label "low"
	input:
	path (csv)
	path(krona_reports_raw)
	path(mappedreads)
	//path(kraken_cons)
	path(abricate)
	path(rmdfile)
	path(mlst)
	output:
	path("Bovreproseq_results_report.html")
	script:
	"""
	
	cp ${csv} samples.csv
	cp ${krona_reports_raw} rawreads.html
	# handle empty mapped reads files
	for i in *mappedreads.txt
	do
	 	if [ \$(wc -l < "\${i}" ) -eq 0 ]
		 then
	 		echo "Amplicon_Name Size Reads" >> \${i}
			echo "NA NA NA" >> \${i}
	 	fi
	done
	# handle empty kraken consensus files
	#for k in *_cons_kraken.csv
	#do
	#	#if [ \$(wc -l < "\${k}" ) -eq 0 ]
	#	#then
	#		echo "C	NO READS FOUND	NA" >> \${k}	
	 #	fi
	#done

	cp ${rmdfile} report.Rmd
	

	Rscript -e 'rmarkdown::render(input="report.Rmd",params=list(csv="samples.csv",krona="rawreads.html"),output_file = "Bovreproseq_results_report.html")'
	"""

}
// performs remote blast of the consensus sequences
process blast_cons {
	publishDir "${params.outdir}/blast/",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path(consensus)
	path (taxdb)
	output:
	path("${SampleName}_report_blast.csv")
	
	script:
	"""

	cp ${taxdb}/* ./
	blastn -db nt -query ${consensus} -out ${SampleName}_blast.csv -outfmt "7 qseqid sseqid length qcovs pident evalue" -max_target_seqs 1 -remote
	grep -v "#" ${SampleName}_blast.csv|sort|uniq > ${SampleName}_report_blast.csv
	sed -i '1i queryid\tsubject_id\talignment length\tquery_coverage\t%identity\tevalue\tscinames' ${SampleName}_report_blast.csv
	
	"""

}


// uses custome database to predict the presence of the amplicons
process abricate{
	publishDir "${params.outdir}/abricate/",mode:"copy"
	label "medium"
	input:
	tuple val(SampleName),path(consensus)
	path(dbdir)
	output:
	path("${SampleName}_abricate.csv")
	script:
	"""
	abricate --datadir ${dbdir} --db Bovreproseq -minid 60  -mincov 60 --quiet ${consensus} 1> ${SampleName}_abricate.csv
	"""
	
}

process mlst {
	publishDir "${params.outdir}/mlst/",mode:"copy"
	label "medium"
	input:
	tuple val(SampleName),path(consensus)
	output:
	path("${SampleName}_MLST_results.csv")
	script:
	"""
	mlst --legacy --scheme campylobacter_nonjejuni_9 ${consensus} > ${SampleName}_MLST.csv
	shopt -s extglob
	ST="\$(tail -n +2 ${SampleName}_MLST.csv | cut -f 3)"
	if [  \$ST -eq "4" -o \$ST -eq "7" -o \$ST -eq "12" ];then
		echo "ORGANISM" > temp.csv
		echo "Campylobacter fetus. venerealis" >> temp.csv
		paste ${SampleName}_MLST.csv temp.csv > ${SampleName}_MLST_results.csv
	
	elif [ \$ST -eq "1" -o \$ST -eq "2" -o \$ST -eq "3" -o \$ST -eq "5" -o \$ST -eq "6" -o \$ST -eq "8" -o \$ST -eq "9" -o \$ST -eq "10" -o \$ST -eq "11" -o \$ST -eq "13" -o \$ST -eq "14" ];then
		echo "ORGANISM" > temp.csv
		echo "Campylobacter fetus. fetus" >> temp.csv
		paste ${SampleName}_MLST.csv temp.csv > ${SampleName}_MLST_results.csv
		
	else 
		
		echo "ORGANISM" > temp.csv
		echo "NA" >> temp.csv
		paste ${SampleName}_MLST.csv temp.csv > ${SampleName}_MLST_results.csv
		
	fi
	"""
}



workflow {
	data=Channel
	.fromPath(params.input)
	merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)})
	reference=file("${baseDir}/Bovreproseq_reference.fasta")
	primerbed=file("${baseDir}/Bovreproseq_primer.bed")
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

	// create consensus
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
	// qc report using split bam out put
	stats=splitbam.out.unfilt_stats
	idxstats=splitbam.out.idxstats
	multiqc(stats.mix(idxstats).collect())
	
	// abricate 
	dbdir=("${baseDir}/Bovreproseq_db")
	abricate(splitbam.out.consensus,dbdir)
	
	//tax=("${baseDir}/taxdb")
	//blast_cons(splitbam.out.consensus,tax,db1)

	mlst(splitbam.out.consensus)
	//generate report
	rmd_file=file("${baseDir}/Bovreproseq_tabbed.Rmd")
	if (params.kraken_db){
		make_report(make_csv.out,krona_kraken.out.raw,splitbam.out.mapped.collect(),abricate.out.collect(),rmd_file,mlst.out.collect())
	}
	if (params.centri_db){
		make_report(make_csv.out,krona_centrifuge.out.raw,splitbam.out.mapped.collect(),abricate.out.collect(),rmd_file,mlst.out.collect())
	}
	
	
}
