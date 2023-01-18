#configfile: "config.yaml"	
Barcodes = ["AM-RNA-496", "AM-RNA-497"]

rule all:
	input: 
		expand("{BC}_Aligned.out.genome1.bam", BC=Barcodes)

# BC1 und BC2 get shortened
rule get_BCs_surviving_10X_pipeline:
	input:
		"barcodes/{BC}/barcodes.tsv.gz"
	output:
		"barcodes/WT_E85_{BC}_10X_barcodes.tsv"
	shell:
		"""
		if file --mime-type {input} | grep -q gzip; then
		zcat {input} | cut -f1 -d'-' > {output}
		else
		cut -f1 -d'-' {input} > {output}
		fi
		"""

rule assignment_of_BC_2_read_reduce_to_10X_BCs_part1: 
	input:
		"data/test_{BC}_R1.fastq.gz"
	output:
		temp("barcodes/WT_E85_{BC}_read_tmp.BC.tsv")
	shell: 
		"""
		zcat {input} | sed 's/ .*//' | paste - - - - | cut -f1,2 | sed 's/^@//' | perl -ane '$F[1]=substr($F[1], 0, 16); print "$F[0]\t$F[1]\n"' > {output} 
		"""

rule assignment_of_BC_2_read_reduce_to_10X_BCs_part2: 
	input:
		bc="barcodes/WT_E85_{BC}_10X_barcodes.tsv",
		bc_read_tmp="barcodes/WT_E85_{BC}_read_tmp.BC.tsv"
	output:
		"barcodes/WT_E85_{BC}_read.BC.tsv"
	shell: 
		"""
		fgrep -f {input.bc} {input.bc_read_tmp} | sort > {output}
		"""
		
#rule join_fastq_files:
#	input: 
#		"data/test_{BC}_R2.fastq.gz"
#	output:
#		"{BC}_R2.fastq.gz"
#	shell:
#		"cat {input} > {output}"

rule alignment_using_STAR:
	input:
		"data/test_{BC}_R2.fastq.gz"
	output:
		temp("{BC}_Aligned.out.bam"),
		temp("{BC}_Log.final.out"),
		temp("{BC}_Log.out"),
		temp("{BC}_Log.progress.out"),
		temp("{BC}_SJ.out.tab"),
		temp(directory("{BC}__STARtmp/"))
	shell:
		"""
		/project/bioinf_meissner/src/STAR/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 20 --genomeDir /project/bioinf_meissner/src/SNPsplit/SNPsplit_v0.3.2 --readFilesIn {input} --readFilesCommand zcat --alignEndsType EndToEnd --outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted --outFileNamePrefix {wildcards.BC}_
		"""

rule assignment_of_reads_to_genome1:
	input:
		"{BC}_Aligned.out.bam"
	output:
		temp("{BC}_Aligned.out.SNPsplit_sort.txt"), 
		temp("{BC}_Aligned.out.SNPsplit_report.txt"),
		temp("{BC}_Aligned.out.unassigned.bam"), 
		"{BC}_Aligned.out.allele_flagged.bam",
		"{BC}_Aligned.out.genome1.bam",
		"{BC}_Aligned.out.genome2.bam"
	shell:
		"""
		perl /project/bioinf_meissner/src/SNPsplit/SNPsplit_v0.3.2/SNPsplit --snp_file /project/bioinf_meissner/src/SNPsplit/SNPsplit_v0.3.2/all_SNPs_CAST_EiJ_GRCm38.txt.gz {input} --samtools_path /project/bioinf_meissner/src/samtools/samtools-1.6/
		"""

#rule clean-up:
#	input:
#	output:
#	shell:
#		"""
#		rm -r ${SAMPLE}__STARtmp
#		"""
