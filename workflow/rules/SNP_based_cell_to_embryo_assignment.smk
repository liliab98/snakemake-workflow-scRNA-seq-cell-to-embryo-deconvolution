rule get_BCs_surviving_10X_pipeline:
	input:
		lambda wildcards: f"{config['Barcodes'][wildcards.sample]}barcodes.tsv.gz"
	output:
		"results/barcodes/WT_E85_{sample}_10X_barcodes.tsv"
	log:
		"logs/get_BCs_surviving_10X_pipeline/{sample}.log"
	benchmark:
		"benchmarks/get_BCs_surviving_10X_pipeline/{sample}.txt"
	shell:
		"""
		if file --mime-type {input} | grep -q gzip$; then
		zcat {input} | cut -f1 -d'-' > {output} 2> {log}
		else
		cut -f1 -d'-' {input} > {output} 2> {log}
		fi
		"""

rule assignment_of_BC_2_read_reduce_to_10X_BCs_part1: 
	input:
		lambda wildcards: f"{config['Samples'][wildcards.sample]}_R1_001.fastq.gz"
	output:
		temp("results/barcodes/WT_E85_{sample}_read.BC.tmp")
	log:
		"logs/assignment_of_BC_2_read_reduce_to_10X_BCs_part1/{sample}.log"
	benchmark:
		"benchmarks/assignment_of_BC_2_read_reduce_to_10X_BCs_part1/{sample}.txt"
	shell: 
		"""
		zcat {input} | 
		sed 's/ .*//' | 
		paste - - - - | 
		cut -f1,2 | 
		sed 's/^@//' | 
		perl -ane '$F[1]=substr($F[1], 0, 16); print "$F[0]\t$F[1]\n"' > {output} 2> {log}
		"""

rule assignment_of_BC_2_read_reduce_to_10X_BCs_part2: 
	input:
		bc="results/barcodes/WT_E85_{sample}_10X_barcodes.tsv",
		bc_read_tmp="results/barcodes/WT_E85_{sample}_read.BC.tmp"
	output:
		"results/barcodes/WT_E85_{sample}_read.BC.tsv"
	log:
		"logs/assignment_of_BC_2_read_reduce_to_10X_BCs_part2/{sample}.log"
	benchmark:
		"benchmarks/assignment_of_BC_2_read_reduce_to_10X_BCs_part2/{sample}.txt"
	shell: 
		"""
		fgrep -f {input.bc} {input.bc_read_tmp} | 
		sort -k1,1 -k2,2 > {output} 2> {log}
		"""

#kann eig raus	
#rule join_fastq_files:
#	input: 
#		"data/test_{BC}_R2.fastq.gz"
#	output:
#		"{BC}_R2.fastq.gz"
#	shell:
#		"cat {input} > {output}"

rule alignment_using_STAR:
	input:
		lambda wildcards: f"{config['Samples'][wildcards.sample]}_R2_001.fastq.gz"
	output:
		"results/alignment/{sample}_Aligned.out.bam", #temp? 
		temp("results/alignment/{sample}_Log.final.out"),
		temp("results/alignment/{sample}_Log.out"),
		temp("results/alignment/{sample}_Log.progress.out"),
		temp("results/alignment/{sample}_SJ.out.tab"),
		temp(directory("results/alignment/{sample}__STARtmp/"))
	params:
		SNPsplit_dir = config["SNPsplit_dir_path"]
	log: 
		"logs/alignment_using_STAR/{sample}.log"
	benchmark:
		"benchmarks/alignment_using_STAR/{sample}.log"
	threads: 20
	conda:
		"../envs/alignment.yaml"
	shell:
		"""
		STAR \
		--runThreadN {threads} \
		--outBAMsortingThreadN {threads} \
		--genomeDir {params.SNPsplit_dir} \
		--readFilesIn {input} \
		--readFilesCommand zcat \
		--alignEndsType EndToEnd \
		--outSAMattributes NH HI NM MD \
		--outSAMtype BAM Unsorted \
		--outFileNamePrefix results/alignment/{wildcards.sample}_ 2>&1 {log}
		"""
#{params.star}
rule assignment_of_reads_to_genome1:
	input:
		"results/alignment/{sample}_Aligned.out.bam"
	output:
		temp("results/alignment/{sample}_Aligned.out.SNPsplit_sort.txt"), 
		temp("results/alignment/{sample}_Aligned.out.SNPsplit_report.txt"),
		temp("results/alignment/{sample}_Aligned.out.unassigned.bam"), 
		"results/alignment/{sample}_Aligned.out.allele_flagged.bam",
		"results/alignment/{sample}_Aligned.out.genome1.bam",	#tmp?
		"results/alignment/{sample}_Aligned.out.genome2.bam"	#tmp?
	params:
		SNPsplit = config["SNPsplit_path"],
		SNPfile_gz = config["SNPfile_gz_path"],
		samtools_dir = config["samtools_dir_path"]	
	log:
		"logs/assignment_of_reads_to_genome1/{sample}.log"
	benchmark:
		"benchmarks/assignment_of_reads_to_genome1/{sample}.log"
	threads: 8
	conda:
		"../envs/snp.yaml"
	shell:
		"""
		perl {params.SNPsplit} \
		--snp_file {params.SNPfile_gz} {input} \
		--samtools_path {params.samtools_dir} 2> {log} 
		"""
#why is 2>&1 not working?