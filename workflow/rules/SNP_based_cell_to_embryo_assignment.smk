rule get_BCs_surviving_10X_pipeline:
    input:
        lambda wildcards: f"{config['Barcodes'][wildcards.sample]}{config['ID']}_barcodes.tsv.gz",
    output:
        temp("results/barcodes/{id}_{sample}_10X_barcodes.tsv"),
    log:
        "logs/get_BCs_surviving_10X_pipeline/{id}_{sample}.log",
    benchmark:
        "benchmarks/get_BCs_surviving_10X_pipeline/{id}_{sample}.txt"
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
        lambda wildcards: f"{config['Samples'][wildcards.sample]}_R1_001.fastq.gz",
    output:
        temp("results/barcodes/{id}_{sample}_read.BC.tmp"),
    log:
        "logs/assignment_of_BC_2_read_reduce_to_10X_BCs_part1/{id}_{sample}.log",
    benchmark:
        "benchmarks/assignment_of_BC_2_read_reduce_to_10X_BCs_part1/{id}_{sample}.txt"
    conda:
        "../envs/tools.yaml"
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
        bc="results/barcodes/{id}_{sample}_10X_barcodes.tsv",
        bc_read_tmp="results/barcodes/{id}_{sample}_read.BC.tmp",
    output:
        "results/barcodes/{id}_{sample}_read.BC.tsv",
    log:
        "logs/assignment_of_BC_2_read_reduce_to_10X_BCs_part2/{id}_{sample}.log",
    benchmark:
        "benchmarks/assignment_of_BC_2_read_reduce_to_10X_BCs_part2/{id}_{sample}.txt"
    shell:
        """
        fgrep -f {input.bc} {input.bc_read_tmp} | 
        sort -k1,1 -k2,2 > {output} 2> {log}
        """


# genome indices are already generated
rule alignment_using_STAR:
    input:
        lambda wildcards: f"{config['Samples'][wildcards.sample]}_R2_001.fastq.gz",
    output:
        "results/alignment/{id}_{sample}_Aligned.out.bam",
        temp("results/alignment/{id}_{sample}_Log.final.out"),
        temp("results/alignment/{id}_{sample}_Log.out"),
        temp("results/alignment/{id}_{sample}_Log.progress.out"),
        temp("results/alignment/{id}_{sample}_SJ.out.tab"),
        temp(directory("results/alignment/{id}_{sample}__STARtmp/")),
    params:
        genome_dir=config["genome_dir_path"],
    log:
        "logs/alignment_using_STAR/{id}_{sample}.log",
    benchmark:
        "benchmarks/alignment_using_STAR/{id}_{sample}.log"
    threads: 20
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --outBAMsortingThreadN {threads} \
        --genomeDir {params.genome_dir} \
        --readFilesIn {input} \
        --readFilesCommand zcat \
        --alignEndsType EndToEnd \
        --outSAMattributes NH HI NM MD \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix results/alignment/{wildcards.id}_{wildcards.sample}_ 2>&1 {log}
        """


rule assignment_of_reads_to_genome1:
    input:
        "results/alignment/{id}_{sample}_Aligned.out.bam",
    output:
        temp("results/alignment/{id}_{sample}_Aligned.out.SNPsplit_sort.txt"),
        temp("results/alignment/{id}_{sample}_Aligned.out.SNPsplit_report.txt"),
        temp("results/alignment/{id}_{sample}_Aligned.out.unassigned.bam"),
        temp("results/alignment/{id}_{sample}_Aligned.out.allele_flagged.bam"),
        "results/alignment/{id}_{sample}_Aligned.out.genome1.bam",
        "results/alignment/{id}_{sample}_Aligned.out.genome2.bam",
    params:
        all_snps=config["all_snps_path"],
    log:
        "logs/assignment_of_reads_to_genome1/{id}_{sample}.log",
    benchmark:
        "benchmarks/assignment_of_reads_to_genome1/{id}_{sample}.log"
    conda:
        "../envs/snpsplit.yaml"
    shell:
        """
        SNPsplit --snp_file {params.all_snps} {input} 2> {log}
        mv {wildcards.id}_{wildcards.sample}* results/alignment/
        """
