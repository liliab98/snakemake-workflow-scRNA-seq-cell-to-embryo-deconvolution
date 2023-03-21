rule filter_for_unambigous_alignments1:
    input:
        genome1="results/alignment/{id}_{sample}_Aligned.out.genome1.bam",
        genome2="results/alignment/{id}_{sample}_Aligned.out.genome2.bam",
    output:
        temp("results/snp/{id}_{sample}_splitread.bed"),
    log:
        "logs/filter_for_unambigous_alignments1/{id}_{sample}.log",
    benchmark:
        "benchmarks/filter_for_unambigous_alignments1/{id}_{sample}.txt"
    conda:
        "../envs/tools.yaml"
    shell:
        """
        samtools merge -n - {input.genome1} {input.genome2} | 
        samtools view -h - | 
        sed 's/XX:Z:./XX:i:/' | 
        perl -ane 'if($_=~m/NH:i:/){{if($_=~m/NH:i:1\t/){{print $_}}}}else{{print $_}}' | 
        samtools view -bS - | 
        bamToBed -split -tag XX -i - | 
        perl -ane 'print "chr$F[0]\t$F[2]\t$F[3]\t$F[4]\tG$F[5]\t$F[6]\n"' | 
        sort -k4,4 > {output} 2> {log}
        """


rule filter_for_unambigous_alignments2:
    input:
        "results/snp/{id}_{sample}_splitread.bed",
    output:
        temp("results/snp/{id}_{sample}_unambiguous.tmp"),
    log:
        "logs/filter_for_unambigous_alignments2/{id}_{sample}.log",
    benchmark:
        "benchmarks/filter_for_unambigous_alignments2/{id}_{sample}.txt"
    conda:
        "../envs/tools.yaml"
    shell:
        """
        less {input} | 
        bedtools groupby -g 4 -c 1,2,3,5 -o distinct,collapse,collapse,distinct | 
        perl -ane 'if($F[1]!~m/,/ && $F[4]!~m/,/){{@s=split(/,/,$F[2]); @e=split(/,/,$F[3]); for($i=0; $i<scalar(@s); $i++){{print "$F[1]\t$s[$i]\t$e[$i]\t$F[0]\t$F[4]\n"}}}}' > {output} 2> {log}
        """


rule filter_for_unambigous_alignments3:
    input:
        "results/snp/{id}_{sample}_unambiguous.tmp",
    output:
        temp("results/snp/{id}_{sample}_unambiguous.bed"),
    log:
        "logs/filter_for_unambigous_alignments3/{id}_{sample}.log",
    benchmark:
        "benchmarks/filter_for_unambigous_alignments3/{id}_{sample}.txt"
    conda:
        "../envs/tools.yaml"
    shell:
        """
        bedtools sort -i {input} > {output} 2> {log}
        """


rule assignment_of_reads_to_SNP:
    input:
        "results/snp/{id}_{sample}_unambiguous.bed",
    output:
        temp("results/snp/{id}_{sample}_SNP.tsv"),
    params:
        all_snps_bed=config["all_snps_bed_path"],
    log:
        "logs/assignment_of_reads_to_SNP/{id}_{sample}.log",
    benchmark:
        "benchmarks/assignment_of_reads_to_SNP/{id}_{sample}.txt"
    conda:
        "../envs/tools.yaml"
    shell:
        """
        bedtools intersect -sorted -wa -wb -a {input} -b {params.all_snps_bed} | 
        perl -ane 'print "$F[8]\t$F[5]\t$F[3]\t$F[4]\n"' > {output} 2> {log}
        """


rule get_only_white_list_SNPs:
    input:
        "results/snp/{id}_{sample}_SNP.tsv",
    output:
        whitelistSNP="results/snp/{id}_{sample}_whitelistSNP.tsv",
        allSNP="results/snp/{id}_{sample}_allSNP.tsv",
    params:
        SNPwhite_list=config["SNPwhite_list_path"],
    log:
        "logs/get_only_white_list_SNPs/{id}_{sample}.log",
    benchmark:
        "benchmarks/get_only_white_list_SNPs/{id}_{sample}.txt"
    conda:
        "../envs/tools.yaml"
    shell:
        """
        fgrep -w -f {params.SNPwhite_list} {input} | 
        perl -ane 'print "$F[1].$F[3]\t$F[2]\t$F[0]\n"' | 
        sort -k2,2 | 
        uniq > {output.whitelistSNP} 2> {log}
        cat {params.SNPwhite_list} {input} | 
        perl -ane 'print "$F[1].$F[3]\t$F[2]\t$F[0]\n"' | 
        sort -k2,2 -k3,3 | 
        uniq > {output.allSNP} 2> {log}
        """


rule convert_read_to_BC_information:
    input:
        whitelistSNP="results/snp/{id}_{sample}_whitelistSNP.tsv",
        allSNP="results/snp/{id}_{sample}_allSNP.tsv",
        readBC="results/barcodes/{id}_{sample}_read.BC.tsv",
    output:
        SNPcount="results/snp/{id}_{sample}_SNPcount.tsv",
        allSNPcount="results/snp/{id}_{sample}_allSNPcount.tsv",
    log:
        "logs/convert_read_to_BC_information/{id}_{sample}.log",
    benchmark:
        "benchmarks/convert_read_to_BC_information/{id}_{sample}.txt"
    conda:
        "../envs/tools.yaml"
    shell:
        """
        join -1 2 -2 1 {input.whitelistSNP} {input.readBC} | 
        perl -ane 'print "$F[3]\t$F[1]\t$F[2]\n"' | 
        sort | 
        bedtools groupby -g 1,2 -c 3 -o  count_distinct > {output.SNPcount} 2> {log}
        join -1 2 -2 1 {input.allSNP} {input.readBC} | 
        perl -ane 'print "$F[3]\t$F[1]\t$F[2]\n"' | 
        sort | 
        bedtools groupby -g 1,2 -c 3 -o  count_distinct > {output.allSNPcount} 2> {log}
        """
