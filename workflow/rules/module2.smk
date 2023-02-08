rule filter_for_unambigous_alignments1:
    input:
        genome1 = "results/alignment/{sample}_Aligned.out.genome1.bam",
        genome2 = "results/alignment/{sample}_Aligned.out.genome2.bam"
    output:
        "results/snp/{sample}_splitread.bed"
    log:
        "logs/filter_for_unambigous_alignments1/{sample}.log"
    shell:
        """
        /project/bioinf_meissner/src/samtools/samtools-1.6/samtools merge -n - {input.genome1} {input.genome2} | 
        /project/bioinf_meissner/src/samtools/samtools-1.6/samtools view -h - | 
        sed 's/XX:Z:./XX:i:/' | 
        perl -ane 'if($_=~m/NH:i:/){{if($_=~m/NH:i:1\t/){{print $_}}}}else{{print $_}}' | 
        /project/bioinf_meissner/src/samtools/samtools-1.6/samtools view -bS - | 
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bamToBed -split -tag XX -i - | 
        perl -ane 'print "chr$F[0]\t$F[2]\t$F[3]\t$F[4]\tG$F[5]\t$F[6]\n"' | 
        sort -k4,4 > {output} 2> {log}
        """

rule filter_for_unambigous_alignments2:
    input:
        "results/snp/{sample}_splitread.bed"
    output:
        temp("results/snp/{sample}_unambiguous.tmp")
    log:
        "logs/filter_for_unambigous_alignments2/{sample}.log"
    shell:
        """
        less {input} | 
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools groupby -g 4 -c 1,2,3,5 -o distinct,collapse,collapse,distinct | 
        perl -ane 'if($F[1]!~m/,/ && $F[4]!~m/,/){{@s=split(/,/,$F[2]); @e=split(/,/,$F[3]); for($i=0; $i<scalar(@s); $i++){{print "$F[1]\t$s[$i]\t$e[$i]\t$F[0]\t$F[4]\n"}}}}' > {output} 2> {log}
        """

rule filter_for_unambigous_alignments3:
    input:
        "results/snp/{sample}_unambiguous.tmp"
    output:
        temp("results/snp/{sample}_unambiguous.bed")
    log:
        "logs/filter_for_unambigous_alignments3/{sample}.log"
    shell:
        """
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools sort -i {input} > {output} 2> {log}
        """

rule assignemnt_of_reads_to_SNP:
    input:
        "results/snp/{sample}_unambiguous.bed"
    output:
        temp("results/snp/{sample}_SNP.tsv")
    log:
        "logs/assignemnt_of_reads_to_SNP/{sample}.log"
    shell:
        """
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools intersect -sorted -wa -wb -a {input} -b /project/bioinf_meissner/src/SNPsplit/SNPsplit_v0.3.2/all_SNPs_CAST_EiJ_GRCm38.bed | 
        perl -ane 'print "$F[8]\t$F[5]\t$F[3]\t$F[4]\n"' > {output} 2> {log}
        """

rule get_only_white_list_SNPs:
    input:
        "results/snp/{sample}_SNP.tsv"
    output:
        whitelistSNP = "results/snp/{sample}_whitelistSNP.tsv",
        allSNP = "results/snp/{sample}_allSNP.tsv"
    params: 
        scrna= config["whitelist"]
    log:
        "logs/get_only_white_list_SNPs/{sample}.log"
    shell:
        """
        fgrep -w -f {params.scrna} {input} | 
        perl -ane 'print "$F[1].$F[3]\t$F[2]\t$F[0]\n"' | 
        sort -k2,2 | 
        uniq > {output.whitelistSNP} 2> {log}
        cat {params.scrna} {input} | 
        perl -ane 'print "$F[1].$F[3]\t$F[2]\t$F[0]\n"' | 
        sort -k2,2 -k3,3 | 
        uniq > {output.allSNP} 2> {log}
        """

rule convert_read_to_BC_information:
    input:
        whitelistSNP = "results/snp/{sample}_whitelistSNP.tsv",
        allSNP = "results/snp/{sample}_allSNP.tsv",
        readBC = "results/barcodes/WT_E85_{sample}_read.BC.tsv"
    output:
        SNPcount = "results/snp-count/{sample}_SNPcount.tsv",
        allSNPcount = "results/snp-count/{sample}_allSNPcount.tsv"
    log:
        "logs/convert_read_to_BC_information/{sample}.log"
    shell:
        """
        join -1 2 -2 1 {input.whitelistSNP} {input.readBC}| 
        perl -ane 'print "$F[3]\t$F[1]\t$F[2]\n"' | 
        sort | 
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools groupby -g 1,2 -c 3 -o  count_distinct > {output.SNPcount} 2> {log}
        join -1 2 -2 1 {input.allSNP} {input.readBC} | 
        perl -ane 'print "$F[3]\t$F[1]\t$F[2]\n"' | 
        sort | 
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools groupby -g 1,2 -c 3 -o  count_distinct > {output.allSNPcount} 2> {log}
        """