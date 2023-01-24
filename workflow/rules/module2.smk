
rule filter_for_unambigous_alignments1:
    input:
        genome1 = "results/alignment/{BC}_Aligned.out.genome1.bam",
        genome2 = "results/alignment/{BC}_Aligned.out.genome2.bam"
    output:
        "results/splitreads/{BC}_splitread.bed"
    shell:
        """
        /project/bioinf_meissner/src/samtools/samtools-1.6/samtools merge -n - {input.genome1} {input.genome2} | 
        /project/bioinf_meissner/src/samtools/samtools-1.6/samtools view -h - | 
        sed 's/XX:Z:./XX:i:/' | 
        perl -ane 'if($_=~m/NH:i:/){{if($_=~m/NH:i:1\t/){{print $_}}}}else{{print $_}}' | 
        /project/bioinf_meissner/src/samtools/samtools-1.6/samtools view -bS - | 
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bamToBed -split -tag XX -i - | 
        perl -ane 'print "chr$F[0]\t$F[2]\t$F[3]\t$F[4]\tG$F[5]\t$F[6]\n"' | 
        sort -k4,4 > {output}
        """

rule filter_for_unambigous_alignments2:
    input:
        "results/splitreads/{BC}_splitread.bed"
    output:
        temp("results/splitreads/{BC}_unambiguous.tmp")
    shell:
        """
        less {input} | 
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools groupby -g 4 -c 1,2,3,5 -o distinct,collapse,collapse,distinct | 
        perl -ane 'if($F[1]!~m/,/ && $F[4]!~m/,/){{@s=split(/,/,$F[2]); @e=split(/,/,$F[3]); for($i=0; $i<scalar(@s); $i++){{print "$F[1]\t$s[$i]\t$e[$i]\t$F[0]\t$F[4]\n"}}}}' > {output}
        """

rule filter_for_unambigous_alignments3:
    input:
        "results/splitreads/{BC}_unambiguous.tmp"
    output:
        temp("results/splitreads/{BC}_unambiguous.bed")
    shell:
        """
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools sort -i {input} > {output}
        """

rule assignemnt_of_reads_to_SNP:
    input:
        "results/splitreads/{BC}_unambiguous.bed"
    output:
        temp("results/splitreads/{BC}_SNP.tsv")
    shell:
        """
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools intersect -sorted -wa -wb -a {input} -b /project/bioinf_meissner/src/SNPsplit/SNPsplit_v0.3.2/all_SNPs_CAST_EiJ_GRCm38.bed | 
        perl -ane 'print "$F[8]\t$F[5]\t$F[3]\t$F[4]\n"' > {output}
        """

rule get_only_white_list_SNPs:
    input:
        "results/splitreads/{BC}_SNP.tsv"
    output:
        whitelistSNP = "results/splitreads/{BC}_whitelistSNP.tsv",
        allSNP = "results/splitreads/{BC}_allSNP.tsv"
    shell:
        """
        fgrep -w -f data/scRNA/WT_SNP_white.list.tsv {input} | 
        perl -ane 'print "$F[1].$F[3]\t$F[2]\t$F[0]\n"' | 
        sort -k2,2 | 
        uniq > {output.whitelistSNP}
        cat data/scRNA/WT_SNP_white.list.tsv {input} | 
        perl -ane 'print "$F[1].$F[3]\t$F[2]\t$F[0]\n"' | 
        sort -k2,2 | 
        uniq > {output.allSNP}
        """

rule convert_read_to_BC_information:
    input:
        whitelistSNP = "results/splitreads/{BC}_whitelistSNP.tsv",
        allSNP = "results/splitreads/{BC}_allSNP.tsv",
        readBC = "results/barcodes/WT_E85_{BC}_read.BC.tsv"
    output:
        SNPcount = "results/splitreads/{BC}_SNPcount.tsv",
        allSNPcount = "results/splitreads/{BC}_allSNPcount.tsv"
    shell:
        """
        join -1 2 -2 1 {input.whitelistSNP} {input.readBC}| 
        perl -ane 'print "$F[3]\t$F[1]\t$F[2]\n"' | 
        sort | 
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools groupby -g 1,2 -c 3 -o  count_distinct > {output.SNPcount}
        join -1 2 -2 1 {input.allSNP} {input.readBC} | 
        perl -ane 'print "$F[3]\t$F[1]\t$F[2]\n"' | 
        sort | 
        /project/bioinf_meissner/src/bedtools/bedtools/bin/bedtools groupby -g 1,2 -c 3 -o  count_distinct > {output.allSNPcount}
        """