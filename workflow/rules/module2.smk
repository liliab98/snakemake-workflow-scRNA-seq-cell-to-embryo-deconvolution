
rule filter_for_unambigous_alignments:
    input:
        genome1 = "{BC}_Aligned.out.genome1.bam",
        genome2 = "{BC}_Aligned.out.genome2.bam"
    output:
        "{BC}_splitread.bed"
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