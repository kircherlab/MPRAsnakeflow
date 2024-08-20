# rule assignment_bbmap_ref:
#     """
#     Create a ref for bbmap
#     """
#     conda:
#         "../../envs/bbmap_samtools_picard_htslib.yaml"
#     input:
#         reference="results/assignment/{assignment}/reference/reference.fa",
#         check="results/assignment/{assignment}/design_check.done",
#     output:
#         reference="results/assignment/{assignment}/reference/reference.fa.blablablu",
#     log:
#         temp("results/logs/assignment/bbmap_ref.{assignment}.log"),
#     shell:
#         """
#         OUTDIR=`dirname {output.reference}`;
#         bbmap.sh ref={input.reference} path=$OUTDIR> {log};
#         """


rule assignment_mapping_bbmap:
    """
    Map the reads to the reference and sort unsing bwa mem
    """
    conda:
        "../../envs/bbmap_samtools_htslib.yaml"
    input:
        reads="results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz",
        check="results/assignment/{assignment}/design_check.done",
        reference="results/assignment/{assignment}/reference/reference.fa",
    output:
        bam=temp(
            "results/assignment/{assignment}/bbmap/merge_split{split}.unsorted.bam"
        ),
        sorted_bam=temp(
            "results/assignment/{assignment}/bbmap/merge_split{split}.mapped.bam"
        ),
    log:
        temp("results/logs/assignment/mapping.bbmap.{assignment}.{split}.log"),
    shell:
        """
        bbmap.sh in={input.reads} ref={input.reference} nodisk -t={threads} out={output.bam} &> {log};
        samtools sort -l 0 -@ {threads} {output.bam} > {output.sorted_bam} &>> {log};
        """


rule assignment_mapping_bbmap_getBCs:
    """
    Get the barcodes.

    BAM/SAM fields:
    - bc_string = $2;
    - ref_name = $3;
    - alignement_start = $4;
    - mapping_quality = $5;
    - cigar = $6;
    - aligned_sequence = $10;
    - edit_distance = $12;
    - alignment_score = $13;
    """
    conda:
        "../../envs/bbmap_samtools_htslib.yaml"
    input:
        "results/assignment/{assignment}/bbmap/merge_split{split}.mapped.bam",
    output:
        temp("results/assignment/{assignment}/BCs/barcodes_bbmap.{split}.tsv"),
    params:
        mapping_quality_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["min_mapping_quality"],
    log:
        temp("results/logs/assignment/mapping.bbmap.getBCs.{assignment}.log"),
    shell:
        """
        export LC_ALL=C # speed up sorting
        samtools view -F 1792 {input} | \
        awk -F"\\t" -v "OFS=\\t" '{{
            split($1,a," ");
            split(a[2],a,":");
            split(a[3],a,",");
            if (a[1] !~ /N/) {{
                if (($5 >= {params.mapping_quality_min}) && ($4 >= 1)) {{
                    print a[1],$3,$4";"$6";"$12";"$13";"$5 
                }} else {{
                    print a[1],"other","NA" 
                }}
            }}
        }}' | sort -k1,1 -k2,2 -k3,3 -S 7G > {output} 2> {log}
        """


# bbmap output bam: M06205:82:000000000-LGYKG:1:1101:16472:10372 XI:Z:AACGACGATCTCTAT,YI:Z:>AAA1AD>1A>BGDG  0       CB_404305_000000        16      42      78=1D59=1X42=1X18=      *       0       0       TGCAAGGATGCAGAGGAAGTTAAGAGGGAAAGTTGCTTTGAGAGGAGGACACTGGGAGGGGTTGGGAGTGGCTCCTGAGGCGGTGATAGGCAGGCAGGCCTGACTTGTCCACAGCTCACCGGAGGCCACCTTGGCAGAACCTGTAGGAAGGGCATGTCTGGCCTCCACACCAGCCCCCTCACTCTTCACCATTTCCCCT AAAAABA?FFFBGFEFCGGAEFHFFHACGCHFH2ADGDHFHC0FEHCEHIGEIC=GHHGHGEFHHGGDHHHHG:GH@IHHGH-G=B>F>(CFH=<EGGGFC=GCF>E>CH6GC=G#<FGC<E=G=<=CCGHIHIGE=#C#=DFGDEFHH#E1110FF//F/A/1/B/00/0E00EB0BB13AFA111B@1BAA>11A>1 NM:i:3  AM:i:42
