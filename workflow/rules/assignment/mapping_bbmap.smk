# rule assignment_bwa_ref:
#     """
#     Create mapping reference for BWA from design file.
#     """
#     conda:
#         "../../envs/bbmap_samtools_picard_htslib.yaml"
#     input:
#         ref="results/assignment/{assignment}/reference/reference.fa",
#         check="results/assignment/{assignment}/design_check.done",
#     output:
#         bwa=expand(
#             "results/assignment/{{assignment}}/reference/reference.fa.{ext}",
#             ext=["fai"] + assignment_bwa_dicts,
#         ),
#         d="results/assignment/{assignment}/reference/reference.fa.dict",
#     log:
#         temp("results/logs/assignment/bwa_ref.{assignment}.log"),
#     shell:
#         """
#         bwa index -a bwtsw {input.ref} &> {log};
#         samtools faidx {input.ref} &>> {log};
#         picard CreateSequenceDictionary REFERENCE={input.ref} OUTPUT={output.d} &>> {log}
#         """


rule assignment_mapping_bbmap:
    """
    Map the reads to the reference and sort unsing bwa mem
    """
    conda:
        "../../envs/bbmap_samtools_picard_htslib.yaml"
    input:
        reads="results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz",
        reference="results/assignment/{assignment}/reference/reference.fa",
    params:
        reference_path = lambda wc, input: os.path.dirname(input.reference),
        bam = lambda wc: "results/assignment/{assignment}/bam/merge_split{split}.unsorted.bam".format(
            assignment=wc.assignment, split=wc.split
        ),
    output:
        sorted_bam=temp("results/assignment/{assignment}/bam/merge_split{split}.mapped.bam"),
    log:
        temp("results/logs/assignment/mapping.{assignment}.{split}.log"),
    shell:
        """
        bbmap.sh in={input.reads} ref={input.reference} -t={threads} path={params.reference_path} out={params.bam} 2> {log};
        samtools sort -l 0 -@ {threads} {params.bam} > {output.sorted_bam} 2>> {log};
        rm -f {params.bam};
        """


rule assignment_getBCs:
    """
    Get the barcodes.
    """
    conda:
        "../../envs/bbmap_samtools_picard_htslib.yaml"
    input:
        "results/assignment/{assignment}/bam/merge_split{split}.mapped.bam",
    output:
        temp("results/assignment/{assignment}/BCs/barcodes_incl_other.{split}.tsv"),
    params:
        alignment_start_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["alignment_start"]["min"],
        alignment_start_max=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["alignment_start"]["max"],
        sequence_length_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["sequence_length"]["min"],
        sequence_length_max=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["sequence_length"]["max"],
        mapping_quality_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["min_mapping_quality"],
    log:
        temp("results/logs/assignment/getBCs.{assignment}.{split}.log"),
    shell:
        """
        export LC_ALL=C # speed up sorting
        samtools view -F 1792 {input} | \
        awk -v "OFS=\\t" '{{
            split($(NF),a,":");
            split(a[3],a,",");
            if (a[1] !~ /N/) {{
                if (($5 >= {params.mapping_quality_min}) && ($4 >= {params.alignment_start_min}) && ($4 <= {params.alignment_start_max}) && (length($10) >= {params.sequence_length_min}) && (length($10) <= {params.sequence_length_max})) {{
                    print a[1],$3,$4";"$6";"$12";"$13";"$5
                }} else {{
                    print a[1],"other","NA"
                }}
            }}
        }}' | sort -k1,1 -k2,2 -k3,3 -S 7G > {output} 2> {log}
        """


rule assignment_getBCs_all_bbmap:
    """
    Get the barcodes.
    """
    conda:
        "../../envs/bbmap_samtools_picard_htslib.yaml"
    input:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    output:
        temp("results/assignment/{assignment}/barcodes_incl_other.tsv.gz"),
    params:
        alignment_start_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["alignment_start"]["min"],
        alignment_start_max=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["alignment_start"]["max"],
        sequence_length_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["sequence_length"]["min"],
        sequence_length_max=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["sequence_length"]["max"],
        mapping_quality_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["min_mapping_quality"],
    log:
        temp("results/logs/assignment/getBCs.{assignment}.log"),
    shell:
        """
        export LC_ALL=C # speed up sorting
        samtools view -F 1792 {input} | \
        awk -v "OFS=\\t" -F"\t" '{{
            read_name = $1;
            bc_string = $2;
            ref_name = $3;
            alignement_start = $4;
            mapping_quality = $5;
            cigar = $6;
            aligned_sequence = $10;
            edit_distance = $12;
            alignment_score = $13;
            split(read_name,barcode_string_start,"XI:Z:"); # get barcode information
            split(barcode_string_start[2], barcode,",");
            if (barcode[1] !~ /N/) {{
                if (($mapping_quality >= {params.mapping_quality_min}) && ($alignement_start >= {params.alignment_start_min}) && ($alignement_start <= {params.alignment_start_max}) && (length($aligned_sequence) >= {params.sequence_length_min}) && (length($aligned_sequence) <= {params.sequence_length_max})) {{
                    print barcode[1],$ref_name,$alignement_start";"$cigar";"$edit_distance";"$alignment_score";"$mapping_quality
                }} else {{
                    print barcode[1],"other","NA"
                }}
            }}
        }}' | sort -k1,1 -k2,2 -k3,3 -S 7G | gzip -c > {output} 2> {log};
        """
# bbmap output bam: M06205:82:000000000-LGYKG:1:1101:16472:10372 XI:Z:AACGACGATCTCTAT,YI:Z:>AAA1AD>1A>BGDG  0       CB_404305_000000        16      42      78=1D59=1X42=1X18=      *       0       0       TGCAAGGATGCAGAGGAAGTTAAGAGGGAAAGTTGCTTTGAGAGGAGGACACTGGGAGGGGTTGGGAGTGGCTCCTGAGGCGGTGATAGGCAGGCAGGCCTGACTTGTCCACAGCTCACCGGAGGCCACCTTGGCAGAACCTGTAGGAAGGGCATGTCTGGCCTCCACACCAGCCCCCTCACTCTTCACCATTTCCCCT AAAAABA?FFFBGFEFCGGAEFHFFHACGCHFH2ADGDHFHC0FEHCEHIGEIC=GHHGHGEFHHGGDHHHHG:GH@IHHGH-G=B>F>(CFH=<EGGGFC=GCF>E>CH6GC=G#<FGC<E=G=<=CCGHIHIGE=#C#=DFGDEFHH#E1110FF//F/A/1/B/00/0E00EB0BB13AFA111B@1BAA>11A>1 NM:i:3  AM:i:42

rule assignment_collect:
    """
    Collect mapped reads.
    """
    conda:
        "../../envs/bbmap_samtools_picard_htslib.yaml"
    input:
        bams=expand(
            "results/assignment/{{assignment}}/bam/merge_split{split}.mapped.bam",
            split=range(0, getSplitNumber()),
        ),
    output:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    log:
        temp("results/logs/assignment/collect.{assignment}.log"),
    shell:
        """
        samtools merge -@ {threads} {output} {input.bams} 2> {log}
        """


rule assignment_idx_bam:
    """
    Index the BAM file
    """
    conda:
        "../../envs/bbmap_samtools_picard_htslib.yaml"
    input:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    output:
        "results/assignment/{assignment}/aligned_merged_reads.bam.bai",
    log:
        "results/logs/assignment/{assignment}/assignment_idx_bam.log",
    shell:
        """
        samtools index {input} 2> {log}
        """


rule assignment_flagstat:
    """
    Run samtools flagstat
    """
    conda:
        "../../envs/bbmap_samtools_picard_htslib.yaml"
    input:
        bam="results/assignment/{assignment}/aligned_merged_reads.bam",
        idx="results/assignment/{assignment}/aligned_merged_reads.bam.bai",
    output:
        "results/assignment/{assignment}/statistic/assignment/bam_stats.txt",
    log:
        temp("results/logs/assignment/flagstat.{assignment}.log"),
    shell:
        """
        samtools flagstat {input.bam} > {output} 2> {log}
        """
