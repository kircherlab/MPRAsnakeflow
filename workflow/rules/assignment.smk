SPLIT_FILES_NUMBER = 100


include: "assignment/statistic.smk"


rule assignment_getInputs:
    conda:
        "../envs/default.yaml"
    input:
        lambda wc: config["assignments"][wc.assignment][wc.R],
    output:
        R1=temp("results/assignment/{assignment}/fastq/{R}.fastq.gz"),
    log:
        temp("results/logs/assignment/getInputs.{assignment}.{R}.log"),
    shell:
        """
        if [[ "$(ls {input} | wc -l)" -eq 1 ]]; then 
            ln -s {input} {output}; 
        else
            zcat {input} | gzip -c > {output};
        fi &> {log}
        """


rule assignment_fastq_split:
    input:
        "results/assignment/{assignment}/fastq/{R}.fastq.gz",
    output:
        temp(
            expand(
                "results/assignment/{{assignment}}/fastq/splits/{{R}}.split{split}.fastq.gz",
                split=range(0, getSplitNumber()),
            ),
        ),
    conda:
        "../envs/fastqsplitter.yaml"
    log:
        temp("results/logs/assignment/fastq_split.{assignment}.{R}.log"),
    params:
        files=lambda wc: " ".join(
            [
                "-o %s" % i
                for i in expand(
                    "results/assignment/{assignment}/fastq/splits/{R}.split{split}.fastq.gz",
                    assignment=wc.assignment,
                    R=wc.R,
                    split=range(0, getSplitNumber()),
                )
            ]
        ),
    shell:
        """
        fastqsplitter -i {input} -t 10 {params.files} &> {log}
        """


rule assignment_merge:
    input:
        R1="results/assignment/{assignment}/fastq/splits/R1.split{split}.fastq.gz",
        R2="results/assignment/{assignment}/fastq/splits/R2.split{split}.fastq.gz",
        R3="results/assignment/{assignment}/fastq/splits/R3.split{split}.fastq.gz",
        script_FastQ2doubleIndexBAM=getScript("count/FastQ2doubleIndexBAM.py"),
        script_MergeTrimReadsBAM=getScript("count/MergeTrimReadsBAM.py"),
    output:
        bam=temp("results/assignment/{assignment}/bam/merge_split{split}.bam"),
    conda:
        "../envs/python27.yaml"
    log:
        temp("results/logs/assignment/merge.{assignment}.{split}.log"),
    shell:
        """
        paste <( zcat {input.R1} ) <( zcat {input.R2} ) <( zcat {input.R3} ) | \
        awk '{{ 
            count+=1; 
            if ((count == 1) || (count == 3)) {{ 
                print $1 
            }} else {{
                print $1$2$3 
            }}; 
            if (count == 4) {{
                count=0 
            }} 
        }}' | \
        python {input.script_FastQ2doubleIndexBAM} -p -s 162 -l 15 -m 0 | \
        python {input.script_MergeTrimReadsBAM} -c '' -f CATTGCGTGAACCGACAATTCGTCGAGGGACCTAATAAC -s AGTTGATCCGGTCCTAGGTCTAGAGCGGGCCCTGGCAGA --mergeoverlap -p > {output} 2> {log}
        """


assignment_bwa_dicts = ["bwt", "sa", "pac", "ann", "amb"]


rule assignment_bwa_ref:
    input:
        lambda wc: config["assignments"][wc.assignment]["reference"],
    output:
        ref="results/assignment/{assignment}/reference/reference.fa",
        bwa=expand(
            "results/assignment/{{assignment}}/reference/reference.fa.{ext}",
            ext=["fai"] + assignment_bwa_dicts,
        ),
        d="results/assignment/{assignment}/reference/reference.fa.dict",
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    log:
        temp("results/logs/assignment/bwa_ref.{assignment}.log"),
    shell:
        """
        cat {input} | awk '{{gsub(/[\\]\\[]/,"_")}}$0' > {output.ref};
        bwa index -a bwtsw {output.ref} &> {log};
        samtools faidx {output.ref} &>> {log};
        picard CreateSequenceDictionary REFERENCE={output.ref} OUTPUT={output.d} &>> {log}
        """


rule assignment_mapping:
    input:
        bams=expand(
            "results/assignment/{{assignment}}/bam/merge_split{split}.bam",
            split=range(0, getSplitNumber()),
        ),
        reference="results/assignment/{assignment}/reference/reference.fa",
        bwa_index=expand(
            "results/assignment/{{assignment}}/reference/reference.fa.{ext}",
            ext=["fai", "dict"] + assignment_bwa_dicts,
        ),
    output:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    log:
        temp("results/logs/assignment/mapping.{assignment}.log"),
    shell:
        """
        bwa mem -t 30 -L 80 -M -C {input.reference} <(
            samtools cat {input.bams} | \
            samtools view -F 514 | \
            awk 'BEGIN{{ OFS="\\n"; FS="\\t" }}{{ print "@"$1" "$12","$13,$10,"+",$11 }}';
        ) | \
        samtools view -Su | \
        samtools sort > {output} 2> {log}
        """


rule assignment_idx_bam:
    input:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    output:
        "results/assignment/{assignment}/aligned_merged_reads.bam.bai",
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    log:
        "results/logs/assignment/{assignment}/assignment_idx_bam.log",
    shell:
        """
        samtools index {input} 2> {log}
        """


rule assignment_flagstat:
    input:
        bam="results/assignment/{assignment}/aligned_merged_reads.bam",
        idx="results/assignment/{assignment}/aligned_merged_reads.bam.bai",
    output:
        "results/assignment/{assignment}/stats/assignment/bam_stats.txt",
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    log:
        temp("results/logs/assignment/flagstat.{assignment}.log"),
    shell:
        """
        samtools flagstat {input.bam} > {output} 2> {log}
        """


rule assignment_getBCs:
    input:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    output:
        "results/assignment/{assignment}/barcodes_incl_other.sorted.tsv.gz",
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    params:
        alignment_start_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_start"
        ]["min"],
        alignment_start_max=lambda wc: config["assignments"][wc.assignment][
            "alignment_start"
        ]["max"],
        sequence_length_min=lambda wc: config["assignments"][wc.assignment][
            "sequence_length"
        ]["min"],
        sequence_length_max=lambda wc: config["assignments"][wc.assignment][
            "sequence_length"
        ]["max"],
    log:
        temp("results/logs/assignment/getBCs.{assignment}.log"),
    shell:
        """
        samtools view -F 1792 {input} | \
        awk -v "OFS=\\t" '{{
            split($(NF),a,":");
            split(a[3],a,",");
            if (a[1] !~ /N/) {{
                if (($5 > 0) && ($4 >= {params.alignment_start_min}) && ($4 <= {params.alignment_start_max}) && (length($10) >= {params.sequence_length_min}) && (length($10) <= {params.sequence_length_max})) {{
                    print a[1],$3,$4";"$6";"$12";"$13";"$5 
                }} else {{
                    print a[1],"other","NA" 
                }}
            }}
        }}' | sort -k1,1 -k2,2 -k3,3 | gzip -c > {output} 2> {log}
        """


rule assignment_filter:
    input:
        assignment="results/assignment/{assignment}/barcodes_incl_other.sorted.tsv.gz",
        script=getScript("assignment/filterAssignmentTsv.py"),
    output:
        "results/assignment/{assignment}/assignment_barcodes.{assignment_config}.sorted.tsv.gz",
    conda:
        "../envs/python3.yaml"
    log:
        temp("results/logs/assignment/filter.{assignment}.{assignment_config}.log"),
    params:
        min_support=lambda wc: config["assignments"][wc.assignment]["configs"][
            wc.assignment_config
        ]["min_support"],
        fraction=lambda wc: config["assignments"][wc.assignment]["configs"][
            wc.assignment_config
        ]["fraction"],
        unknown_other=lambda wc: "-o"
        if config["assignments"][wc.assignment]["configs"][wc.assignment_config][
            "unknown_other"
        ]
        else "",
        ambiguous=lambda wc: "-a"
        if config["assignments"][wc.assignment]["configs"][wc.assignment_config][
            "ambiguous"
        ]
        else "",
    shell:
        """
        zcat  {input.assignment} | \
        python {input.script} \
        -m {params.min_support} -f {params.fraction} {params.unknown_other} {params.ambiguous} | \
        gzip -c > {output} 2> {log}
        """
