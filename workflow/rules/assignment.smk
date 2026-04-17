"""
Assignment workflow

This workflow will asssign barcodes to the designed reference/insert/oligos.
The output is a tabular file that matched barcodes with oligos.
"""


include: "assignment/common.smk"
include: "assignment/preprocessing.smk"
include: "assignment/hybridFWDRead.smk"
include: "assignment/statistic.smk"


rule assignment_check_design:
    """
Check if the design file is correct and no duplicated sequences are present (FWD and reverse).
Also check if no duplicated headers and no illegal characters in header.
"""
    input:
        design=lambda wc: config["assignments"][wc.assignment]["design_file"],
        script=getScript("assignment/check_design_file.py"),
    output:
        temp("results/assignment/{assignment}/reference/reference.fa.fxi"),
        touch("results/assignment/{assignment}/design_check.done"),
        ref="results/assignment/{assignment}/reference/reference.fa",
        ref_tmp=temp("results/assignment/{assignment}/reference/reference.tmp.fa"),
    log:
        log=temp("results/logs/assignment/check_design.{assignment}.log"),
        err="results/assignment/{assignment}/design_check.err",
    conda:
        getCondaEnv("python3.yaml")
    params:
        start=lambda wc: (
            config["assignments"][wc.assignment]["alignment_tool"]["configs"]["alignment_start"]["max"]
            if config["assignments"][wc.assignment]["alignment_tool"]["tool"]
            in [
                "bwa",
                "bwa-additional-filtering",
            ]
            else config["assignments"][wc.assignment]["alignment_tool"]["configs"]["alignment_start"]
        ),
        length=lambda wc: (
            config["assignments"][wc.assignment]["alignment_tool"]["configs"]["sequence_length"]["min"]
            if config["assignments"][wc.assignment]["alignment_tool"]["tool"]
            in [
                "bwa",
                "bwa-additional-filtering",
            ]
            else config["assignments"][wc.assignment]["alignment_tool"]["configs"]["sequence_length"]
        ),
        fast_check=lambda wc: (
            "--fast-dict" if config["assignments"][wc.assignment]["design_check"]["fast"] else "--slow-string-search"
        ),
        sequence_collisions=lambda wc: (
            "sense_antisense" if config["assignments"][wc.assignment]["design_check"]["sequence_collisions"] else "skip"
        ),
        attach_sequence=lambda wc: (
            "--attach-sequence %s %s"
            % (
                config["assignments"][wc.assignment]["strand_sensitive"]["forward_adapter"],
                config["assignments"][wc.assignment]["strand_sensitive"]["reverse_adapter"],
            )
            if config["assignments"][wc.assignment]["strand_sensitive"]["enable"]
            else ""
        ),
    shell:
        """
        trap "cat {log.err}" ERR
        cp {input.design} {output.ref_tmp}
        python {input.script} --input {output.ref_tmp} \
        --output {output.ref} \
        --start {params.start} --length {params.length} \
        {params.fast_check} --sequence-check {params.sequence_collisions} \
        {params.attach_sequence} > {log.log} 2> {log.err};
        """


rule assignment_fastq_split:
    """
Split the fastq files into n files for parallelisation.
n is given by split_read in the configuration file.

Runs only if the design file is correct.
"""
    input:
        fastq=lambda wc: getAssignmentRead(
            wc.assignment,
            wc.read,
        ),
        check="results/assignment/{assignment}/design_check.done",
    output:
        temp(
            expand(
                "results/assignment/{{assignment}}/fastq/splits/{{read}}.split{split}.fastq.gz",
                split=range(
                    0,
                    getAssignmentSplitNumber(),
                ),
            ),
        ),
    log:
        temp("results/logs/assignment/fastq_split.{assignment}.{read}.log"),
    conda:
        getCondaEnv("fastqsplitter.yaml")
    params:
        files=lambda wc: " ".join(
            [
                "-o %s" % i
                for i in expand(
                    "results/assignment/{assignment}/fastq/splits/{read}.split{split}.fastq.gz",
                    assignment=wc.assignment,
                    read=wc.read,
                    split=range(
                        0,
                        getAssignmentSplitNumber(),
                    ),
                )
            ]
        ),
    shell:
        """
        fastqsplitter -i <(zcat {input.fastq}) -t 1 {params.files} &> {log}
        """


rule assignment_attach_idx:
    """
Extract the index sequence and add it to the header.
"""
    input:
        read="results/assignment/{assignment}/fastq/splits/{read}.split{split}.fastq.gz",
        BC="results/assignment/{assignment}/fastq/splits/BC.split{split}.fastq.gz",
        script=getScript("attachBCToFastQ.py"),
        libs=getScript("common.py"),
    output:
        read=temp("results/assignment/{assignment}/fastq/splits/{read}.split{split}.BCattached.fastq.gz"),
    log:
        temp("results/logs/assignment/attach_idx.{assignment}.{split}.{read}.log"),
    conda:
        getCondaEnv("NGmerge.yaml")
    params:
        BC_rev_comp=lambda wc: ("--reverse-complement" if config["assignments"][wc.assignment]["BC_rev_comp"] else ""),
        attach_sequence=lambda wc: (
            "--attach-sequence left %s"
            % (
                config["assignments"][wc.assignment]["strand_sensitive"]["forward_adapter"]
                if wc.read == "FWD"
                else reverse_complement(config["assignments"][wc.assignment]["strand_sensitive"]["reverse_adapter"])
            )
            if config["assignments"][wc.assignment]["strand_sensitive"]["enable"]
            else ""
        ),
    shell:
        """
        python {input.script} -r {input.read} -b {input.BC} {params.BC_rev_comp} {params.attach_sequence} | bgzip -c > {output.read} 2> {log}
        """


rule assignment_merge_NGmerge:
    """
Merge the FWD, REV and BC fastq files into one using NGmerge.
"""
    input:
        FWD="results/assignment/{assignment}/fastq/splits/FWD.split{split}.BCattached.fastq.gz",
        REV="results/assignment/{assignment}/fastq/splits/REV.split{split}.BCattached.fastq.gz",
    output:
        un=temp("results/assignment/{assignment}/fastq/merge_split{split}.un.NGmerge.fastq.gz"),
        join=temp("results/assignment/{assignment}/fastq/merge_split{split}.join.NGmerge.fastq.gz"),
    log:
        "results/logs/assignment/merge_NGmerge.{assignment}.{split}.log",
    conda:
        getCondaEnv("NGmerge.yaml")
    params:
        min_overlap=lambda wc: config["assignments"][wc.assignment]["NGmerge"]["min_overlap"],
        frac_mismatches_allowed=lambda wc: config["assignments"][wc.assignment]["NGmerge"]["frac_mismatches_allowed"],
        min_dovetailed_overlap=lambda wc: config["assignments"][wc.assignment]["NGmerge"]["min_dovetailed_overlap"],
    shell:
        """
        NGmerge \
        -1 {input.FWD} \
        -2 {input.REV} \
        -m {params.min_overlap} -p {params.frac_mismatches_allowed} \
        -d \
        -e {params.min_dovetailed_overlap} \
        -z \
        -o  {output.join} \
        -i -f {output.un} &> {log}
        """


rule assignment_merge_fastqjoin:
    """
Merge the FWD, REV and BC fastq files into one using fastq-join.
"""
    input:
        FWD="results/assignment/{assignment}/fastq/splits/FWD.split{split}.BCattached.fastq.gz",
        REV="results/assignment/{assignment}/fastq/splits/REV.split{split}.BCattached.fastq.gz",
    output:
        un1=temp("results/assignment/{assignment}/fastq/merge_split{split}.un1.fastqjoin.fastq.gz"),
        un2=temp("results/assignment/{assignment}/fastq/merge_split{split}.un2.fastqjoin.fastq.gz"),
        join=temp("results/assignment/{assignment}/fastq/merge_split{split}.join.fastqjoin.fastq.gz"),
    log:
        "results/logs/assignment/merge_fastqjoin.{assignment}.{split}.log",
    conda:
        getCondaEnv("fastq-join.yaml")
    params:
        min_overlap=lambda wc: config["assignments"][wc.assignment]["fastq-join"]["min_overlap"],
        max_pct_mismatch=lambda wc: config["assignments"][wc.assignment]["fastq-join"]["max_pct_mismatch"],
    shell:
        """
        fastq-join -p {params.min_overlap} -m {params.max_pct_mismatch} {input.FWD} {input.REV} \
        -o {output.un1} -o {output.un2} -o {output.join} &> {log}
        """


include: "assignment/mapping_exact.smk"
include: "assignment/mapping_bwa.smk"
include: "assignment/mapping_bbmap.smk"
include: "assignment/mapping_pbmm2.smk"


rule assignment_collectBCs:
    """
Get the barcodes.
"""
    input:
        lambda wc: (
            "results/assignment/{assignment}/BCs/barcodes.pbmm2.tsv"
            if config["assignments"][wc.assignment]["alignment_tool"]["tool"] == "pbmm2"
            else expand(
                "results/assignment/{{assignment}}/BCs/barcodes.{mapper}.{split}.tsv",
                split=range(
                    0,
                    getAssignmentSplitNumber(),
                ),
                mapper=config["assignments"][wc.assignment]["alignment_tool"]["tool"],
            )
        ),
    output:
        "results/assignment/{assignment}/barcodes_incl_other.tsv.gz",
    log:
        temp("results/logs/assignment/collectBCs.{assignment}.log"),
    conda:
        getCondaEnv("default.yaml")
    params:
        batch_size=("--batch-size=%d" % getAssignmentSplitNumber() if getAssignmentSplitNumber() > 1 else ""),
    shell:
        """
        export LC_ALL=C # speed up sort
        sort -S 7G {params.batch_size} --parallel={threads} -k1,1 -k2,2 -k3,3 -m {input} | \
        gzip -c > {output} 2> {log}
        """


rule assignment_filter:
    """
Filter the barcodes file based on the config given in the config-file.
FIXME: Limitation is that oligos cannot have a name ambiguous or other.
"""
    input:
        assignment="results/assignment/{assignment}/barcodes_incl_other.tsv.gz",
        script=getScript("assignment/filterAssignmentTsv.py"),
    output:
        final="results/assignment/{assignment}/assignment_barcodes.{assignment_config}.tsv.gz",
        ambiguous="results/assignment/{assignment}/assignment_barcodes_with_ambiguous.{assignment_config}.tsv.gz",
    log:
        log=temp("results/logs/assignment/filter.{assignment}.{assignment_config}.log"),
        err=temp("results/logs/assignment/filter.{assignment}.{assignment_config}.err"),
    conda:
        getCondaEnv("python3.yaml")
    params:
        min_support=lambda wc: config["assignments"][wc.assignment]["configs"][wc.assignment_config]["min_support"],
        fraction=lambda wc: config["assignments"][wc.assignment]["configs"][wc.assignment_config]["fraction"],
        unknown_other="-o",
        ambiguous="-a",
        bc_length=lambda wc: config["assignments"][wc.assignment]["bc_length"],
    shell:
        """
        trap "cat {log.err}" ERR
        zcat  {input.assignment} | \
        awk -v "OFS=\\t" -F"\\t" '{{if (length($1)=={params.bc_length}){{print $0 }}}}' | \
        python {input.script} \
        -m {params.min_support} -f {params.fraction} {params.unknown_other} {params.ambiguous} | \
        tee >(gzip -c > {output.ambiguous}) | \
        awk -v "OFS=\\t"  -F"\\t" '{{ if (($2 != \"ambiguous\") && ($2 != \"other\")) {{ print $1,$2 }} }}' | \
        gzip -c > {output.final} 2> {log.err};
        gzip -l {output.final} | awk 'NR==2 {{exit($2==0)}}' || {{ echo "Error: Empty barcode file {output.final}. No barcodes detected!" >> {log.err}; exit 1; }}
        """
