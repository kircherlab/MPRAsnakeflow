"""
Assignment workflow

This workflow will asssign barcodes to the designed reference/insert/oligos.
The output is a tabular file that matched barcodes with oligos.
"""


include: "assignment/common.smk"
include: "assignment/hybridFWRead.smk"
include: "assignment/statistic.smk"


rule assignment_check_design:
    """
    Check if the design file is correct and no duplicated sequences are present (FW and reverse).
    Also check if no duplicated headers and no illegal characters in header.
    """
    conda:
        "../envs/python3.yaml"
    input:
        design=lambda wc: config["assignments"][wc.assignment]["design_file"],
        script=getScript("assignment/check_design_file.py"),
    output:
        temp("results/assignment/{assignment}/reference/reference.fa.fxi"),
        touch("results/assignment/{assignment}/design_check.done"),
        ref="results/assignment/{assignment}/reference/reference.fa",
    params:
        start=lambda wc: (
            config["assignments"][wc.assignment]["alignment_tool"]["configs"][
                "alignment_start"
            ]
            if config["assignments"][wc.assignment]["alignment_tool"]["tool"]
            == "exact"
            else config["assignments"][wc.assignment]["alignment_tool"]["configs"][
                "alignment_start"
            ]["max"]
        ),
        length=lambda wc: (
            config["assignments"][wc.assignment]["alignment_tool"]["configs"][
                "sequence_length"
            ]
            if config["assignments"][wc.assignment]["alignment_tool"]["tool"]
            == "exact"
            else config["assignments"][wc.assignment]["alignment_tool"]["configs"][
                "sequence_length"
            ]["min"]
        ),
        fast_check=lambda wc: (
            "--fast-dict"
            if config["assignments"][wc.assignment]["design_check"]["fast"]
            or config["assignments"][wc.assignment]["alignment_tool"]["tool"]
            == "bbmap"
            else "--slow-string-search"
        ),
        check_sequence_collitions=lambda wc: (
            "--perform-sequence-check"
            if config["assignments"][wc.assignment]["design_check"][
                "sequence_collitions"
            ]
            else "--skip-sequence-check"
        ),
    log:
        log=temp("results/logs/assignment/check_design.{assignment}.log"),
        err="results/assignment/{assignment}/design_check.err",
    shell:
        """
        trap "cat {log.err}" ERR
        cp {input.design} {output.ref}
        python {input.script} --input {output.ref} \
        --start {params.start} --length {params.length} \
        {params.fast_check} {params.check_sequence_collitions} > {log.log} 2> {log.err};
        """


rule assignment_fastq_split:
    """
    Split the fastq files into n files for parallelisation.
    n is given by split_read in the configuration file.

    Runs only if the design file is correct.
    """
    input:
        fastq=lambda wc: getAssignmentRead(wc.assignment, wc.read),
        check="results/assignment/{assignment}/design_check.done",
    output:
        temp(
            expand(
                "results/assignment/{{assignment}}/fastq/splits/{{read}}.split{split}.fastq.gz",
                split=range(0, getSplitNumber()),
            ),
        ),
    conda:
        "../envs/fastqsplitter.yaml"
    log:
        temp("results/logs/assignment/fastq_split.{assignment}.{read}.log"),
    params:
        files=lambda wc: " ".join(
            [
                "-o %s" % i
                for i in expand(
                    "results/assignment/{assignment}/fastq/splits/{read}.split{split}.fastq.gz",
                    assignment=wc.assignment,
                    read=wc.read,
                    split=range(0, getSplitNumber()),
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
    conda:
        "../envs/NGmerge.yaml"
    input:
        read="results/assignment/{assignment}/fastq/splits/{read}.split{split}.fastq.gz",
        BC="results/assignment/{assignment}/fastq/splits/BC.split{split}.fastq.gz",
        script=getScript("attachBCToFastQ.py"),
        libs=getScript("common.py"),
    output:
        read=temp(
            "results/assignment/{assignment}/fastq/splits/{read}.split{split}.BCattached.fastq.gz"
        ),
    params:
        BC_rev_comp=lambda wc: (
            "--reverse-complement"
            if config["assignments"][wc.assignment]["BC_rev_comp"]
            else ""
        ),
    log:
        temp("results/logs/assignment/attach_idx.{assignment}.{split}.{read}.log"),
    shell:
        """
        python {input.script} -r {input.read} -b {input.BC} {params.BC_rev_comp} | bgzip -c > {output.read} 2> {log}
        """


rule assignment_merge:
    """
    Merge the FW,REV and BC fastq files into one.
    Extract the index sequence and add it to the header.
    """
    conda:
        "../envs/NGmerge.yaml"
    input:
        FW="results/assignment/{assignment}/fastq/splits/FW.split{split}.BCattached.fastq.gz",
        REV="results/assignment/{assignment}/fastq/splits/REV.split{split}.BCattached.fastq.gz",
    output:
        un=temp("results/assignment/{assignment}/fastq/merge_split{split}.un.fastq.gz"),
        join=temp(
            "results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz"
        ),
    params:
        min_overlap=lambda wc: config["assignments"][wc.assignment]["NGmerge"][
            "min_overlap"
        ],
        frac_mismatches_allowed=lambda wc: config["assignments"][wc.assignment][
            "NGmerge"
        ]["frac_mismatches_allowed"],
        min_dovetailed_overlap=lambda wc: config["assignments"][wc.assignment][
            "NGmerge"
        ]["min_dovetailed_overlap"],
    log:
        temp("results/logs/assignment/merge.{assignment}.{split}.log.gz"),
    shell:
        """
        NGmerge \
        -1 {input.FW} \
        -2 {input.REV} \
        -m {params.min_overlap} -p {params.frac_mismatches_allowed} \
        -d \
        -e {params.min_dovetailed_overlap} \
        -z \
        -o  {output.join} \
        -i -f {output.un} \
        -l >(gzip -c - > {log})
        """


include: "assignment/mapping_exact.smk"
include: "assignment/mapping_bwa.smk"
include: "assignment/mapping_bbmap.smk"


rule assignment_collectBCs:
    """
    Get the barcodes.
    """
    input:
        lambda wc: expand(
            "results/assignment/{{assignment}}/BCs/barcodes_{mapper}.{split}.tsv",
            split=range(0, getSplitNumber()),
            mapper=config["assignments"][wc.assignment]["alignment_tool"]["tool"],
        ),
    output:
        "results/assignment/{assignment}/barcodes_incl_other.tsv.gz",
    params:
        batch_size="--batch-size=%d" % getSplitNumber() if getSplitNumber() > 1 else "",
    conda:
        "../envs/default.yaml"
    log:
        temp("results/logs/assignment/collectBCs.{assignment}.log"),
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
        ambigous="results/assignment/{assignment}/assignment_barcodes_with_ambigous.{assignment_config}.tsv.gz",
    conda:
        "../envs/python3.yaml"
    log:
        log=temp("results/logs/assignment/filter.{assignment}.{assignment_config}.log"),
        err=temp("results/logs/assignment/filter.{assignment}.{assignment_config}.err"),
    params:
        min_support=lambda wc: config["assignments"][wc.assignment]["configs"][
            wc.assignment_config
        ]["min_support"],
        fraction=lambda wc: config["assignments"][wc.assignment]["configs"][
            wc.assignment_config
        ]["fraction"],
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
        tee >(gzip -c > {output.ambigous}) | \
        awk -v "OFS=\\t"  -F"\\t" '{{ if (($2 != \"ambiguous\") && ($2 != \"other\")) {{ print $0 }} }}' | \
        gzip -c > {output.final} 2> {log.err};
        gzip -l {output.final} | awk 'NR==2 {{exit($2==0)}}' || {{ echo "Error: Empty barcode file {output.final}. No barcodes detected!" >> {log.err}; exit 1; }}
        """
