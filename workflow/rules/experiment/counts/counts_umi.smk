######################################
### Everything before assigning BC ###
######################################


### Create_BAM_umi without demultiplexing ###


rule experiment_counts_umi_create_BAM:
    """
Create a BAM file from FASTQ input, merge FWD and REV read and save UMI in XI flag.
"""
    input:
        fwd_fastq=lambda wc: getFWD(wc.project, wc.condition, wc.replicate, wc.type, check_splitting=True, check_trimming=True),
        rev_fastq=lambda wc: getREV(wc.project, wc.condition, wc.replicate, wc.type, check_splitting=True, check_trimming=True),
        umi_fastq=lambda wc: getUMI(wc.project, wc.condition, wc.replicate, wc.type, check_splitting=True, check_trimming=True),
        script_FastQ2doubleIndexBAM=getScript("count/FastQ2doubleIndexBAM_python3.py"),
        module_FastQ2doubleIndexBAM=getScript("count/library_python3.py"),
        script_MergeTrimReadsBAM=getScript("count/MergeTrimReadsBAM_python3.py"),
        module_MergeTrimReadsBAM=getScript("count/MergeTrimReads_python3.py"),
    output:
        "results/experiments/{project}/counts/useUMI.{condition}.{replicate}.{type}.{split}.bam",
    log:
        temp("results/logs/experiment/counts/umi/create_BAM.{project}.{condition}.{replicate}.{type}.{split}.log"),
    conda:
        getCondaEnv("python3.yaml")
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
        umi_length=lambda wc: config["experiments"][wc.project]["umi_length"],
        datasetID="{condition}.{replicate}.{type}",
    shell:
        """
        set +o pipefail;

        fwd_length=`zcat {input.fwd_fastq} | head -2 | tail -1 | wc -c`;
        fwd_length=$(expr $(($fwd_length-1)));

        rev_start=$(expr $(($fwd_length+1)));

        minoverlap=`echo ${{fwd_length}} ${{fwd_length}} {params.bc_length} | awk '{{print ($1+$2-$3-1 < 11) ? $1+$2-$3-1 : 11}}'`;

        echo $rev_start >> {log}
        echo $minoverlap >> {log}

        paste <( zcat {input.fwd_fastq} ) <( zcat {input.rev_fastq}  ) <( zcat {input.umi_fastq} ) | \
        awk '{{if (NR % 4 == 2 || NR % 4 == 0) {{print $1$2$3}} else {{print $1}}}}' | \
        python {input.script_FastQ2doubleIndexBAM} -p -s $rev_start -l 0 -m {params.umi_length} --RG {params.datasetID} | \
        python {input.script_MergeTrimReadsBAM} --FirstReadChimeraFilter '' --adapterFirstRead '' --adapterSecondRead '' -p --mergeoverlap --minoverlap $minoverlap > {output} 2>> {log}
        """


rule experiment_counts_umi_attach_idx:
    """
Attach UMI read to FWD/REV headers before NGmerge.
"""
    input:
        read=lambda wc: getExperimentReads(
            wc.read,
            wc.project,
            wc.condition,
            wc.replicate,
            wc.type,
            check_splitting=True,
            check_trimming=True,
        ),
        umi_fastq=lambda wc: getUMI(wc.project, wc.condition, wc.replicate, wc.type, check_splitting=True, check_trimming=True),
        script=getScript("attachBCToFastQ.py"),
        libs=getScript("common.py"),
    output:
        read=temp(
            "results/experiments/{project}/counts/useUMI.{condition}.{replicate}.{type}.{read}.{split}.UMIattached.fastq.gz"
        ),
    log:
        temp("results/logs/experiment/counts/umi/attach_idx.{project}.{condition}.{replicate}.{type}.{read}.{split}.log"),
    conda:
        getCondaEnv("NGmerge.yaml")
    shell:
        """
        python {input.script} -r {input.read} -b {input.umi_fastq} | bgzip -c > {output.read} 2> {log}
        """


use rule experiment_counts_merge_NGmerge_template as experiment_counts_umi_merge_NGmerge with:
    input:
        FWD="results/experiments/{project}/counts/useUMI.{condition}.{replicate}.{type}.FWD.{split}.UMIattached.fastq.gz",
        REV="results/experiments/{project}/counts/useUMI.{condition}.{replicate}.{type}.REV.{split}.UMIattached.fastq.gz",
    output:
        un=temp("results/experiments/{project}/counts/useUMI.{condition}.{replicate}.{type}.{split}.un.NGmerge.fastq.gz"),
        join="results/experiments/{project}/counts/useUMI.{condition}.{replicate}.{type}.{split}.join.NGmerge.fastq.gz",
    log:
        temp("results/logs/experiment/counts/umi/merge_NGmerge.{project}.{condition}.{replicate}.{type}.{split}.log"),


### START COUNTING ####


rule experiment_counts_umi_raw_counts:
    """
Counting BCsxUMIs from the BAM files.
"""
    input:
        lambda wc: (
            expand(
                getUMIBamFile(wc.project),
                project=wc.project,
                condition=wc.condition,
                replicate=wc.replicate,
                type=wc.type,
                split=range(getMaxExperimentSplitNumber()),
            )
            if config["experiments"][wc.project].get("merge_tool", "NGmerge") == "custom"
            else expand(
                "results/experiments/{{project}}/counts/useUMI.{{condition}}.{{replicate}}.{{type}}.{split}.join.NGmerge.fastq.gz",
                split=range(getMaxExperimentSplitNumber()),
            )
        ),
    output:
        "results/experiments/{project}/counts/useUMI.{condition}.{replicate}.{type}.raw_counts.tsv.gz",
    log:
        temp("results/logs/experiment/counts/umi/raw_counts.{project}.{condition}.{replicate}.{type}.log"),
    conda:
        getCondaEnv("bwa_samtools_picard_htslib.yaml")
    params:
        umi_length=lambda wc: config["experiments"][wc.project]["umi_length"],
        datasetID="{condition}.{replicate}.{type}",
        merge_tool=lambda wc: config["experiments"][wc.project].get("merge_tool", "NGmerge"),
    shell:
        """
        if [[ "{params.merge_tool}" == "custom" ]]; then
            samtools merge -c -o - {input} | samtools view -F 1 -r {params.datasetID} | \
            awk -v 'OFS=\\t' '{{ for (i=12; i<=NF; i++) {{
              if ($i ~ /^XJ:Z:/) print $10,substr($i,6,{params.umi_length})
            }}}}' | \
            sort | uniq -c | \
            awk -v 'OFS=\\t' '{{ print $2,$3,$1 }}' | \
            gzip -c > {output} 2> {log}
        else
            zcat {input} | \
            awk -v 'OFS=\\t' -v umi_len={params.umi_length} '
                NR%4==1 {{
                    umi="";
                    if (match($0, /XI:Z:([^,[:space:]]+)/, m)) umi=substr(m[1], 1, umi_len)
                }}
                NR%4==2 {{if (umi != "") print $1, umi}}
            ' | \
            sort | uniq -c | \
            awk -v 'OFS=\\t' '{{ print $2,$3,$1 }}' | \
            gzip -c > {output} 2> {log}
        fi
        """
