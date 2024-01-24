######################################
### Everything before assigning BC ###
######################################


### Create_BAM_umi without demultiplexing ###


rule counts_umi_create_BAM:
    """
    Create a BAM file from FASTQ input, merge FW and REV read and save UMI in XI flag.
    """
    input:
        fw_fastq=lambda wc: getFW(wc.project, wc.condition, wc.replicate, wc.type),
        rev_fastq=lambda wc: getRev(wc.project, wc.condition, wc.replicate, wc.type),
        umi_fastq=lambda wc: getUMI(wc.project, wc.condition, wc.replicate, wc.type),
        script_FastQ2doubleIndexBAM=getScript("count/FastQ2doubleIndexBAM_python3.py"),
        script_MergeTrimReadsBAM=getScript("count/MergeTrimReadsBAM_python3.py"),
    output:
        "results/experiments/{project}/counts/useUMI.{condition}_{replicate}_{type}.bam",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
        umi_length=lambda wc: config["experiments"][wc.project]["umi_length"],
        datasetID="{condition}_{replicate}_{type}",
    conda:
        "../../envs/python3.yaml"
    log:
        temp(
            "results/logs/counts/umi/create_BAM.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        set +o pipefail;

        fwd_length=`zcat {input.fw_fastq} | head -2 | tail -1 | wc -c`;
        fwd_length=$(expr $(($fwd_length-1)));

        rev_start=$(expr $(($fwd_length+1)));

        minoverlap=`echo ${{fwd_length}} ${{fwd_length}} {params.bc_length} | awk '{{print ($1+$2-$3-1 < 11) ? $1+$2-$3-1 : 11}}'`;

        echo $rev_start >> {log}
        echo $minoverlap >> {log}

        paste <( zcat {input.fw_fastq} ) <( zcat {input.rev_fastq}  ) <( zcat {input.umi_fastq} ) | \
        awk '{{if (NR % 4 == 2 || NR % 4 == 0) {{print $1$2$3}} else {{print $1}}}}' | \
        python {input.script_FastQ2doubleIndexBAM} -p -s $rev_start -l 0 -m {params.umi_length} --RG {params.datasetID} | \
        python {input.script_MergeTrimReadsBAM} --FirstReadChimeraFilter '' --adapterFirstRead '' --adapterSecondRead '' -p --mergeoverlap --minoverlap $minoverlap > {output} 2>> {log}
        """


### START COUNTING ####


rule counts_umi_raw_counts:
    """
    Counting BCsxUMIs from the BAM files.
    """
    conda:
        "../../envs/bwa_samtools_picard_htslib.yaml"
    input:
        lambda wc: getUMIBamFile(wc.project, wc.condition, wc.replicate, wc.type),
    output:
        "results/experiments/{project}/counts/useUMI.{condition}_{replicate}_{type}_raw_counts.tsv.gz",
    params:
        umi_length=lambda wc: config["experiments"][wc.project]["umi_length"],
        datasetID="{condition}_{replicate}_{type}",
    log:
        temp(
            "results/logs/counts/umi/raw_counts.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        samtools view -F 1 -r {params.datasetID} {input} | \
        awk -v 'OFS=\\t' '{{ for (i=12; i<=NF; i++) {{
          if ($i ~ /^XJ:Z:/) print $10,substr($i,6,{params.umi_length})
        }}}}' | \
        sort | uniq -c | \
        awk -v 'OFS=\\t' '{{ print $2,$3,$1 }}' | \
        gzip -c > {output} 2> {log}
        """
