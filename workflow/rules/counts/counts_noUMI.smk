######################################
### Everything before assigning BC ###
######################################


### Create_BAM_umi without demultiplexing ###


rule counts_noUMI_create_BAM:
    """
    Create a BAM file from FASTQ input, merge FW and REV read and save UMI in XI flag.
    """
    input:
        fw_fastq=lambda wc: getFW(wc.project, wc.condition, wc.replicate, wc.type),
        rev_fastq=lambda wc: getRev(wc.project, wc.condition, wc.replicate, wc.type),
        script_FastQ2doubleIndexBAM=getScript("count/FastQ2doubleIndexBAM.py"),
        script_MergeTrimReadsBAM=getScript("count/MergeTrimReadsBAM.py"),
    output:
        "results/experiments/{project}/counts/noUMI.{condition}_{replicate}_{type}.bam",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
        datasetID="{condition}_{replicate}_{type}",
    conda:
        "../../envs/python27.yaml"
    log:
        temp(
            "results/logs/counts/noUMI/create_BAM.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        set +o pipefail;

        fwd_length=`zcat {input.fw_fastq} | head -2 | tail -1 | wc -c`;
        fwd_length=$(expr $(($fwd_length-1)));

        rev_start=$(expr $(($fwd_length+1)));

        minoverlap=`echo ${{fwd_length}} ${{fwd_length}} {params.bc_length} | awk '{{print ($1+$2-$3-1 < 11) ? $1+$2-$3-1 : 11}}'`;

        echo $rev_start
        echo $minoverlap

        paste <( zcat {input.fw_fastq} ) <( zcat {input.rev_fastq}  ) | \
        awk '{{if (NR % 4 == 2 || NR % 4 == 0) {{
                print $1$2
            }} else {{
                print $1
            }}}}' | \
        python {input.script_FastQ2doubleIndexBAM} -p -s $rev_start -l 0 -m 0 --RG {params.datasetID} | \
        python {input.script_MergeTrimReadsBAM} --FirstReadChimeraFilter '' --adapterFirstRead '' --adapterSecondRead '' -p --mergeoverlap --minoverlap $minoverlap > {output} 2> {log}
        """


### START COUNTING ####


rule counts_noUMI_raw_counts:
    """
    Counting BCsxUMIs from the BAM files.
    """
    conda:
        "../../envs/bwa_samtools_picard_htslib.yaml"
    input:
        "results/experiments/{project}/counts/noUMI.{condition}_{replicate}_{type}.bam",
    output:
        "results/experiments/{project}/counts/noUMI.{condition}_{replicate}_{type}_raw_counts.tsv.gz",
    params:
        datasetID="{condition}_{replicate}_{type}",
    log:
        temp(
            "results/logs/counts/noUMI/raw_counts_umi.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        samtools view -F 1 -r {params.datasetID} {input} | \
        awk -v 'OFS=\\t' '{{ print $10 }}' | \
        sort | \
        gzip -c > {output} 2> {log}
        """
