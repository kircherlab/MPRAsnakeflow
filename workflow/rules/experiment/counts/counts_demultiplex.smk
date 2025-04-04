######################################
### Everything before assigning BC ###
######################################

### Create_BAM_umi with demultiplexing ###


rule counts_demultiplex_create_index:
    conda:
        getCondaEnv("python3.yaml")
    input:
        experiment_file=lambda wc: config["experiments"][wc.project]["experiment_file"],
        script=getScript("count/create_demultiplexed_index.py"),
    output:
        "results/experiments/{project}/counts/demultiplex_index.tsv",
    log:
        temp("results/logs/counts/create_demultiplexed_index.{project}.log"),
    shell:
        """
        python {input.script} \
        --experiment {input.experiment_file} \
        --output {output} &> {log}
        """


checkpoint counts_demultiplex_BAM_umi:
    input:
        fw_fastq=lambda wc: getFWWithIndex(wc.project),
        rev_fastq=lambda wc: getRevWithIndex(wc.project),
        umi_fastq=lambda wc: getUMIWithIndex(wc.project),
        index_fastq=lambda wc: getIndexWithIndex(wc.project),
        index_list="results/experiments/{project}/counts/demultiplex_index.tsv",
        script=getScript("count/SplitFastQdoubleIndexBAM.py"),
    output:
        "results/experiments/{project}/counts/demultiplex.{name}.bam",
    params:
        outdir=lambda w, output: os.path.split(output[0])[0],
    conda:
        getCondaEnv("python27.yaml")
    log:
        temp("results/logs/counts/demultiplex_BAM_umi.{project}.{name}.log"),
    shell:
        """
            set +o pipefail;

            umi_length=`zcat {input.umi_fastq} | head -2 | tail -1 | wc -c`;
            umi_length=$(expr $(($umi_length-1)));

            idx_length=`zcat {input.index_fastq} | head -2 | tail -1 | wc -c`;
            idx_length=$(expr $(($idx_length-1)));

            fwd_length=`zcat {input.fw_fastq} | head -2 | tail -1 | wc -c`;
            fwd_length=$(expr $(($fwd_length-1)));

            rev_start=$(expr $(($fwd_length+$idx_length+1)));

            echo $rev_start
            echo $idx_length
            echo $umi_length

            python {input.script} -s $rev_start -l $idx_length -m $umi_length -i {input.index_list} --outdir {params.outdir} --remove --summary --separate_files \
            <(\
        paste <( zcat {input.fw_fastq} ) <( zcat {input.index_fastq} ) <( zcat {input.rev_fastq} ) <( zcat {input.umi_fastq} ) | \
            awk '{{ count+=1; if ((count == 1) || (count == 3)) {{ print $1 }} else {{ print $1$2$3$4 }}; if (count == 4) {{ count=0 }} }}'\
            ) &> {log}
        """


rule counts_demultiplex_aggregate:
    input:
        lambda wc: counts_aggregate_demultiplex_input(wc.project),
    output:
        touch("results/experiments/{project}/counts/demultiplex.done"),


rule counts_demultiplex_mergeTrimReads_BAM_umi:
    input:
        demultiplex="results/experiments/{project}/counts/demultiplex.done",
        script=getScript("count/MergeTrimReadsBAM.py"),
    output:
        "results/experiments/{project}/counts/merged_demultiplex.{condition}_{replicate}_{type}.bam",
    conda:
        getCondaEnv("python27.yaml")
    params:
        bam="results/experiments/{project}/counts/demultiplex.{condition}_{replicate}_{type}.bam",
    log:
        temp(
            "results/logs/counts/mergeTrimReads_demultiplex_BAM_umi.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        samtools view -h {params.bam} | \
        python {input.script} -p --mergeoverlap -f ACCGGTCGCCACCATGGTGAGCAAGGGCGAGGA -s CTTAGCTTTCGCTTAGCGATGTGTTCACTTTGC \
        > {output} 2> {log}
        """
