######################################
### Everything before assigning BC ###
######################################

### Create_BAM_umi with demultiplexing ###


rule create_demultiplexed_index:
    conda:
        "../envs/python3.yaml"
    input:
        experiment_file=lambda wc: config["experiments"][wc.project]["experiment_file"],
        script=getScript("count/create_demultiplexed_index.py"),
    output:
        "results/experiments/{project}/counts/demultiplex_index.tsv",
    log:
        "logs/experiments/{project}/counts/create_demultiplexed_index.log",
    shell:
        """
        python {input.script} \
        --experiment {input.experiment_file} \
        --output {output} > {log}
        """


checkpoint create_demultiplexed_BAM_umi:
    input:
        fw_fastq=lambda wc: getFWWithIndex(wc.project),
        rev_fastq=lambda wc: getRevWithIndex(wc.project),
        umi_fastq=lambda wc: getUMIWithIndex(wc.project),
        index_fastq=lambda wc: getIndexWithIndex(wc.project),
        index_list="results/experiments/{project}/counts/demultiplex_index.tsv",
        script=getScript("count/SplitFastQdoubleIndexBAM.py"),
    output:
        "results/experiments/{project}/counts/demultiplex_{name}.bam",
    params:
        outdir=lambda w, output: os.path.split(output[0])[0],
    conda:
        "../envs/python27.yaml"
    log:
        "logs/experiments/{project}/counts/create_demultiplexed_BAM_umi.{name}.log",
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
            ) > {log}
        """



rule aggregate_demultiplex:
    input:
        lambda wc: aggregate_demultiplex_input(wc.project),
    output:
        touch("results/experiments/{project}/counts/demultiplex.done"),


rule mergeTrimReads_demultiplexed_BAM_umi:
    input:
        demultiplex="results/experiments/{project}/counts/demultiplex.done",
        script=getScript("count/MergeTrimReadsBAM.py"),
    output:
        "results/experiments/{project}/counts/merged_demultiplex_{condition}_{replicate}_{type}.bam",
    conda:
        "../envs/python27.yaml"
    params:
        bam="results/experiments/{project}/counts/demultiplex_{condition}_{replicate}_{type}.bam",
    log:
        "logs/experiments/{project}/counts/mergeTrimReads_demultiplexed_BAM_umi.{condition}_{replicate}_{type}.log",
    shell:
        """
        samtools view -h {params.bam} | \
        python {input.script} -p --mergeoverlap -f ACCGGTCGCCACCATGGTGAGCAAGGGCGAGGA -s CTTAGCTTTCGCTTAGCGATGTGTTCACTTTGC \
        > {output} 2> {log}
        """


### Create_BAM_umi without demultiplexing ###


rule create_BAM_umi:
    input:
        fw_fastq=lambda wc: getFW(wc.project, wc.condition, wc.replicate, wc.type),
        rev_fastq=lambda wc: getRev(wc.project, wc.condition, wc.replicate, wc.type),
        umi_fastq=lambda wc: getUMI(wc.project, wc.condition, wc.replicate, wc.type),
        script_FastQ2doubleIndexBAM=getScript("count/FastQ2doubleIndexBAM.py"),
        script_MergeTrimReadsBAM=getScript("count/MergeTrimReadsBAM.py"),
    output:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}.bam",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
        datasetID="{condition}_{replicate}_{type}",
    conda:
        "../envs/python27.yaml"
    log:
        "logs/experiments/{project}/counts/create_BAM_umi.{condition}_{replicate}_{type}.log",
    shell:
        """
        set +o pipefail;

        umi_length=`zcat {input.umi_fastq} | head -2 | tail -1 | wc -c`;
        umi_length=$(expr $(($umi_length-1)));

        fwd_length=`zcat {input.fw_fastq} | head -2 | tail -1 | wc -c`;
        fwd_length=$(expr $(($fwd_length-1)));

        rev_start=$(expr $(($fwd_length+1)));

        minoverlap=`echo ${{fwd_length}} ${{fwd_length}} {params.bc_length} | awk '{{print ($1+$2-$3-1 < 11) ? $1+$2-$3-1 : 11}}'`;

        echo $rev_start
        echo $umi_length
        echo $minoverlap

        paste <( zcat {input.fw_fastq} ) <( zcat {input.rev_fastq}  ) <( zcat {input.umi_fastq} ) | \
        awk '{{if (NR % 4 == 2 || NR % 4 == 0) {{print $1$2$3}} else {{print $1}}}}' | \
        python {input.script_FastQ2doubleIndexBAM} -p -s $rev_start -l 0 -m $umi_length --RG {params.datasetID} | \
        python {input.script_MergeTrimReadsBAM} --FirstReadChimeraFilter '' --adapterFirstRead '' --adapterSecondRead '' -p --mergeoverlap --minoverlap $minoverlap > {output} 2> {log}
        """


### START COUNTING ####


rule raw_counts_umi:
    """
    Counting BCsxUMIs from the BAM files.
    """
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    input:
        lambda wc: getBamFile(wc.project, wc.condition, wc.replicate, wc.type),
    output:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_raw_counts.tsv.gz",
    params:
        umi_length=lambda wc: config["experiments"][wc.project]["umi_length"],
        datasetID="{condition}_{replicate}_{type}",
    log:
        "logs/experiments/{project}/counts/raw_counts_umi.{condition}_{replicate}_{type}.log",
    shell:
        """
        samtools view -F 1 -r {params.datasetID} {input} | \
        awk -v 'OFS=\\t' '{{ for (i=12; i<=NF; i++) {{
          if ($i ~ /^XJ:Z:/) print $10,substr($i,6,{params.umi_length})
        }}}}' | \
        sort | uniq -c | \
        awk -v 'OFS=\\t' '{{ print $2,$3,$1 }}' | \
        gzip -c > {output}
        """


rule filter_counts:
    """
    Filter the counts to BCs only of the correct length (defined in the config file)
    """
    conda:
        "../envs/default.yaml"
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_raw_counts.tsv.gz",
    output:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
    log:
        "logs/experiments/{project}/counts/filter_counts.{condition}_{replicate}_{type}",
    shell:
        """
        bc={params.bc_length};
        echo $bc;
        zcat {input} | grep -v "N" | \
        awk -v var="$bc" -v 'OFS=\\t' '{{ if (length($1) == var) {{ print }} }}' | \
        sort | \
        gzip -c > {output}
        """


rule final_counts_umi:
    """
    Discarding PCR duplicates (taking BCxUMI only one time)
    """
    conda:
        "../envs/default.yaml"
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    output:
        counts="results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
    log:
        "logs/experiments/{project}/counts/final_counts_umi.{condition}_{replicate}_{type}.log",
    shell:
        """
        zcat {input} | awk '{{print $1}}' | \
        uniq -c | \
        awk -v 'OFS=\\t' '{{ print $2,$1 }}' | \
        gzip -c > {output.counts}
        """


rule final_counts_umi_samplerer:
    """
    Creates full + new distribution DNA files
    """
    input:
        counts="results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
        script=getScript("count/samplerer.py"),
    output:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts.sampling.{config}.tsv.gz",
    conda:
        "../envs/python3.yaml"
    params:
        samplingprop=lambda wc: counts_getSamplingConfig(
            wc.project, wc.config, wc.type, "prop"
        ),
        downsampling=lambda wc: counts_getSamplingConfig(
            wc.project, wc.config, wc.type, "threshold"
        ),
        samplingtotal=lambda wc: counts_getSamplingConfig(
            wc.project, wc.config, wc.type, "total"
        ),
        seed=lambda wc: counts_getSamplingConfig(wc.project, wc.config, wc.type, "seed"),
    log:
        "logs/experiments/{project}/counts/final_counts_umi_samplerer.{condition}_{replicate}_{type}_{config}.log",
    shell:
        """
        python {input.script} --input {input.counts} \
        {params.samplingprop} \
        {params.downsampling} \
        {params.samplingtotal} \
        {params.seed} \
        --output {output} > {log}
        """


rule dna_rna_merge_counts:
    """
    Merge DNA and RNA counts together.
    Is done in two ways. First no not allow zeros in DNA or RNA BCs (withoutZeros).
    Second with zeros, so a BC can be defined only in the DNA or RNA (withZeros)
    """
    conda:
        "../envs/default.yaml"
    input:
        dna=lambda wc: getFinalCounts(wc.project, wc.config, "DNA", wc.raw_or_assigned),
        rna=lambda wc: getFinalCounts(wc.project, wc.config, "RNA", wc.raw_or_assigned),
    output:
        "results/experiments/{project}/{raw_or_assigned}/{condition}_{replicate}.merged.config.{config}.tsv.gz",
    params:
        zero=lambda wc: "false" if withoutZeros(wc.project, wc.config) else "true",
    log:
        "logs/experiments/{project}/{raw_or_assigned}/dna_rna_merge_counts.{condition}_{replicate}_{config}.log",
    shell:
        """
        zero={params.zero};
        if [[ -z "${{zero//false}}" ]]
        then
            echo "Using no zeros";
            join -1 1 -2 1 -t"$(echo -e '\\t')" \
            <( zcat  {input.dna} | sort ) \
            <( zcat {input.rna} | sort) | \
            gzip -c > {output};
        else
            echo "Allowing zeros!"
            join -e 0 -a1 -a2 -t"$(echo -e '\\t')" -o 0 1.2 2.2 \
            <( zcat  {input.dna} | sort ) \
            <( zcat {input.rna}  | sort) | \
            gzip -c > {output}
        fi > {log}
        """
