######################################
### Everything before assigning BC ###
######################################

### Create_BAM_umi with demultiplexing ###


rule create_demultiplexed_index:
    output:
        "results/experiments/{project}/counts/demultiplex_index.tsv",
    run:
        import csv
        import pandas as pd

        exp = getExperiments(wildcards.project).astype(str)

        exp_dna = pd.concat(
            [
                exp[["BC_DNA"]],
                exp[["Condition", "Replicate"]].agg("_".join, axis=1) + "_DNA",
            ],
            axis=1,
        ).rename(columns={"BC_DNA": "BC", 0: "Sample"})

        exp_rna = pd.concat(
            [
                exp[["BC_RNA"]],
                exp[["Condition", "Replicate"]].agg("_".join, axis=1) + "_RNA",
            ],
            axis=1,
        ).rename(columns={"BC_RNA": "BC", 0: "Sample"})

        exp_dna.append(exp_rna).to_csv(output[0], sep="\t", header=False, index=False)


checkpoint create_demultiplexed_BAM_umi:
    input:
        fw_fastq=lambda wc: getFWWithIndex(wc.project),
        rev_fastq=lambda wc: getRevWithIndex(wc.project),
        umi_fastq=lambda wc: getUMIWithIndex(wc.project),
        index_fastq=lambda wc: getIndexWithIndex(wc.project),
        index_list="results/experiments/{project}/counts/demultiplex_index.tsv",
    output:
        "results/experiments/{project}/counts/demultiplex_{name}.bam",
    params:
        outdir="results/experiments/{project}/counts/",
    conda:
        "../envs/python27.yaml"
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

            python {SCRIPTS_DIR}/count/SplitFastQdoubleIndexBAM.py -s $rev_start -l $idx_length -m $umi_length -i {input.index_list} --outdir {params.outdir} --remove --summary --separate_files \
            <(\
        paste <( zcat {input.fw_fastq} ) <( zcat {input.index_fastq} ) <( zcat {input.rev_fastq} ) <( zcat {input.umi_fastq} ) | \
            awk '{{ count+=1; if ((count == 1) || (count == 3)) {{ print $1 }} else {{ print $1$2$3$4 }}; if (count == 4) {{ count=0 }} }}'\
            )
        """


def aggregate_input(project):
    output = []
    conditions = getConditions(project)
    for condition in conditions:
        replicates = getReplicatesOfCondition(project, condition)
        names = expand(
            "{condition}_{replicate}_{type}",
            condition=condition,
            replicate=replicates,
            type=["DNA", "RNA"],
        )
        for name in names:
            with checkpoints.create_demultiplexed_BAM_umi.get(
                project=project, name=name
            ).output[0].open() as f:
                output += [f.name]
    return output


rule aggregate_demultiplex:
    input:
        lambda wc: aggregate_input(wc.project),
    output:
        touch("results/experiments/{project}/counts/demultiplex.done"),


rule mergeTrimReads_demultiplexed_BAM_umi:
    input:
        demultiplex="results/experiments/{project}/counts/demultiplex.done",
    output:
        "results/experiments/{project}/counts/merged_demultiplex_{condition}_{replicate}_{type}.bam",
    conda:
        "../envs/python27.yaml"
    params:
        bam="results/experiments/{project}/counts/demultiplex_{condition}_{replicate}_{type}.bam",
    shell:
        """
        samtools view -h {params.bam} | \
        python {SCRIPTS_DIR}/count/MergeTrimReadsBAM.py -p --mergeoverlap -f ACCGGTCGCCACCATGGTGAGCAAGGGCGAGGA -s CTTAGCTTTCGCTTAGCGATGTGTTCACTTTGC \
        > {output}
        """


### Create_BAM_umi without demultiplexing ###


rule create_BAM_umi:
    input:
        fw_fastq=lambda wc: getFW(wc.project, wc.condition, wc.replicate, wc.type),
        rev_fastq=lambda wc: getRev(wc.project, wc.condition, wc.replicate, wc.type),
        umi_fastq=lambda wc: getUMI(wc.project, wc.condition, wc.replicate, wc.type),
    output:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}.bam",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
        datasetID="{condition}_{replicate}_{type}",
    conda:
        "../envs/python27.yaml"
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
        python {SCRIPTS_DIR}/count/FastQ2doubleIndexBAM.py -p -s $rev_start -l 0 -m $umi_length --RG {params.datasetID} | \
        python {SCRIPTS_DIR}/count/MergeTrimReadsBAM.py --FirstReadChimeraFilter '' --adapterFirstRead '' --adapterSecondRead '' -p --mergeoverlap --minoverlap $minoverlap > {output}
        """


### START COUNTING ####


def getBam(project, condition, replicate, type):
    """
    gelper to get the correct BAM file (demultiplexed or not)
    """
    if config["experiments"][project]["demultiplex"]:
        return "results/%s/counts/merged_demultiplex_%s_%s_%s.bam" % (
            project,
            condition,
            replicate,
            type,
        )
    else:
        return "results/experiments/%s/counts/%s_%s_%s.bam" % (
            project,
            condition,
            replicate,
            type,
        )


rule raw_counts_umi:
    """
    Counting BCsxUMIs from the BAM files.
    """
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    input:
        lambda wc: getBam(wc.project, wc.condition, wc.replicate, wc.type),
    output:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_raw_counts.tsv.gz",
    params:
        umi_length=lambda wc: config["experiments"][wc.project]["umi_length"],
        datasetID="{condition}_{replicate}_{type}",
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
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_raw_counts.tsv.gz",
    output:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
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
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    output:
        counts="results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz"
    shell:
        """
        zcat {input} | awk '{{print $1}}' | \
        uniq -c | \
        awk -v 'OFS=\\t' '{{ print $2,$1 }}' | \
        gzip -c > {output.counts}
        """


def useSampling(project, conf, dna_or_rna):
    return "sampling" in config["experiments"][project]["configs"][conf] and dna_or_rna in config["experiments"][project]["configs"][conf]["sampling"]


def counts_getSamplingConfig(project, conf, dna_or_rna, command):
    if useSampling(project, conf, dna_or_rna):
        if dna_or_rna in config["experiments"][project]["configs"][conf]["sampling"]:
            if command in config["experiments"][project]["configs"][conf]["sampling"][dna_or_rna]:
                return "--%s %f" % (
                    command,
                    config["experiments"][project]["configs"][conf]["sampling"][dna_or_rna][command],
                )

    return ""


rule final_counts_umi_samplerer:
    """
    Creates full + new distribution DNA files
    """
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
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
        seed=lambda wc: counts_getSamplingConfig(
            wc.project, wc.config, wc.type, "seed"
        ),
    wildcard_constraints:
        downsampling="^full",
    shell:
        """
        python {SCRIPTS_DIR}/count/samplerer.py --input {input} \
        {params.samplingprop} \
        {params.downsampling} \
        {params.seed} \
        --output {output}
        """


def getFinalCounts(project, conf, rna_or_dna, raw_or_assigned):
    output = ""
    if raw_or_assigned == "counts":
        if useSampling(project, conf, rna_or_dna):
            output = (
                "results/experiments/{project}/%s/{condition}_{replicate}_%s_final_counts.sampling.{config}.tsv.gz"
                % (raw_or_assigned, rna_or_dna)
            )

        else:
            output = (
                "results/experiments/{project}/%s/{condition}_{replicate}_%s_final_counts.tsv.gz"
                % (raw_or_assigned, rna_or_dna)
            )
    else:
        output = (
            "results/experiments/{project}/%s/{condition}_{replicate}_%s_final_counts.config.{config}.tsv.gz"
            % (raw_or_assigned, rna_or_dna)
        )
    return output


def withoutZeros(project, conf):
    return (
        config["experiments"][project]["configs"][conf]["minDNACounts"] > 0
        and config["experiments"][project]["configs"][conf]["minRNACounts"] > 0
    )


rule dna_rna_merge_counts:
    """
    Merge DNA and RNA counts together.
    Is done in two ways. First no not allow zeros in DNA or RNA BCs (withoutZeros).
    Second with zeros, so a BC can be defined only in the DNA or RNA (withZeros)
    """
    input:
        dna=lambda wc: getFinalCounts(wc.project, wc.config, "DNA", wc.raw_or_assigned),
        rna=lambda wc: getFinalCounts(wc.project, wc.config, "RNA", wc.raw_or_assigned),
    output:
        "results/experiments/{project}/{raw_or_assigned}/{condition}_{replicate}.merged.config.{config}.tsv.gz",
    params:
        zero=lambda wc: "false" if withoutZeros(wc.project, wc.config) else "true",
    shell:
        """
        zero={params.zero};
        if [[ $zero=false ]]
        then
            join -1 1 -2 1 -t"$(echo -e '\\t')" \
            <( zcat  {input.dna} | sort ) \
            <( zcat {input.rna} | sort) | \
            gzip -c > {output}
        else
            join -e 0 -a1 -a2 -t"$(echo -e '\\t')" -o 0 1.2 2.2 \
            <( zcat  {input.dna} | sort ) \
            <( zcat {input.rna}  | sort) | \
            gzip -c > {output}
        fi
        """
