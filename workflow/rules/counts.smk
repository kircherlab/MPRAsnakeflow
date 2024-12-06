######################################
### Everything before assigning BC ###
######################################


include: "counts/counts_demultiplex.smk"
include: "counts/counts_umi.smk"
include: "counts/counts_noUMI.smk"
include: "counts/counts_onlyFW.smk"


rule counts_filter_counts:
    """
    Filter the counts to BCs only of the correct length (defined in the config file)
    """
    conda:
        "../envs/default.yaml"
    input:
        lambda wc: getRawCounts(wc.project, wc.type),
    output:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
    log:
        temp(
            "results/logs/counts/filter_counts.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        bc={params.bc_length};
        echo $bc;
        zcat {input} | grep -v "N" | \
        awk -v var="$bc" -v 'OFS=\\t' '{{ if (length($1) == var) {{ print }} }}' | \
        sort | \
        gzip -c > {output}
        """


rule counts_final_counts:
    """
    Counting BCs.
    Discarding PCR duplicates (taking BCxUMI only one time)
    """
    conda:
        "../envs/default.yaml"
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    output:
        counts="results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
    log:
        temp(
            "results/logs/counts/final_counts_umi.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        zcat {input} | awk '{{print $1}}' | \
        uniq -c | \
        awk -v 'OFS=\\t' '{{ print $2,$1 }}' | \
        gzip -c > {output.counts} 2> {log}
        """


rule counts_final_counts_samplerer:
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
        filtermincounts=lambda wc: counts_getFilterConfig(
            wc.project, wc.config, wc.type, "min_counts"
        ),
    log:
        temp(
            "results/logs/counts/final_counts_umi_samplerer.{project}.{condition}.{replicate}.{type}.{config}.log"
        ),
    shell:
        """
        python {input.script} --input {input.counts} \
        {params.samplingprop} \
        {params.downsampling} \
        {params.samplingtotal} \
        {params.seed} \
        {params.filtermincounts} \
        --output {output} &> {log}
        """


rule counts_dna_rna_merge_counts:
    """
    Merge DNA and RNA counts together.
    Is done in two ways. First no not allow zeros in DNA or RNA BCs (RNA and DNA min_counts not zero).
    Second with zeros, so a BC can be defined only in the DNA or RNA (RNA or DNA min_counts zero)
    """
    conda:
        "../envs/default.yaml"
    input:
        dna=lambda wc: getFinalCounts(wc.project, wc.config, wc.condition, "DNA", wc.raw_or_assigned),
        rna=lambda wc: getFinalCounts(wc.project, wc.config, wc.condition, "RNA", wc.raw_or_assigned),
    output:
        "results/experiments/{project}/{raw_or_assigned}/{condition}_{replicate}.merged.config.{config}.tsv.gz",
    params:
        zero=lambda wc: "false" if withoutZeros(wc.project, wc.config) else "true",
        minRNACounts=lambda wc: counts_getFilterConfig(
            wc.project, wc.config, "RNA", "min_counts"
        ),
        minDNACounts=lambda wc: counts_getFilterConfig(
            wc.project, wc.config, "DNA", "min_counts"
        ),
    log:
        temp(
            "results/logs/{raw_or_assigned}/dna_rna_merge_counts.{project}.{condition}.{replicate}.{config}.log"
        ),
    shell:
        """
        zero={params.zero};
        if [[ -z "${{zero//false}}" ]]
        then
            join -1 1 -2 1 -t"$(echo -e '\\t')" \
            <( zcat  {input.dna} | sort ) \
            <( zcat {input.rna} | sort);
        else
            join -e 0 -a1 -a2 -t"$(echo -e '\\t')" -o 0 1.2 2.2 \
            <( zcat  {input.dna} | sort ) \
            <( zcat {input.rna}  | sort);
        fi  | \
        awk -v 'OFS=\\t' '{{if ($2 >= {params.minDNACounts} && $3 >= {params.minRNACounts}) {{print $0}}}}' | \
        gzip -c > {output} 2> {log}
        """
