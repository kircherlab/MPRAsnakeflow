include: "counts_common.smk"


#################################
## count 10 most frequent UMIs ##
#################################


rule experiment_statistic_counts_frequent_umis:
    """
    Count the 10 most frequent UMIs per condition, replicate and DNA/RNA.
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    output:
        report(
            "results/experiments/{project}/statistic/counts.freqUMIs.{condition}_{replicate}_{type}.txt",
            caption="../../../report/counts/frequent_umis.rst",
            category="{project}",
            subcategory="UMIs",
            labels={
                "condition": "{condition}",
                "replicate": "{replicate}",
                "DNA/RNA": "{type}",
            },
        ),
    log:
        temp(
            "results/logs/experiment/statistic/counts/frequent_umis.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        set +o pipefail;
        zcat {input} | cut -f 2 | sort | uniq -c | sort -nr | head > {output} 2> {log}
        """


##############################
## Barcode base composition ##
##############################


rule experiment_statistic_counts_barcode_base_composition:
    """
    Count the nucleotide composition of the barcodes per condition, replicate and DNA/RNA.
    """
    conda:
        getCondaEnv("python3.yaml")
    input:
        counts="results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
        script=getScript("count/nucleotideCountPerPosition.py"),
    output:
        bc=temp(
            "results/experiments/{project}/counts/{condition}_{replicate}_{type}_final.BC.tsv.gz"
        ),
        stats=report(
            "results/experiments/{project}/statistic/counts/BCNucleotideComposition.{condition}_{replicate}_{type}.tsv.gz",
            caption="../../../report/counts/barcode_base_composition.rst",
            category="{project}",
            subcategory="Barcode nucleotides",
            labels={
                "Analysis": "Nucleotide Composition",
                "Condition": "{condition}",
                "Replicate": "{replicate}",
                "DNA/RNA": "{type}",
            },
        ),
    params:
        name="{condition}_{replicate}_{type}",
    log:
        temp(
            "results/logs/experiment/statistic/counts/barcode_base_composition.{project}.{condition}_{replicate}_{type}.log"
        ),
    shell:
        """
        zcat {input.counts} | awk '{{print $1}}' | gzip -c > {output.bc};
        python {input.script} \
        --column 1 \
        --chunksize 100000 \
        --input {output.bc} \
        --output {output.stats} &> {log}
        """


#################################
## Count statistic of barcodes ##
#################################


rule experiment_statistic_counts_table:
    """
    Count statistic of barcodes and UMIs per condition, replicate and DNA/RNA.
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        lambda wc: (
            "results/experiments/{project}/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz"
            if wc.countType != "raw"
            else getRawCounts(wc.project, wc.type)
        ),
    output:
        temp(
            "results/experiments/{project}/statistic/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz"
        ),
    params:
        cond="{condition}",
        rep="{replicate}",
        type="{type}",
    log:
        temp(
            "results/logs/experiment/statistic/counts/table.{project}.{condition}_{replicate}_{type}_{countType}.log"
        ),
    shell:
        """
        paste <( echo "{params.cond}") <( echo "{params.rep}") <( echo "{params.type}") \
        <( 
            zcat {input} | \
            awk -v OFS='\\t' 'BEGIN{{
                pbar="NA"
            }}{{ 
                count += $NF; umi_sum+=$3; if (pbar != $1) {{ barcodes+=1 }}; pbar=$1 
            }}END{{ 
                if (NR > 0) {{
                    print umi_sum/NR,count,NR,barcodes
                }} else {{
                    print 0,0,0,0
                }}
            }}' 
        ) \
        <( 
            zcat {input} | cut -f 2 | sort -u | wc -l 
        ) | \
        gzip -c > {output} 2> {log}
        """


rule experiment_statistic_counts_stats_merge:
    """
    Merge the count statistic of all replicates and conditions into one table.
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        lambda wc: getCountStats(wc.project, wc.countType),
    output:
        temp("results/experiments/{project}/statistic/counts/count_{countType}.tsv"),
    log:
        temp(
            "results/logs/experiment/statistic/counts/stats_merge.{project}.{countType}.log"
        ),
    shell:
        """
        zcat {input} | sort -k1,1 -k3,3 -k2,2 > {output} 2> {log}
        """


rule experiment_statistic_counts_BC_in_RNA_DNA:
    """
    Count the number of barcodes shared between RNA and DNA per condition and replicate.
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        dna=lambda wc: statistic_counts_BC_in_RNA_DNA_helper(
            project, wc.condition, "DNA", wc.countType
        ),
        rna=lambda wc: statistic_counts_BC_in_RNA_DNA_helper(
            project, wc.condition, "RNA", wc.countType
        ),
    output:
        temp(
            "results/experiments/{project}/statistic/counts/{condition}_{replicate}_{countType}_BC_in_RNA_DNA.tsv.gz"
        ),
    params:
        cond="{condition}",
        rep="{replicate}",
    log:
        temp(
            "results/logs/experiment/statistic/counts/BC_in_RNA_DNA.{project}.{condition}_{replicate}_{countType}.log"
        ),
    shell:
        """
        paste <( echo "{params.cond}") <( echo "{params.rep}") \
        <( join <( zcat {input.dna} | cut -f 1 | sort | uniq ) \
        <( zcat {input.rna} | cut -f 1 | sort | uniq ) | wc -l ) | \
        gzip -c > {output} 2> {log}
        """


rule experiment_statistic_counts_BC_in_RNA_DNA_merge:
    """
    Merge the shared barcodes statistic of all replicates and conditions into one table.
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        getBCinRNADNAStats,
    output:
        temp(
            "results/experiments/{project}/statistic/counts/BC_in_RNA_DNA_{countType}.tsv"
        ),
    log:
        temp(
            "results/logs/experiment/statistic/counts/BC_in_RNA_DNA_merge.{project}.{countType}.log"
        ),
    shell:
        """
        zcat {input} | sort -k1,1 -k2,2 > {output} 2> {log}
        """


rule experiment_statistic_counts_final:
    """
    Combine the final count statistic of all replicates and conditions into one table.
    """
    conda:
        getCondaEnv("r.yaml")
    input:
        counts="results/experiments/{project}/statistic/counts/count_{countType}.tsv",
        shared="results/experiments/{project}/statistic/counts/BC_in_RNA_DNA_{countType}.tsv",
        script=getScript("count/combine_count_stats.R"),
    output:
        stats=report(
            "results/experiments/{project}/statistic/counts.{countType}.tsv",
            caption="../../../report/counts/counts.rst",
            category="{project}",
            subcategory="Counts",
            labels={
                "analysis": "{countType} count statistic",
            },
        ),
    log:
        temp("results/logs/experiment/statistic/counts/final.{project}.{countType}.log"),
    shell:
        """
        Rscript {input.script} --count {input.counts} --shared {input.shared} --output {output} > {log}
        """
