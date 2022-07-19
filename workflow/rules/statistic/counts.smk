#################################
## count 10 most frequent UMIs ##
#################################


# count frequent UMIs per condition, replicate and DNA/RNA
rule statistic_counts_frequent_umis:
    conda:
        "../../envs/default.yaml"
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    output:
        report(
            "results/experiments/{project}/statistic/counts.freqUMIs.{condition}_{replicate}_{type}.txt",
            caption="../../report/counts/frequent_umis.rst",
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
            "results/logs/statistic/counts/frequent_umis.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        set +o pipefail;
        zcat {input} | cut -f 2 | sort | uniq -c | sort -nr | head > {output} 2> {log}
        """


##############################
## Barcode base composition ##
##############################


rule statistic_counts_barcode_base_composition:
    conda:
        "../../envs/python3.yaml"
    input:
        counts="results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
        script=getScript("count/nucleotideCountPerPosition.py"),
    output:
        bc=temp(
            "results/experiments/{project}/counts/{condition}_{replicate}_{type}_final.BC.tsv.gz"
        ),
        stats=report(
            "results/experiments/{project}/statistic/counts/BCNucleotideComposition.{condition}_{replicate}_{type}.tsv.gz",
            caption="../../report/counts/barcode_base_composition.rst",
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
            "results/logs/statistic/counts/barcode_base_composition.{project}.{condition}_{replicate}_{type}.log"
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


# count Reads, Barcodes per UMI, Barcodes and Unique UMIs
rule statistic_counts_table:
    conda:
        "../../envs/default.yaml"
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz",
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
            "results/logs/statistic/counts/table.{project}.{condition}_{replicate}_{type}_{countType}.log"
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


# concat DNA, RNA-counts (rule statistic_counts) for all experiments, and replicates
rule statistic_counts_stats_merge:
    conda:
        "../../envs/default.yaml"
    input:
        lambda wc: getCountStats(wc.project, wc.countType),
    output:
        temp("results/experiments/{project}/statistic/counts/count_{countType}.tsv"),
    log:
        temp("results/logs/statistic/counts/stats_merge.{project}.{countType}.log"),
    shell:
        """
        zcat {input} | sort -k1,1 -k3,3 -k2,2 > {output} 2> {log}
        """


# Statistic of barcodes shared between RNA&DNA per condition and replicate
rule statistic_counts_BC_in_RNA_DNA:
    conda:
        "../../envs/default.yaml"
    input:
        dna="results/experiments/{project}/counts/{condition}_{replicate}_DNA_{countType}_counts.tsv.gz",
        rna="results/experiments/{project}/counts/{condition}_{replicate}_RNA_{countType}_counts.tsv.gz",
    output:
        temp(
            "results/experiments/{project}/statistic/counts/{condition}_{replicate}_{countType}_BC_in_RNA_DNA.tsv.gz"
        ),
    params:
        cond="{condition}",
        rep="{replicate}",
    log:
        temp(
            "results/logs/statistic/counts/BC_in_RNA_DNA.{project}.{condition}_{replicate}_{countType}.log"
        ),
    shell:
        """
        paste <( echo "{params.cond}") <( echo "{params.rep}") \
        <( join <( zcat {input.dna} | cut -f 1 | sort | uniq ) \
        <( zcat {input.rna} | cut -f 1 | sort | uniq ) | wc -l ) | \
        gzip -c > {output} 2> {log}
        """


# concat shared barcodes (rule statistic_BC_in_RNA_DNA) for all experiments, and replicates
rule statistic_counts_BC_in_RNA_DNA_merge:
    conda:
        "../../envs/default.yaml"
    input:
        getBCinRNADNAStats,
    output:
        temp(
            "results/experiments/{project}/statistic/counts/BC_in_RNA_DNA_{countType}.tsv"
        ),
    log:
        temp(
            "results/logs/statistic/counts/BC_in_RNA_DNA_merge.{project}.{countType}.log"
        ),
    shell:
        """
        zcat {input} | sort -k1,1 -k2,2 > {output} 2> {log}
        """


# making final count statistics
rule statistic_counts_final:
    conda:
        "../../envs/r.yaml"
    input:
        counts="results/experiments/{project}/statistic/counts/count_{countType}.tsv",
        shared="results/experiments/{project}/statistic/counts/BC_in_RNA_DNA_{countType}.tsv",
        script=getScript("count/combine_count_stats.R"),
    output:
        stats=report(
            "results/experiments/{project}/statistic/counts.{countType}.tsv",
            caption="../../report/counts/counts.rst",
            category="{project}",
            subcategory="Counts",
            labels={
                "analysis": "{countType} count statistic",
            },
        ),
    log:
        temp("results/logs/statistic/counts/final.{project}.{countType}.log"),
    shell:
        """
        Rscript {input.script} --count {input.counts} --shared {input.shared} --output {output} > {log}
        """
