#################################
## count 10 most frequent UMIs ##
#################################


# count frequent UMIs per condition, replicate and DNA/RNA
rule statistic_frequent_umis:
    conda:
        "../../envs/default.yaml"
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    output:
        freqUMIs=(
            "results/experiments/{project}/stats/counts/freqUMIs_{condition}_{replicate}_{type}.txt"
        ),
    log:
        "logs/experiments/{project}/stats/counts/statistic_frequent_umis.{condition}_{replicate}_{type}.log"
    shell:
        """
        set +o pipefail;
        zcat {input} | cut -f 2 | sort | uniq -c | sort -nr | head > {output.freqUMIs}
        """


##############################
## Barcode base composition ##
##############################


rule statistic_barcode_base_composition:
    conda:
        "../../envs/python3.yaml"
    input:
        counts="results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
        script=getScript("count/nucleotideCountPerPosition.py"),
    output:
        bc=temp(
            "results/experiments/{project}/counts/{condition}_{replicate}_{type}_final.BC.tsv.gz"
        ),
        stats="results/experiments/{project}/stats/counts/BCNucleotideComposition/{condition}_{replicate}_{type}.tsv.gz",
    params:
        name="{condition}_{replicate}_{type}",
    log:
       "logs/experiments/{project}/stats/counts/BCNucleotideComposition/statistic_barcode_base_composition.{condition}_{replicate}_{type}.log"
    shell:
        """
        zcat {input.counts} | awk '{{print $1}}' | gzip -c > {output.bc};
        python {input.script} \
        --column 1 \
        --chunksize 100000 \
        --input {output.bc} \
        --output {output.stats} > {log}
        """


#################################
## Count statistic of barcodes ##
#################################


# count Reads, Barcodes per UMI, Barcodes and Unique UMIs
rule statistic_counts:
    conda:
        "../../envs/default.yaml"
    input:
        "results/experiments/{project}/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz",
    output:
        "results/experiments/{project}/stats/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz",
    params:
        cond="{condition}",
        rep="{replicate}",
        type="{type}",
    log:
        "logs/experiments/{project}/stats/counts/statistic_counts.{condition}_{replicate}_{type}_{countType}.log",
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
        gzip -c > {output}
        """


# concat DNA, RNA-counts (rule statistic_counts) for all experiments, and replicates
rule statistic_count_stats_merge:
    conda:
        "../../envs/default.yaml"
    input:
        getCountStats,
    output:
        "results/experiments/{project}/stats/counts/count_{countType}.tsv",
    log:
        "logs/experiments/{project}/stats/counts/statistic_count_stats_merge.{countType}.log",
    shell:
        """
        zcat {input} | sort -k1,1 -k3,3 -k2,2 > {output} 2> {log}
        """


# Statistic of barcodes shared between RNA&DNA per condition and replicate
rule statistic_BC_in_RNA_DNA:
    conda:
        "../../envs/default.yaml"
    input:
        dna="results/experiments/{project}/counts/{condition}_{replicate}_DNA_{countType}_counts.tsv.gz",
        rna="results/experiments/{project}/counts/{condition}_{replicate}_RNA_{countType}_counts.tsv.gz",
    output:
        "results/experiments/{project}/stats/counts/{condition}_{replicate}_{countType}_BC_in_RNA_DNA.tsv.gz",
    params:
        cond="{condition}",
        rep="{replicate}",
    log:
        "logs/experiments/{project}/stats/counts/statistic_BC_in_RNA_DNA.{condition}_{replicate}_{countType}.log",
    shell:
        """
        paste <( echo "{params.cond}") <( echo "{params.rep}") \
        <( join <( zcat {input.dna} | cut -f 1 | sort | uniq ) \
        <( zcat {input.rna} | cut -f 1 | sort | uniq ) | wc -l ) | \
        gzip -c > {output}
        """



# concat shared barcodes (rule statistic_BC_in_RNA_DNA) for all experiments, and replicates
rule statistic_BC_in_RNA_DNA_merge:
    conda:
        "../../envs/default.yaml"
    input:
        getBCinRNADNAStats,
    output:
        "results/experiments/{project}/stats/counts/BC_in_RNA_DNA_{countType}.tsv",
    log:
        "logs/experiments/{project}/stats/counts/statistic_BC_in_RNA_DNA_merge.{countType}.log"
    shell:
        """
        zcat {input} | sort -k1,1 -k2,2 > {output}
        """


# making final count statistics
rule statistic_counts_final:
    conda:
        "../../envs/r.yaml"
    input:
        counts="results/experiments/{project}/stats/counts/count_{countType}.tsv",
        shared="results/experiments/{project}/stats/counts/BC_in_RNA_DNA_{countType}.tsv",
        script=getScript("count/combine_count_stats.R"),
    output:
        "results/experiments/{project}/stats/statistic_count_{countType}.tsv",
    log:
        "logs/experiments/{project}/stats/statistic_counts_final.{countType}.log"
    shell:
        """
        Rscript {input.script} --count {input.counts} --shared {input.shared} --output {output} > {log}
        """
