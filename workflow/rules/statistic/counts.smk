#################################
## count 10 most frequent UMIs ##
#################################


# count frequent UMIs per condition, replicate and DNA/RNA
rule frequent_umis:
    input:
        "results/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    output:
        freqUMIs=(
            "results/{project}/stats/counts/freqUMIs_{condition}_{replicate}_{type}.txt"
        ),
    shell:
        """
        set +o pipefail;
        zcat {input} | cut -f 2 | sort | uniq -c | sort -nr | head > {output.freqUMIs}
        """


##############################
## Barcode base composition ##
##############################


rule barcode_base_composition:
    conda:
        "../../envs/mpraflow_pandas.yaml"
    input:
        counts="results/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
    output:
        bc=temp(
            "results/{project}/counts/{condition}_{replicate}_{type}_final.BC.tsv.gz"
        ),
        stats="results/{project}/stats/counts/BCNucleotideComposition/{condition}_{replicate}_{type}.tsv.gz",
    params:
        name="{condition}_{replicate}_{type}",
    shell:
        """
        zcat {input.counts} | awk '{{print $2}}' | gzip -c > {output.bc};
        python {SCRIPTS_DIR}/count/nucleotideCountPerPosition.py \
        --column 1 \
        --chunksize 100000 \
        --input {output.bc} \
        --output {output.stats}
        """


#################################
## Count statistic of barcodes ##
#################################


# count Reads, Barcodes per UMI, Barcodes and Unique UMIs
rule statistic_counts:
    input:
        "results/{project}/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz",
    output:
        "results/{project}/stats/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz",
    params:
        cond="{condition}",
        rep="{replicate}",
        type="{type}",
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


# get all counts of experiment (rule statistic_counts)
def getCountStats(wc):
    exp = getExperiments(wc.project)
    output = []
    for index, row in exp.iterrows():
        output += expand(
            "results/{project}/stats/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz",
            project=wc.project,
            condition=row["Condition"],
            replicate=row["Replicate"],
            type=["DNA", "RNA"],
            countType=wc.countType,
        )
    return output


# concat DNA, RNA-counts (rule statistic_counts) for all experiments, and replicates
rule count_stats_merge:
    input:
        getCountStats,
    output:
        "results/{project}/stats/counts/count_{countType}.tsv",
    shell:
        """
        zcat {input} | sort -k1,1 -k3,3 -k2,2 > {output}
        """


# Statistic of barcodes shared between RNA&DNA per condition and replicate
rule statistic_BC_in_RNA_DNA:
    input:
        dna="results/{project}/counts/{condition}_{replicate}_DNA_{countType}_counts.tsv.gz",
        rna="results/{project}/counts/{condition}_{replicate}_RNA_{countType}_counts.tsv.gz",
    output:
        "results/{project}/stats/counts/{condition}_{replicate}_{countType}_BC_in_RNA_DNA.tsv.gz",
    params:
        cond="{condition}",
        rep="{replicate}",
    shell:
        """
        paste <( echo "{params.cond}") <( echo "{params.rep}") \
        <( join <( zcat {input.dna} | cut -f 1 | sort | uniq ) \
        <( zcat {input.rna} | cut -f 1 | sort | uniq ) | wc -l ) | \
        gzip -c > {output}
        """


# get all barcodes of experiment (rule statistic_BC_in_RNA_DNA)
def getBCinRNADNAStats(wc):
    exp = getExperiments(wc.project)
    output = []
    for index, row in exp.iterrows():
        output += expand(
            "results/{project}/stats/counts/{condition}_{replicate}_{countType}_BC_in_RNA_DNA.tsv.gz",
            project=wc.project,
            condition=row["Condition"],
            replicate=row["Replicate"],
            countType=wc.countType,
        )
    return output


# concat shared barcodes (rule statistic_BC_in_RNA_DNA) for all experiments, and replicates
rule stats_BC_in_RNA_DNA_merge:
    input:
        getBCinRNADNAStats,
    output:
        "results/{project}/stats/counts/BC_in_RNA_DNA_{countType}.tsv",
    shell:
        """
        zcat {input} | sort -k1,1 -k2,2 > {output}
        """


# making final count statistics
rule stats_final:
    conda:
        "../../envs/mpraflow_r.yaml"
    input:
        counts="results/{project}/stats/counts/count_{countType}.tsv",
        shared="results/{project}/stats/counts/BC_in_RNA_DNA_{countType}.tsv",
    output:
        "results/{project}/stats/statistic_count_{countType}.tsv",
    shell:
        """
        Rscript {SCRIPTS_DIR}/count/combine_count_stats.R --count {input.counts} --shared {input.shared} --output {output}
        """
