rule assignment_mapping_exact_reference:
    """
    Create reference to map the exact design
    """
    conda:
        "../../envs/default.yaml"
    input:
        lambda wc: config["assignments"][wc.assignment]["reference"],
    output:
        "results/assignment/{assignment}/reference/reference_exact.fa",
    log:
        temp("results/logs/assignment/mapping_exact_reference.{assignment}.log"),
    shell:
        """
        paste <(
            cat {input} | awk '{{if ($1 ~ /^>/) {{ gsub(/[\\]\\[]/,"_"); print substr($1,2)}}}}';
            cat {input} | awk '{{if ($1 ~ /^>/) {{ gsub(/[\\]\\[]/,"_"); print substr($1,2)}}}}';
        ) <(
            cat {input} | awk '{{if ($1 ~ /^[^>]/) {{ seq=seq$1}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}';
            cat {input} | awk '{{if ($1 ~ /^[^>]/) {{ seq=seq$1}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}' | tr ACGTacgt TGCAtgca | rev;
        ) > {output}
        """


rule assignment_mapping_exact:
    """
    Map the reads to the reference and sort using exact match.
    """
    conda:
        "../../envs/default.yaml"
    input:
        reads="results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz",
        reference="results/assignment/{assignment}/reference/reference_exact.fa",
    output:
        temp("results/assignment/{assignment}/BCs/barcodes_exact.{split}.tsv"),
    log:
        temp("results/logs/assignment/mapping_exact.{assignment}.{split}.log"),
    params:
        alignment_start=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["alignment_start"],
        sequence_length=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["sequence_length"],
    shell:
        """
        # Look up exact matches in design file
        export LC_ALL=C # speed up sort

        awk -v "OFS=\\t" 'NR==FNR {{a[$2] = $1; next}} {{if ($3 in a) print $2,a[$3],"{params.sequence_length}M"; else print $2,"other","NA"}}' \
        <(
            cat {input.reference} | awk -v "OFS=\\t" '{{print $1,substr($2, {params.alignment_start},{params.sequence_length})}}'
        ) \
        <(
            zcat {input.reads} | awk 'NR%4==2 || NR%4==1' | paste - -
        ) | \
        sed 's/,YI:Z[^\\t]*//g' | sed 's/XI:Z://g' | \
        sort -k1,1 -k2,2 -k3,3 -S 7G > {output} 2> {log}
        """
