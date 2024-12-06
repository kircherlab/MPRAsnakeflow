def statistic_counts_BC_in_RNA_DNA_helper(project, condition, dna_or_rna, countType):

    replicates = getReplicatesOfConditionType(project, condition, dna_or_rna)

    if countType == "raw":
        output = getRawCounts(project, dna_or_rna)
    else:
        output = (
            "results/experiments/{project}/counts/{condition}_{replicate}_%s_{countType}_counts.tsv.gz"
            % dna_or_rna
        )

    if len(replicates) == 1:
        output = output.replace("{replicate}", replicates[0])

    return output


# get all counts of experiment (rule statistic_counts)
def getCountStats(project, countType):
    exp = getExperiments(project)
    output = []
    for index, row in exp.iterrows():
        condition=row["Condition"]
        for dna_or_rna in ["DNA", "RNA"]:
            replicates = getReplicatesOfConditionType(project, condition, dna_or_rna)
            if len(replicates) == 1:
                replicate = replicates[0]
            else:
                replicate = row["Replicate"]
            output += expand(
                "results/experiments/{{project}}/statistic/counts/{condition}_{replicate}_{type}_{{countType}}_counts.tsv.gz",
                condition=condition,
                replicate=replicate,
                type=dna_or_rna,
            )
    return output
