# preprocessing.smk specific functions


def getExperimentCutadaptAdapters(project, read):
    output = []
    if "adapters" in config["experiments"][project] and read in config["experiments"][project]["adapters"]:
        adapters_config = config["experiments"][project]["adapters"][read]
        if isinstance(adapters_config, list) and isinstance(adapters_config[0], int):
            output = ["-u %d" % u for u in adapters_config]
        else:

            if "three_prime" in adapters_config:
                for adapter in adapters_config["three_prime"]:
                    output.append("-a %s" % adapter)
            if "five_prime" in adapters_config:
                for adapter in adapters_config["five_prime"]:
                    output.append("-g %s" % adapter)

            return " ".join(output)
    return " ".join(output)


def getMaxExperimentSplitNumber() -> int:
    splits = [1]

    for project in getProjects():
        splits += [config["experiments"][project]["split_number"]]

    return max(splits)


# count.smk specific functions


def useUMI(project: str, type: str = "DNA") -> bool:
    """
    Helper to check if UMI should be used
    """
    return "UMI" in experiments[project] or f"{type}_UMI" in experiments[project]


def onlyFWD(project: str, type: str = "DNA") -> bool:
    """
    Helper to check if only forward reads should be used
    """
    return f"{type}_BC_R" not in experiments[project]


def noUMI(project: str, type: str = "DNA") -> bool:
    """
    Helper to check if UMI should not be used
    """
    return (
        "UMI" not in experiments[project]
        and f"{type}_UMI" not in experiments[project]
        and f"{type}_BC_R" in experiments[project]
    )


def useTrimming(project: str, read_type: str) -> bool:
    """
    Helper to check if trimming should be used for a specific read type.
    """
    if "adapters" in config["experiments"][project]:
        if read_type in config["experiments"][project]["read_type"]:
            return True
    return False


def useSplitting(project: str, rnaDna_type: str) -> bool:
    """
    Helper to check if splitting should be used. Will only apply for merging FWD and REV reads (creating BAM file) and split_number is > 1.
    """
    return not onlyFWD(project, rnaDna_type) and config["experiments"][project]["split_number"] > 1


def getExperimentReads(
    read_type: str,
    project: str,
    condition: str,
    replicate: str,
    rnaDna_type: str,
    check_splitting: bool,
    check_trimming: bool,
):
    if read_type == "FWD":
        return getFWD(project, condition, replicate, rnaDna_type, check_splitting, check_trimming)
    elif read_type == "REV":
        return getREV(project, condition, replicate, rnaDna_type, check_splitting, check_trimming)
    elif read_type == "UMI":
        return getUMI(project, condition, replicate, rnaDna_type, check_splitting, check_trimming)
    else:
        raise ValueError("read_type must be one of FWD, REV or UMI")


def getFWD(
    project: str, condition: str, replicate: str, rnaDna_type: str, check_splitting: bool, check_trimming: bool
) -> list[str]:
    if check_trimming and useTrimming(project, "FWD"):
        if check_splitting and useSplitting(project, rnaDna_type):
            return ["results/experiments/{project}/fastq/FWD.trimmed.{condition}.{replicate}.{type}.{split}.fastq.gz"]
        else:
            return ["results/experiments/{project}/fastq/FWD.trimmed.{condition}.{replicate}.{type}.0.fastq.gz"]

    elif check_splitting and useSplitting(project, rnaDna_type):
        return ["results/experiments/{project}/fastq/FWD.split.{condition}.{replicate}.{type}.{split}.fastq.gz"]
    else:
        exp = getExperiments(project)
        exp = exp[exp.Condition == condition]
        exp = exp[exp.Replicate.astype(str) == replicate]
        return [
            "%s/%s" % (config["experiments"][project]["data_folder"], f)
            for f in exp["%s_BC_F" % rnaDna_type].iloc[0].split(";")
        ]


def getFWDWithIndex(project: str) -> list[str]:
    return [
        "%s/%s" % (config["experiments"][project]["data_folder"], f) for f in getExperiments(project).BC_F.iloc[0].split(";")
    ]


def getREV(
    project: str, condition: str, replicate: str, rnaDna_type: str, check_splitting: bool, check_trimming: bool
) -> list[str]:
    if check_trimming and useTrimming(project, "REV"):
        if check_splitting and useSplitting(project, rnaDna_type):
            return ["results/experiments/{project}/fastq/REV.trimmed.{condition}.{replicate}.{type}.{split}.fastq.gz"]
        else:
            return ["results/experiments/{project}/fastq/REV.trimmed.{condition}.{replicate}.{type}.0.fastq.gz"]
    elif check_splitting and useSplitting(project, rnaDna_type):
        return ["results/experiments/{project}/fastq/REV.split.{condition}.{replicate}.{type}.{split}.fastq.gz"]
    else:
        exp = getExperiments(project)
        exp = exp[exp.Condition == condition]
        exp = exp[exp.Replicate.astype(str) == replicate]
        return [
            "%s/%s" % (config["experiments"][project]["data_folder"], f)
            for f in exp["%s_BC_R" % rnaDna_type].iloc[0].split(";")
        ]


def getREVWithIndex(project: str) -> list[str]:
    return [
        "%s%s" % (config["experiments"][project]["data_folder"], f) for f in getExperiments(project).BC_R.iloc[0].split(";")
    ]


def getUMI(
    project: str, condition: str, replicate: str, rnaDna_type: str, check_splitting: bool, check_trimming: bool
) -> list[str]:
    if check_trimming and useTrimming(project, "UMI"):
        if check_splitting and useSplitting(project, rnaDna_type):
            return ["results/experiments/{project}/fastq/UMI.trimmed.{condition}.{replicate}.{type}.{split}.fastq.gz"]
        else:
            return ["results/experiments/{project}/fastq/UMI.trimmed.{condition}.{replicate}.{type}.0.fastq.gz"]
    elif check_splitting and useSplitting(project, rnaDna_type):
        return ["results/experiments/{project}/fastq/UMI.split.{condition}.{replicate}.{type}.{split}.fastq.gz"]
    else:
        exp = getExperiments(project)
        exp = exp[exp.Condition == condition]
        exp = exp[exp.Replicate.astype(str) == replicate]
        return [
            "%s/%s" % (config["experiments"][project]["data_folder"], f)
            for f in exp["%s_UMI" % rnaDna_type].iloc[0].split(";")
        ]


def getUMIWithIndex(project: str) -> list[str]:
    return [config["experiments"][project]["data_folder"] + f for f in getExperiments(project).UMI.iloc[0].split(";")]


def getIndexWithIndex(project: str) -> list[str]:
    return [config["experiments"][project]["data_folder"] + f for f in getExperiments(project).INDEX.iloc[0].split(";")]


def getUMIBamFile(project: str) -> str:
    """
    Helper to get the correct BAM file (demultiplexed or not)
    """

    if config["experiments"][project]["demultiplex"]:
        return "results/experiments/{project}/counts/merged_demultiplex.{condition}.{replicate}.{type}.bam"
    else:
        return "results/experiments/{project}/counts/useUMI.{condition}.{replicate}.{type}.{split}.bam"


def getRawCounts(project: str, dnaRNA_type: str) -> str:
    """
    Helper to get the correct raw counts file (umi/noUMI or just FWD read)
    """
    if useUMI(project, dnaRNA_type):
        if onlyFWD(project, dnaRNA_type):
            return (
                "results/experiments/{project}/counts/onlyFWDUMI.{condition}.{replicate}.%s.raw_counts.tsv.gz" % dnaRNA_type,
            )
        else:
            return "results/experiments/{project}/counts/useUMI.{condition}.{replicate}.%s.raw_counts.tsv.gz" % dnaRNA_type
    elif noUMI(project, dnaRNA_type):
        return "results/experiments/{project}/counts/noUMI.{condition}.{replicate}.%s.raw_counts.tsv.gz" % dnaRNA_type
    elif onlyFWD(project, dnaRNA_type):
        return "results/experiments/{project}/counts/onlyFWD.{condition}.{replicate}.%s.raw_counts.tsv.gz" % dnaRNA_type
    else:
        raise RuntimeError("Error in getRawCounts: no valid option for %s and %s found" % (project, dnaRNA_type))


def counts_aggregate_demultiplex_input(project):
    output = []
    conditions = getConditions(project)
    for condition in conditions:
        replicates = getReplicatesOfCondition(project, condition)
        names = expand(
            "{condition}.{replicate}.{type}",
            condition=condition,
            replicate=replicates,
            type=["DNA", "RNA"],
        )
        for name in names:
            with checkpoints.experiment_counts_demultiplex_BAM_umi.get(project=project, name=name).output[0].open() as f:
                output += [f.name]
    return output


def counts_getFilterConfig(project, conf, dna_or_rna, command):
    value = config["experiments"][project]["configs"][conf]["filter"]["min_%s_counts" % dna_or_rna.lower()]
    filterMap = {"min_counts": "minCounts"}
    if isinstance(value, int):
        return "--%s %d" % (filterMap.get(command, command), value)
    else:
        return "--%s %f" % (filterMap.get(command, command), value)


def counts_getSamplingConfig(project, conf, dna_or_rna, command):
    if useSampling(project, conf, dna_or_rna):
        if dna_or_rna in config["experiments"][project]["configs"][conf]["sampling"]:
            if command in config["experiments"][project]["configs"][conf]["sampling"][dna_or_rna]:
                value = config["experiments"][project]["configs"][conf]["sampling"][dna_or_rna][command]
                if isinstance(value, int):
                    return "--%s %d" % (command, value)
                else:
                    return "--%s %f" % (command, value)

    return ""


def getFinalCounts(project, conf, condition, rna_or_dna, raw_or_assigned):
    output = ""

    replicates = getReplicatesOfConditionType(project, condition, rna_or_dna)
    if len(replicates) > 1:
        replicate = "{replicate}"
    else:
        replicate = replicates[0]

    if raw_or_assigned == "counts":
        if useSampling(project, conf, rna_or_dna):
            output = "results/experiments/{project}/%s/{condition}.%s.%s.final_counts.sampling.{config}.tsv.gz" % (
                raw_or_assigned,
                replicate,
                rna_or_dna,
            )

        else:
            output = "results/experiments/{project}/%s/{condition}.%s.%s.final_counts.tsv.gz" % (
                raw_or_assigned,
                replicate,
                rna_or_dna,
            )
    else:
        output = "results/experiments/{project}/%s/{condition}.%s.%s.final_counts.config.{config}.tsv.gz" % (
            raw_or_assigned,
            replicate,
            rna_or_dna,
        )
    return output


# assigned_counts.smk specific functions


def getAssignmentFile(project, assignment):
    if config["experiments"][project]["assignments"][assignment]["type"] == "file":
        return config["experiments"][project]["assignments"][assignment]["assignment_file"]
    if config["experiments"][project]["assignments"][assignment]["type"] == "config":
        conf = config["experiments"][project]["assignments"][assignment]["assignment_config"]
        name = config["experiments"][project]["assignments"][assignment]["assignment_name"]
        return expand(
            "results/assignment/{assignment}/assignment_barcodes.{config}.tsv.gz",
            assignment=name,
            config=conf,
        )


def assignedCounts_getAssignmentSamplingConfig(project, assignment, command):
    if "sampling" in config["experiments"][project]["assignments"][assignment]:
        if command in config["experiments"][project]["assignments"][assignment]["sampling"]:
            value = config["experiments"][project]["assignments"][assignment]["sampling"][command]
            if isinstance(value, int):
                return "--%s %d" % (command, value)
            else:
                return "--%s %f" % (command, value)

    return ""


# statistic.smk specific functions


# get all barcodes of experiment (rule statistic_BC_in_RNA_DNA)
def getBCinRNADNAStats(wc):
    exp = getExperiments(wc.project)
    output = []
    for index, row in exp.iterrows():
        output += expand(
            "results/experiments/{project}/statistic/counts/{condition}.{replicate}.{countType}_BC_in_RNA_DNA.tsv.gz",
            project=wc.project,
            condition=row["Condition"],
            replicate=row["Replicate"],
            countType=wc.countType,
        )
    return output


def getAssignedCountsStatistic(project, assignment, conf, condition):
    exp = getExperiments(project)
    exp = exp[exp.Condition == condition]
    output = []
    for index, row in exp.iterrows():
        output += [
            "--statistic %s results/experiments/%s/statistic/assigned_counts/%s/%s/%s.%s.merged_assigned_counts.statistic.tsv.gz"
            % (
                str(row["Replicate"]),
                project,
                assignment,
                conf,
                condition,
                str(row["Replicate"]),
            )
        ]
    return output


# get all barcodes of experiment (rule dna_rna_merge_counts_withoutZeros or rule dna_rna_merge_counts_withZeros)
def getMergedCounts(project, raw_or_assigned, condition, conf):
    exp = getExperiments(project)
    exp = exp[exp.Condition == condition]
    files = []
    replicates = []
    for index, row in exp.iterrows():
        files += expand(
            "results/experiments/{project}/{raw_or_assigned}/{condition}.{replicate}.merged.config.{config}.tsv.gz",
            raw_or_assigned=raw_or_assigned,
            project=project,
            condition=condition,
            replicate=row["Replicate"],
            config=conf,
        )
        replicates += [str(row["Replicate"])]
    return [files, replicates]


# variants.smk specific functions


def getVariantsBCThreshold(project):
    return getVariants(project)["min_barcodes"]
