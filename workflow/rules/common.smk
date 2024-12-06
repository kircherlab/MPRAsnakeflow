################################
#### Global functions       ####
################################

SCRIPTS_DIR = "../scripts"


def getScript(name):
    return workflow.source_path("%s/%s" % (SCRIPTS_DIR, name))


##### load config and sample sheets #####
from snakemake.utils import validate
import pandas as pd

validate(config, schema="../schemas/config.schema.yaml")

# load sample sheets
experiments = {}
if "experiments" in config:
    for project in config["experiments"]:
        experiment = pd.read_csv(config["experiments"][project]["experiment_file"])
        validate(experiment, schema="../schemas/experiment_file.schema.yaml")
        experiments[project] = experiment

# validate version of config with MPRAsnakeflow version

import re

# Regular expression to match the first two digits with the dot in the middle
pattern_major_version = r"^(\d+)"
pattern_development_version = r"^(0(\.\d+)?)"


def check_version(pattern, version, config_version):
    # Search for the pattern in the string
    match_version = re.search(pattern, version)

    match_config = re.search(pattern, config_version)

    # Check if a match is found and print the result
    if match_version and match_config:
        if match_version.group(1) != match_config.group(1):
            raise ValueError(
                f"\033[38;2;255;165;0mVersion mismatch: MPRAsnakeflow version is {version}, but config version is {config_version}\033[0m"
            )


if not config["skip_version_check"]:
    check_version(pattern_development_version, version, config["version"])
    check_version(pattern_major_version, version, config["version"])


################################
#### HELPERS AND EXCEPTIONS ####
################################


##### Exceptions #####
class MissingAssignmentInConfigException(Exception):
    """
    Exception raised for if no assignment file is set in the config.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        config_name (string): name of the configuration which assignment is missing.
    """

    def __init__(self, config_name):
        self.config_name = config_name

    def __str__(self):
        return "Config %s has no assignment file defined!" % (self.config_name)


class MissingVariantInConfigException(Exception):
    """
    Exception raised for if no variants config.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        config_name (string): name of the configuration which assignment is missing.
    """

    def __init__(self, config_name):
        self.config_name = config_name

    def __str__(self):
        return "Config %s has no variants defined!" % (self.config_name)


##### get helpers for different things (like Conditions etc) #####
def getAssignments(match_methods=None):
    if "assignments" in config:
        if match_methods:
            output = []
            for assignment in config["assignments"]:
                if (
                    config["assignments"][assignment]["alignment_tool"]["tool"]
                    in match_methods
                ):
                    output.append(assignment)
            return output
        else:
            return list(config["assignments"].keys())
    else:
        return []


def getAssignmentFile(project, assignment):
    if config["experiments"][project]["assignments"][assignment]["type"] == "file":
        return config["experiments"][project]["assignments"][assignment][
            "assignment_file"
        ]
    if config["experiments"][project]["assignments"][assignment]["type"] == "config":
        conf = config["experiments"][project]["assignments"][assignment][
            "assignment_config"
        ]
        name = config["experiments"][project]["assignments"][assignment][
            "assignment_name"
        ]
        return expand(
            "results/assignment/{assignment}/assignment_barcodes.{config}.tsv.gz",
            assignment=name,
            config=conf,
        )


def getProjects():
    if "experiments" in config:
        return list(config["experiments"].keys())
    else:
        return []


def getExperiments(project):
    return experiments[project]


def getConditions(project):
    exp = getExperiments(project)
    return list(exp.Condition.unique())


def getProjectAssignments(project):
    if (
        "assignments" in config["experiments"][project]
        and len(config["experiments"][project]["assignments"]) > 0
    ):
        return list(config["experiments"][project]["assignments"].keys())
    else:
        raise MissingAssignmentInConfigException(project)


def getVariants(project):
    if "variants" in config["experiments"][project]:
        return config["experiments"][project]["variants"]
    else:
        raise MissingVariantInConfigException(project)


def getReplicatesOfCondition(project, condition):
    exp = getExperiments(project)
    exp = exp[exp.Condition == condition]
    return list(exp.Replicate.astype(str))


def getVariantsBCThreshold(project):
    return getVariants(project)["min_barcodes"]


def getFW(project, condition, replicate, rnaDna_type):
    exp = getExperiments(project)
    exp = exp[exp.Condition == condition]
    exp = exp[exp.Replicate.astype(str) == replicate]
    return [
        "%s/%s" % (config["experiments"][project]["data_folder"], f)
        for f in exp["%s_BC_F" % rnaDna_type].iloc[0].split(";")
    ]


def getFWWithIndex(project):
    return [
        "%s/%s" % (config["experiments"][project]["data_folder"], f)
        for f in getExperiments(project).BC_F.iloc[0].split(";")
    ]


def getRev(project, condition, replicate, rnaDna_type):
    exp = getExperiments(project)
    exp = exp[exp.Condition == condition]
    exp = exp[exp.Replicate.astype(str) == replicate]
    return [
        "%s/%s" % (config["experiments"][project]["data_folder"], f)
        for f in exp["%s_BC_R" % rnaDna_type].iloc[0].split(";")
    ]


def getRevWithIndex(project):
    return [
        "%s%s" % (config["experiments"][project]["data_folder"], f)
        for f in getExperiments(project).BC_R.iloc[0].split(";")
    ]


def getUMI(project, condition, replicate, rnaDna_type):
    exp = getExperiments(project)
    exp = exp[exp.Condition == condition]
    exp = exp[exp.Replicate.astype(str) == replicate]
    return [
        "%s/%s" % (config["experiments"][project]["data_folder"], f)
        for f in exp["%s_UMI" % rnaDna_type].iloc[0].split(";")
    ]


def getUMIWithIndex(project):
    return [
        config["experiments"][project]["data_folder"] + f
        for f in getExperiments(project).UMI.iloc[0].split(";")
    ]


def getIndexWithIndex(project):
    return [
        config["experiments"][project]["data_folder"] + f
        for f in getExperiments(project).INDEX.iloc[0].split(";")
    ]


def hasReplicates(project, condition=None):
    if condition == None:
        conditions = getConditions(project)
        for condition in conditions:
            if len(getReplicatesOfCondition(project, condition)) <= 1:
                return False
    else:
        return len(getReplicatesOfCondition(project, condition)) > 1
    return True


def getConfigs(project):
    return list(config["experiments"][project]["configs"].keys())


##### Helper to create output files #####
def getOutputProject_helper(files, betweenReplicates=False):
    """
    Inserts {project} from config into given file.
    When betweenReplicates is True skips projects without replicates in one condition.
    """
    output = []
    projects = getProjects()
    for project in projects:
        if not betweenReplicates or hasReplicates(project):
            for file in files:
                output += expand(
                    file,
                    project=project,
                )
    return output


def getOutputConditionReplicateType_helper(files, project, skip={}):
    """
    Inserts {condition}, {replicate} and {type} from config into given file.
    Can skip projects with the given config set by skip.
    """
    output = []

    for key, value in skip.items():
        if config["experiments"][project][key] == value:
            return []
    conditions = getConditions(project)
    for condition in conditions:
        for file in files:
            for type in ["DNA", "RNA"]:
                replicates = getReplicatesOfConditionType(project, condition, type)
                output += expand(
                    file,
                    project=project,
                    condition=condition,
                    replicate=replicates,
                    type=type,
                )
    return output


def getOutputProjectConditionReplicateType_helper(files, skip={}):
    """
    Inserts {project}, {condition}, {replicate} and {type} from config into given file.
    Can skip projects with the given config set by skip.
    """
    output = []
    projects = getProjects()
    for project in projects:
        # skip projects with the following config
        for file in files:
            fs = (
                expand(
                    file,
                    project=project,
                    condition="{condition}",
                    replicate="{replicate}",
                    type="{type}",
                ),
            )
            output += getOutputConditionReplicateType_helper(
                fs,
                project,
                skip,
            )
    return output


def getOutputProjectConditionConfigType_helper(files):
    """
    Inserts {project}, {condition} and {type} from config into given file.
    """
    output = []
    projects = getProjects()
    for project in projects:
        conditions = getConditions(project)
        for condition in conditions:
            for file in files:
                output += expand(
                    file,
                    project=project,
                    condition=condition,
                    config=getConfigs(project),
                    type=["DNA", "RNA"],
                )
    return output


def getOutputProjectConditionType_helper(file):
    """
    Inserts {project}, {condition} and {type} from config into given file.
    """
    output = []
    projects = getProjects()
    for project in projects:
        conditions = getConditions(project)
        for condition in conditions:
            for file in files:
                output += expand(
                    file,
                    project=project,
                    condition=condition,
                    type=["DNA", "RNA"],
                )
    return output


def getOutputProjectConditionAssignmentConfigType_helper(files):
    """
    Inserts {project}, {condition}, {assignment} and {config} (from configs of project) from config into given file.
    """
    output = []
    projects = getProjects()
    for project in projects:
        try:
            conditions = getConditions(project)
            for condition in conditions:
                for file in files:
                    output += expand(
                        file,
                        project=project,
                        condition=condition,
                        assignment=getProjectAssignments(project),
                        config=getConfigs(project),
                        type=["DNA", "RNA"],
                    )
        except MissingAssignmentInConfigException:
            continue
    return output


def getOutputProjectConditionAssignmentConfig_helper(files):
    """
    Inserts {project}, {condition}, {assignment} and {config} (from configs of project) from config into given file.
    """
    output = []
    projects = getProjects()
    for project in projects:
        try:
            conditions = getConditions(project)
            for condition in conditions:
                for file in files:
                    output += expand(
                        file,
                        project=project,
                        condition=condition,
                        assignment=getProjectAssignments(project),
                        config=getConfigs(project),
                    )
        except MissingAssignmentInConfigException:
            continue
    return output


def getOutputProjectConditionAssignmentConfigThreshold_helper(files):
    """
    Inserts {project}, {condition}, {assignment} {config} (from configs of project) and Threshold from config into given file.
    """
    output = []
    projects = getProjects()
    for project in projects:
        try:
            conditions = getConditions(project)
            for condition in conditions:
                for conf in getConfigs(project):
                    threshold = config["experiments"][project]["configs"][conf][
                        "filter"
                    ]["bc_threshold"]
                    for file in files:
                        output += expand(
                            file,
                            project=project,
                            condition=condition,
                            assignment=getProjectAssignments(project),
                            config=conf,
                            threshold=threshold,
                        )
        except MissingAssignmentInConfigException:
            continue
    return output


def getOutputProjectAssignmentConfig_helper(files, betweenReplicates=False):
    """
    Inserts {project}, {assignment} and {config} (from configs of project) from config into given file.
    When betweenReplicates is True skips projects without replicates in one condition.
    """
    output = []
    projects = getProjects()
    for project in projects:
        if not betweenReplicates or hasReplicates(project):
            for file in files:
                try:
                    output += expand(
                        file,
                        project=project,
                        assignment=getProjectAssignments(project),
                        config=getConfigs(project),
                    )
                except MissingAssignmentInConfigException:
                    continue
    return output


def getOutputProjectConfig_helper(files, betweenReplicates=False):
    """
    Inserts {project}, {config} from config into given file.
    When betweenReplicates is True skips projects without replicates in one condition.
    """
    output = []
    projects = getProjects()
    for project in projects:
        if not betweenReplicates or hasReplicates(project):
            for file in files:
                output += expand(
                    file,
                    project=project,
                    config=getConfigs(project),
                )
    return output


def getOutputVariants_helper(files, betweenReplicates=False):
    """
    Only when variants are set in config file
    Inserts {project}, {condition}, {assignment} and {config} (from configs of project) from config into given file.
    When betweenReplicates is True skips project/condition without replicates in a condition.
    """
    output = []
    projects = getProjects()
    for project in projects:
        conditions = getConditions(project)
        for condition in conditions:
            if "variants" in config["experiments"][project]:
                if hasReplicates(project, condition):
                    for file in files:
                        output += expand(
                            file,
                            project=project,
                            condition=condition,
                            assignment=getProjectAssignments(project),
                            config=list(
                                config["experiments"][project]["configs"].keys()
                            ),
                        )
    return output


def getAssignment_helper(files, match_methods=None):
    output = []
    for file in files:
        output += expand(file, assignment=getAssignments(match_methods))
    return output


def getAssignmentConfig_helper(files):
    output = []
    for assignment in getAssignments():
        for file in files:
            output += expand(
                file,
                assignment=assignment,
                config=config["assignments"][assignment]["configs"].keys(),
            )
    return output


# config functions


def useSampling(project, conf, dna_or_rna):
    return (
        "sampling" in config["experiments"][project]["configs"][conf]
        and dna_or_rna in config["experiments"][project]["configs"][conf]["sampling"]
    )


def withoutZeros(project, conf):
    return (
        config["experiments"][project]["configs"][conf]["filter"]["min_dna_counts"] > 0
        and config["experiments"][project]["configs"][conf]["filter"]["min_rna_counts"]
        > 0
    )


# assignment.smk specific functions


def getSplitNumber():
    splits = [1]

    for assignment in getAssignments():
        splits += [config["assignments"][assignment]["alignment_tool"]["split_number"]]

    return max(splits)


# count.smk specific functions


def getUMIBamFile(project, condition, replicate, type):
    """
    gelper to get the correct BAM file (demultiplexed or not)
    """
    if config["experiments"][project]["demultiplex"]:
        return "results/%s/counts/merged_demultiplex.%s_%s_%s.bam" % (
            project,
            condition,
            replicate,
            type,
        )
    else:
        return "results/experiments/%s/counts/useUMI.%s_%s_%s.bam" % (
            project,
            condition,
            replicate,
            type,
        )


def useUMI(project):
    """
    helper to check if UMI should be used
    """
    return "UMI" in experiments[project] or "DNA_UMI" in experiments[project]


def noUMI(project):
    """
    helper to check if UMI should not be used
    """
    return (
        "UMI" not in experiments[project]
        and "DNA_UMI" not in experiments[project]
        and "DNA_BC_R" in experiments[project]
    )


def onlyFWByLength(project):
    """
    helper to check if only forward reads should be used (length option)
    """
    return (
        "UMI" not in experiments[project]
        and "DNA_BC_R" not in experiments[project]
        and "adapter" not in config["experiments"][project]
    )


def onlyFWbyCutadapt(project):
    """
    helper to check if only forward reads should be used (cutadapt option)
    """
    return (
        "UMI" not in experiments[project]
        and "DNA_BC_R" not in experiments[project]
        and "adapter" in config["experiments"][project]
    )


def getRawCounts(project, type):
    """
    Helper to get the correct raw counts file (umi/noUMI or just FW read)
    """
    if useUMI(project):
        return (
            "results/experiments/{project}/counts/useUMI.{condition}_{replicate}_%s_raw_counts.tsv.gz"
            % type
        )
    elif noUMI(project):
        return (
            "results/experiments/{project}/counts/noUMI.{condition}_{replicate}_%s_raw_counts.tsv.gz"
            % type
        )
    elif onlyFWByLength(project):
        return (
            "results/experiments/{project}/counts/onlyFWByLength.{condition}_{replicate}_%s_raw_counts.tsv.gz"
            % type
        )
    elif onlyFWbyCutadapt(project):
        return (
            "results/experiments/{project}/counts/onlyFWByCutadapt.{condition}_{replicate}_%s_raw_counts.tsv.gz"
            % type
        )
    else:
        raise RuntimeError(
            "Error in getRawCounts: no valid option for %s and %s found"
            % (project, type)
        )


def counts_aggregate_demultiplex_input(project):
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
            with checkpoints.counts_demultiplex_BAM_umi.get(
                project=project, name=name
            ).output[0].open() as f:
                output += [f.name]
    return output


def counts_getFilterConfig(project, conf, dna_or_rna, command):
    value = config["experiments"][project]["configs"][conf]["filter"][
        "min_%s_counts" % dna_or_rna.lower()
    ]
    filterMap = {"min_counts": "minCounts"}
    if isinstance(value, int):
        return "--%s %d" % (filterMap.get(command, command), value)
    else:
        return "--%s %f" % (filterMap.get(command, command), value)


def counts_getSamplingConfig(project, conf, dna_or_rna, command):
    if useSampling(project, conf, dna_or_rna):
        if dna_or_rna in config["experiments"][project]["configs"][conf]["sampling"]:
            if (
                command
                in config["experiments"][project]["configs"][conf]["sampling"][
                    dna_or_rna
                ]
            ):
                value = config["experiments"][project]["configs"][conf]["sampling"][
                    dna_or_rna
                ][command]
                if isinstance(value, int):
                    return "--%s %d" % (command, value)
                else:
                    return "--%s %f" % (command, value)

    return ""


def getReplicatesOfConditionType(project, condition, rna_or_dna):
    exp = getExperiments(project)

    replicates = getReplicatesOfCondition(project, condition)

    if f"{rna_or_dna}_BC_F" in exp.columns:

        exp_filter = exp[exp.Condition == condition]

        if len(replicates) > 1 and exp_filter[f"{rna_or_dna}_BC_F"].nunique() == 1:
            return [replicates[0]]

    return replicates


def getFinalCounts(project, conf, condition, rna_or_dna, raw_or_assigned):
    output = ""

    replicates = getReplicatesOfConditionType(project, condition, rna_or_dna)
    if len(replicates) > 1:
        replicate = "{replicate}"
    else:
        replicate = replicates[0]

    if raw_or_assigned == "counts":
        if useSampling(project, conf, rna_or_dna):
            output = (
                "results/experiments/{project}/%s/{condition}_%s_%s_final_counts.sampling.{config}.tsv.gz"
                % (raw_or_assigned, replicate, rna_or_dna)
            )

        else:
            output = (
                "results/experiments/{project}/%s/{condition}_%s_%s_final_counts.tsv.gz"
                % (raw_or_assigned, replicate, rna_or_dna)
            )
    else:
        output = (
            "results/experiments/{project}/%s/{condition}_%s_%s_final_counts.config.{config}.tsv.gz"
            % (raw_or_assigned, replicate, rna_or_dna)
        )
    return output


# assigned_counts.smk specific functions


def assignedCounts_getAssignmentSamplingConfig(project, assignment, command):
    if "sampling" in config["experiments"][project]["assignments"][assignment]:
        if (
            command
            in config["experiments"][project]["assignments"][assignment]["sampling"]
        ):
            value = config["experiments"][project]["assignments"][assignment][
                "sampling"
            ][command]
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
            "results/experiments/{project}/statistic/counts/{condition}_{replicate}_{countType}_BC_in_RNA_DNA.tsv.gz",
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
            "--statistic %s results/experiments/%s/statistic/assigned_counts/%s/%s/%s_%s_merged_assigned_counts.statistic.tsv.gz"
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
            "results/experiments/{project}/{raw_or_assigned}/{condition}_{replicate}.merged.config.{config}.tsv.gz",
            raw_or_assigned=raw_or_assigned,
            project=project,
            condition=condition,
            replicate=row["Replicate"],
            config=conf,
        )
        replicates += [str(row["Replicate"])]
    return [files, replicates]
