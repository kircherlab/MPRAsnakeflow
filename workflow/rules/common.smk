################################
#### Global functions       ####
################################

SCRIPTS_DIR = "../scripts"
ENVS_DIR = "../envs"


def getWorkflowFile(dir_name, name):
    return workflow.source_path("%s/%s" % (dir_name, name))


def getScript(name):
    return getWorkflowFile(SCRIPTS_DIR, name)


def getCondaEnv(name):
    return getWorkflowFile(ENVS_DIR, name)


##### load config and sample sheets #####
from snakemake.utils import validate
import pandas as pd

# Workaround: validate() is broken from Snakemake 9.5.1 to snakemake 9.14.7 in remote jobs
if version.parse(snakemake.__version__) >= version.parse("9.5.1") and version.parse(
    snakemake.__version__
) <= version.parse("9.14.7"):
    from snakemake_interface_executor_plugins.settings import ExecMode

    # Use the global 'workflow' variable directly as recommended by Snakemake
    if workflow.remote_exec:
        old_exec_mode = workflow.exec_mode
        workflow.workflow_settings.exec_mode = ExecMode.DEFAULT
        validate(config, schema="../schemas/config.schema.yaml")
        workflow.workflow_settings.exec_mode = old_exec_mode
    else:
        validate(config, schema="../schemas/config.schema.yaml")
else:
    validate(config, schema="../schemas/config.schema.yaml")

# load sample sheets
experiments = {}
if "experiments" in config:
    for project in config["experiments"]:
        experiment = pd.read_csv(config["experiments"][project]["experiment_file"])
        validate(experiment, schema="../schemas/experiment_file.schema.yaml")
        experiments[project] = experiment

# validate version of config with MPRAsnakeflow version


def check_version(config_version: version.Version):
    version_check_fail = False
    if __version__.major == 0 and config_version.major == 0:
        # for major version 0, only check minor version
        if __version__.minor != config_version.minor:
            version_check_fail = True
    else:
        if __version__.major != config_version.major:
            version_check_fail = True
    if version_check_fail:
        raise ValueError(
            f"\033[38;2;255;165;0mVersion mismatch: MPRAsnakeflow version is {__version__.public}, but config version is {config_version.public}\033[0m"
        )


if not config["skip_version_check"]:
    check_version(version.parse(config["version"]))


# modify config based on certain rules


def modify_config(config):
    # update sequence length with when adapters are added via strand sensitive is enabled
    if "assignments" in config:
        for assignment in config["assignments"].keys():
            if config["assignments"][assignment]["strand_sensitive"]["enable"]:
                add_length = len(
                    config["assignments"][assignment]["strand_sensitive"][
                        "forward_adapter"
                    ]
                ) + len(
                    config["assignments"][assignment]["strand_sensitive"][
                        "reverse_adapter"
                    ]
                )
                if config["assignments"][assignment]["alignment_tool"]["tool"] == "bwa":
                    config["assignments"][assignment]["alignment_tool"]["configs"][
                        "sequence_length"
                    ]["min"] += add_length
                    config["assignments"][assignment]["alignment_tool"]["configs"][
                        "sequence_length"
                    ]["max"] += add_length
                else:
                    config["assignments"][assignment]["alignment_tool"]["configs"][
                        "sequence_length"
                    ] += add_length
    return config


config = modify_config(config)
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


##### get helpers for different things (like conditions etc) #####


def getAssignments(match_methods=None):
    """
    Retrieve assignment configurations, optionally filtered by alignment tool method.

    Args:
        match_methods (list, optional): List of alignment tool names to filter by.
            If provided, only assignments using these tools are returned.
            Defaults to None.

    Returns:
        list: List of assignment keys/names. If match_methods is specified,
            returns only assignments whose alignment tool is in match_methods.
            Returns empty list if "assignments" key is not in config.

    Examples:
        - getAssignments() -> returns all assignment keys from config
        - getAssignments(["bowtie2", "bwa"]) -> returns assignments using bowtie2 or bwa
    """
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


def getReplicatesOfConditionType(project, condition, rna_or_dna):
    exp = getExperiments(project)

    replicates = getReplicatesOfCondition(project, condition)

    if f"{rna_or_dna}_BC_F" in exp.columns:

        exp_filter = exp[exp.Condition == condition]

        if len(replicates) > 1 and exp_filter[f"{rna_or_dna}_BC_F"].nunique() == 1:
            return [replicates[0]]

    return replicates


def hasReplicates(project, condition=None):
    if condition == None:
        conditions = getConditions(project)
        for condition in conditions:
            if len(getReplicatesOfCondition(project, condition)) <= 1:
                return False
    else:
        return len(getReplicatesOfCondition(project, condition)) > 1
    return True


def getReplicatesOfConditionType(project, condition, rna_or_dna):
    exp = getExperiments(project)

    replicates = getReplicatesOfCondition(project, condition)

    if f"{rna_or_dna}_BC_F" in exp.columns:

        exp_filter = exp[exp.Condition == condition]

        if len(replicates) > 1 and exp_filter[f"{rna_or_dna}_BC_F"].nunique() == 1:
            return [replicates[0]]

    return replicates


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


def reverse_complement(seq):
    complementary = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(reversed([complementary[i] for i in seq]))


def getSplitNumber():
    splits = [1]

    for assignment in getAssignments():
        splits += [config["assignments"][assignment]["alignment_tool"]["split_number"]]

    return max(splits)
