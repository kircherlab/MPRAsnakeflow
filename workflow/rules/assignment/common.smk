"""
Common functions used by the assignment workflow.
"""

assignment_bwa_dicts = ["bwt", "sa", "pac", "ann", "amb"]


def hasBCRead(assignment):
    """
    Return True if the assignment has a BC read.
    """
    return "BC" in config["assignments"][assignment]


def hasLinker(assignment):
    """
    Return True if the assignment contains a linker.
    """
    return "linker" in config["assignments"][assignment]


def hasLinkerLength(assignment):
    """
    Return True if the assignment contains a linker length.
    """
    return "linker_length" in config["assignments"][assignment]


def hasOnlyForwardRead(assignment):
    """
    Return True if the assignment contains only a forward read.
    """
    return "REV" not in config["assignments"][assignment]


def useAssignmentAdapterTrimming(assignment, read):
    """
    Return True if adapter trimming should be used for the given read in the assignment.
    """
    return (
        "adapters" in config["assignments"][assignment]
        and read in config["assignments"][assignment]["adapters"]
    )


def getAssignmentRead(assignment, read):
    """
    Return the correct assignment read.
    """
    if hasBCRead(assignment) or read == "REV":
        return (
            "results/assignment/{assignment}/fastq/{read}.trimmed.fastq.gz"
            if useAssignmentAdapterTrimming(assignment, read)
            else config["assignments"][assignment][read]
        )
    elif hasLinker(assignment):
        return "results/assignment/{assignment}/fastq/{read}.byCutadapt.fastq.gz"
    elif hasLinkerLength(assignment):
        return "results/assignment/{assignment}/fastq/{read}.byLength.fastq.gz"
    else:
        raise RuntimeError(
            "Wrong assignment configuration. Cannot find corerct combinations of reads for assignment %s"
            % assignment
        )


def getMappingRead(assignment: str) -> str:
    """
    Return the start read for the assignment. Can be a joined read or only the forward read.
    """
    if hasOnlyForwardRead(assignment):
        return "results/assignment/{assignment}/fastq/splits/FWD.split{split}.BCattached.fastq.gz"
    else:
        return "results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz"


def getAssignmentCutadaptAdapters(assignment, read):
    if (
        "adapters" in config["assignments"][assignment]
        and read in config["assignments"][assignment]["adapters"]
    ):
        adapters_config = config["assignments"][assignment]["adapters"][read]
        if isinstance(adapters_config is list) and isinstance(adapters_config[0], int):
            return " ".join(["-u %d" % u for u in adapters_config])
        else:
            return " ".join(
                [
                    (
                        "-g %s" % adapter[0]
                        if adapter[1] == "5prime"
                        else "-a %s" % adapter[0]
                    )
                    for adapter in adapters_config
                ]
            )
    else:
        return ""
