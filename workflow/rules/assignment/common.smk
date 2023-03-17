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


def getAssignmentRead(assignment, read):
    """
    Return the correct assignment read.
    """
    if hasBCRead(assignment) or read == "REV":
        return config["assignments"][assignment][read]
    elif hasLinker(assignment):
        return "results/assignment/{assignment}/fastq/{read}.byCutadapt.fastq.gz"
    elif hasLinkerLength(assignment):
        return "results/assignment/{assignment}/fastq/{read}.byLength.fastq.gz"
    else:
        raise RuntimeError(
            "Wrong assignment configuration. Cannot find corerct combinations of reads for assignment %s"
            % assignment
        )
