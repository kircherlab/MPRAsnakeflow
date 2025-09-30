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


def hasAdapters(assignment):
    """
    Return True if the assignment contains a 3 or 5 adapter for removal.
    """
    return has3PrimeAdapters(assignment) or has5PrimeAdapters(assignment)


def has3PrimeAdapters(assignment):
    """
    Return True if the assignment contains a 3' adapter for removal.
    """
    return (
        "adapters" in config["assignments"][assignment]
        and "3prime" in config["assignments"][assignment]["adapters"]
    )


def has5PrimeAdapters(assignment):
    """
    Return True if the assignment contains a 5' adapter for removal.
    """
    return (
        "adapters" in config["assignments"][assignment]
        and "5prime" in config["assignments"][assignment]["adapters"]
    )


def hasOnlyForwardRead(assignment):
    """
    Return True if the assignment contains only a forward read.
    """
    return "REV" not in config["assignments"][assignment]


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


def getMappingRead(assignment: str) -> str:
    """
    Return the final reads for mapping after joining, maybe after adapter removal.
    """
    if has5PrimeAdapters(assignment):
        return (
            "results/assignment/{assignment}/fastq/merge_split{split}.5prime.fastq.gz"
        )
    elif has3PrimeAdapters(assignment):
        return (
            "results/assignment/{assignment}/fastq/merge_split{split}.3prime.fastq.gz"
        )
    elif hasOnlyForwardRead(assignment):
        return "results/assignment/{assignment}/fastq/splits/FW.split{split}.BCattached.fastq.gz"
    else:
        return "results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz"


def getAdapterRemovalReads(assignment: str, five_prime: bool) -> str:
    """
    Return the final reads for mapping after joining, maybe after adapter removal.
    """
    if five_prime and has3PrimeAdapters(assignment):
        return (
            "results/assignment/{assignment}/fastq/merge_split{split}.3prime.fastq.gz"
        )
    elif hasOnlyForwardRead(assignment):
        return "results/assignment/{assignment}/fastq/splits/FW.split{split}.BCattached.fastq.gz"
    else:
        return "results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz"
