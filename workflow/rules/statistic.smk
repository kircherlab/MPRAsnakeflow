####################################
# Statistic of barcodes and oligos #
####################################


##################
## subworkflows ##
##################


# statistic on BC counts (not assigned)
include: "statistic/counts.smk"
# statistic on assigned counts (BCs and oligos)
include: "statistic/assignment.smk"
# statistic on correlation of BCs/coligos
include: "statistic/correlation.smk"
# BC overlap between replicates
include: "statistic/BCoverlapBetweenReplicates.smk"
