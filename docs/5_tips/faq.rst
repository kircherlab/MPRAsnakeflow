.. _FAQ:

==========================
Frequently Asked Questions
==========================

If you have more questions, please write us a ticket on `GitHub <https://github.com/kircherlab/MPRAsnakeflow/issues>`_.

Is it possible to differentiate between sense and antisense?
------------------------------------------------------------
Usually not, because reads will map to both sequence strands equally. Then, the assignment of the barcode becomes ambiguous and is discarded. However, we have a workaround that adds unique sequence adapters to both ends of the oligos for the reference FASTA and the FASTQs. Now, all mapping strategies should be able to differentiate between sense and antisense. To enable this, use the config option: :code:`strand_sensitive: {enable: true}`.

The design/reference file check failed. Why?
--------------------------------------------
The design file must meet the following requirements:

* **Unique headers**: Each sequence must have a unique sequence ID starting from :code:`>` to the first whitespace or newline.
* **No special characters in headers**: Mapping tools create a reference dictionary and cannot handle all characters. Additionally, most databases (like SRA) have a restricted character set for headers.
* **Unique sequences**: Sequences must be different in both sense and antisense directions. Otherwise, the mapper places the read to both IDs, and the barcode becomes ambiguous and is discarded. 

When you allow min/max start/lengths for sequences (e.g., in BWA mapping), ensure that the smallest substring is unique across all other (sub)sequences. If you have antisense collisions and want to keep strand sensitivity, enable it using the option :code:`strand_sensitive: {enable: true}` in the config file (see the previous question).

MPRAsnakeflow is not able to create a Conda environment
--------------------------------------------------------
If you encounter an error like:

    Caused by: json.decoder.JSONDecodeError: Extra data: line 1 column 2785 (char 2784)#

Try the following steps:

1. Remove the incomplete metadata:

    .. code-block:: bash

        rm -r .snakemake/metadata .snakemake/incomplete

2. Retry running MPRAsnakeflow. If the error persists, delete the entire :code:`.snakemake` folder and rerun:

    .. code-block:: bash

        rm -r .snakemake

Can I use STARR-seq with MPRAsnakeflow?
---------------------------------------
No, not yet ;-)

The pipeline is giving an error **"BUG: Out of jobs ready to be started, but not all files built yet."** How can I fix this?
----------------------------------------------------------------------------------------------------------------------------
This error is likely caused by internal issues in Snakemake. Please update Snakemake to the latest version to resolve this issue.
