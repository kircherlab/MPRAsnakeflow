.. _FAQ:

==========================
Frequently Asked Questions
==========================

If you have more question please write us a ticket on `github <https://github.com/kircherlab/MPRAsnakeflow/issues>`_.

MPRAsnakeflow is not able to create a Conda environment
    If you get a message like::

        Caused by: json.decoder.JSONDecodeError: Extra data: line 1 column 2785 (char 2784)#

    Try to do the following steps ::

        rm -r .snakemake/metadata .snakemake/incomplete

    Afterwards try MPRAsnakeflow again. If the above error still occurs, rerun after deleting the entire ``.snakemake`` folder.



Can I use STARR-seq with MPRAsnakeflow?
    No! Not yet ;-)


The pipeline is giving an error **"BUG: Out of jobs ready to be started, but not all files built yet."** and won't run. How can I fix this?
    Please update snakemake, as this error is highly likely to have occured from snakemake internal issues. 
