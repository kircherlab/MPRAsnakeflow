.. _FAQ:

==========================
Frequently Asked Questions
==========================

If you have more question please write us a ticket on `github <https://github.com/kircherlab/MPRAsnakeflow/issues>`_.


Is it possible to differentiate between sense and antisense?
    Usually not because reads will map to both sequence strands equally. Then assignment of the barcode becomes ambiguous and is discarded. But we have a workaround that will add unique sequence adapters to both ends to the oligos, for the reference fasta and the fastqs. Now all mapping strategies should be able to differentiate between sense and antisense. To enable use the config :code:`strand_sensitive: {enable: true}`.

The design/reference file check failed, why?
    The design file has to have:
    
    * Unique headers. Each sequence has to have a unique sequence/id starting from :code:`>` to the first whitespace or newline.
    * No special characters within the headers. This is because mapping tools create a reference dictionary and cannot handle all characters. In addition, most databases (like SRA) have their restricted character set for the header.
    * Unique sequences. They have to be different in sense and antisense directions. Otherwise, the mapper places the read to both IDs and the barcode gets ambiguous and is discarded. When you allow min/max start/lengths for sequences (e.g. in BWA mapping) be aware that the smallest substring has to be unique across all other (sub) sequences. If you have antisense collisions and want to keep the strand sensitivity you can enable it by using the option :code:`strand_sensitive: {enable: true}` in the config file (see question before).

MPRAsnakeflow is not able to create a Conda environment
    If you get a message like:

        Caused by: json.decoder.JSONDecodeError: Extra data: line 1 column 2785 (char 2784)#

    Try to do the following steps:

        rm -r .snakemake/metadata .snakemake/incomplete

    Afterwards try MPRAsnakeflow again. If the above error still occurs, rerun after deleting the entire :code:`.snakemake` folder.

Can I use STARR-seq with MPRAsnakeflow?
    No! Not yet ;-)

The pipeline is giving an error **"BUG: Out of jobs ready to be started, but not all files built yet."** and won't run. How can I fix this?
    Please update snakemake, as this error is highly likely to have occurred from snakemake internal issues.
