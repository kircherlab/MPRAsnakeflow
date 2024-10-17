.. _FAQ:

==========================
Frequently Asked Questions
==========================

If you have more question please write us a ticket on `github <https://github.com/kircherlab/MPRAsnakeflow/issues>`_.


Is it possible to differntiate beteween sense and antisense?
    No! Or not directly. The reason why we are not able to do this is that reads will map to both sequence strands equally. Then assignment of the barcode becomes ambigous and is discarded. But when dsigning oligos you can add short sequence fragment on the start and on the end of the sequence that ar edifferent sense and antisense. These sequences should not be trimmed away during demultiplexing and have to be in the design file. For the lentiMPRA dsign we have 15bp adpaters on both ends for integration of the sequence. They can be used for that purpose.

The design/reference file check faild, why?
    The design file has to have:
        * Unique headers. Each sequence has to have a unique sequence/id strating from :code:`>` to the first whitespace or newline.
        * No special characters within the headers. This is because mapping tools create a reference dictionary and cannot handle all characters. In addition most databases (like SRA) have their restricted character set for the header.
        * Unique sequences. They have to be different. Otherwise mapper place the read to both IDs and the barcode get ambigous and is discarded. Wenn you allow min/max start/lengths for sequences (e.g. in BWA mapping) be aware that the smalles substring has to be unqiue across all other (sub) sequences.


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
