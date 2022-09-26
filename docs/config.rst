.. _Config:

=====================
Config File
=====================

The config file is a yaml file that contains the configuration. Different runs can be configured. We recommend to use one config file per MPRA experiment or MPRA roject. But in theory many different experiments can be configured in only one file. It is divided into :code:`global` (generell settings), :code:`assignments` (assigment workflow), and :code:`experiments` (count workflow including variants). This is a full example file with all possible configurations. 

.. include:: ../config/example_config.yaml
   :code: yaml
   :linenos:


Note that teh config file is conrolled by jscon schema. This means that the config file is validated against the schema. If the config file is not valid, the program will exit with an error message. The schema is located in :download:`../workflow/schemas/config.schema.yaml`.

----------------
General settings
----------------

The general settings are located in the :code:`global` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_global
   :end-before: start_assignments

split_number
  To parallize mapping for assignment the reads are split into :code:`split_number` files. E.g. setting to 300 this means that the reads are split into 300 files and each file is mapped in parallel. This is only usefull when using on a cluster. Running the workflow only on one machine the default value shopuld be used. Default is set to 1. 

--------------------
Assignment workflow
--------------------

The assignment workflow is configured in the :code:`assignments` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_assignments
   :end-before: start_experiments

Each asignment you want to process you have to giv him a name like :code:`example_assignment`. The name is used to name the output files.

sequence_length
    Defines the :code:`min` and :code:`max` of a :code:`sequence_length` specify . :code:`sequence_length` is basically the length of a sequence alignment to an oligo in the reference file. Because there can be insertion and deletions we recommend to vary it a bit around the exact length (e.g. +-5). In theory this option enables designs with multiple sequence lengths.
alignment_start
    Defines the :code:`min` and :code:`max` of the start of the alignment in an oligo. When using adapters you have to set basically the length of the adapter. Otherwise 1 will be the choice for most cases. We also recommend to vary this value a bit because the start might not be exact after the adapter. E.g. by +-1.
R1
    List of forward read files in gzipped fastq format. The full or relative path to the files should be used. Same order in R1, R2, and R3 is important.
R2
    List of index read files in gzipped fastq format. The full or relative path to the files should be used. Same order in R1, R2, and R3 is important.
R3
    list of reverse read files in gzipped fastq format. The full or relative path to the files should be used. Same order in R1, R2, and R3 is important.
reference
    Design file (full or relative path) in fasta format. The design file should contain the oligos in fasta format. The header should contain the oligo name and should be unique. The sequence should be the sequence of the oligo and must also be uniq. When having multiple oligo name swith the same sequence please merge them into one fatsa entry. The oligo name later used to link barcode to oligo. The sequence is used to map the reads to the oligos. Adapters can be in the seuqence and therefore :code:`alignment_start` has to be adjusted.
configs
    After mapping the reads to the design file and extracting the barcodes per oligo the configuration (using different names) can be used to generate multiple filtering and configuration settings of the final maq oligo to barcode. Each configuration is a dictionary with the following keys:
        min_support
            Minimum number of same BC that map to teh same oligo. Larger value gives more evidence to be correct. But can remove lot's of BCs (depedning on the complexity, sequencing depth and quality of sequencing). Recommended option is :code:`3`.
        fraction
            Minumum fraction of same BC that map to teh same oligo. E.g. :code:`0.7` means that at least 70% of the BC map to the same oligo. Larger value gives more evidence to be correct. But can remove lot's of BCs (depedning on the complexity, sequencing depth and quality of sequencing). Recommended option is :code:`0.7`.
        unknown_other
            (Optional) Shows not mapped BCs in the final output map. Not recommended to use as mapping file fore the experiment workflow. But can be usefull for debugging. Default is :code:`false`. 
        ambigous
            (Optional) Shows ambigous BCs in the final output map. Not recommended to use as mapping file fore the experiment workflow. But can be usefull for debugging. Default is :code:`false`.



--------------------------------------
Experiment workflow (including counts)
--------------------------------------

The experiment workflow is configured in the :code:`experiments` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_experiments
   :end-before: end_experiments


Experiment

.. literalinclude:: ../workflow/schemas/samples.schema.yaml
   :language: yaml
