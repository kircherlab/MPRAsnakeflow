.. _Config:

=====================
Config File
=====================

The config file is a yaml file that contains the configuration. Different runs can be configured. We recommend to use one config file per MPRA experiment or MPRA project. But in theory many different experiments can be configured in only one file. It is divided into :code:`global` (generell settings), :code:`assignments` (assigment workflow), and :code:`experiments` (count workflow including variants). This is a full example file with default configurations. :download:`config/example_config.yaml <../config/example_config.yaml>`.

.. literalinclude:: ../config/example_config.yaml
   :language: yaml
   :linenos:


Note that the config file is conrolled by json schema. This means that the config file is validated against the schema. If the config file is not valid, the program will exit with an error message. The schema is located in :download:`workflow/schemas/config.schema.yaml <../workflow/schemas/config.schema.yaml>`.

----------------
General settings
----------------

The general settings are located in the :code:`global` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_global
   :end-before: start_assignments

:assignments:
    Global parameters that hold for the assignment workflow.

    :split_number:
        To parallize mapping for assignment the reads are split into :code:`split_number` files. E.g. setting to 300 this means that the reads are split into 300 files and each file is mapped in parallel. This is only usefull when using on a cluster. Running the workflow only on one machine the default value shopuld be used. Default is set to 1. 

--------------------
Assignment workflow
--------------------

The assignment workflow is configured in the :code:`assignments` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_assignments
   :end-before: start_experiments

Each assignment you want to process you have to giv him a name like :code:`example_assignment`. The name is used to name the output files.

:alignment_tool:
    Alignment tool configuration that is used to map the reads to the oligos.
    
    :tool:
        Alignment tool that is used. Currently :code:`bwa` and :code:`exact` is supported.
    :configs:
        Configurations of the alignment tool selected.

        :sequence_length (bwa):
            Defines the :code:`min` and :code:`max` of a :code:`sequence_length` specify . :code:`sequence_length` is basically the length of a sequence alignment to an oligo in the design file. Because there can be insertion and deletions we recommend to vary it a bit around the exact length (e.g. +-5). In theory this option enables designs with multiple sequence lengths.
        :alignment_start (bwa):
            Defines the :code:`min` and :code:`max` of the start of the alignment in an oligo. When using adapters you have to set basically the length of the adapter. Otherwise 1 will be the choice for most cases. We also recommend to vary this value a bit because the start might not be exact after the adapter. E.g. by +-1.
        :min_mapping_quality (bwa):
            (Optinal) Defines the minimum mapping quality (MAPQ) of the alinment to an oligo. When using oligos with only 1bp difference it is recommended to set it to 1. For regions only with larger edit distances 30 or 40 might be a good choice. Default :code:`1`. 
        :sequence_length (exact):
            Defines the :code:`sequence_length` which is the length of a sequence alignment to an oligo in the design file. Only one length design is supported.
        :alignment_start (exact):
            Defines the start of the alignment in an oligo. When using adapters you have to set basically the length of the adapter. Otherwise 1 will be the choice for most cases.

:bc_length:
    Length of the barcode. Must match with the length of :code:`BC`.
:BC_rev_comp:
    (Optional) If set to :code:`true` the barcode of is reverse complemented. Default is :code:`false`.
:linker_length:
    (Optional) Length of the linker. Only needed if you don't have a barcode read and the barcode is in the FW read with the structure: BC+Linker+Insert. The fixed length is used for the linker after a fixed length of BC. The recommended option is :code:`linker` by defining the exact linker sequence and using cutadapt for trimming.
:linker:
    (Optional) Length of the linker. Only needed if you don't have a barcode read and the barcode is in the FW read with the structure: BC+Linker+Insert. Uses cutadapt to trim the linker to get the barcode as well as the starting of the insert.
:FW:
    List of forward read files in gzipped fastq format. The full or relative path to the files should be used. Same order in FW, BC, and REV is important.
:REV:
    list of reverse read files in gzipped fastq format. The full or relative path to the files should be used. Same order in FW, BC, and REV is important.
:BC:
    List of index read files in gzipped fastq format. The full or relative path to the files should be used. Same order in FW, BC, and REV is important.
:NGmerge:
    (Optional) Options for NGmerge. NGmerge is used merge FW and REV reads. The following options are possible (we recommend to use the default values):

    :min_overlap:
        (Optional) Minimum overlap of the reads. Default :code:`20`.
    :frac_mismatches_allowed:
        (Optional) Fraction of mismatches allowed in the overlap. Default :code:`0.1`.
    :min_dovetailed_overlap:
        (Optional) Minimum dovetailed overlap. Default :code:`10`.

:design_file:
    Design file (full or relative path) in fasta format. The design file should contain the oligos in fasta format. The header should contain the oligo name and should be unique. The sequence should be the sequence of the oligo and must also be unique. When having multiple oligo names with the same sequence please merge them into one fasta entry. The oligo name later used to link barcode to oligo. The sequence is used to map the reads to the oligos. Adapters can be in the seuqence and therefore :code:`alignment_start` has to be adjusted.
:design_check:
    (Optional) Options for checking your design fasta file. Design  file cannot have :code:`[` or :code:`]`, duplicated headers and for best performance sequences should not be identical.

    :fast:
        (Optional) Using a simple dictionary to find identical sequences. This is faster but uses only the whole (or center part depending on start/length) of the design file. Cannot find substrings as part of any sequence. Set to false for more correct, but slower, search. Default :code:`true`.
    :sequence_collitions:
        (Optional) Check if there are identical sequences in the design file. Default :code:`true`.

:configs:
    After mapping the reads to the design file and extracting the barcodes per oligo the configuration (using different names) can be used to generate multiple filtering and configuration settings of the final maq oligo to barcode. Use `<your_config_name>: {}` to use the default values for the keys. Each configuration is a dictionary with the following keys:
    
    :min_support:
        Minimum number of same BC that map to the same oligo. Larger value gives more evidence to be correct. But can remove lot's of BCs (depedning on the complexity, sequencing depth and quality of sequencing). Recommended option is :code:`3`.
    :fraction:
        Minumum fraction of same BC that map to the same oligo. E.g. :code:`0.7` means that at least 70% of the BC map to the same oligo. Larger value gives more evidence to be correct. But can remove lot's of BCs (depedning on the complexity, sequencing depth and quality of sequencing). Recommended option is :code:`0.7`.
    :unknown_other:
            (Optional) Shows not mapped BCs in the final output map. Not recommended to use as mapping file fore the experiment workflow. But can be usefull for debugging. Default is :code:`false`. 
    :ambigous:
            (Optional) Shows ambigous BCs in the final output map. Not recommended to use as mapping file fore the experiment workflow. But can be usefull for debugging. Default is :code:`false`.

--------------------------------------
Experiment workflow (including counts)
--------------------------------------

The experiment workflow is configured in the :code:`experiments` section. Each experiment run (contains one experiment file with all replicates of an experiment). The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_experiments
   :end-before: end_experiments

:bc_length:
    Length of the barcode. This is used to extract the barcode from the index read. The barcode is extracted from the first :code:`bc_length` bases of the index read. When no reverse read is given and :code:`adapter` is not set teh exact length is used to extract the DNA BC from the FW read.
:umi_length:
    (Optional) Length of the UMI. This is used to extract the UMI from the index read. The UMI is extracted from the last :code:`umi_length` bases of the index read. Please provide if you use UMIs.
:adapter:
    (Optional) Adapter sequence in the FW read when no reverse read is given. This is used to trim the sequence and retrieve the BC using cutadapt.
:data_folder:
    Folder where the fastq files are located. Files are defined in the :code:`experiment_file`. The full or relative path to the folder should be used.
:experiment_file:
    Path to the experiment file. The full or relative path to the file should be used. The experiment file is a comma separated file and is decribed in the `Experiment file`_ section.
:demultiplex:
    (Optional) If set to :code:`true` the reads are demultiplexed. This means that the reads are split into different files for each barcode. This is usefull for further analysis. Default is :code:`false`.
:label_file:
    (Optional) Path to the label file. The full or relative path to the file should be used. The label file is a tab separated file and contais the oligo name and the label of it. The oligo name should be the same as in the design file. The label is used to group the oligos in the final output, e.g. for plotting. 
    
    .. code-block:: text

        insert1_name label1
        insert2_name label1
        insert3_name label2

:assignments:
    Per experiments multiple assignments can be defined (naming them differently). Everey assignment name contains the following configurations:

    :type:
        Can be :code:`file` or :code:`config`. :code:`file` means that you use a mapping file which is tab separated and gzipped. It contains in the first column the barcode and in the second column the oligo name. This file can be generated by the `Assignment workflow`_. When using :code:`config`this means that you are referring to a assignment that is specified in this config file.
    :assignment_file:
        When using :code:`file` please insert the path to the assignment file (tsv.gz). When using :code:`config` please set the name of the config previously described the assignment that should be used.
    :assignment_name:
        When using :code:`config` please insert the name of the assignment specified in the config file.
    :assignment_config:
        When using :code:`config` please insert the name config of the :code:`assignment_name` you want to use.
    :sampling:
        (Optional) Options Randomly removing barcodes in the assignment. Just for debug reasons.

        :prop:
            Sample down the BCs in the assignment file to this proporion.
        :total:
            Sample down the BCs in the assignment file to this number.


:configs:
    Each experiment run can have multiple configurations including filter and sampling options.

    :filter:
        (Optional) Filter options. These options are available

        :bc_threshold:
            Minimum number of different BCs required per oligo. A higher value normally increases the correlation betwene replicates but also reduces the number of final oligos. Default option is :code:`10`.
        :DNA:
            Settings for DNA

            :min_counts:
                Mimimum number of DNA counts per barcode. When set to :code:`0` a pseudo count is added. Default option is :code:`1`.
        :RNA:
            Settings for DNA

            :min_counts:
                Mimimum number of RNA counts per barcode. When set to :code:`0` a pseudo count is added. Default option is :code:`1`.
    :sampling:
        (Optional) Options for sampling counts and barcodes. Just for debug reasons.

        :DNA:
            Settings for sampling DNA counts.

            :threshold:
                Maximum threshold for DNA counts assigned to a BC.
            :prop:
                Sample down the DNA counts to this proporion.
            :total:
                Sample down the DNA counts to this number.
            :seed:
                Seed for the random DNA sampling.
        
        :RNA:
            Settings for sampling RNA counts.

            :threshold:
                Maximum threshold for RNA counts assigned to a BC.
            :prop:
                Sample down the RNA counts to this proporion.
            :total:
                Sample down the RNA counts to this number.
            :seed:
                Seed for the random RNA sampling.

=====================
Experiment file
=====================

Here we have 4 different options:

------------------------------
Forward, reverse, and UMI read
------------------------------

Experiment file has a header with :code:`Condition`, :code:`Replicate`, :code:`DNA_BC_F`, :code:`DNA_UMI`, :code:`DNA_BC_R`, :code:`RNA_BC_F`, :code:`RNA_UMI`, and :code:`RNA_BC_R`. Condition together with replicate have to be a uniqe name. Both field entries are not allowed to have :code:`_` and :code:`.`. Multiple file names are allowd seperating them via :code:`;`. An example experiment file can be found here: :download:`resources/example_experiment.csv <../resources/example_experiment.csv>`.

.. literalinclude:: ../resources/example_experiment.csv
   :language: text

------------------------
Forward and reverse read
------------------------

Experiment file has a header with :code:`Condition`, :code:`Replicate`, :code:`DNA_BC_F`, :code:`DNA_BC_R`, :code:`RNA_BC_F`, and :code:`RNA_BC_R`. Condition together with replicate have to be a uniqe name. Both field entries are not allowed to have :code:`_` and :code:`.`. Multiple file names are allowd seperating them via :code:`;`.

------------------
Only forward read
------------------

Experiment file has a header with :code:`Condition`, :code:`Replicate`, :code:`DNA_BC_F`, and :code:`RNA_BC_F`. Condition together with replicate have to be a uniqe name. Both field entries are not allowed to have :code:`_` and :code:`.`. Multiple file names are allowd seperating them via :code:`;`.


-------------------------------------------------------
Forward, reverse, and UMI read using demultiplex option
-------------------------------------------------------

Experiment file has a header with :code:`Condition`, :code:`Replicate`, :code:`BC_DNA`, :code:`BC_RNA`, :code:`BC_F`, :code:`BC_R`, :code:`UMI`, and :code:`INDEX`. Condition together with replicate have to be a uniqe name. Both field entries are not allowed to have :code:`_` and :code:`.`. Multiple file names are allowd seperating them via :code:`;`.
