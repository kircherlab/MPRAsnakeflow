.. _Config:

=====================
Config File
=====================

The config file is a yaml file that contains the configuration. Different runs can be configured. We recommend to use one config file per MPRA experiment or MPRA roject. But in theory many different experiments can be configured in only one file. It is divided into :code:`global` (generell settings), :code:`assignments` (assigment workflow), and :code:`experiments` (count workflow including variants). This is a full example file with all possible configurations. 

.. include:: ../config/example_config.yaml
   :code: yaml


Note that teh config file is conrolled by jscon schema. This means that the config file is validated against the schema. If the config file is not valid, the program will exit with an error message. The schema is located in :download:`../workflow/schemas/config.schema.yaml`.

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :linenos:

----------------
General settings
----------------

The general settings are located in the :code:`global` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :start-after: start_global
   :end-before: start_assignments

split_number
  To parallize mapping for assignment the reads are split into :code:`split_number` files. E.g. setting to 300 this means that the reads are split into 300 files and each file is mapped in parallel.

--------------------
Assignment workflow
--------------------

The assignment workflow is configured in the :code:`assignments` section. The following settings are possible:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
   :language: yaml
   :lines: 91-198

Each asignment you want to process you have to giv him a name like :code:`example_assignment`. The name is used to name the output files.

sequence_length
  Defines the :code:`min` and :code:`max` of a :code:`sequence_length` specify . :code:`sequence_length` is basically the length of a sequence alignment to an oligo in the reference file. Because there can be insertion and deletions we recommend to vary it a bit around the exact length (e.g. +-5). In theory this option enables designs with multiple sequence lengths.
alignment_start
    Defines the :code:`min` and :code:`max` of the start of the alignment in an oligo. When using adapters you have to set basically the length of the adapter. Otherwise 1 will be the choice for most cases. We also recommend to vary this value a bit because the start might not be exact after the adapter. E.g. by +-1. 


--------------------------------------
Experiment workflow (including counts)
--------------------------------------

The experiment workflow is configured in the :code:`experiments` section. The following settings are possible:

.. include:: ../workflow/schemas/config.schema.yaml
   :code: yaml
   :start-after: experiments:


Experiment

.. include:: ../workflow/schemas/samples.schema.yaml
   :code: yaml
