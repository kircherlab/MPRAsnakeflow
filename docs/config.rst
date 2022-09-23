.. _Config:

=====================
Config File
=====================

The config file is a yaml file that contains the configuration. Different runs can be configured. We recommend to use one config file per MPRA experiment or MPRA roject. But in theory many different experiments can be configured in only one file. It is divided into :code:`global` (generell settings), :code:`assignments` (assigment workflow), and :code:`experiments` (count workflow including variants). This is a full example file with all possible configurations. 

.. include:: ../config/example_config.yaml
   :code: yaml


Note that teh config file is conrolled by jscon schema. This means that the config file is validated against the schema. If the config file is not valid, the program will exit with an error message. The schema is located in :code:`../workflow/schemas/config.schema.yaml`.

.. include:: ../workflow/schemas/config.schema.yaml
   :code: yaml

----------------
General settings
----------------

The general settings are located in the :code:`global` section. The following settings are possible:

.. include:: ../workflow/schemas/config.schema.yaml
   :code: yaml
   :start-after: start_global
   :end-before: start_assignments

:split_number
  To parallize mapping for assignment the reads are split into :code:`split_number` files. E.g. setting to 300 this means that the reads are split into 300 files and each file is mapped in parallel.

--------------------
Assignment workflow
--------------------

The assignment workflow is configured in the :code:`assignments` section. The following settings are possible:

.. include:: ../workflow/schemas/config.schema.yaml
   :code: yaml
   :start-after: split_number
   :end-before: experiments

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
