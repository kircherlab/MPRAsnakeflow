.. _Adapter trimming:

=================
Adapter trimming
=================

It is important to state that MPRAsnakeflow requires the core information as input. This means that reads should not have unnecessary sequences left and right of barcodes, oligos or UMIs. There is no difference between the assignment and experiment workflows. There are only two exceptions.

1. Within the assignment workflow, it is possible that the barcode and the oligo can be within one read separated by a spacer with a fixed sequence or length. This can be defined within the config file; see :ref:`Config`. Use it without a barcode read and the option :code:`linker_length` or :code:`linker`.
2. Within the experiment workflow, when single-end reads are used (only FWD or FWD + UMI), it retrieves the barcode from the first N positions of the read, no matter how long the read is. It can also be the last N positions by using the option :code:`bc_extraction: end` (default is :code:`start`; see :ref:`Config`).

Otherwise, you have to trim your reads to run the workflow. Since MPRAsnakeflow version 0.6.0, we provide an optional adapter trimming step using `Cutadapt <https://cutadapt.readthedocs.io/>`_ (before it was only supported partially within the assignment step). You can enable it by setting the config option :code:`adapters: {<READ_TYPE>: ...}` in both the experiment and assignment workflows. By default, no adapters are trimmed.

.. note::

    Adapter trimming can only be used when all reads of a certain type (FWD, REV, UMI, BC) have the same adapters. For example, when you sequenced FWD and REV twice, once each with 150bp and once with 151bp, you might need to remove adapters of different lengths between the two runs. This is currently not supported but can be done manually before running MPRAsnakeflow. See :ref:`Complex example` for an example where this is exactly the case within the assignment workflow.

    Also it is not possible to trim adapter different DNA and RNA counts. for FWD (or REV) count reads they have to be performed in the same way.

We use the adapter trimming option in two examples: :ref:`Complex example` and :ref:`Plasmid example` if you want to see real examples.

Configure adapter trimming
--------------------------

This is done via the config file. You can specify adapters for each read type (FWD, REV, UMI, BC) separately. You can use one of two options:

1. Specify a sequence length for trimming (5' and/or 3' end)
2. Specify the actual adapter sequences to be trimmed (5' and/or 3' end; multiple adapter sequences allowed)

Both options cannot be used together for one read type.

As an example, we want to trim the first 5 and the last 10 bases from a FWD read:

.. code-block:: yaml

    adapters:
      FWD:
        - 5
        - -10

You see we just name the read type (FWD) and provide a list with two entries. Positive entries indicate trimming from the start (5' end) and negative numbers indicate trimming from the end (3' end). This is equal to the :code:`-u` option in `Cutadapt <https://cutadapt.readthedocs.io/>`_.

When we want to use the exact sequence of an adapter, we can use a config like this. Now we trim the start of the FWD and REV reads as well as the end of the BC read with specific sequences:

.. code-block:: yaml

    adapters:
      FWD:
        five_prime:
            - AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
            - AGATCGGTCGAAADT
      REV:
        five_prime:
            - AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
      BC:
        three_prime:
            - ACACTCTTTCCCTACACGACGCTCTTCCGATCT
