.. _api_de:

run-dea
=======

Run differential expression (DE)

.. code-block:: bash

   usage: run-dea [--help] [--batch]
                  [--lfcThreshold LFCTHRESHOLD] [--alpha ALPHA]
                  [--symbolCol SYMBOLCOL] [--orfCol ORFCOL] [--delim DELIM]
                  [--lfcShrinkMethod LFCSHRINKMETHOD] config

Positional arguments
--------------------
``config``
    Yaml config file (same as used for ``run-htseq-workflow``)

Flags
-----
-h, --help          show this help message and exit
-b, --batch         Use 'batch' from sample table

Optional Arguments
------------------
-l, --lfcThreshold  LFC threshold [default: ``log2(1.2)``]
-a, --alpha         FDR threshold [default: ``0.05``]
-s, --symbolCol     Column for features (symbol/names) - HTSeq workflow only [default: ``2``]
-o, --orfCol        Column for ORF types - HTSeq workflow only
-d, --delim         Field separator for the count tables {TAB,CSV} [default: ``TAB``]
--lfcShrinkMethod   LFC shrinkage method {apeglm,ashr} [default: ``apeglm``]
