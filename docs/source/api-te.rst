.. _api_te:

run-tea
=======

Run translation efficiency (TE)

.. code-block:: bash

   usage: run-tea [--help] [--batch] [--filter]
                  [--method METHOD] [--lfcThreshold LFCTHRESHOLD] [--alpha ALPHA]
                  [--symbolCol SYMBOLCOL] [--orfCol ORFCOL] [--delim DELIM]
                  [--lfcShrinkMethod LFCSHRINKMETHOD] config


Positional Arguments
--------------------
``config``
    Yaml config file (same as used for ``run-htseq-workflow``)

Flags
-----
-h, --help          show this help message and exit
-b, --batch         Use 'batch' from sample table
-f, --filter        Filter features with 0 counts in each assay separately

Optional Arguments
------------------
-m, --method        Method to use [default: ``LRT``]
-l, --lfcThreshold  LFC threshold [default: ``log2(1.2)``]
-a, --alpha         FDR threshold [default: ``0.05``]
-s, --symbolCol     Column for features (symbol/names) - HTSeq workflow only [default: ``2``]
-o, --orfCol        Column for ORF types - HTSeq workflow only
-d, --delim         Field separator for the count tables {TAB,CSV} [default: ``TAB``]
--lfcShrinkMethod   LFC shrinkage method {apeglm,ashr} [default: ``apeglm``]
