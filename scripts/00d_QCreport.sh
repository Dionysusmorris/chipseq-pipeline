#!/bin/bash


# fastqc report
parallel -j 3 --eta --bar \
  '
  fastqc -t 3 "{}" && echo "{/.} complete"
  ' \
::: ./*sample*.fastq.gz

