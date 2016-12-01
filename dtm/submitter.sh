#!/bin/bash

echo ""

n_splits=108
for split in $(seq 1 $n_splits); do
    echo "starting job "${split}
    Rscript ./antibiotic_dtm.R ${split}
done;

echo "Done submitting jobs"
