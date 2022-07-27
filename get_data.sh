#!/bin/bash

## Get DeepSequence data
## This is a copy of the script at https://github.com/debbiemarkslab/DeepSequence/blob/master/examples/download_alignments.sh

curl -o alignments.tar.gz https://marks.hms.harvard.edu/deepsequence/alignments.tar.gz
tar -xvf alignments.tar.gz
rm alignments.tar.gz
