#!/bin/bash

for model_path in ./vae/models/*pt ./ardca/models/*jld; do
    samples_path=${model_path}_samples.h5
    if [ ! -s "$samples_path" ]; then
        continue
    fi
    md_out_file=${model_path}.md
    if [ -s "$md_out_file" ]; then
        :
    else
        python ./mean_dimension_from_samples.py --sample_file "$samples_path" --out_file "$md_out_file" --bootstrap 1000
    fi
done
