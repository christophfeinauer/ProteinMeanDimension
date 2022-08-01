#/bin/bash
set -x
set -e

for model_path in ./models/*pt; do
    out_path=${model_path}_samples.h5
    if [ -s ${out_path} ]; then
        continue
    fi
    python ./sample.py  --model_path $model_path  --latent_dim=5 --num_hidden_units=40 --sample_dist=U --out_file=${out_path} --nsamples=1000 --elbo_samples 5000 --batch_size 5
done
