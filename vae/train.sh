#/bin/bash

for train_fasta in $(ls ../data/*a2m); do

   train_fasta_basename=$(basename $train_fasta)
   out_file_prefix=./models/${train_fasta_basename}

   chk0="${out_file_prefix}"_checkpoint_0.pt
   if [ -s "${chk0}" ]; then
        echo "found $chk0, skipping"
        continue
    fi

    python ./train_on_fasta.py --fasta_train_path $train_fasta --checkpoint --latent_dim 5 --weight_decay 0.01 --num_hidden_units 40 --checkpoint --out_file_prefix ${out_file_prefix}

done
