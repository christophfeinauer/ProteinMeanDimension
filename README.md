# Mean Dimension of Generative Models for Protein Sequences


This repository contains the code for reproducing all plots in the paper.


## Setup

### Clone repository and enter the directory

Make sure `git-lfs` is installed by running

```bash
$ git lfs install --skip-repo
```

on the shell.

Then, execute


```bash
$ git clone git@github.com:christophfeinauer/ProteinMeanDimension.git && cd ProteinMeanDimension
```

### Create Julia Environment


```
$ julia --project=. -e 'using Pkg; Pkg.instantiate(;verbose=true)'
```

## Download Data

Run

```
bash get_data.sh
```

to get the deepsequence alignments. These are used for the

## Calculating the mean dimension on a new model

For calculating the mean dimension for a new model you can use the code in `mean_dimension_from_samples.py`. It expects a [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) file with two datasets in it, one called `samples`, which should contain the protein sequences with amino acids mapped to integer indices, and one called `logp`, which should contain the corresponding log probabilities for the samples. For the models used in the paper (ArDCA and VAE), the scripts in `ardca/` and `vae/` will produce these files.

The layout of these datasets is a bit intricate since the calculation of the log probability and the evaluation of the mean dimension is decoupled in this code, to make it more efficient.

The mean dimension is based on estimating the contribution of single positions to the variance of the log probability under the uniform distribution. The input to the code is therefore the log probability of sequences where single amino acids have been exchanged.

For every of the `N` positions the calculation of the mean dimension is based on `nsamples_per_position` samples. The dataset `samples` should therefore be of size `(nsamples_per_position, N, q, N)`, where `q` is the number of possible amino acids (typically `q=21`) and  `(:, i, :, :)` are the sequences used for estimating the contribution of position `i`, where `(m, i, a, :)` is a single sequence of length `N`. For a given index `m`, the sequences `(m, i, :, :)` should only differ in the position `i` and `(m, i, a, :)` should contain an `a` in position `i`. The code then uses comparisons between `(m, i, a, :)` and `(m, i, b, :)` for calculating the mean dimension.

The dataset `logp` should be of size `(nsamples_per_position, N, q)` and contain at index `(m, i, a)` the log probability in the model of the sequence in `(m, i, a, :)`. 

To make this more clear and assuming that the function `get_logp(seq)` calculates the log probability of a sequence in the model, then the following pseudocode illustrates how to create the datasets:

```
logp = zeros(nsamples_per_position, N ,q)
samples = zeros(nsamples_per_position, N, q ,N)
for i in 1:N
    for m in 1:nsamples_per_position
        seq = rand(1:q, N)
        for a in 1:q
            seq[i] = a
            logp[m, i, a] = get_logp(seq)
            samples[m, i, a, :] = seq[:]
        end
    end
end
```

The code in `mean_dimension_from_samples.py` uses `samples` only for consistency checks. If you are sure that you got everything right you can comment out these and pass only `logp`.


