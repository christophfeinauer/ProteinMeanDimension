import argparse
import h5py
import numpy as np
from tqdm import tqdm


def estimate_mean_dimension(logp, nsamples, q, N):

    contribs = []

    for i in range(N):

        logp_i = logp[:, i, :]

        logp_i_mean = np.mean(logp_i, axis=1)
        logp_i_rand = logp_i[range(nsamples_per_N), np.random.randint(low=0, high=q, size=nsamples_per_N)]

        contribs.append(np.var(logp_i_rand - logp_i_mean))

    var_total = np.var(logp.reshape(nsamples, q)[range(nsamples), np.random.randint(low=0, high=q, size=nsamples)])

    md = sum(contribs) / var_total

    return md, contribs, var_total


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("--sample_file", type=str, required=True)
    parser.add_argument("--out_file", type=str, required=True)
    parser.add_argument("--bootstrap", type=int, default=1000)

    args = parser.parse_args()

    with h5py.File(args.sample_file, 'r') as fid:
        logp = fid["logp"][()]
        samples = fid["samples"][()]

    nsamples, q, N = samples.shape

    assert nsamples * q == logp.shape[0]
    assert nsamples % N == 0

    nsamples_per_N = nsamples // N

    logp = logp.reshape(nsamples_per_N, N, q)
    samples = samples.reshape(nsamples_per_N, N, q, N)

    md, contribs, var_total = estimate_mean_dimension(logp, nsamples, q, N)

    if args.bootstrap > 0:

        md_vec = np.zeros(args.bootstrap)
        contribs_mat = np.zeros((N, args.bootstrap))
        var_total_vec = np.zeros(args.bootstrap)

        print("bootstrapping..")
        for k in tqdm(range(args.bootstrap)):
            logp_bootstrap = logp[np.random.randint(low=0, high=nsamples_per_N, size=nsamples_per_N), :, :]
            md_bootstrap, contribs_bootstrap, var_total_bootstrap = estimate_mean_dimension(logp, nsamples, q, N)
            md_vec[k] = md_bootstrap
            contribs_mat[:, k] = contribs_bootstrap[:]
            var_total_vec[k] = var_total_bootstrap

        md_std = np.std(md_vec)
        contribs_std = np.std(contribs_mat, axis=1)
        var_total_std = np.std(var_total_vec)

    with open(args.out_file, "w") as fid:

        if args.bootstrap == 0:
            print(md, file=fid)
            for i in range(N):
                print(contribs[i], file=fid)
            print(var_total, file=fid)
        else:
            print("{} {}".format(md, md_std), file=fid)
            for i in range(N):
                print("{} {}".format(contribs[i], contribs_std[i]), file=fid)
            print("{} {}".format(var_total, var_total_std), file=fid)
