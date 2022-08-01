import matplotlib.pyplot as plt
import numpy as np
from read_fasta import read_fasta
import seaborn as sns
sns.set_style("white")


def read_contribs(md_file):

    contribs = []

    with open(md_file) as fid:
        for line_number, line in enumerate(fid):

            if line_number < 1:
                continue

            contribs.append(float(line.split()[0]))

    return np.array(contribs[:-1])


if __name__ == '__main__':

    alignment_file = "../data/GAL4_YEAST_1_b0.6.a2m"
    md_file = "./models/GAL4_YEAST_1_b0.6.a2m_checkpoint_10000.pt.md"

    palette = sns.color_palette(palette='Set2', n_colors=2)

    fig, axes = plt.subplots(1,2, figsize=(10,4))

    contribs = read_contribs(md_file)

    seq_msa, _, q = read_fasta(alignment_file)

    M, N = seq_msa.shape

    assert seq_msa.shape[1] == len(contribs)

    f = np.zeros((N, q))

    for m in range(M):
        for n in range(N):
            f[n, seq_msa[m, n]] += 1

    f = f / M

    entropies = []
    for n in range(N):
        f_n = f[n, :]
        f_n_nz = np.array(list(filter(lambda x: x > 0, f_n)))
        entropy = - sum(f_n_nz * np.log(f_n_nz))
        entropies.append(entropy)

    sns.scatterplot(x=entropies, y=contribs, ax=axes[0], marker='^', alpha=0.8, color=palette[0])
    sns.barplot(x=list(range(1,N+1)), y=contribs, color=palette[1], ax=axes[1])

    axes[0].set_xlabel("Site Entropy")
    axes[0].set_ylabel("Influence")

    axes[1].set_xlabel("Position")
    axes[1].set_ylabel("Influence")
    axes[1].set_xticks([0, 9, 19, 29, 39, 49, 59])



    fig.savefig("entropy_plot_vae.pdf")
