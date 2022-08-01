import matplotlib.pyplot as plt
import os
import seaborn as sns
from collections import defaultdict
import numpy as np
import matplotlib.lines as mlines
import argparse
sns.set_style('white')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("--dir", default="./models")

    args = parser.parse_args()

    markers = ['^']

    proteins = ['BRCA1', 'GAL4', 'UBC9', 'SUMO1']

    palette = sns.color_palette(palette='Set2', n_colors=len(proteins))

    md_dict = defaultdict(list)
    sr_dict = defaultdict(list)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    #axins = [None]*2
    #axins[0] = axes[0].inset_axes([0.3, 0.2, 0.47, 0.47])
    #axins[1] = axes[1].inset_axes([0.3, 0.35, 0.47, 0.47])

    for file in os.listdir(args.dir):

        file_path = os.path.join(args.dir, file)

        if not file.endswith("pt"):
            continue

        if "checkpoint" not in file:
            continue

        epoch = int(file.split("_")[-1].split(".")[0])

        plist = [file.startswith(protein) for protein in proteins]
        if not any(plist):
            continue
        else:
            protein_ind = plist.index(True)
            protein = proteins[protein_ind]

        md_file_path = file_path + ".md"
        sr_file_path = file_path + ".sr"

        if os.path.isfile(md_file_path):
            with open(md_file_path) as fid:
                data = np.loadtxt(md_file_path)
                md = data[0][0]
                md_std = data[0][1]
                md_dict[protein].append((epoch, md, md_std))

        if os.path.isfile(sr_file_path):
            with open(sr_file_path) as fid:
                data = np.loadtxt(sr_file_path)
                sr = data.item()
                sr_dict[protein].append((epoch, sr))

    for protein_ind, protein in enumerate(proteins):

        print(protein)

        epochs_mds = np.array([k[0] for k in md_dict[protein]])
        mds = np.array([k[1] for k in md_dict[protein]])
        std_mds = np.array([k[2] for k in md_dict[protein]])

        print(epochs_mds)

        epochs_srs = np.array([k[0] for k in sr_dict[protein]])
        srs = np.array([k[1] for k in sr_dict[protein]])

        sortperm_mds = epochs_mds.argsort()
        mds = mds[sortperm_mds]
        std_mds = std_mds[sortperm_mds]
        epochs_mds = epochs_mds[sortperm_mds]

        sortperm_srs = epochs_srs.argsort()
        srs = srs[sortperm_srs]
        epochs_srs = epochs_srs[sortperm_srs]


        sns.lineplot(x=epochs_srs, y=srs, ax=axes[0], color=palette[protein_ind], marker=markers[0], alpha=0.8, ls='--')
        sns.lineplot(x=epochs_mds, y=mds, ax=axes[1], color=palette[protein_ind], marker=markers[0], alpha=0.8, ls='--')

        #sns.lineplot(x=epochs_srs, y=srs, ax=axins[0], color=palette[protein_ind], marker=markers[0], alpha=0.8, ls='--')
        #sns.lineplot(x=epochs_mds, y=mds, ax=axins[1], color=palette[protein_ind], marker=markers[0], alpha=0.8, ls='--')

        axes[0].set_xlabel("Epoch")
        axes[0].set_ylabel("Spearman Correlation")
        axes[1].set_xlabel("Epoch")
        axes[1].set_ylabel("Mean Dimension")

        axes[0].set_xscale("symlog")
        axes[1].set_xscale("symlog")

        #axins[0].set_xlim(0,250)
        #axins[1].set_xlim(0,250)

        axes[0].set_xlim(0, 10000)
        axes[1].set_xlim(0, 10000)

        axes[0].yaxis.set_ticks_position("both")
        axes[1].yaxis.set_ticks_position("both")

        axes[0].minorticks_on()
        axes[1].minorticks_on()

    # create color legend
    handles = []
    for protein_ind, protein in enumerate(proteins):
        handle = mlines.Line2D([], [], color=palette[protein_ind], marker='^', linestyle='None', markersize=10, label=protein)
        handles.append(handle)

    lgd = fig.legend(handles=handles, ncol=len(proteins), bbox_to_anchor=(0.5, .04), loc='center', frameon=False)

    fig.tight_layout()

    fig.subplots_adjust(bottom=0.2)

    plt.savefig("md_sr_overlay_vae.pdf")
