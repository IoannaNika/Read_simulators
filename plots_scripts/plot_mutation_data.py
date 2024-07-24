import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


mutations = "data/tuples_pacbio_sars_cov_2_rev_compl/dataset/samples_mutations.tsv"

# load the mutations
mutations = pd.read_csv(mutations, sep="\t")

res = {}

grs = mutations["genomic_region"].unique()

for m in mutations["genomic_region"].unique(): 
    # get the mutations in the region
    region = mutations[mutations["genomic_region"] == m]

    pos_mutations_in_region = list(region[region["label"] == "positive"]["mutations"])
    neg_mutations_in_region = list(region[region["label"] != "positive"]["mutations"])
    res[m] = {"positive": pos_mutations_in_region, "negative": neg_mutations_in_region}


# sort genomic regions based on the start position
sorted_regions = sorted(res.keys(), key=lambda x: int(x.split("_")[0]))
sorted_xtick_regions = [x.replace("_", "-") for x in sorted_regions]
# plot the mutations
fig, ax = plt.subplots()
for i, region in enumerate(sorted_regions):
    pos_mutations = res[region]["positive"]
    neg_mutations = res[region]["negative"]
    # plot boxplots next to each other
    # plot positive mutations
    ax.boxplot(pos_mutations, positions=[i*2 -0.5], widths=0.6, patch_artist=True, showfliers=False, boxprops=dict(facecolor="#E1BE6A"))
    # plot negative mutations
    ax.boxplot(neg_mutations, positions=[i*2 +0.5], widths=0.6, patch_artist=True, showfliers=False, boxprops=dict(facecolor="#40B0A6"))

pos_patch = mpatches.Patch(color='#E1BE6A', label='Positive')
neg_patch = mpatches.Patch(color='#40B0A6', label='Negative')
fig.legend(handles=[pos_patch, neg_patch], labels=['Positive', 'Negative'], 
           ncol=2, loc='upper center', bbox_to_anchor=(0.5, 0.95))
        

ax.set_xticks(range(0, len(sorted_regions)*2, 2))
ax.set_xticklabels(sorted_xtick_regions, rotation=90)

plt.savefig("plots_scripts/mutation_distribution_training_set.pdf", format="pdf", bbox_inches="tight", dpi=300)

##################################################################

mutations = "data/lumc_data/lumc_dataset_pairs_ed_n_m.tsv"

# load the mutations
mutations = pd.read_csv(mutations, sep="\t")

res = {}

grs = mutations["start"].unique()

for m in mutations["start"].unique(): 
    # get the mutations in the region
    region = mutations[mutations["start"] == m]

    pos_mutations_in_region = list(region[region["label"] == "positive"]["n_mutations"])
    neg_mutations_in_region = list(region[region["label"] != "positive"]["n_mutations"])
    res[m] = {"positive": pos_mutations_in_region, "negative": neg_mutations_in_region}


# sort genomic regions based on the start position
sorted_regions = sorted(res.keys(), key=lambda x: int(x))

# plot the mutations
fig, ax = plt.subplots()
for i, region in enumerate(sorted_regions):
    pos_mutations = res[region]["positive"]
    neg_mutations = res[region]["negative"]
    # plot boxplots next to each other
    # plot positive mutations
    ax.boxplot(pos_mutations, positions=[i*2-0.5], widths=0.6, patch_artist=True, showfliers=False, boxprops=dict(facecolor="#E1BE6A"))
    # plot negative mutations
    ax.boxplot(neg_mutations, positions=[i*2+0.5], widths=0.6, patch_artist=True, showfliers=False, boxprops=dict(facecolor="#40B0A6"))

pos_patch = mpatches.Patch(color='#E1BE6A', label='Positive')
neg_patch = mpatches.Patch(color='#40B0A6', label='Negative')
fig.legend(handles=[pos_patch, neg_patch], labels=['Positive', 'Negative'], 
           ncol=2, loc='upper center', bbox_to_anchor=(0.5, 0.95))

ax.set_xticks(range(0, len(sorted_regions)*2, 2))
updated = [x for x in sorted_xtick_regions if int(x.split("-")[0]) in sorted_regions]
ax.set_xticklabels(updated, rotation=90)

plt.savefig("plots_scripts/mutation_distribution_lumc_tuples.pdf", format="pdf", bbox_inches="tight", dpi=300)
