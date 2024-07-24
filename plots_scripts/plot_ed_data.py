import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


eds = "data/tuples_pacbio_sars_cov_2_rev_compl/dataset/samples.tsv"

# load the edit distances
eds = pd.read_csv(eds, sep="\t")

res = {}

grs = eds["genomic_region"].unique()

for m in eds["genomic_region"].unique(): 
    region = eds[eds["genomic_region"] == m]
    pos_eds_in_region = list(region[region["label"] == "positive"]["edit_distance"])
    neg_eds_in_region = list(region[region["label"] != "positive"]["edit_distance"])
    res[m] = {"positive": pos_eds_in_region, "negative": neg_eds_in_region}


# sort genomic regions based on the start position
sorted_regions = sorted(res.keys(), key=lambda x: int(x.split("_")[0]))
sorted_xtick_regions = [x.replace("_", "-") for x in sorted_regions]
# plot the mutations
fig, ax = plt.subplots()
for i, region in enumerate(sorted_regions):
    pos_eds = res[region]["positive"]
    neg_eds = res[region]["negative"]
    # plot boxplots next to each other
    # plot positive mutations
    ax.boxplot(pos_eds, positions=[i*2 -0.5], widths=0.6, patch_artist=True, showfliers=False, boxprops=dict(facecolor="#E1BE6A"))
    # plot negative mutations
    ax.boxplot(neg_eds, positions=[i*2 +0.5], widths=0.6, patch_artist=True, showfliers=False, boxprops=dict(facecolor="#40B0A6"))

pos_patch = mpatches.Patch(color='#E1BE6A', label='Positive')
neg_patch = mpatches.Patch(color='#40B0A6', label='Negative')
fig.legend(handles=[pos_patch, neg_patch], labels=['Positive', 'Negative'], 
           ncol=2, loc='upper center', bbox_to_anchor=(0.5, 0.95))
        

ax.set_xticks(range(0, len(sorted_regions)*2, 2))
ax.set_xticklabels(sorted_xtick_regions, rotation=90)

plt.savefig("plots_scripts/eds_distribution_training_set.pdf", format="pdf", bbox_inches="tight", dpi=300)

##################################################################

eds = "data/lumc_data/lumc_dataset_pairs_ed_n_m.tsv"

# load the mutations
eds = pd.read_csv(eds, sep="\t")

res = {}

grs = eds["start"].unique()

for m in eds["start"].unique(): 
    # get the mutations in the region
    region = eds[eds["start"] == m]

    pos_eds_in_region = list(region[region["label"] == "positive"]["edit_distance"])
    neg_eds_in_region = list(region[region["label"] != "positive"]["edit_distance"])
    res[m] = {"positive": pos_eds_in_region, "negative": neg_eds_in_region}


# sort genomic regions based on the start position
sorted_regions = sorted(res.keys(), key=lambda x: int(x))

fig, ax = plt.subplots()
for i, region in enumerate(sorted_regions):
    pos_eds = res[region]["positive"]
    neg_eds = res[region]["negative"]
    # plot boxplots next to each other
    # plot positive mutations
    ax.boxplot(pos_eds, positions=[i*2-0.5], widths=0.6, patch_artist=True, showfliers=False, boxprops=dict(facecolor="#E1BE6A"))
    # plot negative mutations
    ax.boxplot(neg_eds, positions=[i*2+0.5], widths=0.6, patch_artist=True, showfliers=False, boxprops=dict(facecolor="#40B0A6"))

pos_patch = mpatches.Patch(color='#E1BE6A', label='Positive')
neg_patch = mpatches.Patch(color='#40B0A6', label='Negative')
fig.legend(handles=[pos_patch, neg_patch], labels=['Positive', 'Negative'], 
           ncol=2, loc='upper center', bbox_to_anchor=(0.5, 0.95))

ax.set_xticks(range(0, len(sorted_regions)*2, 2))
updated = [x for x in sorted_xtick_regions if int(x.split("-")[0]) in sorted_regions]
ax.set_xticklabels(updated, rotation=90)
plt.savefig("plots_scripts/eds_distribution_lumc_tuples.pdf", format="pdf", bbox_inches="tight", dpi=300)
