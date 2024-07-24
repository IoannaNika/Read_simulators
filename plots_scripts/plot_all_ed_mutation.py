import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


mutations_train = "data/tuples_pacbio_sars_cov_2_rev_compl/dataset/samples_mutations.tsv"
ed_train = "data/tuples_pacbio_sars_cov_2_rev_compl/dataset/samples.tsv"
# load the mutations
mutations_train = pd.read_csv(mutations_train, sep="\t")
ed_train = pd.read_csv(ed_train, sep="\t")

lumc_test = "data/lumc_data/lumc_dataset_pairs_ed_n_m.tsv"

lumc_test = pd.read_csv(lumc_test, sep="\t")

# plot boxplot for training set mutations and edit distance across positive and negative samples and for the test set across positive and negative samples
# plot all boxplots in one figure and save it as a pdf
fig, ax = plt.subplots(1, 1, figsize=(10, 5))

# plot the boxplot for the training set
mutations_train_pos = mutations_train[mutations_train["label"] == "positive"]["mutations"]
mutations_train_neg = mutations_train[mutations_train["label"] != "positive"]["mutations"]
ed_train_pos = ed_train[ed_train["label"] == "positive"]["edit_distance"]
ed_train_neg = ed_train[ed_train["label"] != "positive"]["edit_distance"]
lumc_test_pos_ed = lumc_test[lumc_test["label"] == "positive"]["edit_distance"]
lumc_test_neg_ed = lumc_test[lumc_test["label"] != "positive"]["edit_distance"]
lumc_test_pos_mutations = lumc_test[lumc_test["label"] == "positive"]["n_mutations"]
lumc_test_neg_mutations = lumc_test[lumc_test["label"] != "positive"]["n_mutations"]

# plot the boxplot for the training set
# plot the boxplot for the training set
ax.boxplot([mutations_train_pos, mutations_train_neg, ed_train_pos, ed_train_neg, lumc_test_pos_ed, lumc_test_neg_ed, lumc_test_pos_mutations, lumc_test_neg_mutations],
            labels=["Mutation\nTraining\nPos", "Mutations\nTraining\nNeg", "Edit Distance\nTraining\nPos", "Edit Distance\nTraining\nNeg", "LUMC\nEdit Distance\nPos", "LUMC\nEdit Distance\nNeg", "LUMC\nMutations\nPos", "LUMC\nMutations\nNeg"],
            patch_artist=True, showfliers=False)

plt.grid()

plt.savefig("plots_scripts/all_ed_mutations.pdf")

