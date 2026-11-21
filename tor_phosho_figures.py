import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt
from upsetplot import generate_counts, plot, UpSet, from_contents

############################################# Phospho-proteome data

############ Time course
# Time course phospho-proteome dataset
path_to_tm_phospho_proteome = r"./re_analysis/phospho_proteome_time_course/"

# Control(Afb2) time course
ph_afb2_30m_tm = pd.read_csv(path_to_tm_phospho_proteome + 'afb2_time_course_phospho/afb2_30min_signi_time_course.tsv',
                             sep='\t')
ph_afb2_1h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'afb2_time_course_phospho/afb2_1hr_signi_time_course.tsv',
                            sep='\t')
ph_afb2_2h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'afb2_time_course_phospho/afb2_2hr_signi_time_course.tsv',
                            sep='\t')
ph_afb2_3h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'afb2_time_course_phospho/afb2_3hr_signi_time_course.tsv',
                            sep='\t')
ph_afb2_5h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'afb2_time_course_phospho/afb2_5hr_signi_time_course.tsv',
                            sep='\t')

# TORC1 time course
ph_torc1_30m_tm = pd.read_csv(
    path_to_tm_phospho_proteome + 'torc1_time_course_phospho/TORC1_30_min_signi_time_course.tsv', sep='\t')
ph_torc1_1h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'torc1_time_course_phospho/TORC1_1_hr_signi_time_course.tsv',
                             sep='\t')
ph_torc1_2h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'torc1_time_course_phospho/TORC1_2_hr_signi_time_course.tsv',
                             sep='\t')
ph_torc1_3h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'torc1_time_course_phospho/TORC1_3_hr_signi_time_course.tsv',
                             sep='\t')
ph_torc1_5h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'torc1_time_course_phospho/TORC1_5_hr_signi_time_course.tsv',
                             sep='\t')

# TORC2 time course
ph_torc2_30m_tm = pd.read_csv(
    path_to_tm_phospho_proteome + 'torc2_time_course_phospho/TORC2_30_min_signi_time_course.tsv', sep='\t')
ph_torc2_1h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'torc2_time_course_phospho/TORC2_1_hr_signi_time_course.tsv',
                             sep='\t')
ph_torc2_2h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'torc2_time_course_phospho/TORC2_2_hr_signi_time_course.tsv',
                             sep='\t')
ph_torc2_3h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'torc2_time_course_phospho/TORC2_3_hr_signi_time_course.tsv',
                             sep='\t')
ph_torc2_5h_tm = pd.read_csv(path_to_tm_phospho_proteome + 'torc2_time_course_phospho/TORC2_5_hr_signi_time_course.tsv',
                             sep='\t')

######### Pair Wise
path_to_pr_phospho_proteome = r"./re_analysis/phospho_proteome_pair_wise/"
# TORC1 pair wise
ph_torc1_0h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc1_pairwise_phospho/torc1_zero_hr_norm.tsv', sep='\t')
ph_torc1_30m_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc1_pairwise_phospho/torc1_30_min_norm.tsv', sep='\t')
ph_torc1_1h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc1_pairwise_phospho/torc1_1hr_norm.tsv', sep='\t')
ph_torc1_2h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc1_pairwise_phospho/torc1_2hr_norm.tsv', sep='\t')
ph_torc1_3h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc1_pairwise_phospho/torc1_3hr_norm.tsv', sep='\t')
ph_torc1_5h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc1_pairwise_phospho/torc1_5hr_norm.tsv', sep='\t')

# TORC2 pair wise
ph_torc2_0h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc2_pairwise_phospho/torc2_zero_hr_norm.tsv', sep='\t')
ph_torc2_30m_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc2_pairwise_phospho/torc2_30_min_norm.tsv', sep='\t')
ph_torc2_1h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc2_pairwise_phospho/torc2_1hr_norm.tsv', sep='\t')
ph_torc2_2h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc2_pairwise_phospho/torc2_2hr_norm.tsv', sep='\t')
ph_torc2_3h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc2_pairwise_phospho/torc2_3hr_norm.tsv', sep='\t')
ph_torc2_5h_pr = pd.read_csv(path_to_pr_phospho_proteome + 'torc2_pairwise_phospho/torc2_5hr_norm.tsv', sep='\t')

########################################################################################################################

# getting all hypo-phosphorylated from TORC1
ph_torc1_30m_pr_hypo = ph_torc1_30m_pr[ph_torc1_30m_pr['log2ratio_Condition2'] < 0]
ph_torc1_1h_pr_hypo = ph_torc1_1h_pr[ph_torc1_1h_pr['log2ratio_Condition2'] < 0]
ph_torc1_2h_pr_hypo = ph_torc1_2h_pr[ph_torc1_2h_pr['log2ratio_Condition2'] < 0]
ph_torc1_3h_pr_hypo = ph_torc1_3h_pr[ph_torc1_3h_pr['log2ratio_Condition2'] < 0]
ph_torc1_5h_pr_hypo = ph_torc1_5h_pr[ph_torc1_5h_pr['log2ratio_Condition2'] < 0]

# getting all hyper-phosphorylated from TORC1
ph_torc1_30m_pr_hyper = ph_torc1_30m_pr[ph_torc1_30m_pr['log2ratio_Condition2'] > 0]
ph_torc1_1h_pr_hyper = ph_torc1_1h_pr[ph_torc1_1h_pr['log2ratio_Condition2'] > 0]
ph_torc1_2h_pr_hyper = ph_torc1_2h_pr[ph_torc1_2h_pr['log2ratio_Condition2'] > 0]
ph_torc1_3h_pr_hyper = ph_torc1_3h_pr[ph_torc1_3h_pr['log2ratio_Condition2'] > 0]
ph_torc1_5h_pr_hyper = ph_torc1_5h_pr[ph_torc1_5h_pr['log2ratio_Condition2'] > 0]

# getting all hypo-phosphorylated from TORC2
ph_torc2_30m_pr_hypo = ph_torc2_30m_pr[ph_torc2_30m_pr['log2ratio_Condition2'] < 0]
ph_torc2_1h_pr_hypo = ph_torc2_1h_pr[ph_torc2_1h_pr['log2ratio_Condition2'] < 0]
ph_torc2_2h_pr_hypo = ph_torc2_2h_pr[ph_torc2_2h_pr['log2ratio_Condition2'] < 0]
ph_torc2_3h_pr_hypo = ph_torc2_3h_pr[ph_torc2_3h_pr['log2ratio_Condition2'] < 0]
ph_torc2_5h_pr_hypo = ph_torc2_5h_pr[ph_torc2_5h_pr['log2ratio_Condition2'] < 0]

# getting all hyper-phosphorylated from TORC2
ph_torc2_30m_pr_hyper = ph_torc2_30m_pr[ph_torc2_30m_pr['log2ratio_Condition2'] > 0]
ph_torc2_1h_pr_hyper = ph_torc2_1h_pr[ph_torc2_1h_pr['log2ratio_Condition2'] > 0]
ph_torc2_2h_pr_hyper = ph_torc2_2h_pr[ph_torc2_2h_pr['log2ratio_Condition2'] > 0]
ph_torc2_3h_pr_hyper = ph_torc2_3h_pr[ph_torc2_3h_pr['log2ratio_Condition2'] > 0]
ph_torc2_5h_pr_hyper = ph_torc2_5h_pr[ph_torc2_5h_pr['log2ratio_Condition2'] > 0]

path_to_save_phospho_graphs = r"./re_analysis/phospho_proteome_graphs/"

# Upset plot for TORC1 and TORC2

################################################### FOR TORC1
# Building upset plot for all-time points for deregulated phosphosites
torc1_deregulated = from_contents({"0.5 hr": ph_torc1_30m_pr['phospho_site'],
                                   "1 hr": ph_torc1_1h_pr['phospho_site'],
                                   "2 hr": ph_torc1_2h_pr['phospho_site'],
                                   "3 hr": ph_torc1_3h_pr["phospho_site"],
                                   "5 hr": ph_torc1_5h_pr["phospho_site"]})
torc1_deregulated_table = torc1_deregulated.reset_index()
torc1_deregulated_table.to_excel(path_to_save_phospho_graphs + "/TORC1_deregulated_upset.xlsx", header=True,
                                 index=False)
torc1_deregulated_plot = UpSet(torc1_deregulated, subset_size="count", show_counts=True, sort_by='degree',
                               sort_categories_by='-input')
torc1_deregulated_plot.plot()
plt.suptitle("TORC1 deregulated phospho-peptides")
plt.savefig(path_to_save_phospho_graphs + "TORC1_deregulated_upset_plot.tif", dpi=600)
plt.close()

# Building upset plot for all-time points for hypo phosphosites
torc1_hyper = from_contents({"0.5 hr hyper": ph_torc1_30m_pr_hyper['phospho_site'],
                             "1 hr hyper": ph_torc1_1h_pr_hyper['phospho_site'],
                             "2 hr hyper": ph_torc1_2h_pr_hyper['phospho_site'],
                             "3 hr hyper": ph_torc1_3h_pr_hyper["phospho_site"],
                             "5 hr hyper": ph_torc1_5h_pr_hyper["phospho_site"]})
torc1_hyper_table = torc1_hyper.reset_index()
torc1_hyper_table.to_excel(path_to_save_phospho_graphs + "/TORC1_hyper_upset.xlsx", header=True, index=False)
torc1_hyper_plot = UpSet(torc1_hyper, subset_size="count", show_counts=True, sort_by='degree',
                         sort_categories_by='-input')
torc1_hyper_plot.plot()
plt.suptitle("TORC1 Hyper-phosphorylated peptides")
plt.savefig(path_to_save_phospho_graphs + "TORC1_hyper_upset_plot.tif", dpi=600)
plt.close()

# Building upset plot for all-time points for hypo phosphosites
torc1_hypo = from_contents({"0.5 hr hypo": ph_torc1_30m_pr_hypo['phospho_site'],
                            "1 hr hypo": ph_torc1_1h_pr_hypo['phospho_site'],
                            "2 hr hypo": ph_torc1_2h_pr_hypo['phospho_site'],
                            "3 hr hypo": ph_torc1_3h_pr_hypo["phospho_site"],
                            "5 hr hypo": ph_torc1_5h_pr_hypo["phospho_site"]})
torc1_hypo_table = torc1_hypo.reset_index()
torc1_hypo_table.to_excel(path_to_save_phospho_graphs + "/TORC1_hypo_upset.xlsx", header=True, index=False)
torc1_hypo_plot = UpSet(torc1_hypo, subset_size="count", show_counts=True, sort_by='degree',
                        sort_categories_by='-input')
torc1_hypo_plot.plot()
plt.suptitle("TORC1 Hypo-phosphorylated peptides")
plt.savefig(path_to_save_phospho_graphs + "TORC1_hypo_upset_plot.tif", dpi=600)
plt.close()

################################################### FOR TORC2

# Building upset plot for all-time points for deregulated phosphosites
torc2_deregulated = from_contents({"0.5 hr": ph_torc2_30m_pr['phospho_site'],
                                   "1 hr": ph_torc2_1h_pr['phospho_site'],
                                   "2 hr": ph_torc2_2h_pr['phospho_site'],
                                   "3 hr": ph_torc2_3h_pr["phospho_site"],
                                   "5 hr": ph_torc2_5h_pr["phospho_site"]})
torc2_deregulated_table = torc2_deregulated.reset_index()
torc2_deregulated_table.to_excel(path_to_save_phospho_graphs + "/TORC2_deregulated_upset.xlsx", header=True,
                                 index=False)
torc2_deregulated_plot = UpSet(torc2_deregulated, subset_size="count", show_counts=True, sort_by='degree',
                               sort_categories_by='-input')
torc2_deregulated_plot.plot()
plt.suptitle("TORC2 deregulated phospho-peptides")
plt.savefig(path_to_save_phospho_graphs + "TORC2_deregulated_upset_plot.tif", dpi=600)
plt.close()

# Building upset plot for all-time points for hypo phosphosites
torc2_hyper = from_contents({"0.5 hr hyper": ph_torc2_30m_pr_hyper['phospho_site'],
                             "1 hr hyper": ph_torc2_1h_pr_hyper['phospho_site'],
                             "2 hr hyper": ph_torc2_2h_pr_hyper['phospho_site'],
                             "3 hr hyper": ph_torc2_3h_pr_hyper["phospho_site"],
                             "5 hr hyper": ph_torc2_5h_pr_hyper["phospho_site"]})
torc2_hyper_table = torc2_hyper.reset_index()
torc2_hyper_table.to_excel(path_to_save_phospho_graphs + "/TORC2_hyper_upset.xlsx", header=True, index=False)
torc2_hyper_plot = UpSet(torc2_hyper, subset_size="count", show_counts=True, sort_by='degree',
                         sort_categories_by='-input')
torc2_hyper_plot.plot()
plt.suptitle("TORC2 Hyper-phosphorylated peptides")
plt.savefig(path_to_save_phospho_graphs + "TORC2_hyper_upset_plot.tif", dpi=600)
plt.close()

# Building upset plot for all-time points for hypo phosphosites
torc2_hypo = from_contents({"0.5 hr hypo": ph_torc2_30m_pr_hypo['phospho_site'],
                            "1 hr hypo": ph_torc2_1h_pr_hypo['phospho_site'],
                            "2 hr hypo": ph_torc2_2h_pr_hypo['phospho_site'],
                            "3 hr hypo": ph_torc2_3h_pr_hypo["phospho_site"],
                            "5 hr hypo": ph_torc2_5h_pr_hypo["phospho_site"]})
torc2_hypo_table = torc2_hypo.reset_index()
torc2_hypo_table.to_excel(path_to_save_phospho_graphs + "/TORC2_hypo_upset.xlsx", header=True, index=False)
torc2_hypo_plot = UpSet(torc2_hypo, subset_size="count", show_counts=True, sort_by='degree',
                        sort_categories_by='-input')
torc2_hypo_plot.plot()
plt.suptitle("TORC2 Hypo-phosphorylated peptides")
plt.savefig(path_to_save_phospho_graphs + "TORC2_hypo_upset_plot.tif", dpi=600)
plt.close()

# Venn diagram for 30min-1hr common and 1hr-2hr common
################################################### FOR TORC1

# For 30 min and 1 hr
ph_torc1_30m_1h_com = ph_torc1_30m_pr.merge(ph_torc1_1h_pr, on='phospho_site', how='inner', suffixes=('_30min', '_1hr'))
ph_torc1_30m_1h_com.to_csv(path_to_save_phospho_graphs + 'TORC1_30m_1hr_common_phospho.tsv', index=False, header=True,
                           sep='\t')

ph_torc1_30m_only_comp_1h = set(ph_torc1_30m_pr['phospho_site']) - set(ph_torc1_1h_pr['phospho_site'])
ph_torc1_1h_only_comp_30m = set(ph_torc1_1h_pr['phospho_site']) - set(ph_torc1_30m_pr['phospho_site'])

# Plotting venn diagram
a = venn2(subsets=(len(ph_torc1_30m_only_comp_1h), len(ph_torc1_1h_only_comp_30m), ph_torc1_30m_1h_com.shape[0]),
          set_labels=('0.5 hr', '1 hr'), set_colors=('orange', 'blue'), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(25)
    text.set_weight('bold')
    text.set_fontname('arial')
for text in a.set_labels:
    text.set_fontsize(20)
    text.set_weight('bold')
    text.set_fontname('arial')
plt.title("TORC1 phospho common between 30 min and 1 hr")
plt.savefig(path_to_save_phospho_graphs + 'TORC1_30m_1h_common_venn.tif', dpi=600)
plt.close()

# For 1 hr and 2 hr

ph_torc1_1h_2h_com = ph_torc1_1h_pr.merge(ph_torc1_2h_pr, on='phospho_site', how='inner', suffixes=('_1hr', '_2hr'))
ph_torc1_1h_2h_com.to_csv(path_to_save_phospho_graphs + 'TORC1_1h_2h_common_phospho.tsv', index=False, header=True,
                          sep='\t')

ph_torc1_1h_only_comp_2h = set(ph_torc1_1h_pr['phospho_site']) - set(ph_torc1_2h_pr['phospho_site'])
ph_torc1_2h_only_comp_1h = set(ph_torc1_2h_pr['phospho_site']) - set(ph_torc1_1h_pr['phospho_site'])

# Plotting venn diagram
a = venn2(subsets=(len(ph_torc1_1h_only_comp_2h), len(ph_torc1_2h_only_comp_1h), ph_torc1_1h_2h_com.shape[0]),
          set_labels=('1 hr', '2 hr'), set_colors=('orange', 'blue'), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(25)
    text.set_weight('bold')
    text.set_fontname('arial')
for text in a.set_labels:
    text.set_fontsize(20)
    text.set_weight('bold')
    text.set_fontname('arial')
plt.title("TORC1 phospho common between 1 hr & 2 hr")
plt.savefig(path_to_save_phospho_graphs + 'TORC1_1h_2h_common_venn.tif', dpi=600)
plt.close()

######################################### FOR TORC2

# For 30 min and 1 hr
ph_torc2_30m_1h_com = ph_torc2_30m_pr.merge(ph_torc2_1h_pr, on='phospho_site', how='inner', suffixes=('_30min', '_1hr'))
ph_torc2_30m_1h_com.to_csv(path_to_save_phospho_graphs + 'TORC2_30m_1hr_common_phospho.tsv', index=False, header=True,
                           sep='\t')

ph_torc2_30m_only_comp_1h = set(ph_torc2_30m_pr['phospho_site']) - set(ph_torc2_1h_pr['phospho_site'])
ph_torc2_1h_only_comp_30m = set(ph_torc2_1h_pr['phospho_site']) - set(ph_torc2_30m_pr['phospho_site'])

# Plotting venn diagram
a = venn2(subsets=(len(ph_torc2_30m_only_comp_1h), len(ph_torc2_1h_only_comp_30m), ph_torc2_30m_1h_com.shape[0]),
          set_labels=('0.5 hr', '1 hr'), set_colors=('orange', 'blue'), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(25)
    text.set_weight('bold')
    text.set_fontname('arial')
for text in a.set_labels:
    text.set_fontsize(20)
    text.set_weight('bold')
    text.set_fontname('arial')
plt.title("TORC2 phospho common between 30 min and 1 hr")
plt.savefig(path_to_save_phospho_graphs + 'TORC2_30m_1h_common_venn.tif', dpi=600)
plt.close()

# For 1 hr and 2 hr

ph_torc2_1h_2h_com = ph_torc2_1h_pr.merge(ph_torc2_2h_pr, on='phospho_site', how='inner', suffixes=('_1hr', '_2hr'))
ph_torc2_1h_2h_com.to_csv(path_to_save_phospho_graphs + 'TORC2_1h_2h_common_phospho.tsv', index=False, header=True,
                          sep='\t')

ph_torc2_1h_only_comp_2h = set(ph_torc2_1h_pr['phospho_site']) - set(ph_torc2_2h_pr['phospho_site'])
ph_torc2_2h_only_comp_1h = set(ph_torc2_2h_pr['phospho_site']) - set(ph_torc2_1h_pr['phospho_site'])

# Plotting venn diagram
a = venn2(subsets=(len(ph_torc2_1h_only_comp_2h), len(ph_torc2_2h_only_comp_1h), ph_torc2_1h_2h_com.shape[0]),
          set_labels=('1 hr', '2 hr'), set_colors=('orange', 'blue'), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(25)
    text.set_weight('bold')
    text.set_fontname('arial')
for text in a.set_labels:
    text.set_fontsize(20)
    text.set_weight('bold')
    text.set_fontname('arial')
plt.title("TORC2 phospho common between 1 hr & 2 hr")
plt.savefig(path_to_save_phospho_graphs + 'TORC2_1h_2h_common_venn.tif', dpi=600)
plt.close()


################################ Heatmap for TORC1 and TORC2


def df_for_heatmap_prt(data, compare_data, value_name):
    df1 = pd.pivot_table(data=data, index='phospho_site', columns='time_stamp', values=value_name)
    # df1.fillna(0, inplace=True)
    df2 = df1.loc[:, ['30_min', '1_hr', '2_hr', '3_hr', '5_hr']]
    df3 = df2[df2.index.isin(compare_data['phospho_site'].to_list())]
    df4 = df3.sort_values('30_min', ascending=True)
    return df4


def add_cross_at_na(data):
    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            if pd.isna(data.iloc[y, x]):
                ax.plot([x, x + 1], [y, y + 1], color='gray', lw=0.5)
                ax.plot([x, x + 1], [y + 1, y], color='gray', lw=0.5)


########## For TORC1

torc1_all_time_data = pd.concat([
    ph_torc1_30m_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
    ph_torc1_1h_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
    ph_torc1_2h_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
    ph_torc1_3h_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
    ph_torc1_5h_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
])

# Comparing TORC1 30m-1hr common with rest time points
ph_torc1_30m_1h_com_fold_change = df_for_heatmap_prt(torc1_all_time_data, ph_torc1_30m_1h_com,
                                                     'log2ratio_Condition2').add_suffix("_log2ratio")
ph_torc1_30m_1h_com_pvalue = df_for_heatmap_prt(torc1_all_time_data, ph_torc1_30m_1h_com,
                                                'pValue_Condition2').add_suffix("_pValue")
ph_torc1_30m_1h_com_all_data = pd.concat([ph_torc1_30m_1h_com_fold_change, ph_torc1_30m_1h_com_pvalue], axis=1)

ph_torc1_30m_1h_com_all_data_down = ph_torc1_30m_1h_com_all_data[
    ph_torc1_30m_1h_com_all_data['30_min_log2ratio'] < 0].reset_index()
ph_torc1_30m_1h_com_all_data_down.to_csv(path_to_save_phospho_graphs + "ph_torc1_30m_1h_com_all_data_down.tsv",
                                         header=True, index=False, sep='\t')
ph_torc1_30m_1h_com_all_data_up = ph_torc1_30m_1h_com_all_data[
    ph_torc1_30m_1h_com_all_data['30_min_log2ratio'] > 0].reset_index()
ph_torc1_30m_1h_com_all_data_up.to_csv(path_to_save_phospho_graphs + "ph_torc1_30m_1h_com_all_data_up.tsv", header=True,
                                       index=False, sep='\t')

# Heatmap for TORC1 down
ph_torc1_30m_1h_com_all_data_down_heatmap = ph_torc1_30m_1h_com_all_data_down[['phospho_site', '30_min_log2ratio',
                                                                               '1_hr_log2ratio', '2_hr_log2ratio',
                                                                               '3_hr_log2ratio',
                                                                               '5_hr_log2ratio']].set_index(
    'phospho_site')
plt.figure(figsize=(25, 3))
ax = sns.heatmap(ph_torc1_30m_1h_com_all_data_down_heatmap.T,
                 xticklabels=True, yticklabels=['0.5 hr', '1 hr', '2 hr', '3 hr', '5 hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=False, linewidths=0.5)
# ax.xaxis.tick_top()
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")
ax.set_facecolor('white')
add_cross_at_na(ph_torc1_30m_1h_com_all_data_down_heatmap.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=2, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC1 hypo phosphorylated sites common from 30min & 1hr')
plt.subplots_adjust(bottom=0.2)
plt.savefig(path_to_save_phospho_graphs + "ph_torc1_30m_1h_com_all_data_down_heatmap.tif", dpi=600)
plt.close()

# Heatmap for TORC1 UP
ph_torc1_30m_1h_com_all_data_up_heatmap = ph_torc1_30m_1h_com_all_data_up[['phospho_site', '30_min_log2ratio',
                                                                           '1_hr_log2ratio', '2_hr_log2ratio',
                                                                           '3_hr_log2ratio',
                                                                           '5_hr_log2ratio']].set_index('phospho_site')
plt.figure(figsize=(25, 3))
ax = sns.heatmap(ph_torc1_30m_1h_com_all_data_up_heatmap.T,
                 xticklabels=True, yticklabels=['0.5 hr', '1 hr', '2 hr', '3 hr', '5 hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=False, linewidths=0.5)
# ax.xaxis.tick_top()
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")

ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.set_facecolor('white')
add_cross_at_na(ph_torc1_30m_1h_com_all_data_up_heatmap.T)
ax.tick_params(axis='x', which='major', labelsize=2, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC1 hyper phosphorylated sites common from 30min & 1hr')
plt.subplots_adjust(bottom=0.24)
plt.savefig(path_to_save_phospho_graphs + "ph_torc1_30m_1h_com_all_data_up_heatmap.tif", dpi=600)
plt.close()

# Heatmap for TORC1 commonly deregulated in all time points

ph_torc1_all_time_deregulated = \
torc1_deregulated_table[torc1_deregulated_table['0.5 hr'] & torc1_deregulated_table['1 hr']
                        & torc1_deregulated_table['2 hr'] & torc1_deregulated_table['3 hr']
                        & torc1_deregulated_table['5 hr']]['id'].reset_index(drop=True)
ph_torc1_all_time_deregulated_df = ph_torc1_all_time_deregulated.to_frame('phospho_site')
torc1_all_time_deregulated = df_for_heatmap_prt(torc1_all_time_data, ph_torc1_all_time_deregulated_df,
                                                'log2ratio_Condition2').add_suffix("_log2ratio")

plt.figure(figsize=(15, 5))
ax = sns.heatmap(torc1_all_time_deregulated.T,
                 xticklabels=True, yticklabels=['0.5 hr', '1 hr', '2 hr', '3 hr', '5 hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5)
ax.xaxis.tick_top()
ax.set_xlabel("phospho sites")
ax.set_ylabel("Time")
ax.set_facecolor('white')
add_cross_at_na(torc1_all_time_deregulated.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=6, labelbottom=False, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=7)
ax.set_title('TORC1 all time deregulated')
plt.savefig(path_to_save_phospho_graphs + "ph_torc1_all_time_deregulated_heatmap_1.tif", dpi=600)
plt.close()

########## For TORC2

torc2_all_time_data = pd.concat([
    ph_torc2_30m_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
    ph_torc2_1h_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
    ph_torc2_2h_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
    ph_torc2_3h_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
    ph_torc2_5h_pr[['time_stamp', 'log2ratio_Condition2', 'phospho_site', 'pValue_Condition2']],
])

# Comparing TORC1 30m-1hr common with rest time points
ph_torc2_30m_1h_com_fold_change = df_for_heatmap_prt(torc2_all_time_data, ph_torc2_30m_1h_com,
                                                     'log2ratio_Condition2').add_suffix("_log2ratio")
ph_torc2_30m_1h_com_pvalue = df_for_heatmap_prt(torc2_all_time_data, ph_torc2_30m_1h_com,
                                                'pValue_Condition2').add_suffix("_pValue")
ph_torc2_30m_1h_com_all_data = pd.concat([ph_torc2_30m_1h_com_fold_change, ph_torc2_30m_1h_com_pvalue], axis=1)
ph_torc2_30m_1h_com_all_data.to_csv(path_to_save_phospho_graphs + "ph_torc2_30m_1h_com_all_data.tsv", sep='\t')

ph_torc2_30m_1h_com_all_data_down = ph_torc2_30m_1h_com_all_data[
    ph_torc2_30m_1h_com_all_data['30_min_log2ratio'] < 0].reset_index()
ph_torc2_30m_1h_com_all_data_down.to_csv(path_to_save_phospho_graphs + "ph_torc2_30m_1h_com_all_data_down.tsv",
                                         header=True, index=False, sep='\t')
ph_torc2_30m_1h_com_all_data_up = ph_torc2_30m_1h_com_all_data[
    ph_torc2_30m_1h_com_all_data['30_min_log2ratio'] > 0].reset_index()
ph_torc2_30m_1h_com_all_data_up.to_csv(path_to_save_phospho_graphs + "ph_torc2_30m_1h_com_all_data_up.tsv", header=True,
                                       index=False, sep='\t')

# Heatmap for TORC1 down
ph_torc2_30m_1h_com_all_data_down_heatmap = ph_torc2_30m_1h_com_all_data_down[['phospho_site', '30_min_log2ratio',
                                                                               '1_hr_log2ratio', '2_hr_log2ratio',
                                                                               '3_hr_log2ratio',
                                                                               '5_hr_log2ratio']].set_index(
    'phospho_site')
plt.figure(figsize=(25, 3))
ax = sns.heatmap(ph_torc2_30m_1h_com_all_data_down_heatmap.T,
                 xticklabels=True, yticklabels=['0.5 hr', '1 hr', '2 hr', '3 hr', '5 hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=False, linewidths=0.5)
# ax.xaxis.tick_top()
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")
ax.set_facecolor('white')
add_cross_at_na(ph_torc2_30m_1h_com_all_data_down_heatmap.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=5, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC2 hypo phosphorylated sites common from 30min & 1hr')
plt.subplots_adjust(bottom=0.3)
plt.savefig(path_to_save_phospho_graphs + "ph_torc2_30m_1h_com_all_data_down_heatmap.tif", dpi=600)
plt.close()

# Heatmap for TORC1 UP
ph_torc2_30m_1h_com_all_data_up_heatmap = ph_torc2_30m_1h_com_all_data_up[['phospho_site', '30_min_log2ratio',
                                                                           '1_hr_log2ratio', '2_hr_log2ratio',
                                                                           '3_hr_log2ratio',
                                                                           '5_hr_log2ratio']].set_index('phospho_site')
plt.figure(figsize=(20, 3))
ax = sns.heatmap(ph_torc2_30m_1h_com_all_data_up_heatmap.T,
                 xticklabels=True, yticklabels=['0.5 hr', '1 hr', '2 hr', '3 hr', '5 hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=False, linewidths=0.5)
# ax.xaxis.tick_top()
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")
ax.set_facecolor('white')
add_cross_at_na(ph_torc2_30m_1h_com_all_data_up_heatmap.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=2, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC2 hyper phosphorylated sites common from 30min & 1hr')
plt.subplots_adjust(bottom=0.23)
plt.savefig(path_to_save_phospho_graphs + "ph_torc2_30m_1h_com_all_data_up_heatmap.tif", dpi=600)
plt.close()

# Heatmap for TORC2 commonly deregulated in all time points

ph_torc2_all_time_deregulated = \
torc2_deregulated_table[torc2_deregulated_table['0.5 hr'] & torc2_deregulated_table['1 hr']
                        & torc2_deregulated_table['2 hr'] & torc2_deregulated_table['3 hr']
                        & torc2_deregulated_table['5 hr']]['id'].reset_index(drop=True)
ph_torc2_all_time_deregulated_df = ph_torc2_all_time_deregulated.to_frame('phospho_site')
torc2_all_time_deregulated = df_for_heatmap_prt(torc2_all_time_data, ph_torc2_all_time_deregulated_df,
                                                'log2ratio_Condition2').add_suffix("_log2ratio")

plt.figure(figsize=(15, 5))
ax = sns.heatmap(torc2_all_time_deregulated.T,
                 xticklabels=True, yticklabels=['0.5 hr', '1 hr', '2 hr', '3 hr', '5 hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5)
ax.xaxis.tick_top()
ax.set_xlabel("phospho sites")
ax.set_ylabel("Time")
ax.set_facecolor('white')
add_cross_at_na(torc2_all_time_deregulated.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=7, labelbottom=False, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=7)
ax.set_title('TORC2 all time deregulated')
plt.savefig(path_to_save_phospho_graphs + "ph_torc2_all_time_deregulated_heatmap_1.tif", dpi=600)
plt.close()

########################################## TORC1 & TORC2 common

ph_torc1_torc2_30m_1h_com = ph_torc1_30m_1h_com.merge(ph_torc2_30m_1h_com, how='inner', on="phospho_site",
                                                      suffixes=("_TORC1", "_TORC2"))
ph_torc1_torc2_30m_1h_com.to_csv(path_to_save_phospho_graphs + "TORC1_TORC2_30m_1hr_common.tsv", index=False,
                                 header=True, sep='\t')

ph_torc1_only_compare_torc2_30m_1h = set(ph_torc1_30m_1h_com["phospho_site"]) - set(ph_torc2_30m_1h_com["phospho_site"])
ph_torc2_only_compare_torc1_30m_1h = set(ph_torc2_30m_1h_com["phospho_site"]) - set(ph_torc1_30m_1h_com["phospho_site"])

# Plotting venn diagram
c = venn2(subsets=(
len(ph_torc1_only_compare_torc2_30m_1h), len(ph_torc2_only_compare_torc1_30m_1h), ph_torc1_torc2_30m_1h_com.shape[0]),
          set_labels=('TORC1', 'TORC2'), set_colors=('orange', 'blue'), alpha=0.6)
for text in c.subset_labels:
    text.set_fontsize(25)
    text.set_fontname('arial')
for text in c.set_labels:
    text.set_fontsize(20)
    text.set_fontname('arial')
plt.title("TORC1 and TORC2 common phospho_sites(30m-1hr)")
plt.savefig(path_to_save_phospho_graphs + 'TORC1_TORC2_common_phospho_sites_30m-1hr_venn.tif', dpi=600)
plt.close()

# TORC1-TORC2 common in TORC1 dataset heatmap
ph_torc1_torc2_common_in_torc1 = df_for_heatmap_prt(torc1_all_time_data, ph_torc1_torc2_30m_1h_com,
                                                    'log2ratio_Condition2').add_suffix("_log2ratio").reset_index()

ph_torc1_torc2_common_in_torc1_heatmap = ph_torc1_torc2_common_in_torc1[['phospho_site', '30_min_log2ratio',
                                                                         '1_hr_log2ratio', '2_hr_log2ratio',
                                                                         '3_hr_log2ratio', '5_hr_log2ratio']].set_index(
    'phospho_site')
plt.figure(figsize=(25, 3))
ax = sns.heatmap(ph_torc1_torc2_common_in_torc1_heatmap.T,
                 xticklabels=True, yticklabels=['0.5 hr', '1 hr', '2 hr', '3 hr', '5 hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5)
ax.xaxis.tick_top()
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")
ax.set_facecolor('white')
add_cross_at_na(ph_torc1_torc2_common_in_torc1_heatmap.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=5, labelbottom=False, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=5)
ax.set_title('TORC1 & TORC2 common phosphorylated sites in TORC1 dataset')
plt.savefig(path_to_save_phospho_graphs + "ph_torc1_torc2_common_in_torc1_heatmap.tif", dpi=600)
plt.close()

# TORC1-TORC2 common in TORC2 dataset heatmap
ph_torc1_torc2_common_in_torc2 = df_for_heatmap_prt(torc2_all_time_data, ph_torc1_torc2_30m_1h_com,
                                                    'log2ratio_Condition2').add_suffix("_log2ratio").reset_index()

ph_torc1_torc2_common_in_torc2_heatmap = ph_torc1_torc2_common_in_torc2[['phospho_site', '30_min_log2ratio',
                                                                         '1_hr_log2ratio', '2_hr_log2ratio',
                                                                         '3_hr_log2ratio', '5_hr_log2ratio']].set_index(
    'phospho_site')
plt.figure(figsize=(25, 3))
ax = sns.heatmap(ph_torc1_torc2_common_in_torc2_heatmap.T,
                 xticklabels=True, yticklabels=['0.5 hr', '1 hr', '2 hr', '3 hr', '5 hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5)
ax.xaxis.tick_top()
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")
ax.set_facecolor('white')
add_cross_at_na(ph_torc1_torc2_common_in_torc2_heatmap.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=5, labelbottom=False, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=5)
ax.set_title('TORC1 & TORC2 common phosphorylated sites in TORC2 dataset')
plt.savefig(path_to_save_phospho_graphs + "ph_torc1_torc2_common_in_torc2_heatmap.tif", dpi=600)
plt.close()

#########################################################################
###### Loading proteome data
# TORC1 proteome data
path_to_proteome_data = r"./re_analysis\proteome_pair_wise/"
prt_torc1_30min = pd.read_csv(path_to_proteome_data + '/torc1_proteome/prt_TORC1_30m_norm.tsv', sep='\t')
prt_torc1_1hr = pd.read_csv(path_to_proteome_data + '/torc1_proteome/prt_TORC1_1h_norm.tsv', sep='\t')
prt_torc1_2hr = pd.read_csv(path_to_proteome_data + '/torc1_proteome/prt_TORC1_2h_norm.tsv', sep='\t')
prt_torc1_3hr = pd.read_csv(path_to_proteome_data + '/torc1_proteome/prt_TORC1_3h_norm.tsv', sep='\t')
prt_torc1_5hr = pd.read_csv(path_to_proteome_data + '/torc1_proteome/prt_TORC1_5h_norm.tsv', sep='\t')

prt_torc1_all_time_data = pd.concat([
    prt_torc1_30min[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc1_1hr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc1_2hr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc1_3hr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc1_5hr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']]
])

prt_torc1_all_time_data_pivot = pd.pivot_table(data=prt_torc1_all_time_data, index='geneName', columns='time_stamp',
                                               values='log2ratio_Condition2')
prt_torc1_all_time_data_pivot.fillna(0, inplace=True)
prt_torc1_all_time_data_pivot.sort_values('30_min', inplace=True)
prt_torc1_all_time_data_pivot1 = prt_torc1_all_time_data_pivot.loc[:,
                                 ['30_min', '1_hr', '2_hr', '3_hr', '5_hr']].reset_index()

# TORC2 proteome data
prt_torc2_30min = pd.read_csv(path_to_proteome_data + '/torc2_proteome/prt_TORC2_30m_norm.tsv', sep='\t')
prt_torc2_1hr = pd.read_csv(path_to_proteome_data + '/torc2_proteome/prt_TORC2_1h_norm.tsv', sep='\t')
prt_torc2_2hr = pd.read_csv(path_to_proteome_data + '/torc2_proteome/prt_TORC2_2h_norm.tsv', sep='\t')
prt_torc2_3hr = pd.read_csv(path_to_proteome_data + '/torc2_proteome/prt_TORC2_3h_norm.tsv', sep='\t')
prt_torc2_5hr = pd.read_csv(path_to_proteome_data + '/torc2_proteome/prt_TORC2_5h_norm.tsv', sep='\t')

prt_torc2_all_time_data = pd.concat([
    prt_torc2_30min[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc2_1hr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc2_2hr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc2_3hr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc2_5hr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']]
])

prt_torc2_all_time_data_pivot = pd.pivot_table(data=prt_torc2_all_time_data, index='geneName', columns='time_stamp',
                                               values='log2ratio_Condition2')
prt_torc2_all_time_data_pivot.fillna(0, inplace=True)
prt_torc2_all_time_data_pivot.sort_values('30_min', inplace=True)
prt_torc2_all_time_data_pivot1 = prt_torc2_all_time_data_pivot.loc[:, ['30_min', '1_hr', '2_hr', '3_hr', '5_hr']]
prt_torc2_all_time_data_pivot1.to_csv(path_to_save_phospho_graphs + 'protein_torc2_all_time_data.tsv', sep='\t')
# Looking at TORC1 30m-1hr phosphosites distribution in alltime points in phospho and total proteome data

# torc1 phospho hyper and corresponding proteome change
ph_torc1_30m_1h_com_all_data_up1 = ph_torc1_30m_1h_com_all_data_up.reset_index()
ph_torc1_30m_1h_com_all_data_up1['geneName'] = ph_torc1_30m_1h_com_all_data_up1['phospho_site'].str.split('_').str[0]

torc1_phopho_proteome_up = ph_torc1_30m_1h_com_all_data_up1.merge(prt_torc1_all_time_data_pivot1, how='left',
                                                                  on='geneName', suffixes=('_phospho', '_proteome'))
torc1_phopho_proteome_up2 = torc1_phopho_proteome_up.sort_values('phospho_site')

# Heatmap
torc1_phopho_proteome_heatmap = torc1_phopho_proteome_up2[
    ['phospho_site', '30_min_log2ratio', '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio',
     '30_min', '1_hr', '2_hr', '3_hr', '5_hr']].set_index('phospho_site')

# plt.figure(figsize=(15, 25))
plt.figure(figsize=(25, 5))
ax = sns.heatmap(torc1_phopho_proteome_heatmap.T,
                 xticklabels=True,
                 yticklabels=['0.5hr_phospho', '1hr_phospho', '2hr_phospho', '3hr_phospho', '5hr_phospho',
                              '0.5hr_prot', '1hr_prot', '2hr_prot', '3hr_prot', '5hr_prot'],
                 cmap=sns.color_palette("PRGn", as_cmap=True),
                 cbar=True, cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False, linewidths=0.5
                 )
# ax.xaxis.tick_top()
ax.set_facecolor('white')
add_cross_at_na(torc1_phopho_proteome_heatmap.T)
ax.hlines([5], *ax.get_xlim())
ax.set_xlabel("phospho sites")
ax.set_ylabel("Time")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=2, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC1 hyper phosphorylated sites common from 30min & 1hr with its corresponding proteome change')
plt.savefig(path_to_save_phospho_graphs + "torc1_30m_1h_hyper_phospho_protein_heatmap_horizontal.tif", dpi=600)
plt.close()

# # Heatmap for torc1 phospho hyper and corresponding proteome change


# torc1 phospho hyper and corresponding proteome change
ph_torc1_30m_1h_com_all_data_down1 = ph_torc1_30m_1h_com_all_data_down.reset_index()
ph_torc1_30m_1h_com_all_data_down1['geneName'] = ph_torc1_30m_1h_com_all_data_down1['phospho_site'].str.split('_').str[
    0]

torc1_phopho_proteome_down = ph_torc1_30m_1h_com_all_data_down1.merge(prt_torc1_all_time_data_pivot1, how='left',
                                                                      on='geneName', suffixes=('_phospho', '_proteome'))
torc1_phopho_proteome_down2 = torc1_phopho_proteome_down.sort_values('phospho_site')

# Heatmap
torc1_phopho_proteome_heatmap_down = torc1_phopho_proteome_down2[
    ['phospho_site', '30_min_log2ratio', '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio',
     '30_min', '1_hr', '2_hr', '3_hr', '5_hr']].set_index('phospho_site')

plt.figure(figsize=(25, 5))
# plt.figure(figsize=(25, 5))
ax = sns.heatmap(torc1_phopho_proteome_heatmap_down.T,
                 xticklabels=True,
                 yticklabels=['0.5hr_phospho', '1hr_phospho', '2hr_phospho', '3hr_phospho', '5hr_phospho',
                              '0.5hr_prot', '1hr_prot', '2hr_prot', '3hr_prot', '5hr_prot'],
                 cmap=sns.color_palette("PRGn", as_cmap=True),
                 cbar=True, cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False, linewidths=0.5
                 )
# ax.xaxis.tick_top()
ax.hlines([5], *ax.get_xlim())
ax.set_xlabel("phospho_site")
ax.set_ylabel("Time")
ax.set_facecolor('white')
add_cross_at_na(torc1_phopho_proteome_heatmap_down.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=2, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC1 hypo phosphorylated sites common from 30min & 1hr with its corresponding proteome change')
plt.subplots_adjust(bottom=0.25)
plt.savefig(path_to_save_phospho_graphs + "torc1_30m_1h_hypo_phospho_protein_heatmap_horizontal.tif", dpi=600)
plt.close()

##########################################################################
## TORC2

# Looking at TORC2 30m-1hr phosphosites distribution in alltime points in phospho and total proteome data

# torc2 phospho hyper and corresponding proteome change
ph_torc2_30m_1h_com_all_data_up1 = ph_torc2_30m_1h_com_all_data_up.reset_index()
ph_torc2_30m_1h_com_all_data_up1['geneName'] = ph_torc2_30m_1h_com_all_data_up1['phospho_site'].str.split('_').str[0]

torc2_phopho_proteome_up = ph_torc2_30m_1h_com_all_data_up1.merge(prt_torc2_all_time_data_pivot1, how='left',
                                                                  on='geneName', suffixes=('_phospho', '_proteome'))
torc2_phopho_proteome_up2 = torc2_phopho_proteome_up.sort_values('phospho_site')

# Heatmap
torc2_phopho_proteome_heatmap = torc2_phopho_proteome_up2[
    ['phospho_site', '30_min_log2ratio', '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio',
     '30_min', '1_hr', '2_hr', '3_hr', '5_hr']].set_index('phospho_site')

# plt.figure(figsize=(15, 25))
plt.figure(figsize=(25, 5))
ax = sns.heatmap(torc2_phopho_proteome_heatmap.T,
                 xticklabels=True,
                 yticklabels=['0.5hr_phospho', '1hr_phospho', '2hr_phospho', '3hr_phospho', '5hr_phospho',
                              '0.5hr_prot', '1hr_prot', '2hr_prot', '3hr_prot', '5hr_prot'],
                 cmap=sns.color_palette("PRGn", as_cmap=True),
                 cbar=True, cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False, linewidths=0.5
                 )
# ax.xaxis.tick_top()
ax.hlines([5], *ax.get_xlim())
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")
ax.set_facecolor('white')
add_cross_at_na(torc2_phopho_proteome_heatmap.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=2, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC2 hyper phosphorylated sites common from 30min & 1hr with its corresponding proteome change')
plt.subplots_adjust(bottom=0.2)
plt.savefig(path_to_save_phospho_graphs + "torc2_30m_1h_hyper_phospho_protein_heatmap_horizontal.tif", dpi=600)
plt.close()

# torc2 phospho hypo and corresponding proteome change
ph_torc2_30m_1h_com_all_data_down1 = ph_torc2_30m_1h_com_all_data_down.reset_index()
ph_torc2_30m_1h_com_all_data_down1['geneName'] = ph_torc2_30m_1h_com_all_data_down1['phospho_site'].str.split('_').str[
    0]

torc2_phopho_proteome_down = ph_torc2_30m_1h_com_all_data_down1.merge(prt_torc2_all_time_data_pivot1, how='left',
                                                                      on='geneName', suffixes=('_phospho', '_proteome'))
torc2_phopho_proteome_down2 = torc2_phopho_proteome_down.sort_values('phospho_site')

# Heatmap
torc2_phopho_proteome_heatmap_down = torc2_phopho_proteome_down2[
    ['phospho_site', '30_min_log2ratio', '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio',
     '30_min', '1_hr', '2_hr', '3_hr', '5_hr']].set_index('phospho_site')

# plt.figure(figsize=(15, 25))
plt.figure(figsize=(25, 5))
ax = sns.heatmap(torc2_phopho_proteome_heatmap_down.T,
                 xticklabels=True,
                 yticklabels=['0.5hr_phospho', '1hr_phospho', '2hr_phospho', '3hr_phospho', '5hr_phospho',
                              '0.5hr_prot', '1hr_prot', '2hr_prot', '3hr_prot', '5hr_prot'],
                 cmap=sns.color_palette("PRGn", as_cmap=True),
                 cbar=False, cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False, linewidths=0.5
                 )
# ax.xaxis.tick_top()
ax.hlines([5], *ax.get_xlim())
ax.set_ylabel("Time")
ax.set_xlabel("phospho_site")
ax.set_facecolor('white')
add_cross_at_na(torc2_phopho_proteome_heatmap_down.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=8, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC2 hypo phosphorylated sites common from 30min & 1hr with its corresponding proteome change')
plt.subplots_adjust(bottom=0.25)
plt.savefig(path_to_save_phospho_graphs + "torc2_30m_1h_hypo_phospho_protein_heatmap_horizontal.tif", dpi=600)
plt.close()

# Comparing TORC1 all time deregulated with its corresponding proteome

torc1_all_time_deregulated_1 = torc1_all_time_deregulated.reset_index()
torc1_all_time_deregulated_1['geneName'] = torc1_all_time_deregulated_1['phospho_site'].str.split('_').str[0]

torc1_all_time_deregulated_with_proteome = torc1_all_time_deregulated_1.merge(prt_torc1_all_time_data_pivot1,
                                                                              how='left', on='geneName',
                                                                              suffixes=('_phospho', '_proteome'))

# Heatmap
torc1_all_time_deregulated_with_proteome_heatmap = torc1_all_time_deregulated_with_proteome[
    ['phospho_site', '30_min_log2ratio', '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio',
     '30_min', '1_hr', '2_hr', '3_hr', '5_hr']].set_index('phospho_site')

# plt.figure(figsize=(15, 25))
plt.figure(figsize=(15, 5))
ax = sns.heatmap(torc1_all_time_deregulated_with_proteome_heatmap.T,
                 xticklabels=True,
                 yticklabels=['0.5hr_phospho', '1hr_phospho', '2hr_phospho', '3hr_phospho', '5hr_phospho',
                              '0.5hr_prot', '1hr_prot', '2hr_prot', '3hr_prot', '5hr_prot'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar=True, cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False, linewidths=0.5
                 )
# ax.xaxis.tick_top()
ax.hlines([5], *ax.get_xlim())
ax.set_ylabel("Time")
ax.set_xlabel("phospho_site")
ax.set_facecolor('white')
add_cross_at_na(torc1_all_time_deregulated_with_proteome_heatmap.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=8, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC1 all time deregulated phosphosites with its corresponding proteome change')
plt.subplots_adjust(bottom=0.25)
plt.savefig(path_to_save_phospho_graphs + "torc1_all_time_deregulated_phosphosites_protein_heatmap_horizontal.tif",
            dpi=600)
plt.close()

# Comparing TORC1 all time deregulated with its corresponding proteome

torc2_all_time_deregulated_1 = torc2_all_time_deregulated.reset_index()
torc2_all_time_deregulated_1['geneName'] = torc2_all_time_deregulated_1['phospho_site'].str.split('_').str[0]

torc2_all_time_deregulated_with_proteome = torc2_all_time_deregulated_1.merge(prt_torc2_all_time_data_pivot1,
                                                                              how='left', on='geneName',
                                                                              suffixes=('_phospho', '_proteome'))

# Heatmap
torc2_all_time_deregulated_with_proteome_heatmap = torc2_all_time_deregulated_with_proteome[
    ['phospho_site', '30_min_log2ratio', '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio',
     '30_min', '1_hr', '2_hr', '3_hr', '5_hr']].set_index('phospho_site')

# plt.figure(figsize=(15, 25))
plt.figure(figsize=(15, 5))
ax = sns.heatmap(torc2_all_time_deregulated_with_proteome_heatmap.T,
                 xticklabels=True,
                 yticklabels=['0.5hr_phospho', '1hr_phospho', '2hr_phospho', '3hr_phospho', '5hr_phospho',
                              '0.5hr_prot', '1hr_prot', '2hr_prot', '3hr_prot', '5hr_prot'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar=True, cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False, linewidths=0.5
                 )
# ax.xaxis.tick_top()
ax.hlines([5], *ax.get_xlim())
ax.set_ylabel("Time")
ax.set_xlabel("phospho_site")
ax.set_facecolor('white')
add_cross_at_na(torc2_all_time_deregulated_with_proteome_heatmap.T)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=8, labelbottom=True, bottom=False, top=False)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.set_title('TORC2 all time deregulated phosphosites with its corresponding proteome change')
plt.subplots_adjust(bottom=0.25)
plt.savefig(path_to_save_phospho_graphs + "torc2_all_time_deregulated_phosphosites_protein_heatmap_horizontal.tif",
            dpi=600)
plt.close()
