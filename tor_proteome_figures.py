import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt
from upsetplot import generate_counts, plot, UpSet, from_contents

############################### getting total proteome dataset

# Time course proteome dataset
path_to_tm_proteome = r"./re_analysis/proteome_time_course/"

# Control(Afb2) time course proteome
prt_afb2_30m_tm = pd.read_csv(path_to_tm_proteome+"/proteome_ctr_time_course/proteome_ctr_30min_signi.tsv", sep='\t')
prt_afb2_1h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_ctr_time_course/proteome_ctr_1h_signi.tsv", sep='\t')
prt_afb2_2h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_ctr_time_course/proteome_ctr_2h_signi.tsv", sep='\t')
prt_afb2_3h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_ctr_time_course/proteome_ctr_3h_signi.tsv", sep='\t')
prt_afb2_5h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_ctr_time_course/proteome_ctr_5h_signi.tsv", sep='\t')

# TORC1 time course proteome
prt_torc1_30m_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc1_time_course/proteome_torc1_30min_signi.tsv", sep='\t')
prt_torc1_1h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc1_time_course/proteome_torc1_1hr_signi.tsv", sep='\t')
prt_torc1_2h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc1_time_course/proteome_torc1_2hr_signi.tsv", sep='\t')
prt_torc1_3h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc1_time_course/proteome_torc1_3hr_signi.tsv", sep='\t')
prt_torc1_5h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc1_time_course/proteome_torc1_5hr_signi.tsv", sep='\t')

# TORC2 time course proteome
prt_torc2_30m_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc2_time_course/proteome_torc2_30min_signi.tsv", sep='\t')
prt_torc2_1h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc2_time_course/proteome_torc2_1hr_signi.tsv", sep='\t')
prt_torc2_2h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc2_time_course/proteome_torc2_2hr_signi.tsv", sep='\t')
prt_torc2_3h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc2_time_course/proteome_torc2_3hr_signi.tsv", sep='\t')
prt_torc2_5h_tm = pd.read_csv(path_to_tm_proteome+"/proteome_torc2_time_course/proteome_torc2_5hr_signi.tsv", sep='\t')

# Pairwise total proteome
path_to_pr_proteome = r"./re_analysis/proteome_pair_wise/"

# TORC1 pairwise data
prt_torc1_0h_pr = pd.read_csv(path_to_pr_proteome+"/torc1_proteome/prt_TORC1_0h_signi.tsv", sep='\t')
prt_torc1_30m_pr = pd.read_csv(path_to_pr_proteome+"/torc1_proteome/prt_TORC1_30m_norm.tsv", sep='\t')
prt_torc1_1h_pr = pd.read_csv(path_to_pr_proteome+"/torc1_proteome/prt_TORC1_1h_norm.tsv", sep='\t')
prt_torc1_2h_pr = pd.read_csv(path_to_pr_proteome+"/torc1_proteome/prt_TORC1_2h_norm.tsv", sep='\t')
prt_torc1_3h_pr = pd.read_csv(path_to_pr_proteome+"/torc1_proteome/prt_TORC1_3h_norm.tsv", sep='\t')
prt_torc1_5h_pr = pd.read_csv(path_to_pr_proteome+"/torc1_proteome/prt_TORC1_5h_norm.tsv", sep='\t')

# TORC2 pairwise data
prt_torc2_0h_pr = pd.read_csv(path_to_pr_proteome+"/torc2_proteome/prt_TORC2_0h_signi.tsv", sep='\t')
prt_torc2_30m_pr = pd.read_csv(path_to_pr_proteome+"/torc2_proteome/prt_TORC2_30m_norm.tsv", sep='\t')
prt_torc2_1h_pr = pd.read_csv(path_to_pr_proteome+"/torc2_proteome/prt_TORC2_1h_norm.tsv", sep='\t')
prt_torc2_2h_pr = pd.read_csv(path_to_pr_proteome+"/torc2_proteome/prt_TORC2_2h_norm.tsv", sep='\t')
prt_torc2_3h_pr = pd.read_csv(path_to_pr_proteome+"/torc2_proteome/prt_TORC2_3h_norm.tsv", sep='\t')
prt_torc2_5h_pr = pd.read_csv(path_to_pr_proteome+"/torc2_proteome/prt_TORC2_5h_norm.tsv", sep='\t')

path_to_save_proteome_result = r"./re_analysis/proteome_graphs_results/"
# To check how many are common between time course and pairwise data at protein level
# For TORC1
prt_torc1_tm_pr_comp = [
    [len(set(prt_torc1_30m_tm['geneName'])), len(set(prt_torc1_30m_pr['geneName'])), len(set(prt_torc1_30m_tm['geneName']) & set(prt_torc1_30m_pr['geneName']))],
    [len(set(prt_torc1_1h_tm['geneName'])), len(set(prt_torc1_1h_pr['geneName'])), len(set(prt_torc1_1h_tm['geneName']) & set(prt_torc1_1h_pr['geneName']))],
    [len(set(prt_torc1_2h_tm['geneName'])), len(set(prt_torc1_2h_pr['geneName'])), len(set(prt_torc1_2h_tm['geneName']) & set(prt_torc1_2h_pr['geneName']))],
    [len(set(prt_torc1_3h_tm['geneName'])), len(set(prt_torc1_3h_pr['geneName'])), len(set(prt_torc1_3h_tm['geneName']) & set(prt_torc1_3h_pr['geneName']))],
    [len(set(prt_torc1_5h_tm['geneName'])), len(set(prt_torc1_5h_pr['geneName'])), len(set(prt_torc1_5h_tm['geneName']) & set(prt_torc1_5h_pr['geneName']))],
]
prt_torc1_tm_pr_comp_df = pd.DataFrame(prt_torc1_tm_pr_comp, columns=["Time Course(TC)", "Pair Wise(PW)", "Common between TC & PW"])
prt_torc1_tm_pr_comp_df.to_excel(path_to_save_proteome_result+'torc1_tm_pr_comparison.xlsx', index=False, header=True)

# Comparing 30min data with 1 hr and 1 hr with 2 hr

# TORC1 30min and 1 hr (only pairwise data is considered)
prt_torc1_30m_1hr_com_pr = prt_torc1_30m_pr.merge(prt_torc1_1h_pr, how='inner', on='geneName', suffixes=('_30_min', '_1_hr'))
prt_torc1_30m_1hr_com_pr.to_csv(path_to_save_proteome_result + 'proteome_torc1_30min_1hr_common.tsv', index=False, sep='\t')

prt_torc1_only_30m_compare_to_1hr_pr = set(prt_torc1_30m_pr['geneName']) - set(prt_torc1_1h_pr['geneName'])
prt_torc1_only_1hr_compare_to_30m_pr = set(prt_torc1_1h_pr['geneName']) - set(prt_torc1_30m_pr['geneName'])

# Plotting venn diagram for TORC1 proteome 30 min and 1 hr
a = venn2(subsets=(len(prt_torc1_only_30m_compare_to_1hr_pr), len(prt_torc1_only_1hr_compare_to_30m_pr), prt_torc1_30m_1hr_com_pr.shape[0]),
      set_labels=('0.5 hour', '1 hour'), set_colors=('orange', 'blue'), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(25)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in a.set_labels:
    text.set_fontsize(20)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("Deregulated proteins upon TORC1 depletion")
plt.savefig(path_to_save_proteome_result+'TORC1_proteome_30m_1h_venn.tif', dpi=600)
plt.close()

# TORC1 1 hr and 2 hr
prt_torc1_1hr_2hr_com_pr = prt_torc1_1h_pr.merge(prt_torc1_2h_pr, how='inner', on='geneName', suffixes=('_1_hr', '_2_hr'))
prt_torc1_1hr_2hr_com_pr.to_csv(path_to_save_proteome_result + 'proteome_torc1_1hr_2hr_common.tsv', index=False, sep='\t')

prt_torc1_only_1hr_compare_to_2hr_pr = set(prt_torc1_1h_pr['geneName']) - set(prt_torc1_2h_pr['geneName'])
prt_torc1_only_2hr_compare_to_1hr_pr = set(prt_torc1_2h_pr['geneName']) - set(prt_torc1_1h_pr['geneName'])

# Plotting venn diagram for TORC1 proteome 1 hr and 2 hr
a = venn2(subsets=(len(prt_torc1_only_1hr_compare_to_2hr_pr), len(prt_torc1_only_2hr_compare_to_1hr_pr), prt_torc1_1hr_2hr_com_pr.shape[0]),
      set_labels=('1 hour', '2 hour'), set_colors=('orange', 'blue'), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(25)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in a.set_labels:
    text.set_fontsize(20)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("Deregulated proteins upon TORC1 depletion")
plt.savefig(path_to_save_proteome_result+'TORC1_proteome_1h_2h_venn.tif', dpi=600)
plt.close()

# Plotting venn diagram for TORC1 30 min, 1hr and 2 hr
a = venn3(subsets=(set(prt_torc1_30m_pr['geneName']), set(prt_torc1_1h_pr['geneName']), set(prt_torc1_2h_pr['geneName'])),
          set_labels=('30 min', '1 hour', '2 hour'))
plt.title("Deregulated proteins upon TORC1 depletion")
plt.savefig(path_to_save_proteome_result+'TORC1_proteome_30m_1h_2h_venn.png', dpi=600)
plt.close()

torc1_all_time_data = pd.concat([
    prt_torc1_30m_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc1_1h_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc1_2h_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc1_3h_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc1_5h_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
])


def df_for_heatmap_prt(data, compare_data, value_name):
    df1 = pd.pivot_table(data=data, index='geneName', columns='time_stamp', values=value_name)
    # df1.fillna(0, inplace=True)
    df2 = df1.loc[:, ['30_min', '1_hr', '2_hr', '3_hr', '5_hr']]
    df3 = df2[df2.index.isin(compare_data['geneName'].to_list())]
    df4 = df3.sort_values('30_min', ascending=True)
    return df4


def add_cross_at_na(data):
    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            if pd.isna(data.iloc[y, x]):
                ax.plot([x, x + 1], [y, y + 1], color='gray', lw=0.5)
                ax.plot([x, x + 1], [y + 1, y], color='gray', lw=0.5)


# comparing TORC1 30min-1hr common with rest time points
prt_torc1_30m_1hr_com_pr_fold_change = df_for_heatmap_prt(torc1_all_time_data, prt_torc1_30m_1hr_com_pr, 'log2ratio_Condition2').add_suffix("_log2ratio")
prt_torc1_30m_1hr_com_pr_p_value = df_for_heatmap_prt(torc1_all_time_data, prt_torc1_30m_1hr_com_pr, 'pValue_Condition2').add_suffix("_pValue")
prt_torc1_30m_1hr_com_pr_all_data = pd.concat([prt_torc1_30m_1hr_com_pr_fold_change, prt_torc1_30m_1hr_com_pr_p_value],  axis=1)

prt_torc1_30m_1hr_com_pr_all_data_down = prt_torc1_30m_1hr_com_pr_all_data[prt_torc1_30m_1hr_com_pr_all_data['30_min_log2ratio'] < 0].reset_index()
prt_torc1_30m_1hr_com_pr_all_data_down.to_csv(path_to_save_proteome_result+'prt_torc1_30m_1hr_com_pr_all_data_down.tsv', index=False, sep='\t')
prt_torc1_30m_1hr_com_pr_all_data_up = prt_torc1_30m_1hr_com_pr_all_data[prt_torc1_30m_1hr_com_pr_all_data['30_min_log2ratio'] > 0].reset_index()
prt_torc1_30m_1hr_com_pr_all_data_up.to_csv(path_to_save_proteome_result+'prt_torc1_30m_1hr_com_pr_all_data_up.tsv', index=False, sep='\t')

# heatmap for torc1 down and up
prt_torc1_30m_1hr_com_pr_all_data_down_heatmap = prt_torc1_30m_1hr_com_pr_all_data_down[['geneName', '30_min_log2ratio',
                                                 '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio']].set_index('geneName')
plt.figure(figsize=(7, 10))
ax = sns.heatmap(prt_torc1_30m_1hr_com_pr_all_data_down_heatmap,
                 yticklabels=True, xticklabels=['0.5hr', '1hr', '2hr', '3hr', '5hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5 )
# ax.xaxis.tick_top()
ax.set_ylabel("Proteins")
ax.set_xlabel("Time")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.set_facecolor('white')
add_cross_at_na(prt_torc1_30m_1hr_com_pr_all_data_down_heatmap)
# ax.tick_params(axis='x', which='major', labelsize=9, labelbottom = True)
ax.set_title('TORC1 down regulated proteins common from 30min & 1hr')
plt.savefig(path_to_save_proteome_result+"prt_torc1_30m_1hr_com_pr_all_data_down_heatmap.tif", dpi=600)
plt.close()

prt_torc1_30m_1hr_com_pr_all_data_up_heatmap = prt_torc1_30m_1hr_com_pr_all_data_up[['geneName', '30_min_log2ratio',
                                                 '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio']].set_index('geneName')
plt.figure(figsize=(10, 17))
ax = sns.heatmap(prt_torc1_30m_1hr_com_pr_all_data_up_heatmap,
                 yticklabels=True, xticklabels=['0.5hr', '1hr', '2hr', '3hr', '5hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5)
# ax.xaxis.tick_top()
ax.set_ylabel("Proteins")
ax.set_xlabel("Time")
ax.set_facecolor('white')
add_cross_at_na(prt_torc1_30m_1hr_com_pr_all_data_up_heatmap)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='y', which='major', labelsize=7)
ax.tick_params(axis='x', which='major', labelsize=7)
ax.set_title('TORC1 up regulated proteins common from 30min & 1hr')
plt.savefig(path_to_save_proteome_result+"prt_torc1_30m_1hr_com_pr_all_data_up_heatmap.tif", dpi=600)
plt.close()

# For TORC2
prt_torc2_tm_pr_comp = [
    [len(set(prt_torc2_30m_tm['geneName'])), len(set(prt_torc2_30m_pr['geneName'])), len(set(prt_torc2_30m_tm['geneName']) & set(prt_torc2_30m_pr['geneName']))],
    [len(set(prt_torc2_1h_tm['geneName'])), len(set(prt_torc2_1h_pr['geneName'])), len(set(prt_torc2_1h_tm['geneName']) & set(prt_torc2_1h_pr['geneName']))],
    [len(set(prt_torc2_2h_tm['geneName'])), len(set(prt_torc2_2h_pr['geneName'])), len(set(prt_torc2_2h_tm['geneName']) & set(prt_torc2_2h_pr['geneName']))],
    [len(set(prt_torc2_3h_tm['geneName'])), len(set(prt_torc2_3h_pr['geneName'])), len(set(prt_torc2_3h_tm['geneName']) & set(prt_torc2_3h_pr['geneName']))],
    [len(set(prt_torc2_5h_tm['geneName'])), len(set(prt_torc2_5h_pr['geneName'])), len(set(prt_torc2_5h_tm['geneName']) & set(prt_torc2_5h_pr['geneName']))],
]
prt_torc2_tm_pr_comp_df = pd.DataFrame(prt_torc2_tm_pr_comp, columns=["Time Course(TC)", "Pair Wise(PW)", "Common between TC & PW"])
prt_torc2_tm_pr_comp_df.to_excel(path_to_save_proteome_result+'torc2_tm_pr_comparison.xlsx', index=False, header=True)

# Comparing 30min data with 1 hr and 1 hr with 2 hr

# TORC2 30min and 1 hr (only pairwise data is considered)
prt_torc2_30m_1hr_com_pr = prt_torc2_30m_pr.merge(prt_torc2_1h_pr, how='inner', on='geneName', suffixes=('_30_min', '_1_hr'))
prt_torc2_30m_1hr_com_pr.to_csv(path_to_save_proteome_result + 'proteome_torc2_30min_1hr_common.tsv', index=False, sep='\t')

prt_torc2_only_30m_compare_to_1hr_pr = set(prt_torc2_30m_pr['geneName']) - set(prt_torc2_1h_pr['geneName'])
prt_torc2_only_1hr_compare_to_30m_pr = set(prt_torc2_1h_pr['geneName']) - set(prt_torc2_30m_pr['geneName'])

# Plotting venn diagram for TORC2 proteome 30 min and 1 hr
b = venn2(subsets=(len(prt_torc2_only_30m_compare_to_1hr_pr), len(prt_torc2_only_1hr_compare_to_30m_pr), prt_torc2_30m_1hr_com_pr.shape[0]),
      set_labels=('0.5 hour', '1 hour'), set_colors=('orange', 'blue'), alpha=0.6)
for text in b.subset_labels:
    text.set_fontsize(25)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in b.set_labels:
    text.set_fontsize(20)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("Deregulated proteins upon TORC2 depletion")
plt.savefig(path_to_save_proteome_result+'TORC2_proteome_30m_1h_venn.tif', dpi=600)
plt.close()

# TORC2 1 hr and 2 hr
prt_torc2_1hr_2hr_com_pr = prt_torc2_1h_pr.merge(prt_torc2_2h_pr, how='inner', on='geneName', suffixes=('_1_hr', '_2_hr'))
prt_torc2_1hr_2hr_com_pr.to_csv(path_to_save_proteome_result + 'proteome_torc2_1hr_2hr_common.tsv', index=False, sep='\t')

prt_torc2_only_1hr_compare_to_2hr_pr = set(prt_torc2_1h_pr['geneName']) - set(prt_torc2_2h_pr['geneName'])
prt_torc2_only_2hr_compare_to_1hr_pr = set(prt_torc2_2h_pr['geneName']) - set(prt_torc2_1h_pr['geneName'])

# Plotting venn diagram for TORC2 proteome 30 min and 1 hr
b = venn2(subsets=(len(prt_torc2_only_1hr_compare_to_2hr_pr), len(prt_torc2_only_2hr_compare_to_1hr_pr), prt_torc2_1hr_2hr_com_pr.shape[0]),
      set_labels=('1 hour', '2 hour'), set_colors=('orange', 'blue'), alpha=0.6)
for text in b.subset_labels:
    text.set_fontsize(25)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in b.set_labels:
    text.set_fontsize(20)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("Deregulated proteins upon TORC2 depletion")
plt.savefig(path_to_save_proteome_result+'TORC2_proteome_1h_2h_venn.tif', dpi=600)
plt.close()

# Plotting venn diagram for TORC1 30 min, 1hr and 2 hr
a = venn3(subsets=(set(prt_torc2_30m_pr['geneName']), set(prt_torc2_1h_pr['geneName']), set(prt_torc2_2h_pr['geneName'])),
          set_labels=('30 min', '1 hour', '2 hour'))
plt.title("Deregulated proteins upon TORC2 depletion")
plt.savefig(path_to_save_proteome_result+'TORC2_proteome_30m_1h_2h_venn.tif', dpi=600)
plt.close()

torc2_all_time_data = pd.concat([
    prt_torc2_30m_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc2_1h_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc2_2h_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc2_3h_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
    prt_torc2_5h_pr[['time_stamp', 'log2ratio_Condition2', 'geneName', 'pValue_Condition2']],
])


# comparing TORC2 30min-1hr common with rest time points
prt_torc2_30m_1hr_com_pr_fold_change = df_for_heatmap_prt(torc2_all_time_data, prt_torc2_30m_1hr_com_pr, 'log2ratio_Condition2').add_suffix("_log2ratio")
prt_torc2_30m_1hr_com_pr_p_value = df_for_heatmap_prt(torc2_all_time_data, prt_torc2_30m_1hr_com_pr, 'pValue_Condition2').add_suffix("_pValue")
prt_torc2_30m_1hr_com_pr_all_data = pd.concat([prt_torc2_30m_1hr_com_pr_fold_change, prt_torc2_30m_1hr_com_pr_p_value],  axis=1)

prt_torc2_30m_1hr_com_pr_all_data_down = prt_torc2_30m_1hr_com_pr_all_data[prt_torc2_30m_1hr_com_pr_all_data['30_min_log2ratio'] < 0].reset_index()
prt_torc2_30m_1hr_com_pr_all_data_down.to_csv(path_to_save_proteome_result+'prt_torc2_30m_1hr_com_pr_all_data_down.tsv', index=False, sep='\t')
prt_torc2_30m_1hr_com_pr_all_data_up = prt_torc2_30m_1hr_com_pr_all_data[prt_torc2_30m_1hr_com_pr_all_data['30_min_log2ratio'] > 0].reset_index()
prt_torc2_30m_1hr_com_pr_all_data_up.to_csv(path_to_save_proteome_result+'prt_torc2_30m_1hr_com_pr_all_data_up.tsv', index=False, sep='\t')

# heatmap for torc1 down and up
prt_torc2_30m_1hr_com_pr_all_data_down_heatmap = prt_torc2_30m_1hr_com_pr_all_data_down[['geneName', '30_min_log2ratio',
                                                 '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio']].set_index('geneName')
plt.figure(figsize=(7, 10))
ax = sns.heatmap(prt_torc2_30m_1hr_com_pr_all_data_down_heatmap,
                 yticklabels=True, xticklabels=['30min', '1hr', '2hr', '3hr', '5hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5 )
# ax.xaxis.tick_top()
ax.set_ylabel("Proteins")
ax.set_xlabel("Time")
ax.set_facecolor('white')
add_cross_at_na(prt_torc2_30m_1hr_com_pr_all_data_down_heatmap)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
# ax.tick_params(axis='x', which='major', labelsize=9, labelbottom = True)
ax.set_title('TORC2 down regulated proteins common from 30min & 1hr')
plt.savefig(path_to_save_proteome_result+"prt_torc2_30m_1hr_com_pr_all_data_down_heatmap.tif", dpi=600)
plt.close()

prt_torc2_30m_1hr_com_pr_all_data_up_heatmap = prt_torc2_30m_1hr_com_pr_all_data_up[['geneName', '30_min_log2ratio',
                                                 '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio']].set_index('geneName')
plt.figure(figsize=(10, 17))
ax = sns.heatmap(prt_torc2_30m_1hr_com_pr_all_data_up_heatmap,
                 yticklabels=True, xticklabels=['30min', '1hr', '2hr', '3hr', '5hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5)
# ax.xaxis.tick_top()
ax.set_ylabel("Proteins")
ax.set_xlabel("Time")
ax.set_facecolor('white')
add_cross_at_na(prt_torc2_30m_1hr_com_pr_all_data_up_heatmap)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.tick_params(axis='x', which='major', labelsize=8)
ax.set_title('TORC2 up regulated proteins common from 30min & 1hr')
plt.savefig(path_to_save_proteome_result+"prt_torc2_30m_1hr_com_pr_all_data_up_heatmap.tif", dpi=600)
plt.close()


# comparing TORC2 1hr-2hr common with rest time points
prt_torc2_1hr_2hr_com_pr_fold_change = df_for_heatmap_prt(torc2_all_time_data, prt_torc2_1hr_2hr_com_pr, 'log2ratio_Condition2').add_suffix("_log2ratio")
prt_torc2_1hr_2hr_com_pr_p_value = df_for_heatmap_prt(torc2_all_time_data, prt_torc2_1hr_2hr_com_pr, 'pValue_Condition2').add_suffix("_pValue")
prt_torc2_1hr_2hr_com_pr_all_data = pd.concat([prt_torc2_1hr_2hr_com_pr_fold_change, prt_torc2_1hr_2hr_com_pr_p_value],  axis=1)

prt_torc2_1hr_2hr_com_pr_all_data_down = prt_torc2_1hr_2hr_com_pr_all_data[prt_torc2_1hr_2hr_com_pr_all_data['1_hr_log2ratio'] < 0].reset_index()
prt_torc2_1hr_2hr_com_pr_all_data_down.to_csv(path_to_save_proteome_result+'prt_torc2_1hr_2hr_com_pr_all_data_down.tsv', index=False, sep='\t')
prt_torc2_1hr_2hr_com_pr_all_data_up = prt_torc2_1hr_2hr_com_pr_all_data[prt_torc2_1hr_2hr_com_pr_all_data['1_hr_log2ratio'] > 0].reset_index()
prt_torc2_1hr_2hr_com_pr_all_data_up.to_csv(path_to_save_proteome_result+'prt_torc2_1hr_2hr_com_pr_all_data_up.tsv', index=False, sep='\t')

# heatmap for torc1 down and up
prt_torc2_1hr_2hr_com_pr_all_data_down_heatmap = prt_torc2_1hr_2hr_com_pr_all_data_down[['geneName', '30_min_log2ratio',
                                                 '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio']].set_index('geneName')
plt.figure(figsize=(7, 10))
ax = sns.heatmap(prt_torc2_1hr_2hr_com_pr_all_data_down_heatmap,
                 yticklabels=True, xticklabels=['0.5hr', '1hr', '2hr', '3hr', '5hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5 )
# ax.xaxis.tick_top()
ax.set_ylabel("Proteins")
ax.set_xlabel("Time")
ax.set_facecolor('white')
add_cross_at_na(prt_torc2_1hr_2hr_com_pr_all_data_down_heatmap)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
# ax.tick_params(axis='x', which='major', labelsize=9, labelbottom = True)
ax.set_title('TORC2 down regulated proteins common from 1hr & 2hr')
plt.savefig(path_to_save_proteome_result+"prt_torc2_1hr_2hr_com_pr_all_data_down_heatmap.tif", dpi=600)
plt.close()

prt_torc2_1hr_2hr_com_pr_all_data_up_heatmap = prt_torc2_1hr_2hr_com_pr_all_data_up[['geneName', '30_min_log2ratio',
                                                 '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio']].set_index('geneName')
plt.figure(figsize=(10, 17))
ax = sns.heatmap(prt_torc2_1hr_2hr_com_pr_all_data_up_heatmap,
                 yticklabels=True, xticklabels=['0.5hr', '1hr', '2hr', '3hr', '5hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'}, center=0, vmin=-4, vmax=4, square=True, linewidths=0.5)
# ax.xaxis.tick_top()
ax.set_ylabel("Proteins")
ax.set_xlabel("Time")
ax.set_facecolor('white')
add_cross_at_na(prt_torc2_1hr_2hr_com_pr_all_data_up_heatmap)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='y', which='major', labelsize=8)
ax.tick_params(axis='x', which='major', labelsize=8)
ax.set_title('TORC2 up regulated proteins common from 1hr & 2hr')
plt.savefig(path_to_save_proteome_result+"prt_torc2_1hr_2hr_com_pr_all_data_up_heatmap.tif", dpi=600)
plt.close()

#################### Common between TORC1 and TORC2 ########################
# Between 30 min and 1 hour
prt_torc1_torc2_30m_1h_com = prt_torc1_30m_1hr_com_pr.merge(prt_torc2_30m_1hr_com_pr, on="geneName", how="inner", suffixes=('_TORC1', '_TORC2'))
prt_torc1_torc2_30m_1h_com.to_csv(path_to_save_proteome_result+"TORC1_TORC2_common_prot_30m_1hr.tsv", index=False, header=True, sep='\t')

prt_torc1_only_compare_torc2_30m_1hr_pr = set(prt_torc1_30m_1hr_com_pr['geneName']) - set(prt_torc2_30m_1hr_com_pr['geneName'])
prt_torc2_only_compare_torc1_30m_1hr_pr = set(prt_torc2_30m_1hr_com_pr['geneName']) - set(prt_torc1_30m_1hr_com_pr['geneName'])

# Plotting venn diagram
c = venn2(subsets=(len(prt_torc1_only_compare_torc2_30m_1hr_pr), len(prt_torc2_only_compare_torc1_30m_1hr_pr), prt_torc1_torc2_30m_1h_com.shape[0]),
      set_labels=('TORC1', 'TORC2'), set_colors=('orange', 'blue'), alpha=0.6)
for text in c.subset_labels:
    text.set_fontsize(25)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in c.set_labels:
    text.set_fontsize(20)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("TORC1 and TORC2 common (30m-1hr)")
plt.savefig(path_to_save_proteome_result+'TORC1_TORC2_common_30m-1hr_venn.tif', dpi=600)
plt.close()

# Between 1 hr and 2 hr

prt_torc1_torc2_1hr_2hr_com = prt_torc1_1hr_2hr_com_pr.merge(prt_torc2_1hr_2hr_com_pr, on="geneName", how="inner", suffixes=('_TORC1', '_TORC2'))
prt_torc1_torc2_1hr_2hr_com.to_csv(path_to_save_proteome_result+"TORC1_TORC2_common_prot_1hr_2hr.tsv", index=False, header=True, sep='\t')

prt_torc1_only_compare_torc2_1hr_2hr_pr = set(prt_torc1_1hr_2hr_com_pr['geneName']) - set(prt_torc2_1hr_2hr_com_pr['geneName'])
prt_torc2_only_compare_torc1_1hr_2hr_pr = set(prt_torc2_1hr_2hr_com_pr['geneName']) - set(prt_torc1_1hr_2hr_com_pr['geneName'])

# Plotting venn diagram
d = venn2(subsets=(len(prt_torc1_only_compare_torc2_1hr_2hr_pr), len(prt_torc2_only_compare_torc1_1hr_2hr_pr), prt_torc1_torc2_1hr_2hr_com.shape[0]),
      set_labels=('TORC1', 'TORC2'), set_colors=('orange', 'blue'), alpha=0.6)
for text in d.subset_labels:
    text.set_fontsize(25)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in d.set_labels:
    text.set_fontsize(20)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("TORC1 and TORC2 common (1hr-2hr)")
plt.savefig(path_to_save_proteome_result+'TORC1_TORC2_common_1hr-2hr_venn.tif', dpi=600)
plt.close()
