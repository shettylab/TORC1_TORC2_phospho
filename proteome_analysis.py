import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from upsetplot import UpSet, from_contents

aux_ctr_path = r"./re_analysis/proteome_time_course/proteome_ctr_time_course/"
torc_prot_path = r"./TOR_AID_Mass_spec/Ctrl_Avo3_Kog1_Proteome_pairwise_Raw/"

# Auxin effect data (Afb2 time course significant data)
prt_auxin_effect_30min = pd.read_csv(aux_ctr_path + r"proteome_ctr_30min_signi.tsv", sep='\t')
prt_auxin_effect_1h = pd.read_csv(aux_ctr_path + r"proteome_ctr_1h_signi.tsv", sep='\t')
prt_auxin_effect_2h = pd.read_csv(aux_ctr_path + r"proteome_ctr_2h_signi.tsv", sep='\t')
prt_auxin_effect_3h = pd.read_csv(aux_ctr_path + r"proteome_ctr_3h_signi.tsv", sep='\t')
prt_auxin_effect_5h = pd.read_csv(aux_ctr_path + r"proteome_ctr_5h_signi.tsv", sep='\t')

# TORC2 data
prt_avo_zero = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C0hvsA0h\SQ_P317_MH_JH_jh002_C0hvsA0h_PROTEIN.tsv",
                           sep="\t")
prt_avo_30m = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C30minvsA30min\SQ_P317_MH_JH_jh002_C30minvsA30min_PROTEIN.tsv",
                          sep='\t')
prt_avo_1hr = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C1hvsA1h\SQ_P317_MH_JH_jh002_C1hvsA1h_PROTEIN.tsv",
                          sep="\t")
prt_avo_2hr = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C2hvsA2h\SQ_P317_MH_JH_jh002_C2hvsA2h_PROTEIN.tsv",
                          sep='\t')
prt_avo_3hr = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C3hvsA3h\SQ_P317_MH_JH_jh002_C3hvsA3h_PROTEIN.tsv",
                          sep='\t')
prt_avo_5hr = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C5hvsA5h\SQ_P317_MH_JH_jh002_C5hvsA5h_PROTEIN.tsv",
                          sep='\t')

# TORC1 data
prt_kog_zero = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C0hvsK0h\SQ_P317_MH_JH_jh002_C0hvsK0h_PROTEIN.tsv",
                           sep="\t")
prt_kog_30m = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C30minvsK30min\SQ_P317_MH_JH_jh002_C30minvsK30min_PROTEIN.tsv",
                          sep="\t")
prt_kog_1hr = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C1hvsK1h\SQ_P317_MH_JH_jh002_C1hvsK1h_PROTEIN.tsv",
                          sep="\t")
prt_kog_2hr = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C2hvsK2h\SQ_P317_MH_JH_jh002_C2hvsK2h_PROTEIN.tsv",
                          sep="\t")
prt_kog_3hr = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C3hvsK3h\SQ_P317_MH_JH_jh002_C3hvsK3h_PROTEIN.tsv",
                          sep="\t")
prt_kog_5hr = pd.read_csv(torc_prot_path + r"SQ_P317_MH_JH_jh002_C5hvsK5h\SQ_P317_MH_JH_jh002_C5hvsK5h_PROTEIN.tsv",
                          sep="\t")

pvalue_limit = 0.05
fold_change = 0


def modified_data(df, time_stamp, auxin_data, var=str()):
    data_1 = df[['ac', 'geneName', 'proteinDescription', 'nbPeptides', 'log2ratio_Condition2',
                 'pValue_Condition2', 'medianInt_Condition2', 'qValue_Condition2']]
    # Get only significant entries
    data_signi = data_1[data_1['pValue_Condition2'] <= pvalue_limit]
    data_2 = data_signi[(data_signi['log2ratio_Condition2'] > fold_change) |
                        (data_signi['log2ratio_Condition2'] < -fold_change)]
    data_2['time_stamp'] = str(time_stamp)
    # Removing auxin effect
    auxin_data = auxin_data.rename(columns={'log2ratio_Condition'+var: 'log2ratio_Condition2'})
    df_auxin_data_merge = data_2.merge(auxin_data, how='inner', on='geneName', suffixes=("_df", "_aux"))
    df_auxin_data_merge['potential_auxin_effect'] = df_auxin_data_merge['log2ratio_Condition2_df'] * df_auxin_data_merge['log2ratio_Condition2_aux']
    aux_effect_prot = df_auxin_data_merge[df_auxin_data_merge['potential_auxin_effect'] > 0]['geneName']
    aux_removed_data = data_2[~(data_2['geneName'].isin(aux_effect_prot))]

    return aux_removed_data


def norm_with_zero_hr(df, zero_hr_df):  # This check tagging effect
    df_zero_hr_df_merge = df.merge(zero_hr_df, how='inner', on='geneName', suffixes=("_df", "_zero_hr"))
    df_zero_hr_df_merge["zero_hr_effect"] = df_zero_hr_df_merge["log2ratio_Condition2_df"] * df_zero_hr_df_merge["log2ratio_Condition2_zero_hr"]
    zero_hr_effect = df_zero_hr_df_merge[df_zero_hr_df_merge["zero_hr_effect"] > 0]
    print("Number of entries in given df:", df.shape[0])
    print("number of proteins affected due to tagging effect (zero_hr_effect):", zero_hr_effect.shape[0])
    df_without_zero_hr_effect = df[~(df['geneName'].isin(zero_hr_effect['geneName']))]
    print("Number of entries after removing zero hr entries:", df_without_zero_hr_effect.shape[0])
    return df_without_zero_hr_effect, zero_hr_effect


# For zero hours data
def modified_data_for_zero(df, time_stamp):
    data_1 = df[['ac', 'geneName', 'proteinDescription', 'nbPeptides', 'log2ratio_Condition2',
                 'pValue_Condition2', 'medianInt_Condition2', 'qValue_Condition2']]
    data_1['time_stamp'] = str(time_stamp)
    # Get only significant entries
    data_signi = data_1[data_1['pValue_Condition2'] < pvalue_limit]
    # Considering differentially expressed phosphosites
    up_down_prt = data_signi[(data_signi['log2ratio_Condition2'] > fold_change) |
                             (data_signi['log2ratio_Condition2'] < -fold_change)]
    return up_down_prt

######################################################  AVO3 ANALYSIS ##################################################


torc2_result_path = r"./re_analysis/proteome_pair_wise/torc2_proteome/"
# Avo dataset
prt_avo_0h_signi_df = modified_data_for_zero(prt_avo_zero, 'zero_hr')
prt_avo_0h_signi_df.to_csv(torc2_result_path + r'prt_TORC2_0h_signi.tsv', index=False, sep="\t")
prt_avo_0h_signi = prt_avo_0h_signi_df['geneName'].to_list()

# Modifying 30 min data
prt_avo_30m_signi = modified_data(prt_avo_30m, '30_min', prt_auxin_effect_30min, var=str(2))
prt_avo_30m_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_avo_30m_norm, prt_avo_30m_tag_effect = norm_with_zero_hr(prt_avo_30m_signi, prt_avo_0h_signi_df)
prt_avo_30m_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_avo_30m_norm_down = prt_avo_30m_norm[prt_avo_30m_norm['log2ratio_Condition2'] < 0]
prt_avo_30m_norm_up = prt_avo_30m_norm[prt_avo_30m_norm['log2ratio_Condition2'] > 0]
prt_avo_30m_norm.to_csv(torc2_result_path + r'prt_TORC2_30m_norm.tsv', index=False, sep="\t")

# Modifying 1 hour data
prt_avo_1h_signi = modified_data(prt_avo_1hr, "1_hr", prt_auxin_effect_1h, var=str(3))
prt_avo_1h_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_avo_1h_norm, prt_avo_1h_tag_effect = norm_with_zero_hr(prt_avo_1h_signi, prt_avo_0h_signi_df)
prt_avo_1h_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_avo_1h_norm_down = prt_avo_1h_norm[prt_avo_1h_norm['log2ratio_Condition2'] < 0]
prt_avo_1h_norm_up = prt_avo_1h_norm[prt_avo_1h_norm['log2ratio_Condition2'] > 0]
prt_avo_1h_norm.to_csv(torc2_result_path + r'prt_TORC2_1h_norm.tsv', index=False, sep="\t")

# Modifying 2 hour data
prt_avo_2h_signi = modified_data(prt_avo_2hr, "2_hr", prt_auxin_effect_2h, var=str(4))
prt_avo_2h_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_avo_2h_norm, prt_avo_2h_tag_effect = norm_with_zero_hr(prt_avo_2h_signi, prt_avo_0h_signi_df)
prt_avo_2h_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_avo_2h_norm_down = prt_avo_2h_norm[prt_avo_2h_norm['log2ratio_Condition2'] < 0]
prt_avo_2h_norm_up = prt_avo_2h_norm[prt_avo_2h_norm['log2ratio_Condition2'] > 0]
prt_avo_2h_norm.to_csv(torc2_result_path + r'prt_TORC2_2h_norm.tsv', index=False, sep='\t')

# Modifying 3 hour data
prt_avo_3h_signi = modified_data(prt_avo_3hr, "3_hr", prt_auxin_effect_3h, var=str(5))
prt_avo_3h_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_avo_3h_norm, prt_avo_3h_tag_effect = norm_with_zero_hr(prt_avo_3h_signi, prt_avo_0h_signi_df)
prt_avo_3h_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_avo_3h_norm_down = prt_avo_3h_norm[prt_avo_3h_norm['log2ratio_Condition2'] < 0]
prt_avo_3h_norm_up = prt_avo_3h_norm[prt_avo_3h_norm['log2ratio_Condition2'] > 0]
prt_avo_3h_norm.to_csv(torc2_result_path + r'prt_TORC2_3h_norm.tsv', index=False, sep='\t')

# Modifying 5 hour data
prt_avo_5h_signi = modified_data(prt_avo_5hr, "5_hr", prt_auxin_effect_5h, var=str(6))
prt_avo_5h_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_avo_5h_norm, prt_avo_5h_tag_effect = norm_with_zero_hr(prt_avo_5h_signi, prt_avo_0h_signi_df)
prt_avo_5h_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_avo_5h_norm_down = prt_avo_5h_norm[prt_avo_5h_norm['log2ratio_Condition2'] < 0]
prt_avo_5h_norm_up = prt_avo_5h_norm[prt_avo_5h_norm['log2ratio_Condition2'] > 0]
prt_avo_5h_norm.to_csv(torc2_result_path + r'prt_TORC2_5h_norm.tsv', index=False, sep='\t')

# Building upset plot for TORC2 down proteins
prt_torc2_all_down = from_contents({"0.5hr": prt_avo_30m_norm_down['geneName'],
                                    "1hr": prt_avo_1h_norm_down['geneName'],
                                    "2hr": prt_avo_2h_norm_down['geneName'],
                                    "3hr": prt_avo_3h_norm_down['geneName'],
                                    "5hr": prt_avo_5h_norm_down['geneName']})

prt_torc2_all_down_table = prt_torc2_all_down.reset_index()
prt_torc2_all_down_table.to_excel(torc2_result_path + "prt_torc2_all_down_upset.xlsx", header=True,
                                  index=False)
prt_torc2_all_down_plot = UpSet(prt_torc2_all_down, subset_size="count", show_counts=True, sort_by='degree',
                                sort_categories_by='-input')
prt_torc2_all_down_plot.plot()
plt.suptitle("TORC2 down regulated proteins")
plt.savefig(torc2_result_path + "prt_torc2_all_down_upset_plot.tif", dpi=600)
plt.close()

# Building upset plot for TORC2 up proteins
prt_torc2_all_up = from_contents({"0.5hr": prt_avo_30m_norm_up['geneName'],
                                  "1hr": prt_avo_1h_norm_up['geneName'],
                                  "2hr": prt_avo_2h_norm_up['geneName'],
                                  "3hr": prt_avo_3h_norm_up['geneName'],
                                  "5hr": prt_avo_5h_norm_up['geneName']})

prt_torc2_all_up_table = prt_torc2_all_up.reset_index()
prt_torc2_all_up_table.to_excel(torc2_result_path + "prt_torc2_all_up_upset_ncf.xlsx", header=True, index=False)
prt_torc2_all_up_plot = UpSet(prt_torc2_all_up, subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
prt_torc2_all_up_plot.plot()
plt.suptitle("TORC2 up regulated proteins")
plt.savefig(torc2_result_path + "prt_torc2_all_up_upset_plot.tif", dpi=600)
plt.close()

# building upset plot for TORC2 differentially expressed proteins
prt_torc2_all_diff = from_contents({"0.5hr": prt_avo_30m_norm['geneName'],
                                    "1hr": prt_avo_1h_norm['geneName'],
                                    "2hr": prt_avo_2h_norm['geneName'],
                                    "3hr": prt_avo_3h_norm['geneName'],
                                    "5hr": prt_avo_5h_norm['geneName']})
prt_torc2_all_diff_table = prt_torc2_all_diff.reset_index()
prt_torc2_all_diff_table.to_excel(torc2_result_path + "prt_torc2_all_diff_upset_ncf.xlsx", header=True, index=False)
prt_torc2_all_diff_plot = UpSet(prt_torc2_all_diff, subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
prt_torc2_all_diff_plot.plot()
plt.suptitle("TORC2 differentially expressed proteins")
plt.savefig(torc2_result_path + "prt_torc2_all_diff_upset_plot.tif", dpi=600)
plt.close()

# Identifying common between 30min and 1hour from TORC2 dataset
prt_30m_1h_torc2_com = prt_avo_30m_norm.merge(prt_avo_1h_norm, how='inner', on='geneName', suffixes=('_30_min', '_1_hr'))
prt_1h_2h_torc2_com = prt_avo_1h_norm.merge(prt_avo_2h_norm, how='inner', on='geneName', suffixes=('_1_hr', '_2_hr'))

################################################ TORC1 proteome analysis ###############################################

torc1_result_path = r"./re_analysis/proteome_pair_wise/torc1_proteome/"
# kog dataset
prt_kog_0h_signi_df = modified_data_for_zero(prt_kog_zero, 'zero_hr')
prt_kog_0h_signi_df.to_csv(torc1_result_path + 'prt_TORC1_0h_signi.tsv', index=False, sep='\t')
prt_kog_0h_signi = prt_kog_0h_signi_df['geneName'].to_list()

# Modifying 30 min data
prt_kog_30m_signi = modified_data(prt_kog_30m, '30_min', prt_auxin_effect_30min, var=str(2))
prt_kog_30m_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_kog_30m_norm, prt_kog_30m_tag_effect = norm_with_zero_hr(prt_kog_30m_signi, prt_kog_0h_signi_df)
prt_kog_30m_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_kog_30m_norm_down = prt_kog_30m_norm[prt_kog_30m_norm['log2ratio_Condition2'] < 0]
prt_kog_30m_norm_up = prt_kog_30m_norm[prt_kog_30m_norm['log2ratio_Condition2'] > 0]
prt_kog_30m_norm.to_csv(torc1_result_path+'prt_TORC1_30m_norm.tsv', index=False, sep='\t')

# Modifying 1 hour data
prt_kog_1h_signi = modified_data(prt_kog_1hr, "1_hr", prt_auxin_effect_1h, var=str(3))
prt_kog_1h_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_kog_1h_norm, prt_kog_1h_tag_effect = norm_with_zero_hr(prt_kog_1h_signi, prt_kog_0h_signi_df)
prt_kog_1h_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_kog_1h_norm_down = prt_kog_1h_norm[prt_kog_1h_norm['log2ratio_Condition2'] < 0]
prt_kog_1h_norm_up = prt_kog_1h_norm[prt_kog_1h_norm['log2ratio_Condition2'] > 0]
prt_kog_1h_norm.to_csv(torc1_result_path + 'prt_TORC1_1h_norm.tsv', index=False, sep='\t')

# Modifying 2 hour data
prt_kog_2h_signi = modified_data(prt_kog_2hr, "2_hr", prt_auxin_effect_2h, var=str(4))
prt_kog_2h_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_kog_2h_norm, prt_kog_2h_tag_effect = norm_with_zero_hr(prt_kog_2h_signi, prt_kog_0h_signi_df)
prt_kog_2h_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_kog_2h_norm_down = prt_kog_2h_norm[prt_kog_2h_norm['log2ratio_Condition2'] < 0]
prt_kog_2h_norm_up = prt_kog_2h_norm[prt_kog_2h_norm['log2ratio_Condition2'] > 0]
prt_kog_2h_norm.to_csv(torc1_result_path+'prt_TORC1_2h_norm.tsv', index=False, sep='\t')

# Modifying 3 hour data
prt_kog_3h_signi = modified_data(prt_kog_3hr, "3_hr", prt_auxin_effect_3h, var=str(5))
prt_kog_3h_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_kog_3h_norm, prt_kog_3h_tag_effect = norm_with_zero_hr(prt_kog_3h_signi, prt_kog_0h_signi_df)
prt_kog_3h_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_kog_3h_norm_down = prt_kog_3h_norm[prt_kog_3h_norm['log2ratio_Condition2'] < 0]
prt_kog_3h_norm_up = prt_kog_3h_norm[prt_kog_3h_norm['log2ratio_Condition2'] > 0]
prt_kog_3h_norm.to_csv(torc1_result_path + 'prt_TORC1_3h_norm.tsv', index=False, sep='\t')

# Modifying 5 hour data
prt_kog_5h_signi = modified_data(prt_kog_5hr, "5_hr", prt_auxin_effect_5h, var=str(6))
prt_kog_5h_signi.sort_values('pValue_Condition2', ignore_index=True, inplace=True)
# Removing zero hours up and down sites
prt_kog_5h_norm, prt_kog_5h_tag_effect = norm_with_zero_hr(prt_kog_5h_signi, prt_kog_0h_signi_df)
prt_kog_5h_norm.sort_values('log2ratio_Condition2', inplace=True)
prt_kog_5h_norm_down = prt_kog_5h_norm[prt_kog_5h_norm['log2ratio_Condition2'] < 0]
prt_kog_5h_norm_up = prt_kog_5h_norm[prt_kog_5h_norm['log2ratio_Condition2'] > 0]
prt_kog_5h_norm.to_csv(torc1_result_path + 'prt_TORC1_5h_norm.tsv', index=False, sep='\t')

# Building upset plot for TORC1 down proteins
prt_torc1_all_down = from_contents({"0.5hr": prt_kog_30m_norm_down['geneName'],
                                    "1hr": prt_kog_1h_norm_down['geneName'],
                                    "2hr": prt_kog_2h_norm_down['geneName'],
                                    "3hr": prt_kog_3h_norm_down['geneName'],
                                    "5hr": prt_kog_5h_norm_down['geneName']})

prt_torc1_all_down_table = prt_torc1_all_down.reset_index()
prt_torc1_all_down_table.to_excel(torc1_result_path + "prt_TORC1_all_down_upset_ncf.xlsx", header=True, index=False)
prt_torc1_all_down_plot = UpSet(prt_torc1_all_down, subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
prt_torc1_all_down_plot.plot()
plt.suptitle("TORC1 down regulated proteins")
plt.savefig(torc1_result_path + "prt_TORC1_all_down_upset_plot.tif", dpi=600)
plt.close()

# Building upset plot for TORC1 up proteins
prt_torc1_all_up = from_contents({"0.5hr": prt_kog_30m_norm_up['geneName'],
                                  "1hr": prt_kog_1h_norm_up['geneName'],
                                  "2hr": prt_kog_2h_norm_up['geneName'],
                                  "3hr": prt_kog_3h_norm_up['geneName'],
                                  "5hr": prt_kog_5h_norm_up['geneName']})

prt_torc1_all_up_table = prt_torc1_all_up.reset_index()
prt_torc1_all_up_table.to_excel(torc1_result_path + "prt_torc1_all_up_upset_ncf.xlsx", header=True, index=False)
prt_torc1_all_up_plot = UpSet(prt_torc1_all_up, subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
prt_torc1_all_up_plot.plot()
plt.suptitle("TORC1 up regulated proteins")
plt.savefig(torc1_result_path + "prt_torc1_all_up_upset_plot_ncf.tif", dpi=600)
plt.close()

# building upset plot for TORC1 differentially expressed proteins
prt_torc1_all_diff = from_contents({"0.5hr_deregulated": prt_kog_30m_norm['geneName'],
                                    "1hr_deregulated": prt_kog_1h_norm['geneName'],
                                    "2hr_deregulated": prt_kog_2h_norm['geneName'],
                                    "3hr_deregulated": prt_kog_3h_norm['geneName'],
                                    "5hr_deregulated": prt_kog_5h_norm['geneName']})
prt_torc1_all_diff_table = prt_torc1_all_diff.reset_index()
prt_torc1_all_diff_table.to_excel(torc1_result_path + r"prt_TORC1_all_diff_upset_ncf.xlsx", header=True, index=False)
prt_torc1_all_diff_plot = UpSet(prt_torc1_all_diff, subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
prt_torc1_all_diff_plot.plot()
plt.suptitle("TORC1 differentially expressed proteins")
plt.savefig(torc1_result_path + "prt_TORC1_all_diff_upset_plot.tif", dpi=600)
plt.close()
#########################

# Identifying common between 30min and 1hour from TORC1 dataset
prt_30m_1h_torc1_com = prt_kog_30m_norm.merge(prt_kog_1h_norm, how='inner', on='geneName',
                                              suffixes=('_30_min', '_1_hr'))
prt_1h_2h_torc1_com = prt_kog_1h_norm.merge(prt_kog_2h_norm, how='inner', on='geneName',
                                              suffixes=('_1_hr', '_2_hr'))
