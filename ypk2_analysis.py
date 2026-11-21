import os
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt
import re

os.chdir(r"D:\PhD\Python_script")
# Loading dataset

# TORC2 AUX and NO-AUX compared samples
path_of_aux_no_aux = r"./TOR_AID_Mass_spec/Ypk2_Avo3AID_phosphoproteome/Raw_files/Noaux_Aux_pairwise/"

torc2_aux_noaux_15min = pd.read_csv(
    path_of_aux_no_aux + "SQ_P317_MH_SS_ss040_1-4VS5-8/SQ_P317_MH_SS_ss040_1-4VS5-8_PEPTIDE.tsv", sep='\t')
torc2_aux_noaux_30min = pd.read_csv(
    path_of_aux_no_aux + "SQ_P317_MH_SS_ss040_17-20VS21-24/SQ_P317_MH_SS_ss040_17-20VS21-24_PEPTIDE.tsv", sep='\t')
torc2_aux_noaux_60min = pd.read_csv(
    path_of_aux_no_aux + "SQ_P317_MH_SS_ss040_33-36VS37-40/SQ_P317_MH_SS_ss040_33-36VS37-40_PEPTIDE.tsv", sep='\t')

# YPK2-CA AUX and NO-AUX compared samples

ypk2_aux_noaux_15min = pd.read_csv(
    path_of_aux_no_aux + "SQ_P317_MH_SS_ss040_9-12VS13-16/SQ_P317_MH_SS_ss040_9-12VS13-16_PEPTIDE.tsv", sep='\t')
ypk2_aux_noaux_30min = pd.read_csv(
    path_of_aux_no_aux + "SQ_P317_MH_SS_ss040_25-28VS29-32/SQ_P317_MH_SS_ss040_15-28VS29-32_PEPTIDE.tsv", sep='\t')
ypk2_aux_noaux_60min = pd.read_csv(
    path_of_aux_no_aux + "SQ_P317_MH_SS_ss040_41-44VS45-48/SQ_P317_MH_SS_ss040_41-44VS45-48_PEPTIDE.tsv", sep='\t')

# TORC2 AUX and YPK2-CA compared samples
path_of_torc2_ypk2ca_aux = r"./TOR_AID_Mass_spec/Ypk2_Avo3AID_phosphoproteome/Raw_files/Avo3AID_VS_YPK2CA/"

torc2_ypk2_aux_15min = pd.read_csv(
    path_of_torc2_ypk2ca_aux + "SQ_P317_MH_SS_ss040_5-8VS13-16/SQ_P317_MH_SS_ss040_5-8VS13-16_PEPTIDE.tsv", sep='\t')
torc2_ypk2_aux_30min = pd.read_csv(
    path_of_torc2_ypk2ca_aux + "SQ_P317_MH_SS_ss040_21-24VS29-32/SQ_P317_MH_SS_ss040_21-24VS29-32_PEPTIDE.tsv",
    sep='\t')
torc2_ypk2_aux_60min = pd.read_csv(
    path_of_torc2_ypk2ca_aux + "SQ_P317_MH_SS_ss040_37-40VS45-48/SQ_P317_MH_SS_ss040_37-40VS45-48_PEPTIDE.tsv",
    sep='\t')

fold_change = 0


def signi_data(df, time_stamp, ctr=None, df_var=str(2), ctr_var=None):  # To check auxin effect
    data_1 = df[['ac', 'geneName', 'modifCoord', 'ptm', 'motifX', 'log2ratio_Condition' + df_var, 'pValue_Condition' + df_var,
                 'medianInt_Condition' + df_var]]
    # sorting data based on median intensity of condition2
    data_1.sort_values('medianInt_Condition' + df_var, ascending=False, ignore_index=True, inplace=True)
    # Removing NA values from coordinates values
    data_1.dropna(subset=['modifCoord'], inplace=True)
    # Converting double entries into single
    data_1['modifCoord'] = data_1['modifCoord'].astype(str)
    # Creating new accession ID by combing Gene name and AA coordinates
    data_1['phospho_site'] = data_1['geneName'] + '_' + data_1['modifCoord']
    data_1['time_stamp'] = str(time_stamp)
    # Get only significant entries
    data_signi = data_1[data_1['pValue_Condition' + df_var] <= 0.05]
    # dropping duplicated phosphosites
    data_signi_1 = data_signi.drop_duplicates('phospho_site', keep='first')
    # Considering differentially expressed phosphosites
    up_down_prt_1 = data_signi_1[(data_signi_1['log2ratio_Condition' + df_var] > fold_change) |
                                 (data_signi_1['log2ratio_Condition' + df_var] < fold_change)]
    print("number of initial signi entries:", up_down_prt_1.shape[0])
    up_down_prt = up_down_prt_1.rename(columns={'log2ratio_Condition' + df_var: 'log2ratio_Condition2',
                                                'pValue_Condition' + df_var: 'pValue_Condition2',
                                                'medianInt_Condition' + df_var: 'medianInt_Condition2'})
    if ctr is not None:  # removing auxin effect
        ctr_data = ctr.rename(columns={'log2ratio_Condition' + str(ctr_var): 'log2ratio_Condition2'})
        df_auxin_data_merge = up_down_prt.merge(ctr_data, how='inner', on='phospho_site', suffixes=("_df", "_aux"))
        df_auxin_data_merge['potential_auxin_effect'] = df_auxin_data_merge["log2ratio_Condition2_df"] * \
                                                        df_auxin_data_merge['log2ratio_Condition2_aux']
        aux_effect_prot = df_auxin_data_merge[df_auxin_data_merge['potential_auxin_effect'] > 0]['phospho_site']
        aux_removed_data = up_down_prt[~(up_down_prt['phospho_site'].isin(aux_effect_prot))]
        print("Number of entries which are auxin effect:", aux_effect_prot.shape[0])
        aux_opposite_entries = df_auxin_data_merge[df_auxin_data_merge['potential_auxin_effect'] < 0]['phospho_site']
        aux_opposite_entries1 = df_auxin_data_merge[df_auxin_data_merge['phospho_site'].isin(aux_opposite_entries)]
        print("phospho sites with opposite fold change:", aux_opposite_entries.shape[0])
        return aux_removed_data, aux_opposite_entries1
    else:
        return up_down_prt_1


# To get motif sequence from data
def motif_selection(dataframe, seq_length):
    dataframe1 = dataframe.copy()
    dataframe1['motifX'] = dataframe1['motifX'].str.split(',', expand=False)
    dataframe1['modifCoord'] = dataframe1['modifCoord'].str.split(',', expand=False)
    dataframe1['ptm'] = dataframe1['ptm'].str.split('|', expand=False)
    dataframe2 = dataframe1.explode(['motifX', 'modifCoord', 'ptm'])
    # dataframe3 = dataframe2[dataframe2['ptm'].str.contains("Phospho")]
    dataframe2['new_phospho_site'] = dataframe2['geneName'] + '_' + dataframe2['modifCoord']
    dataframe2['motifX'] = dataframe2['motifX'].apply(lambda x: ''.join(x.split('*')))
    dataframe3 = dataframe2[dataframe2['motifX'].str.len() > seq_length]
    df_with_motif = dataframe3.drop_duplicates("new_phospho_site")
    return df_with_motif


# Saving data into fasta format
def df_to_fasta(dataframe, df_col_name_for_header, df_col_name_for_seq, fasta_file):
    with open(str(fasta_file), 'w') as file:
        for index, rows in dataframe.iterrows():
            file.write('>' + str(rows[df_col_name_for_header]) + '\n'
                       + str(rows[df_col_name_for_seq]) + '\n')


def add_cross_at_na(data):
    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            if pd.isna(data.iloc[y,x]):
                ax.plot([x, x + 1], [y, y + 1], color='gray', lw=0.5)
                ax.plot([x, x + 1], [y + 1, y], color='gray', lw=0.5)


# Getting significant entries from each data
path_to_save_result = r"./re_analysis/ypk2_results/"

# TORC2 Aux NO-Aux 15 min
torc2_aux_noaux_15min_signi = signi_data(torc2_aux_noaux_15min, '15_min', df_var=str(2))
torc2_aux_noaux_15min_signi.to_csv(path_to_save_result + "TORC2_AUX_NO-AUX_15min.tsv", index=False, header=True,
                                   sep='\t')

torc2_aux_noaux_15min_signi_hypo = torc2_aux_noaux_15min_signi[
    torc2_aux_noaux_15min_signi["log2ratio_Condition2"] < fold_change]
torc2_aux_noaux_15min_signi_hypo_motif = motif_selection(torc2_aux_noaux_15min_signi_hypo, 8)
df_to_fasta(torc2_aux_noaux_15min_signi_hypo_motif, 'phospho_site', 'motifX',
            path_to_save_result + "TORC2_AUX_NO-AUX_15min_hypo_motif_seq.txt")

torc2_aux_noaux_15min_signi_hyper = torc2_aux_noaux_15min_signi[
    torc2_aux_noaux_15min_signi["log2ratio_Condition2"] > fold_change]
torc2_aux_noaux_15min_signi_hyper_motif = motif_selection(torc2_aux_noaux_15min_signi_hyper, 8)
df_to_fasta(torc2_aux_noaux_15min_signi_hyper_motif, 'phospho_site', 'motifX',
            path_to_save_result + "TORC2_AUX_NO-AUX_15min_hyper_motif_seq.txt")

# Auxin effect
afb2_time_course_data_path = r"./re_analysis/phospho_proteome_time_course/afb2_time_course_phospho/"
afb2_30min_tc = pd.read_csv(afb2_time_course_data_path + "afb2_30min_signi_time_course.tsv", sep='\t')

# TORC2 Aux NO-Aux 30 min
torc2_aux_noaux_30min_signi, torc2_aux_noaux_30min_aux_opposite_entry = signi_data(torc2_aux_noaux_30min, '30_min',
                                                                                   afb2_30min_tc, ctr_var=str(2),
                                                                                   df_var=str(6))
torc2_aux_noaux_30min_signi.to_csv(path_to_save_result + "TORC2_AUX_NO-AUX_30min.tsv", index=False, header=True, sep='\t')

torc2_aux_noaux_30min_signi_hypo = torc2_aux_noaux_30min_signi[torc2_aux_noaux_30min_signi["log2ratio_Condition2"] < fold_change]
torc2_aux_noaux_30min_signi_hypo_motif = motif_selection(torc2_aux_noaux_30min_signi_hypo, 8)
df_to_fasta(torc2_aux_noaux_30min_signi_hypo_motif, 'phospho_site', 'motifX', path_to_save_result + "TORC2_AUX_NO-AUX_30min_hypo_motif_seq.txt")

torc2_aux_noaux_30min_signi_hyper = torc2_aux_noaux_30min_signi[torc2_aux_noaux_30min_signi["log2ratio_Condition2"] > fold_change]
torc2_aux_noaux_30min_signi_hyper_motif = motif_selection(torc2_aux_noaux_30min_signi_hyper, 8)
df_to_fasta(torc2_aux_noaux_30min_signi_hyper_motif, 'phospho_site', 'motifX', path_to_save_result + "TORC2_AUX_NO-AUX_30min_hyper_motif_seq.txt")

# TORC2 Aux NO-Aux 60 min
afb2_60min_tc = pd.read_csv(afb2_time_course_data_path + "afb2_1hr_signi_time_course.tsv", sep='\t')

torc2_aux_noaux_60min_signi, torc2_aux_noaux_60min_aux_opposite_entry = \
    signi_data(torc2_aux_noaux_60min, '60_min', afb2_60min_tc, ctr_var=str(3), df_var=str(10))
torc2_aux_noaux_60min_signi.to_csv(path_to_save_result + "TORC2_AUX_NO-AUX_60min.tsv", index=False, header=True, sep='\t')

torc2_aux_noaux_60min_signi_hypo = torc2_aux_noaux_60min_signi[torc2_aux_noaux_60min_signi["log2ratio_Condition2"] < fold_change]
torc2_aux_noaux_60min_signi_hypo_motif = motif_selection(torc2_aux_noaux_60min_signi_hypo, 8)
df_to_fasta(torc2_aux_noaux_60min_signi_hypo_motif, 'phospho_site', 'motifX', path_to_save_result + "TORC2_AUX_NO-AUX_60min_hypo_motif_seq.txt")

torc2_aux_noaux_60min_signi_hyper = torc2_aux_noaux_60min_signi[torc2_aux_noaux_60min_signi["log2ratio_Condition2"] > fold_change]
torc2_aux_noaux_60min_signi_hyper_motif = motif_selection(torc2_aux_noaux_60min_signi_hyper, 8)
df_to_fasta(torc2_aux_noaux_60min_signi_hyper_motif, 'phospho_site', 'motifX',path_to_save_result + "TORC2_AUX_NO-AUX_60min_hyper_motif_seq.txt")

# YPK2 Aux NoAux 15 min
ypk2_aux_noaux_15min_signi = signi_data(ypk2_aux_noaux_15min, '15_min', df_var=str(4))
ypk2_aux_noaux_15min_signi.to_csv(path_to_save_result+'YPK2_AUX_NO-AUX_15min.tsv', index=False, sep='\t')


ypk2_aux_noaux_15min_signi_hypo = ypk2_aux_noaux_15min_signi[ypk2_aux_noaux_15min_signi["log2ratio_Condition4"] < fold_change]
# ypk2_aux_noaux_15min_signi_hypo_motif = motif_selection(ypk2_aux_noaux_15min_signi_hypo, 8)
# df_to_fasta(ypk2_aux_noaux_15min_signi_hypo_motif, 'phospho_site', 'motifX', path_to_save_result+'ypk2_aux_noaux_15min_signi_hypo_motif_seq.txt')

ypk2_aux_noaux_15min_signi_hyper = ypk2_aux_noaux_15min_signi[ypk2_aux_noaux_15min_signi["log2ratio_Condition4"] > fold_change]
# ypk2_aux_noaux_15min_signi_hyper_motif = motif_selection(ypk2_aux_noaux_15min_signi_hyper, 8)
# df_to_fasta(ypk2_aux_noaux_15min_signi_hyper_motif, 'phospho_site', 'motifX', path_to_save_result+'ypk2_aux_noaux_15min_signi_hyper_motif_seq.txt')

# YPK2 Aux NoAux 30 min
ypk2_aux_noaux_30min_signi = signi_data(ypk2_aux_noaux_30min, '30_min', df_var=str(8))
ypk2_aux_noaux_30min_signi.to_csv(path_to_save_result+'YPK2_AUX_NO-AUX_30min.tsv', index=False, sep='\t')

ypk2_aux_noaux_30min_signi_hypo = ypk2_aux_noaux_30min_signi[ypk2_aux_noaux_30min_signi["log2ratio_Condition8"] < fold_change]
ypk2_aux_noaux_30min_signi_hypo_motif = motif_selection(ypk2_aux_noaux_30min_signi_hypo, 8)
df_to_fasta(ypk2_aux_noaux_30min_signi_hypo_motif, 'phospho_site', 'motifX', path_to_save_result+'ypk2_aux_noaux_30min_signi_hypo_motif_seq.txt')

ypk2_aux_noaux_30min_signi_hyper = ypk2_aux_noaux_30min_signi[ypk2_aux_noaux_30min_signi["log2ratio_Condition8"] > fold_change]
ypk2_aux_noaux_30min_signi_hyper_motif = motif_selection(ypk2_aux_noaux_30min_signi_hyper, 8)
df_to_fasta(ypk2_aux_noaux_30min_signi_hyper_motif, 'phospho_site', 'motifX', path_to_save_result+'ypk2_aux_noaux_30min_signi_hyper_motif_seq.txt')

# YPK2 Aux NoAux 60 min
ypk2_aux_noaux_60min_signi = signi_data(ypk2_aux_noaux_60min, '60_min', df_var=str(12))
ypk2_aux_noaux_60min_signi.to_csv(path_to_save_result+'YPK2_AUX_NO-AUX_60min.tsv', index=False, sep='\t')

ypk2_aux_noaux_60min_signi_hypo = ypk2_aux_noaux_60min_signi[ypk2_aux_noaux_60min_signi["log2ratio_Condition12"] < fold_change]
# ypk2_aux_noaux_60min_signi_hypo_motif = motif_selection(ypk2_aux_noaux_60min_signi_hypo, 8)
# df_to_fasta(ypk2_aux_noaux_60min_signi_hypo_motif, 'phospho_site', 'motifX', path_to_save_result+'ypk2_aux_noaux_60min_signi_hypo_motif_seq.txt')

ypk2_aux_noaux_60min_signi_hyper = ypk2_aux_noaux_60min_signi[ypk2_aux_noaux_60min_signi["log2ratio_Condition12"] > fold_change]
# ypk2_aux_noaux_60min_signi_hyper_motif = motif_selection(ypk2_aux_noaux_60min_signi_hyper, 8)
# df_to_fasta(ypk2_aux_noaux_60min_signi_hyper_motif, 'phospho_site', 'motifX', path_to_save_result+'ypk2_aux_noaux_60min_signi_hyper_motif_seq.txt')

# YPK2 independent substrates motif
ypk2_30m_60m_union = ypk2_aux_noaux_30min_signi.merge(ypk2_aux_noaux_60min_signi, how='outer', on="phospho_site", suffixes=("_30min", "_60min"))
ypk2_30m_60m_union_hypo = ypk2_30m_60m_union[ypk2_30m_60m_union['log2ratio_Condition8']<0]
ypk2_30m_60m_union_hypo1 = ypk2_30m_60m_union_hypo[['ac_30min', 'geneName_30min', 'motifX_30min', 'modifCoord_30min', 'ptm_30min', 'log2ratio_Condition8',
                                                    'pValue_Condition8', 'phospho_site', 'log2ratio_Condition12', 'pValue_Condition12']]
ypk2_30m_60m_union_hypo2 = ypk2_30m_60m_union_hypo1.rename(columns={'geneName_30min': 'geneName',
                                                                    'motifX_30min': 'motifX',
                                                                    'modifCoord_30min': 'modifCoord',
                                                                    'ptm_30min': 'ptm'})
ypk2_30m_60m_union_hypo3 = motif_selection(ypk2_30m_60m_union_hypo2, 8)
ypk2_30m_60m_union_hypo3['AGC motif'] = ypk2_30m_60m_union_hypo3['motifX'].str.contains(r"R.R..[ST]", regex=True)
ypk2_30m_60m_union_hypo3['PKC motif'] = ypk2_30m_60m_union_hypo3['motifX'].str.contains(r"(R..[ST].R)|([ST].R)", regex=True)
ypk2_30m_60m_union_hypo3['TOR motif'] = ypk2_30m_60m_union_hypo3['motifX'].str.contains(r"[F].[ALVIPFMW][ED][ALVIPFMW]", regex=True)
ypk2_30m_60m_union_hypo3['seq_length'] = ypk2_30m_60m_union_hypo3['motifX'].str.len()
ypk2_30m_60m_union_hypo3.to_csv(path_to_save_result+'ypk2_30m_60m_union_motif.tsv', sep='\t')
ypk2_30m_60m_union_hypo4 = ypk2_30m_60m_union_hypo3[ypk2_30m_60m_union_hypo3['ptm'].str.contains('Phospho')]
highest_seq_len = ypk2_30m_60m_union_hypo4['seq_length'].value_counts()
potential_torc2_motif = ypk2_30m_60m_union_hypo4[~((ypk2_30m_60m_union_hypo3['AGC motif']) &
                                                 (ypk2_30m_60m_union_hypo3['PKC motif']) &
                                                 (ypk2_30m_60m_union_hypo3['TOR motif']))]
potential_torc2_motif_1 = potential_torc2_motif[potential_torc2_motif['seq_length']==highest_seq_len.index[0]]
df_to_fasta(potential_torc2_motif_1, 'phospho_site', 'motifX', path_to_save_result+'potential_torc2_motif.txt')


############################################ graphs for
# venn diagram for TORC2 all-time points
plt.figure(figsize=(12, 12))
a = venn3([set(torc2_aux_noaux_15min_signi['phospho_site']), set(torc2_aux_noaux_30min_signi['phospho_site']),
           set(torc2_aux_noaux_60min_signi['phospho_site'])],
          set_labels=('15_min'+'\nN='+str(torc2_aux_noaux_15min_signi.shape[0]),
                      '30_min'+'\nN='+str(torc2_aux_noaux_30min_signi.shape[0]),
                      '60_min'+'\nN='+str(torc2_aux_noaux_60min_signi.shape[0])), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(20)
    text.set_fontname('arial')
for text in a.set_labels:
    text.set_fontsize(18)
    text.set_fontname('arial')
plt.title("TORC2 phospho common between 15 min, 30 min and 60 min")
plt.savefig(path_to_save_result+'TORC2_15m_30m_60m_common_venn.tif', dpi=600)
plt.close()

# common between TORC2 30min and 60 min
torc2_30m_60m_com = torc2_aux_noaux_30min_signi.merge(torc2_aux_noaux_60min_signi, how='inner', on='phospho_site', suffixes=('_30min', '_60_min'))
torc2_30m_60m_com.to_csv(path_to_save_result+'TORC2_30m-60m_common_from_ypk2_expt.tsv', index=False, header=True, sep='\t')

# venn diagram for 30m_60m_common
plt.figure(figsize=(12, 12))
a = venn2([set(torc2_aux_noaux_30min_signi['phospho_site']), set(torc2_aux_noaux_60min_signi['phospho_site'])],
          set_labels=('30 min'+'\nN='+ str(torc2_aux_noaux_30min_signi.shape[0]), '60 min'+'\nN='+ str(torc2_aux_noaux_60min_signi.shape[0])),
          set_colors=('#D13444', '#D1C134'), alpha=1)
for text in a.subset_labels:
    text.set_fontsize(30)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in a.set_labels:
    text.set_fontsize(18)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("TORC2(YPK_Data) phospho common between 30 min and 60 min")
plt.savefig(path_to_save_result+'TORC2(YPK)_30m_60m_common_venn.tif', dpi=600)
plt.close()

# Loading TORC2 entries from TORC2 pairwise expt (where TORC2 Aux compared with Afb2 Aux)
ph_torc2_30m_signi = pd.read_csv(r"./re_analysis/phospho_proteome_pair_wise/torc2_pairwise_phospho/torc2_30_min_norm.tsv", sep='\t')
ph_torc2_1h_signi = pd.read_csv(r"./re_analysis/phospho_proteome_pair_wise/torc2_pairwise_phospho/torc2_1hr_norm.tsv", sep='\t')
ph_torc2_30m_1h_com = pd.read_csv(r"./re_analysis/phospho_proteome_graphs/TORC2_30m_1hr_common_phospho.tsv", sep='\t')

# Comparing TORC2 30min-60min common(YPK expt) AND TORC2 30min-1hr common(pairwise expt)

torc2_30m_60m_com_from_both_expt = torc2_30m_60m_com.merge(ph_torc2_30m_1h_com, how='inner', on='phospho_site', suffixes=('_ypk', '_pairwise'))
torc2_30m_60m_com_from_both_expt.to_csv(path_to_save_result+'TORC2_30m-60m_common_from_both_expt.tsv', index=False, header=True, sep='\t')

# Venn diagram for TORC2 30min-60min common(YPK expt) AND TORC2 30min-1hr common(pairwise expt)
plt.figure(figsize=(12, 12))
a = venn2([set(torc2_30m_60m_com['phospho_site']), set(ph_torc2_30m_1h_com['phospho_site'])],
          set_labels=('TORC2 30m-60m (YPK)'+'\nN='+ str(torc2_30m_60m_com.shape[0]), 'TORC2 30min-1hr (PairWise)'+'\nN='+ str(ph_torc2_30m_1h_com.shape[0])),
          set_colors=(('orange'), ('blue')), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(30)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in a.set_labels:
    text.set_fontsize(18)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("TORC2 30m-60m(YPK_Data) & 30min-1hr (PairWise) common")
plt.savefig(path_to_save_result+'TORC2_30m-60m_common_from_both(ypk-pairwise)_expt_venn.tif', dpi=600)
plt.close()

# Getting TORC2 30m-60m(YPK2) specific phospho-site compare to TORC2 30m-1hr(pairwise)
torc2_30m_60m_com_only = torc2_30m_60m_com[~(torc2_30m_60m_com['phospho_site'].isin(ph_torc2_30m_1h_com['phospho_site']))]
torc2_30m_60m_com_only.to_csv(path_to_save_result+'TORC2_30-60m(YPK)_specific_compare_TO_TORC2_30m-1h(pairwise).tsv', index=False, sep='\t')

# Getting TORC2 30m-1hr(pairwise) specific phospho-site compare to TORC2 30m-60m(YPK)
ph_torc2_30m_1h_com_only = ph_torc2_30m_1h_com[~(ph_torc2_30m_1h_com['phospho_site'].isin(torc2_30m_60m_com['phospho_site']))]
ph_torc2_30m_1h_com_only.to_csv(path_to_save_result+'TORC2_30m-1h(pairwise)_specific_compare_TO_TORC2_30m-60m(ypk2).tsv', index=False, sep='\t')


# YPK2-CA data set analysis

# venn diagram for YPK2-CA all-time points
plt.figure(figsize=(12, 12))
a = venn3([set(ypk2_aux_noaux_15min_signi['phospho_site']), set(ypk2_aux_noaux_30min_signi['phospho_site']),
           set(ypk2_aux_noaux_60min_signi['phospho_site'])],
          set_labels=('15_min'+'\nN='+str(ypk2_aux_noaux_15min_signi.shape[0]),
                      '30_min'+'\nN='+str(ypk2_aux_noaux_30min_signi.shape[0]),
                      '60_min'+'\nN='+str(ypk2_aux_noaux_60min_signi.shape[0])), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(25)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in a.set_labels:
    text.set_fontsize(18)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("YPK2-CA phospho common between 15 min, 30 min and 60 min")
plt.savefig(path_to_save_result+'YPK2-CA_15m_30m_60m_common_venn.tif', dpi=600)
plt.close()

# YPK2-CA 30min and 60 min common
ypk2_30m_60m_com = ypk2_aux_noaux_30min_signi.merge(ypk2_aux_noaux_60min_signi, how='inner', on='phospho_site', suffixes=('_30_min', '_60_min'))
ypk2_30m_60m_com.to_csv(path_to_save_result+'YPK2_30m_60m_common.tsv', index=False, sep='\t')

# Venn diagram for YPK2-CA 30min and 60 min common
plt.figure(figsize=(12, 12))
a = venn2([set(ypk2_aux_noaux_30min_signi['phospho_site']), set(ypk2_aux_noaux_60min_signi['phospho_site'])],
          set_labels=('YPK2-CA 30 min'+'\nN='+str(ypk2_aux_noaux_30min_signi.shape[0]), 'YPK2-CA 60 min'+'\nN='+str(ypk2_aux_noaux_60min_signi.shape[0])),
          set_colors=('#2DE087', '#2D87E0'), alpha=0.8)
for text in a.subset_labels:
    text.set_fontsize(30)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in a.set_labels:
    text.set_fontsize(18)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("YPK2-CA phospho common between 30 min and 60 min")
plt.savefig(path_to_save_result+'ypk2-CA_phospho_common_between_30min_60_min_venn.tif', dpi=600)
plt.close()
# common between YPK2-CA 30m-60m common and TORC2 30m-60m common from both expt(YPK2 & Pairwise)

ypk2_independent_torc2_both_expt = ypk2_30m_60m_com.merge(torc2_30m_60m_com_from_both_expt, how='inner', on='phospho_site', suffixes=('_YPK2-CA', '_TORC2'))
ypk2_independent_torc2_both_expt.to_csv(path_to_save_result+'YPK2_independent_TORC2_potential_candidate_compare_both_expt.tsv', index=False, sep='\t')


# Venn diagram for potential YPK2 independent targets
plt.figure(figsize=(12, 12))
a = venn2([set(ypk2_30m_60m_com['phospho_site']), set(torc2_30m_60m_com_from_both_expt['phospho_site'])],
          set_labels=('YPK2-CA 30m-60m'+'\nN='+ str(ypk2_30m_60m_com.shape[0]), 'TORC2 30min-1hr both expt'+'\nN='+ str(torc2_30m_60m_com_from_both_expt.shape[0])),
          set_colors=('orange', 'blue'), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(30)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in a.set_labels:
    text.set_fontsize(18)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("Potential YPK2 independent TORC2 candidates compare both expt")
plt.savefig(path_to_save_result+'Potential_YPK2_independent_TORC2_candidates_compare_both_expt_venn.tif', dpi=600)
plt.close()

# Venn diagram for potential YPK2 independent targets when compare with TORC2 30m-1hr (pairwise) expt

ypk2_independent_torc2_pairwise = ypk2_30m_60m_com.merge(ph_torc2_30m_1h_com, how='inner', on='phospho_site', suffixes=('_YPK2-CA', '_TORC2(pairwise)'))
ypk2_independent_torc2_pairwise.to_csv(path_to_save_result+'YPK2_independent_TORC2_potential_candidate_compare_TORC2_pairwise.tsv', index=False, sep='\t')


# Venn diagram for potential YPK2 independent targets
plt.figure(figsize=(12, 12))
a = venn2([set(ypk2_30m_60m_com['phospho_site']), set(ph_torc2_30m_1h_com['phospho_site'])],
          set_labels=('YPK2-CA 30m-60m'+'\nN='+ str(ypk2_30m_60m_com.shape[0]), 'TORC2 30min-1hr pairwise'+'\nN='+ str(ph_torc2_30m_1h_com.shape[0])),
          set_colors=('orange', 'blue'), alpha=0.6)
for text in a.subset_labels:
    text.set_fontsize(30)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in a.set_labels:
    text.set_fontsize(18)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("Potential YPK2 independent TORC2 candidates compare TORC2 pairwise")
plt.savefig(path_to_save_result+'Potential_YPK2_independent_TORC2_candidates_compare_TORC2_pairwise_venn.tif', dpi=600)
plt.close()

# Venn diagram for potential YPK2 independent targets when compare with TORC2 30m-60m (ypk2) expt

ypk2_independent_torc2 = ypk2_30m_60m_com.merge(torc2_30m_60m_com, how='inner', on='phospho_site', suffixes=('_YPK2-CA', '_TORC2(YPK2-CA)'))
ypk2_independent_torc2.to_csv(path_to_save_result+'YPK2_independent_TORC2_potential_candidate_compare_YPK2-CA.tsv', index=False, sep='\t')

# Venn diagram for potential YPK2 independent targets
plt.figure(figsize=(12, 12))
a = venn2([set(ypk2_30m_60m_com['phospho_site']), set(torc2_30m_60m_com['phospho_site'])],
          set_labels=('YPK2-CA 30m-60m'+'\nN='+ str(ypk2_30m_60m_com.shape[0]), 'TORC2 30min-1hr(YPK2-CA)'+'\nN='+ str(torc2_30m_60m_com.shape[0])),
          set_colors=('#2DE0E0', '#ffab54'), alpha=0.8)
for text in a.subset_labels:
    text.set_fontsize(30)
    text.set_fontname('arial')
    text.set_weight('bold')
for text in a.set_labels:
    text.set_fontsize(18)
    text.set_fontname('arial')
    text.set_weight('bold')
plt.title("Potential YPK2 independent TORC2 candidates compare TORC2 YPK2-CA")
plt.savefig(path_to_save_result+'Potential_YPK2_independent_TORC2_candidates_compare_TORC2_YPK2-CA_venn.tif', dpi=600)
plt.close()

# Getting motif for YPK2 independent candidates
ypk2_independent_torc2_hypo = ypk2_independent_torc2[ypk2_independent_torc2['log2ratio_Condition8'] < 0]
ypk2_independent_torc2_hypo1 = ypk2_independent_torc2_hypo[['ac_30_min', 'geneName_30_min', 'modifCoord_30_min', 'ptm_30_min',
       'motifX_30_min', 'log2ratio_Condition8', 'pValue_Condition8', 'phospho_site', 'log2ratio_Condition12', 'pValue_Condition12']]
ypk2_independent_torc2_hypo2 = ypk2_independent_torc2_hypo1.rename(columns=({'geneName_30_min': 'geneName',
                               'modifCoord_30_min': 'modifCoord', 'ptm_30_min': 'ptm', 'motifX_30_min': 'motifX'}))
ypk2_independent_torc2_hypo3 = motif_selection(ypk2_independent_torc2_hypo2, 8)

ypk2_independent_torc2_hypo3['AGC motif'] = ypk2_independent_torc2_hypo3['motifX'].str.contains(r"R.R..[ST]", regex=True)
ypk2_independent_torc2_hypo3['PKC motif'] = ypk2_independent_torc2_hypo3['motifX'].str.contains(r"(R..[ST].R)|([ST].R)", regex=True)
ypk2_independent_torc2_hypo3['TOR motif'] = ypk2_independent_torc2_hypo3['motifX'].str.contains(r"[F].[ALVIPFMW][ED][ALVIPFMW]", regex=True)
ypk2_independent_torc2_hypo3['seq_length'] = ypk2_independent_torc2_hypo3['motifX'].str.len()
ypk2_independent_torc2_hypo3.to_csv(path_to_save_result+'ypk2_independent_torc2_hypo_motif.tsv', sep='\t')
ypk2_independent_torc2_hypo4 = ypk2_independent_torc2_hypo3[ypk2_independent_torc2_hypo3['ptm'].str.contains('Phospho')]
highest_seq_len = ypk2_independent_torc2_hypo4['seq_length'].value_counts()
ypk2_independent_torc2_potential_torc2_motif = ypk2_independent_torc2_hypo4[~((ypk2_independent_torc2_hypo3['AGC motif']) &
                                                 (ypk2_independent_torc2_hypo3['PKC motif']) &
                                                 (ypk2_independent_torc2_hypo3['TOR motif']))]
ypk2_independent_torc2_potential_torc2_motif_1 = ypk2_independent_torc2_potential_torc2_motif[
    ypk2_independent_torc2_potential_torc2_motif['seq_length'] == highest_seq_len.index[0]]
df_to_fasta(ypk2_independent_torc2_potential_torc2_motif_1, 'phospho_site', 'motifX', path_to_save_result+'ypk2_independent_torc2_potential_torc2_motif.txt')


ypk2_independent_torc2_in_ypk2_30_60_com = ypk2_30m_60m_com[ypk2_30m_60m_com['phospho_site'].isin(ypk2_independent_torc2['phospho_site'])][['phospho_site', 'log2ratio_Condition8', 'log2ratio_Condition12']]
ypk2_independent_torc2_in_ypk2_30_60_com.rename(columns={'log2ratio_Condition8':'30min_YPK2', 'log2ratio_Condition12':'60min_YPK2'}, inplace=True)
ypk2_independent_torc2_in_ypk2_30_60_com['30_60_avg_ypk2'] = ypk2_independent_torc2_in_ypk2_30_60_com[['30min_YPK2', '60min_YPK2']].mean(axis=1)

ypk2_independent_torc2_in_torc2_30_60_com = torc2_30m_60m_com[torc2_30m_60m_com['phospho_site'].isin(ypk2_independent_torc2['phospho_site'])][['phospho_site', 'log2ratio_Condition2_30min', 'log2ratio_Condition2_60_min']]
ypk2_independent_torc2_in_torc2_30_60_com.rename(columns={'log2ratio_Condition2_30min': '30min_TORC2', 'log2ratio_Condition2_60_min': '60min_TORC2'}, inplace=True)
ypk2_independent_torc2_in_torc2_30_60_com['30_60_avg_torc2'] = ypk2_independent_torc2_in_torc2_30_60_com[['30min_TORC2', '60min_TORC2']].mean(axis=1)

ypk2_independent_torc2_in_ypk2_torc2 = ypk2_independent_torc2_in_ypk2_30_60_com.merge(ypk2_independent_torc2_in_torc2_30_60_com, on='phospho_site').sort_values('30min_YPK2')
ypk2_independent_torc2_in_ypk2_torc2['opposite_pattern'] = ypk2_independent_torc2_in_ypk2_torc2['30min_YPK2'] * ypk2_independent_torc2_in_ypk2_torc2['60min_YPK2'] * ypk2_independent_torc2_in_ypk2_torc2['30min_TORC2']*ypk2_independent_torc2_in_ypk2_torc2['60min_TORC2']
ypk2_independent_torc2_in_ypk2_torc2['phospho_diff'] = ypk2_independent_torc2_in_ypk2_torc2['30_60_avg_ypk2']-ypk2_independent_torc2_in_ypk2_torc2['30_60_avg_torc2']
ypk2_independent_torc2_in_ypk2_torc2.to_excel(path_to_save_result + 'ypk2_independent_torc2_analysis.xlsx', header=True)

ypk2_independent_torc2_heatmap = ypk2_independent_torc2_in_ypk2_torc2[['phospho_site', '30min_YPK2', '60min_YPK2','30min_TORC2', '60min_TORC2']]
ypk2_independent_torc2_heatmap.set_index('phospho_site', inplace=True)
ypk2_independent_torc2_heatmap.to_csv(path_to_save_result+'ypk2_independent_torc2_heatmap_data.tsv', sep='\t')
# Heatmap
plt.figure(figsize=(20, 3))
ax = sns.heatmap(ypk2_independent_torc2_heatmap.T,
                 xticklabels=True, yticklabels=['30min_YPK2', '60min_YPK2','30min_TORC2', '60min_TORC2'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False,
                 linewidths=0.5 )
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=4, labelbottom = True, bottom=True, top = False)
ax.tick_params(axis='y', which='major', labelsize=7)
plt.subplots_adjust(bottom=0.25)
ax.set_title('Potential YPK2 independent TORC2 candidates')
plt.savefig(path_to_save_result+"potential_ypk2_independent_torc2_candidates_heatmap_horizontal.tif", dpi=600)
plt.close()

# potential ypk2_independent_torc2_candidates only hyper
ypk2_independent_torc2_heatmap_hyper = ypk2_independent_torc2_heatmap[ypk2_independent_torc2_heatmap['30min_YPK2'] > 0]
plt.figure(figsize=(20, 4))
ax = sns.heatmap(ypk2_independent_torc2_heatmap_hyper.T,
                 xticklabels=True, yticklabels=['30min_YPK2', '60min_YPK2','30min_TORC2', '60min_TORC2'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False,
                 linewidths=0.5 )
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=6, labelbottom = True, bottom=True, top = False)
ax.tick_params(axis='y', which='major', labelsize=7)
plt.subplots_adjust(bottom=0.25)
ax.set_title('Potential YPK2 independent TORC2 only hyper-phospho candidates')
plt.savefig(path_to_save_result+"potential_ypk2_independent_torc2_candidates_heatmap_hyper_horizontal.tif", dpi=600)
plt.close()

# potential ypk2_independent_torc2_candidates only hypo
ypk2_independent_torc2_heatmap_hypo = ypk2_independent_torc2_heatmap[ypk2_independent_torc2_heatmap['30min_YPK2'] < 0]
plt.figure(figsize=(20, 3))
ax = sns.heatmap(ypk2_independent_torc2_heatmap_hypo.T,
                 xticklabels=True, yticklabels=['30min_YPK2', '60min_YPK2','30min_TORC2', '60min_TORC2'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False,
                 linewidths=0.5 )
ax.set_ylabel("Time")
ax.set_xlabel("phospho sites")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=6, labelbottom = True, bottom=True, top = False)
ax.tick_params(axis='y', which='major', labelsize=7)
plt.subplots_adjust(bottom=0.25)
ax.set_title('Potential YPK2 independent TORC2 only hypo-phospho candidates')
plt.savefig(path_to_save_result+"potential_ypk2_independent_torc2_candidates_heatmap_hypo_horizontal.tif", dpi=600)
plt.close()
# 30m-1hr common from (YPK2-TORC2) and 30m-1hr common from TORC2(Pairwise expt)

ypk2_independent_torc2_compare_ph_30m_1h_common = ypk2_independent_torc2.merge(ph_torc2_30m_1h_com, on='phospho_site',
                                                                               how='inner')
ypk2_independent_torc2_compare_ph_30m_1h_common.to_csv(path_to_save_result+'YPK2_independent_TORC2_compare_pairwise_30m_60m_common.tsv',
                                                       sep='\t', index=False)
# Venn diagram for 30m-1hr common from (YPK2-TORC2) and 30m-1hr common from TORC2(Pairwise expt)
# This gives more confident candidates that identified in both data set pairwise TORC2 expt and YPK2 expt
plt.figure(figsize=(12, 12))
a = venn2([set(ypk2_independent_torc2['phospho_site']), set(ph_torc2_30m_1h_com['phospho_site'])],
          set_labels=('YPK2 Independent'+'\nN='+ str(ypk2_independent_torc2.shape[0]), 'TORC2 30min-1hr (Pairwise)'+'\nN='+ str(ph_torc2_30m_1h_com.shape[0])),
          set_colors=('#dbffdf', '#d1abd1'), alpha=0.9)
for text in a.subset_labels:
    text.set_fontsize(30)
    text.set_weight('bold')
    text.set_fontname('Arial')
for text in a.set_labels:
    text.set_fontsize(18)
    text.set_weight('bold')
    text.set_fontname('Arial')
plt.title("Potential YPK2 independent TORC2 candidates compare with TORC2 pairwise")
plt.savefig(path_to_save_result+'Potential_YPK2_independent_TORC2_candidates_compare_torc2_pairwsie_venn.tif', dpi=600)
plt.close()

# Getting TORC2 all time (30min to 5hr) data for heat map
ph_torc2_30m_1h_com_all_time_data = pd.read_csv(r'.\re_analysis\phospho_proteome_graphs\ph_torc2_30m_1h_com_all_data.tsv', sep='\t')

ypk2_independent_torc2_compare_ph_30m_1h_common_in_torc2_all_time = ph_torc2_30m_1h_com_all_time_data[
    ph_torc2_30m_1h_com_all_time_data['phospho_site'].isin(ypk2_independent_torc2_compare_ph_30m_1h_common['phospho_site'])]

ypk2_independent_in_torc2_ph = ypk2_independent_torc2_compare_ph_30m_1h_common_in_torc2_all_time[
    ['phospho_site', '30_min_log2ratio', '1_hr_log2ratio', '2_hr_log2ratio', '3_hr_log2ratio', '5_hr_log2ratio']]

# Heatmap
ypk2_independent_in_torc2_ph.set_index('phospho_site', inplace=True)


plt.figure(figsize=(5, 10))
ax = sns.heatmap(ypk2_independent_in_torc2_ph,
                 yticklabels=True, xticklabels=['0.5_hr', '1_hr', '2_hr', '3_hr', '5_hr'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False,
                 linewidths=0.5 )
ax.set_ylabel("phospho sites")
ax.set_xlabel("Time")
ax.set_facecolor('white')
add_cross_at_na(ypk2_independent_in_torc2_ph)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=10, labelbottom = True, bottom=True, top = False)
ax.tick_params(axis='y', which='major', labelsize=7)
plt.subplots_adjust(left=0.25)
ax.set_title('YPK2 independent TORC2 in TORC2 pairwise')
plt.savefig(path_to_save_result+"YPK2_independent_TORC2_in_TORC2_pairwise_vertical.tif", dpi=600)
plt.close()

# ypk2 independent compared pairwise TORC2 with phospho all time and protein all time

prt_torc2_all_time = pd.read_csv(r"./re_analysis/phospho_proteome_graphs/protein_torc2_all_time_data.tsv", sep='\t')

ypk2_independent_in_torc2_ph1 = ypk2_independent_in_torc2_ph.reset_index()
ypk2_independent_in_torc2_ph1['geneName'] = ypk2_independent_in_torc2_ph1['phospho_site'].str.split('_').str[0]

ypk2_independent_torc2_prot = ypk2_independent_in_torc2_ph1.merge(prt_torc2_all_time, how='left', on='geneName', suffixes=('_phospho', '_proteome'))
ypk2_independent_torc2_prot1 = ypk2_independent_torc2_prot.sort_values('phospho_site')

# Heatmap
ypk2_independent_torc2_prot_heatmap = ypk2_independent_torc2_prot[['phospho_site', '30_min_log2ratio', '1_hr_log2ratio', '2_hr_log2ratio',
       '3_hr_log2ratio', '5_hr_log2ratio', '30_min', '1_hr', '2_hr', '3_hr', '5_hr']].set_index('phospho_site')
ypk2_independent_torc2_prot_heatmap.to_csv(path_to_save_result+'ypk2_independent_torc2_prot_heatmap.tsv', sep='\t')

plt.figure(figsize=(8, 12))
ax = sns.heatmap(ypk2_independent_torc2_prot_heatmap,
                 yticklabels=True, xticklabels=['0.5hr_phospho', '1hr_phospho', '2hr_phospho', '3hr_phospho', '5hr_phospho',
                                                '0.5hr_prot', '1hr_prot', '2hr_prot', '3hr_prot', '5hr_prot'],
                 cmap=sns.color_palette("coolwarm", as_cmap=True),
                 cbar_kws={'label': 'log2 fold change'},
                 center=0, vmin=-4, vmax=4, square=False,
                 linewidths=0.5 )
ax.vlines([5], color='black', linestyle='--', linewidth=2, *ax.get_ylim())
ax.set_ylabel("phospho sites")
ax.set_xlabel("Time")
ax.set_facecolor('white')
add_cross_at_na(ypk2_independent_torc2_prot_heatmap)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.tick_params(axis='x', which='major', labelsize=10, labelbottom = True, bottom=True, top = False)
ax.tick_params(axis='y', which='major', labelsize=7)
plt.subplots_adjust(bottom=0.18)
ax.set_title('YPK2 independent TORC2 in TORC2 pairwise & Proteome')
plt.savefig(path_to_save_result+"YPK2_independent_TORC2_in_TORC2_pairwise_&_proteome_vertical.tif", dpi=600)
plt.close()
