import pandas as pd

# getting control(afb2) time course data
afb2_time_course_data_path = r"./re_analysis/phospho_proteome_time_course/afb2_time_course_phospho/"

afb2_30min_tc = pd.read_csv(afb2_time_course_data_path+"afb2_30min_signi_time_course.tsv", sep='\t')
afb2_1h_tc = pd.read_csv(afb2_time_course_data_path+"afb2_1hr_signi_time_course.tsv", sep='\t')
afb2_2h_tc = pd.read_csv(afb2_time_course_data_path+"afb2_2hr_signi_time_course.tsv", sep='\t')
afb2_3h_tc = pd.read_csv(afb2_time_course_data_path+"afb2_3hr_signi_time_course.tsv", sep='\t')
afb2_5h_tc = pd.read_csv(afb2_time_course_data_path+"afb2_5hr_signi_time_course.tsv", sep='\t')

# Getting control(afb2) pairwise data

# TORC2 (Avo3) pairwise data

torc2_pairwise_data_path = r"./TOR_AID_Mass_spec/Ctrl_Avo3_Kog1_Phosphoproteome_pairwise_Raw/Avo/"

torc2_zero_hr = pd.read_csv(torc2_pairwise_data_path+"SQ_P317_MH_JH_jh001_C0hvsA0h/SQ_P317_MH_JH_jh001_C0hvsA0h_PEPTIDE.tsv", sep="\t")
torc2_30_min = pd.read_csv(torc2_pairwise_data_path+"SQ_P317_MH_JH_jh001_C30minvsA30min/SQ_P317_MH_JH_jh001_C30minvsA30min_PEPTIDE.tsv", sep="\t")
torc2_1hr = pd.read_csv(torc2_pairwise_data_path+"SQ_P317_MH_JH_jh001_C1hvsA1h/SQ_P317_MH_JH_jh001_C1hvsA1h_PEPTIDE.tsv", sep="\t")
torc2_2hr = pd.read_csv(torc2_pairwise_data_path+"SQ_P317_MH_JH_jh001_C2hvsA2h/SQ_P317_MH_JH_jh001_C2hvsA2h_PEPTIDE.tsv", sep="\t")
torc2_3hr = pd.read_csv(torc2_pairwise_data_path+"SQ_P317_MH_JH_jh001_C3hvsA3h/SQ_P317_MH_JH_jh001_C3hvsA3h_PEPTIDE.tsv", sep="\t")
torc2_5hr = pd.read_csv(torc2_pairwise_data_path+"SQ_P317_MH_JH_jh001_C5hvsA5h/SQ_P317_MH_JH_jh001_C5hvsA5h_PEPTIDE.tsv", sep="\t")

# TORC1 (Kog1) pairwise data

torc1_pairwise_data_path = r"./TOR_AID_Mass_spec/Ctrl_Avo3_Kog1_Phosphoproteome_pairwise_Raw/kog/"

torc1_zero_hr = pd.read_csv(torc1_pairwise_data_path+"SQ_P317_MH_JH_jh001_C0hvsK0h/SQ_P317_MH_JH_jh001_C0hvsK0h_PEPTIDE.tsv", sep="\t")
torc1_30_min = pd.read_csv(torc1_pairwise_data_path+"SQ_P317_MH_JH_jh001_C30minvsK30min/SQ_P317_MH_JH_jh001_C30minvsK30min_PEPTIDE.tsv", sep="\t")
torc1_1hr = pd.read_csv(torc1_pairwise_data_path+"SQ_P317_MH_JH_jh001_C1hvsK1h/SQ_P317_MH_JH_jh001_C1hvsK1h_PEPTIDE.tsv", sep="\t")
torc1_2hr = pd.read_csv(torc1_pairwise_data_path+"SQ_P317_MH_JH_jh001_C2hvsK2h/SQ_P317_MH_JH_jh001_C2hvsK2h_PEPTIDE.tsv", sep="\t")
torc1_3hr = pd.read_csv(torc1_pairwise_data_path+"SQ_P317_MH_JH_jh001_C3hvsK3h/SQ_P317_MH_JH_jh001_C3hvsK3h_PEPTIDE.tsv", sep="\t")
torc1_5hr = pd.read_csv(torc1_pairwise_data_path+"SQ_P317_MH_JH_jh001_C5hvsK5h/SQ_P317_MH_JH_jh001_C5hvsK5h_PEPTIDE.tsv", sep="\t")

fold_change = 0


def signi_data(df, time_stamp, ctr=None, df_var=str(2), ctr_var=None):  # To check auxin effect
    data_1 = df[['ac', 'geneName', 'modifCoord','ptm', 'motifX', 'log2ratio_Condition' + df_var, 'pValue_Condition' + df_var,
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
    up_down_prt = data_signi_1[(data_signi_1['log2ratio_Condition' + df_var] > fold_change) |
                               (data_signi_1['log2ratio_Condition' + df_var] < fold_change)]
    print("number of initial signi entries:", up_down_prt.shape[0])
    if ctr is not None:  # removing auxin effect
        ctr_data = ctr.rename(columns={'log2ratio_Condition' + str(ctr_var): 'log2ratio_Condition2'})
        df_auxin_data_merge = up_down_prt.merge(ctr_data, how='inner', on='phospho_site', suffixes=("_df", "_aux"))
        df_auxin_data_merge['potential_auxin_effect'] = df_auxin_data_merge["log2ratio_Condition" + str(df_var) + "_df"] * \
                                                        df_auxin_data_merge['log2ratio_Condition' + str(df_var) + "_aux"]
        aux_effect_prot = df_auxin_data_merge[df_auxin_data_merge['potential_auxin_effect'] > 0]['phospho_site']
        aux_removed_data = up_down_prt[~(up_down_prt['phospho_site'].isin(aux_effect_prot))]
        print("Number of entries which are auxin effect:", aux_effect_prot.shape[0])
        aux_opposite_entries = df_auxin_data_merge[df_auxin_data_merge['potential_auxin_effect'] < 0]['phospho_site']
        aux_opposite_entries1 = df_auxin_data_merge[df_auxin_data_merge['phospho_site'].isin(aux_opposite_entries)]
        print("phospho sites with opposite fold change:", aux_opposite_entries.shape[0])
        return aux_removed_data, aux_opposite_entries1
    else:
        return up_down_prt


def norm_with_zero_hr(df, zero_hr_df):  # This check tagging effect
    df_zero_hr_df_merge = df.merge(zero_hr_df, how='inner', on='phospho_site', suffixes=("_df", "_zero_hr"))
    df_zero_hr_df_merge["zero_hr_effect"] = df_zero_hr_df_merge["log2ratio_Condition2_df"] * df_zero_hr_df_merge["log2ratio_Condition2_zero_hr"]
    zero_hr_effect = df_zero_hr_df_merge[df_zero_hr_df_merge["zero_hr_effect"] > 0]
    print("Number of entries in given df:", df.shape[0])
    print("number of phospho site affected due to tagging effect (zero_hr_effect):", zero_hr_effect.shape[0])
    df_without_zero_hr_effect = df[~(df['phospho_site'].isin(zero_hr_effect['phospho_site']))]
    print("Number of entries after removing zero hr entries:", df_without_zero_hr_effect.shape[0])
    return df_without_zero_hr_effect, zero_hr_effect


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
    df_with_motif['seq_length'] = df_with_motif['motifX'].str.len()
    return df_with_motif


# Saving data into fasta format
def df_to_fasta(dataframe, df_col_name_for_header, df_col_name_for_seq, fasta_file):
    with open(str(fasta_file), 'w') as file:
        for index, rows in dataframe.iterrows():
            file.write('>' + str(rows[df_col_name_for_header]) + '\n'
                       + str(rows[df_col_name_for_seq]) + '\n')


# TORC2 (Avo3) data
path_to_save_torc2_processed_data = r"./re_analysis/phospho_proteome_pair_wise/torc2_pairwise_phospho/"

# Represent the tagging effect (Avo3-AID compared with Afb2 at zero hr)
torc2_zero_hr_signi = signi_data(torc2_zero_hr, "zero_hr")
torc2_zero_hr_signi.to_csv(path_to_save_torc2_processed_data+"torc2_zero_hr_norm.tsv", header=True, index=False, sep='\t')

# TORC2 30 min is compared with afb2 30 min and require only zero min control gene
# Auxin effect removed
torc2_30_min_signi, torc2_30_min_aux_opposite_entry = signi_data(torc2_30_min, '30_min', afb2_30min_tc, df_var=str(2), ctr_var=str(2))
# Tagging effect removed
torc2_30_min_norm, torc2_30_min_zero_hr_effect = norm_with_zero_hr(torc2_30_min_signi, torc2_zero_hr_signi)
# Saving to tsv format
torc2_30_min_norm.to_csv(path_to_save_torc2_processed_data+"torc2_30_min_norm.tsv", header=True, index=False, sep='\t')
#
torc2_30_min_norm_hypo = torc2_30_min_norm[torc2_30_min_norm["log2ratio_Condition2"] < fold_change]
torc2_30_min_norm_hypo_motif = motif_selection(torc2_30_min_norm_hypo, 8)
torc2_30_min_norm_hypo_motif2 = torc2_30_min_norm_hypo_motif[torc2_30_min_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc2_30_min_norm_hypo_motif2_highest_seq_len = torc2_30_min_norm_hypo_motif2['seq_length'].value_counts()
torc2_30_min_norm_hypo_motif3 = torc2_30_min_norm_hypo_motif2[torc2_30_min_norm_hypo_motif2['seq_length']==torc2_30_min_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_30_min_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_30min_norm_hypo_motif_seq.txt")

torc2_30_min_norm_hyper = torc2_30_min_norm[torc2_30_min_norm["log2ratio_Condition2"] > fold_change]
torc2_30_min_norm_hyper_motif = motif_selection(torc2_30_min_norm_hyper, 8)
torc2_30_min_norm_hyper_motif2 = torc2_30_min_norm_hyper_motif[torc2_30_min_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc2_30_min_norm_hyper_motif2_highest_seq_len = torc2_30_min_norm_hyper_motif2['seq_length'].value_counts()
torc2_30_min_norm_hyper_motif3 = torc2_30_min_norm_hyper_motif2[torc2_30_min_norm_hyper_motif2['seq_length'] == torc2_30_min_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_30_min_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_30min_norm_hyper_motif_seq.txt")

# Torc2 1hr
# Auxin effect removed
torc2_1hr_signi, torc2_1hr_aux_opposite_entry = signi_data(torc2_1hr, '1_hr', afb2_1h_tc, df_var=str(2), ctr_var=str(3))
# Tagging effect removed
torc2_1hr_norm, torc2_1hr_zero_hr_effect = norm_with_zero_hr(torc2_1hr_signi, torc2_zero_hr_signi)
# Saving to tsv format
torc2_1hr_norm.to_csv(path_to_save_torc2_processed_data+"torc2_1hr_norm.tsv", header=True, index=False, sep='\t')
#
torc2_1hr_norm_hypo = torc2_1hr_norm[torc2_1hr_norm["log2ratio_Condition2"] < fold_change]
torc2_1hr_norm_hypo_motif = motif_selection(torc2_1hr_norm_hypo, 8)
torc2_1hr_norm_hypo_motif2 = torc2_1hr_norm_hypo_motif[torc2_1hr_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc2_1hr_norm_hypo_motif2_highest_seq_len = torc2_1hr_norm_hypo_motif2['seq_length'].value_counts()
torc2_1hr_norm_hypo_motif3 = torc2_1hr_norm_hypo_motif2[torc2_1hr_norm_hypo_motif2['seq_length']==torc2_1hr_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_1hr_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_1hr_norm_hypo_motif_seq.txt")

torc2_1hr_norm_hyper = torc2_1hr_norm[torc2_1hr_norm["log2ratio_Condition2"] > fold_change]
torc2_1hr_norm_hyper_motif = motif_selection(torc2_1hr_norm_hyper, 8)
torc2_1hr_norm_hyper_motif2 = torc2_1hr_norm_hyper_motif[torc2_1hr_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc2_1hr_norm_hyper_motif2_highest_seq_len = torc2_1hr_norm_hyper_motif2['seq_length'].value_counts()
torc2_1hr_norm_hyper_motif3 = torc2_1hr_norm_hyper_motif2[torc2_1hr_norm_hyper_motif2['seq_length'] == torc2_1hr_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_1hr_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_1hr_norm_hyper_motif_seq.txt")

# Torc2 2hr
# Auxin effect removed
torc2_2hr_signi, torc2_2hr_aux_opposite_entry = signi_data(torc2_2hr, '2_hr', afb2_2h_tc, df_var=str(2), ctr_var=str(4))
# Tagging effect removed
torc2_2hr_norm, torc2_2hr_zero_hr_effect = norm_with_zero_hr(torc2_2hr_signi, torc2_zero_hr_signi)
# Saving to tsv format
torc2_2hr_norm.to_csv(path_to_save_torc2_processed_data+"torc2_2hr_norm.tsv", header=True, index=False, sep='\t')
#
torc2_2hr_norm_hypo = torc2_2hr_norm[torc2_2hr_norm["log2ratio_Condition2"] < fold_change]
torc2_2hr_norm_hypo_motif = motif_selection(torc2_2hr_norm_hypo, 8)
torc2_2hr_norm_hypo_motif2 = torc2_2hr_norm_hypo_motif[torc2_2hr_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc2_2hr_norm_hypo_motif2_highest_seq_len = torc2_2hr_norm_hypo_motif2['seq_length'].value_counts()
torc2_2hr_norm_hypo_motif3 = torc2_2hr_norm_hypo_motif2[torc2_2hr_norm_hypo_motif2['seq_length']==torc2_2hr_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_2hr_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_2hr_norm_hypo_motif_seq.txt")

torc2_2hr_norm_hyper = torc2_2hr_norm[torc2_2hr_norm["log2ratio_Condition2"] > fold_change]
torc2_2hr_norm_hyper_motif = motif_selection(torc2_2hr_norm_hyper, 8)
torc2_2hr_norm_hyper_motif2 = torc2_2hr_norm_hyper_motif[torc2_2hr_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc2_2hr_norm_hyper_motif2_highest_seq_len = torc2_2hr_norm_hyper_motif2['seq_length'].value_counts()
torc2_2hr_norm_hyper_motif3 = torc2_2hr_norm_hyper_motif2[torc2_2hr_norm_hyper_motif2['seq_length'] == torc2_2hr_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_2hr_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_2hr_norm_hyper_motif_seq.txt")

# Torc2 3hr
# Auxin effect removed
torc2_3hr_signi, torc2_3hr_aux_opposite_entry = signi_data(torc2_3hr, '3_hr', afb2_3h_tc, df_var=str(2), ctr_var=str(5))
# Tagging effect removed
torc2_3hr_norm, torc2_3hr_zero_hr_effect = norm_with_zero_hr(torc2_3hr_signi, torc2_zero_hr_signi)
# Saving to tsv format
torc2_3hr_norm.to_csv(path_to_save_torc2_processed_data+"torc2_3hr_norm.tsv", header=True, index=False, sep='\t')
#
torc2_3hr_norm_hypo = torc2_3hr_norm[torc2_3hr_norm["log2ratio_Condition2"] < fold_change]
torc2_3hr_norm_hypo_motif = motif_selection(torc2_3hr_norm_hypo, 8)
torc2_3hr_norm_hypo_motif2 = torc2_3hr_norm_hypo_motif[torc2_3hr_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc2_3hr_norm_hypo_motif2_highest_seq_len = torc2_3hr_norm_hypo_motif2['seq_length'].value_counts()
torc2_3hr_norm_hypo_motif3 = torc2_3hr_norm_hypo_motif2[torc2_3hr_norm_hypo_motif2['seq_length']==torc2_3hr_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_3hr_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_3hr_norm_hypo_motif_seq.txt")

torc2_3hr_norm_hyper = torc2_3hr_norm[torc2_3hr_norm["log2ratio_Condition2"] > fold_change]
torc2_3hr_norm_hyper_motif = motif_selection(torc2_3hr_norm_hyper, 8)
torc2_3hr_norm_hyper_motif2 = torc2_3hr_norm_hyper_motif[torc2_3hr_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc2_3hr_norm_hyper_motif2_highest_seq_len = torc2_3hr_norm_hyper_motif2['seq_length'].value_counts()
torc2_3hr_norm_hyper_motif3 = torc2_3hr_norm_hyper_motif2[torc2_3hr_norm_hyper_motif2['seq_length'] == torc2_3hr_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_3hr_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_3hr_norm_hyper_motif_seq.txt")


# Torc2 5hr
# Auxin effect removed
torc2_5hr_signi, torc2_5hr_aux_opposite_entry = signi_data(torc2_5hr, '5_hr', afb2_5h_tc, df_var=str(2), ctr_var=str(6))
# Tagging effect removed
torc2_5hr_norm, torc2_5hr_zero_hr_effect = norm_with_zero_hr(torc2_5hr_signi, torc2_zero_hr_signi)
# Saving to tsv format
torc2_5hr_norm.to_csv(path_to_save_torc2_processed_data+"torc2_5hr_norm.tsv", header=True, index=False, sep='\t')
#
torc2_5hr_norm_hypo = torc2_5hr_norm[torc2_5hr_norm["log2ratio_Condition2"] < fold_change]
torc2_5hr_norm_hypo_motif = motif_selection(torc2_5hr_norm_hypo, 8)
torc2_5hr_norm_hypo_motif2 = torc2_5hr_norm_hypo_motif[torc2_5hr_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc2_5hr_norm_hypo_motif2_highest_seq_len = torc2_5hr_norm_hypo_motif2['seq_length'].value_counts()
torc2_5hr_norm_hypo_motif3 = torc2_5hr_norm_hypo_motif2[torc2_5hr_norm_hypo_motif2['seq_length']==torc2_5hr_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_5hr_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_1hr_norm_hypo_motif_seq.txt")

torc2_5hr_norm_hyper = torc2_5hr_norm[torc2_5hr_norm["log2ratio_Condition2"] > fold_change]
torc2_5hr_norm_hyper_motif = motif_selection(torc2_5hr_norm_hyper, 8)
torc2_5hr_norm_hyper_motif2 = torc2_5hr_norm_hyper_motif[torc2_5hr_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc2_5hr_norm_hyper_motif2_highest_seq_len = torc2_5hr_norm_hyper_motif2['seq_length'].value_counts()
torc2_5hr_norm_hyper_motif3 = torc2_5hr_norm_hyper_motif2[torc2_5hr_norm_hyper_motif2['seq_length'] == torc2_5hr_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc2_5hr_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc2_processed_data+"torc2_1hr_norm_hyper_motif_seq.txt")


# TORC1 dataset

# TORC1 (Kog1) data
path_to_save_torc1_processed_data = r"./re_analysis/phospho_proteome_pair_wise/torc1_pairwise_phospho/"

# Represent the tagging effect (Kog1-AID compared with Afb2 at zero hr i.e. without auxin treatment)
torc1_zero_hr_signi = signi_data(torc1_zero_hr, "zero_hr")
torc1_zero_hr_signi.to_csv(path_to_save_torc1_processed_data+"torc1_zero_hr_norm.tsv", header=True, index=False, sep='\t')


# TORC1 30 min is compared with afb2 30 min and require only zero min control gene
# Auxin effect removed
torc1_30_min_signi, torc1_30_min_aux_opposite_entry = signi_data(torc1_30_min, '30_min', afb2_30min_tc, df_var=str(2), ctr_var=str(2))
# Tagging effect removed
torc1_30_min_norm, torc1_30_min_zero_hr_effect = norm_with_zero_hr(torc1_30_min_signi, torc1_zero_hr_signi)
# Saving to tsv format
torc1_30_min_norm.to_csv(path_to_save_torc1_processed_data+"torc1_30_min_norm.tsv", header=True, index=False, sep='\t')
#
torc1_30_min_norm_hypo = torc1_30_min_norm[torc1_30_min_norm["log2ratio_Condition2"] < fold_change]
torc1_30_min_norm_hypo_motif = motif_selection(torc1_30_min_norm_hypo, 8)
torc1_30_min_norm_hypo_motif2 = torc1_30_min_norm_hypo_motif[torc1_30_min_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc1_30_min_norm_hypo_motif2_highest_seq_len = torc1_30_min_norm_hypo_motif2['seq_length'].value_counts()
torc1_30_min_norm_hypo_motif3 = torc1_30_min_norm_hypo_motif2[torc1_30_min_norm_hypo_motif2['seq_length']==torc1_30_min_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_30_min_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_30min_norm_hypo_motif_seq.txt")

torc1_30_min_norm_hyper = torc1_30_min_norm[torc1_30_min_norm["log2ratio_Condition2"] > fold_change]
torc1_30_min_norm_hyper_motif = motif_selection(torc1_30_min_norm_hyper, 8)
torc1_30_min_norm_hyper_motif2 = torc1_30_min_norm_hyper_motif[torc1_30_min_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc1_30_min_norm_hyper_motif2_highest_seq_len = torc1_30_min_norm_hyper_motif2['seq_length'].value_counts()
torc1_30_min_norm_hyper_motif3 = torc1_30_min_norm_hyper_motif2[torc1_30_min_norm_hyper_motif2['seq_length'] == torc1_30_min_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_30_min_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_30min_norm_hyper_motif_seq.txt")

# TORC1 1 hr
# Auxin effect removed
torc1_1hr_signi, torc1_1hr_aux_opposite_entry = signi_data(torc1_1hr, '1_hr', afb2_1h_tc, df_var=str(2), ctr_var=str(3))
# Tagging effect removed
torc1_1hr_norm, torc1_1hr_zero_hr_effect = norm_with_zero_hr(torc1_1hr_signi, torc1_zero_hr_signi)
# Saving to tsv format
torc1_1hr_norm.to_csv(path_to_save_torc1_processed_data+"torc1_1hr_norm.tsv", header=True, index=False, sep='\t')
#
torc1_1hr_norm_hypo = torc1_1hr_norm[torc1_1hr_norm["log2ratio_Condition2"] < fold_change]
torc1_1hr_norm_hypo_motif = motif_selection(torc1_1hr_norm_hypo, 8)
torc1_1hr_norm_hypo_motif2 = torc1_1hr_norm_hypo_motif[torc1_1hr_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc1_1hr_norm_hypo_motif2_highest_seq_len = torc1_1hr_norm_hypo_motif2['seq_length'].value_counts()
torc1_1hr_norm_hypo_motif3 = torc1_1hr_norm_hypo_motif2[torc1_1hr_norm_hypo_motif2['seq_length']==torc1_1hr_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_1hr_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_1hr_norm_hypo_motif_seq.txt")

torc1_1hr_norm_hyper = torc1_1hr_norm[torc1_1hr_norm["log2ratio_Condition2"] > fold_change]
torc1_1hr_norm_hyper_motif = motif_selection(torc1_1hr_norm_hyper, 8)
torc1_1hr_norm_hyper_motif2 = torc1_1hr_norm_hyper_motif[torc1_1hr_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc1_1hr_norm_hyper_motif2_highest_seq_len = torc1_1hr_norm_hyper_motif2['seq_length'].value_counts()
torc1_1hr_norm_hyper_motif3 = torc1_1hr_norm_hyper_motif2[torc1_1hr_norm_hyper_motif2['seq_length'] == torc1_1hr_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_1hr_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_1hr_norm_hyper_motif_seq.txt")

# TORC1 2 hr
# Auxin effect removed
torc1_2hr_signi, torc1_2hr_aux_opposite_entry = signi_data(torc1_2hr, '2_hr', afb2_2h_tc, df_var=str(2), ctr_var=str(4))
# Tagging effect removed
torc1_2hr_norm, torc1_2hr_zero_hr_effect = norm_with_zero_hr(torc1_2hr_signi, torc1_zero_hr_signi)
# Saving to tsv format
torc1_2hr_norm.to_csv(path_to_save_torc1_processed_data+"torc1_2hr_norm.tsv", header=True, index=False, sep='\t')
#
torc1_2hr_norm_hypo = torc1_2hr_norm[torc1_2hr_norm["log2ratio_Condition2"] < fold_change]
torc1_2hr_norm_hypo_motif = motif_selection(torc1_2hr_norm_hypo, 8)
torc1_2hr_norm_hypo_motif2 = torc1_2hr_norm_hypo_motif[torc1_2hr_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc1_2hr_norm_hypo_motif2_highest_seq_len = torc1_2hr_norm_hypo_motif2['seq_length'].value_counts()
torc1_2hr_norm_hypo_motif3 = torc1_2hr_norm_hypo_motif2[torc1_2hr_norm_hypo_motif2['seq_length']==torc1_2hr_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_2hr_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_2hr_norm_hypo_motif_seq.txt")

torc1_2hr_norm_hyper = torc1_2hr_norm[torc1_2hr_norm["log2ratio_Condition2"] > fold_change]
torc1_2hr_norm_hyper_motif = motif_selection(torc1_2hr_norm_hyper, 8)
torc1_2hr_norm_hyper_motif2 = torc1_2hr_norm_hyper_motif[torc1_2hr_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc1_2hr_norm_hyper_motif2_highest_seq_len = torc1_2hr_norm_hyper_motif2['seq_length'].value_counts()
torc1_2hr_norm_hyper_motif3 = torc1_2hr_norm_hyper_motif2[torc1_2hr_norm_hyper_motif2['seq_length'] == torc1_2hr_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_2hr_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_2hr_norm_hyper_motif_seq.txt")


# TORC1 3 hr
# Auxin effect removed
torc1_3hr_signi, torc1_3hr_aux_opposite_entry = signi_data(torc1_3hr, '3_hr', afb2_3h_tc, df_var=str(2), ctr_var=str(5))
# Tagging effect removed
torc1_3hr_norm, torc1_3hr_zero_hr_effect = norm_with_zero_hr(torc1_3hr_signi, torc1_zero_hr_signi)
# Saving to tsv format
torc1_3hr_norm.to_csv(path_to_save_torc1_processed_data+"torc1_3hr_norm.tsv", header=True, index=False, sep='\t')
#
torc1_3hr_norm_hypo = torc1_3hr_norm[torc1_3hr_norm["log2ratio_Condition2"] < fold_change]
torc1_3hr_norm_hypo_motif = motif_selection(torc1_3hr_norm_hypo, 8)
torc1_3hr_norm_hypo_motif2 = torc1_3hr_norm_hypo_motif[torc1_3hr_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc1_3hr_norm_hypo_motif2_highest_seq_len = torc1_3hr_norm_hypo_motif2['seq_length'].value_counts()
torc1_3hr_norm_hypo_motif3 = torc1_3hr_norm_hypo_motif2[torc1_3hr_norm_hypo_motif2['seq_length']==torc1_3hr_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_3hr_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_3hr_norm_hypo_motif_seq.txt")

torc1_3hr_norm_hyper = torc1_3hr_norm[torc1_3hr_norm["log2ratio_Condition2"] > fold_change]
torc1_3hr_norm_hyper_motif = motif_selection(torc1_3hr_norm_hyper, 8)
torc1_3hr_norm_hyper_motif2 = torc1_3hr_norm_hyper_motif[torc1_3hr_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc1_3hr_norm_hyper_motif2_highest_seq_len = torc1_3hr_norm_hyper_motif2['seq_length'].value_counts()
torc1_3hr_norm_hyper_motif3 = torc1_3hr_norm_hyper_motif2[torc1_3hr_norm_hyper_motif2['seq_length'] == torc1_3hr_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_3hr_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_3hr_norm_hyper_motif_seq.txt")


# TORC1 5 hr
# Auxin effect removed
torc1_5hr_signi, torc1_5hr_aux_opposite_entry = signi_data(torc1_5hr, '5_hr', afb2_5h_tc, df_var=str(2), ctr_var=str(6))
# Tagging effect removed
torc1_5hr_norm, torc1_5hr_zero_hr_effect = norm_with_zero_hr(torc1_5hr_signi, torc1_zero_hr_signi)
# Saving to tsv format
torc1_5hr_norm.to_csv(path_to_save_torc1_processed_data+"torc1_5hr_norm.tsv", header=True, index=False, sep='\t')
#
torc1_5hr_norm_hypo = torc1_5hr_norm[torc1_5hr_norm["log2ratio_Condition2"] < fold_change]
torc1_5hr_norm_hypo_motif = motif_selection(torc1_5hr_norm_hypo, 8)
torc1_5hr_norm_hypo_motif2 = torc1_5hr_norm_hypo_motif[torc1_5hr_norm_hypo_motif['ptm'].str.contains('Phospho')]
torc1_5hr_norm_hypo_motif2_highest_seq_len = torc1_5hr_norm_hypo_motif2['seq_length'].value_counts()
torc1_5hr_norm_hypo_motif3 = torc1_5hr_norm_hypo_motif2[torc1_5hr_norm_hypo_motif2['seq_length']==torc1_5hr_norm_hypo_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_5hr_norm_hypo_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_5hr_norm_hypo_motif_seq.txt")

torc1_5hr_norm_hyper = torc1_5hr_norm[torc1_5hr_norm["log2ratio_Condition2"] > fold_change]
torc1_5hr_norm_hyper_motif = motif_selection(torc1_5hr_norm_hyper, 8)
torc1_5hr_norm_hyper_motif2 = torc1_5hr_norm_hyper_motif[torc1_5hr_norm_hyper_motif['ptm'].str.contains('Phospho')]
torc1_5hr_norm_hyper_motif2_highest_seq_len = torc1_5hr_norm_hyper_motif2['seq_length'].value_counts()
torc1_5hr_norm_hyper_motif3 = torc1_5hr_norm_hyper_motif2[torc1_5hr_norm_hyper_motif2['seq_length'] == torc1_5hr_norm_hyper_motif2_highest_seq_len.index[0]]
df_to_fasta(torc1_5hr_norm_hyper_motif3, 'phospho_site', 'motifX', path_to_save_torc1_processed_data+"torc1_5hr_norm_hyper_motif_seq.txt")
