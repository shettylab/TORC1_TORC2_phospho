import pandas as pd
import matplotlib.pyplot as plt

# Getting time course proteome control data

prt_tm_ctr = pd.read_csv(r".\TOR_AID_Mass_spec\Ctrl_Avo3_Kog1_Proteome_time_course_Raw\SQ_P317_MH_JH_jh002_timecoures_Control\SQ_P317_MH_JH_jh002_timecoures_Control_PROTEIN.tsv", sep="\t")

prt_tm_ctr_30m = prt_tm_ctr[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition2', 'cv_Condition2', 'log2ratio_Condition2',  'fracNAFeatures.Condition2', 'pValue_Condition2',
                             'FTestPValue', 'FTestQValue']]
prt_tm_ctr_1h = prt_tm_ctr[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition3', 'cv_Condition3', 'log2ratio_Condition3',  'fracNAFeatures.Condition3', 'pValue_Condition3',
                             'FTestPValue', 'FTestQValue']]
prt_tm_ctr_2h = prt_tm_ctr[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition4', 'cv_Condition4', 'log2ratio_Condition4',  'fracNAFeatures.Condition4', 'pValue_Condition4',
                             'FTestPValue', 'FTestQValue']]
prt_tm_ctr_3h = prt_tm_ctr[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition5', 'cv_Condition5', 'log2ratio_Condition5',  'fracNAFeatures.Condition5', 'pValue_Condition5',
                             'FTestPValue', 'FTestQValue']]
prt_tm_ctr_5h = prt_tm_ctr[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition6', 'cv_Condition6', 'log2ratio_Condition6',  'fracNAFeatures.Condition6', 'pValue_Condition6',
                             'FTestPValue', 'FTestQValue']]

pvalue_limit = 0.05
fold_change = 0


def modified_data(df, time_stamp, ctr, var=str()):
    data_1 = df[['ac', 'geneName', 'proteinName', 'nbPeptides', 'medianInt_Condition'+var, 'log2ratio_Condition'+var, 'pValue_Condition'+var,
                 'FTestPValue']]
    # Get only significant entries
    data_signi = data_1[data_1['pValue_Condition'+var] <= pvalue_limit]
    data_signi['time_stamp'] = str(time_stamp)
    # Considering differentially expressed phosphosites
    up_down_prt = data_signi[(data_signi['log2ratio_Condition'+var] > fold_change) |
                               (data_signi['log2ratio_Condition'+var] < -fold_change)]
    if ctr is not None:  # removing auxin effect
        # ctr_data = ctr.rename(columns={'log2ratio_Condition' + var: 'log2ratio_Condition2'})
        df_auxin_data_merge = up_down_prt.merge(ctr, how='inner', on='geneName', suffixes=("_df", "_aux"))
        df_auxin_data_merge['potential_auxin_effect'] = df_auxin_data_merge["log2ratio_Condition"+str(var)+"_df"] * \
                                                        df_auxin_data_merge['log2ratio_Condition'+str(var)+"_aux"]
        aux_effect_prot = df_auxin_data_merge[df_auxin_data_merge['potential_auxin_effect'] > 0]['geneName']
        aux_removed_data = up_down_prt[~(up_down_prt['geneName'].isin(aux_effect_prot))]
        return aux_removed_data
    else:
        return up_down_prt


prt_tm_path = r"./re_analysis/proteome_time_course/"

prt_tm_ctr_30m_signi = modified_data(prt_tm_ctr_30m, "30_min", ctr=None, var=str(2))
prt_tm_ctr_30m_signi.to_csv(prt_tm_path + 'proteome_ctr_time_course/proteome_ctr_30min_signi.tsv', index=False, sep="\t")

prt_tm_ctr_1h_signi = modified_data(prt_tm_ctr_1h, "1_hr", ctr=None, var=str(3))
prt_tm_ctr_1h_signi.to_csv(prt_tm_path + 'proteome_ctr_time_course/proteome_ctr_1h_signi.tsv', index=False, sep="\t")

prt_tm_ctr_2h_signi = modified_data(prt_tm_ctr_2h, "2_hr", ctr=None, var=str(4))
prt_tm_ctr_2h_signi.to_csv(prt_tm_path + 'proteome_ctr_time_course/proteome_ctr_2h_signi.tsv', index=False, sep="\t")

prt_tm_ctr_3h_signi = modified_data(prt_tm_ctr_3h, "3_hr", ctr=None, var=str(5))
prt_tm_ctr_3h_signi.to_csv(prt_tm_path + 'proteome_ctr_time_course/proteome_ctr_3h_signi.tsv', index=False, sep="\t")

prt_tm_ctr_5h_signi = modified_data(prt_tm_ctr_5h, "5_hr", ctr=None, var=str(6))
prt_tm_ctr_5h_signi.to_csv(prt_tm_path + 'proteome_ctr_time_course/proteome_ctr_5h_signi.tsv', index=False, sep="\t")

# TORC1 time course dataset
prt_tm_torc1 = pd.read_csv(r".\TOR_AID_Mass_spec\Ctrl_Avo3_Kog1_Proteome_time_course_Raw\SQ_P317_MH_JH_jh002_timecoures_Kog1\SQ_P317_MH_JH_jh002_timecoures_Kog1_PROTEIN.tsv", sep="\t")

prt_tm_torc1_30m = prt_tm_torc1[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition2', 'cv_Condition2', 'log2ratio_Condition2',  'fracNAFeatures.Condition2', 'pValue_Condition2',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc1_1h = prt_tm_torc1[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition3', 'cv_Condition3', 'log2ratio_Condition3',  'fracNAFeatures.Condition3', 'pValue_Condition3',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc1_2h = prt_tm_torc1[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition4', 'cv_Condition4', 'log2ratio_Condition4',  'fracNAFeatures.Condition4', 'pValue_Condition4',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc1_3h = prt_tm_torc1[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition5', 'cv_Condition5', 'log2ratio_Condition5',  'fracNAFeatures.Condition5', 'pValue_Condition5',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc1_5h = prt_tm_torc1[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition6', 'cv_Condition6', 'log2ratio_Condition6',  'fracNAFeatures.Condition6', 'pValue_Condition6',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc1_30m_signi = modified_data(prt_tm_torc1_30m, "30_min", ctr=prt_tm_ctr_30m_signi, var=str(2))
prt_tm_torc1_30m_signi.to_csv(prt_tm_path + 'proteome_torc1_time_course/proteome_torc1_30min_signi.tsv', index=False, sep="\t")

prt_tm_torc1_1h_signi = modified_data(prt_tm_torc1_1h, "1_hr", ctr=prt_tm_ctr_1h_signi, var=str(3))
prt_tm_torc1_1h_signi.to_csv(prt_tm_path + 'proteome_torc1_time_course/proteome_torc1_1hr_signi.tsv', index=False, sep="\t")

prt_tm_torc1_2h_signi = modified_data(prt_tm_torc1_2h, "2_hr", ctr=prt_tm_ctr_2h_signi, var= str(4))
prt_tm_torc1_2h_signi.to_csv(prt_tm_path + 'proteome_torc1_time_course/proteome_torc1_2hr_signi.tsv', index=False, sep="\t")

prt_tm_torc1_3h_signi = modified_data(prt_tm_torc1_3h, "3_hr", ctr=prt_tm_ctr_3h_signi, var= str(5))
prt_tm_torc1_3h_signi.to_csv(prt_tm_path + 'proteome_torc1_time_course/proteome_torc1_3hr_signi.tsv', index=False, sep="\t")

prt_tm_torc1_5h_signi = modified_data(prt_tm_torc1_5h, "5_hr", ctr=prt_tm_ctr_5h_signi, var= str(6))
prt_tm_torc1_5h_signi.to_csv(prt_tm_path + 'proteome_torc1_time_course/proteome_torc1_5hr_signi.tsv', index=False, sep="\t")


# TORC2 time course dataset
prt_tm_torc2 = pd.read_csv(r".\TOR_AID_Mass_spec\Ctrl_Avo3_Kog1_Proteome_time_course_Raw\SQ_P317_MH_JH_jh002_timecoures_Avo3\SQ_P317_MH_JH_jh002_timecoures_Avo3_PROTEIN.tsv", sep="\t")

prt_tm_torc2_30m = prt_tm_torc2[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition2', 'cv_Condition2', 'log2ratio_Condition2',  'fracNAFeatures.Condition2', 'pValue_Condition2',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc2_1h = prt_tm_torc2[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition3', 'cv_Condition3', 'log2ratio_Condition3',  'fracNAFeatures.Condition3', 'pValue_Condition3',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc2_2h = prt_tm_torc2[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition4', 'cv_Condition4', 'log2ratio_Condition4',  'fracNAFeatures.Condition4', 'pValue_Condition4',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc2_3h = prt_tm_torc2[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition5', 'cv_Condition5', 'log2ratio_Condition5',  'fracNAFeatures.Condition5', 'pValue_Condition5',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc2_5h = prt_tm_torc2[['proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue', 'nbPeptides', 'allAccessions',
                             'medianInt_Condition6', 'cv_Condition6', 'log2ratio_Condition6',  'fracNAFeatures.Condition6', 'pValue_Condition6',
                             'FTestPValue', 'FTestQValue']]

prt_tm_torc2_30m_signi = modified_data(prt_tm_torc2_30m, "30_min", ctr=prt_tm_ctr_30m_signi, var=str(2))
prt_tm_torc2_30m_signi.to_csv(prt_tm_path + 'proteome_torc2_time_course/proteome_torc2_30min_signi.tsv', index=False, sep="\t")

prt_tm_torc2_1h_signi = modified_data(prt_tm_torc2_1h, "1_hr", ctr=prt_tm_ctr_1h_signi, var=str(3))
prt_tm_torc2_1h_signi.to_csv(prt_tm_path + 'proteome_torc2_time_course/proteome_torc2_1hr_signi.tsv', index=False, sep="\t")

prt_tm_torc2_2h_signi = modified_data(prt_tm_torc2_2h, "2_hr", ctr=prt_tm_ctr_2h_signi, var=str(4))
prt_tm_torc2_2h_signi.to_csv(prt_tm_path + 'proteome_torc2_time_course/proteome_torc2_2hr_signi.tsv', index=False, sep="\t")

prt_tm_torc2_3h_signi = modified_data(prt_tm_torc2_3h, "3_hr", ctr=prt_tm_ctr_3h_signi, var=str(5))
prt_tm_torc2_3h_signi.to_csv(prt_tm_path + 'proteome_torc2_time_course/proteome_torc2_3hr_signi.tsv', index=False, sep="\t")

prt_tm_torc2_5h_signi = modified_data(prt_tm_torc2_5h, "5_hr", ctr=prt_tm_ctr_5h_signi, var=str(6))
prt_tm_torc2_5h_signi.to_csv(prt_tm_path + 'proteome_torc2_time_course/proteome_torc2_5hr_signi.tsv', index=False, sep="\t")

