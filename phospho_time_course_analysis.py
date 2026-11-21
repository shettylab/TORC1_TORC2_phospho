"""
To get time point-wise separation of data from avo3 (TORC2), kog1(TORC1), and control (afb2) dataset.
Get significant entries from each section
"""

import pandas as pd
from matplotlib import pyplot as plt
from upsetplot import generate_counts, plot, UpSet, from_contents

phospho_time_course_afb2 = pd.read_csv(
    r".\TOR_AID_Mass_spec\Ctrl_Avo3_Kog1_Phosphoproteome_time_course_Raw\SQ_P317_MH_JH_jh001_timecoures_C\SQ_P317_MH_JH_jh001_timecoures_C_PEPTIDE.tsv",
    sep="\t")

# Separating each time point data from Afb2 (control) data

afb2_zero_hr = phospho_time_course_afb2[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '001_jh-001_PO4r_D21', '019_jh-001_PO4_D21', '037_jh-001_PO4_D21', '055_jh-001_PO4_D21',
     'medianInt_Condition1', 'cv_Condition1', 'fracNAFeatures.Condition1']]

afb2_30_min = phospho_time_course_afb2[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '002_jh-001_PO4r_D21', '020_jh-001_PO4_D21', '038_jh-001_PO4_D21', '056_jh-001_PO4_D21',
     'medianInt_Condition2', 'cv_Condition2', 'log2ratio_Condition2', 'fracNAFeatures.Condition2',
     'pValue_Condition2', 'FTestPValue', 'FTestQValue']]

afb2_one_hr = phospho_time_course_afb2[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '003_jh-001_PO4r_D21', '021_jh-001_PO4_D21', '039_jh-001_PO4_D21', '057_jh-001_PO4_D21',
     'medianInt_Condition3', 'cv_Condition3', 'log2ratio_Condition3', 'fracNAFeatures.Condition3',
     'pValue_Condition3', 'FTestPValue', 'FTestQValue']]

afb2_two_hr = phospho_time_course_afb2[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '004_jh-001_PO4r_D21', '022_jh-001_PO4_D21', '040_jh-001_PO4_D21', '058_jh-001_PO4_D21',
     'medianInt_Condition4', 'cv_Condition4', 'log2ratio_Condition4', 'fracNAFeatures.Condition4',
     'pValue_Condition4', 'FTestPValue', 'FTestQValue']]

afb2_three_hr = phospho_time_course_afb2[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '005_jh-001_PO4r_D21', '023_jh-001_PO4_D21', '041_jh-001_PO4_D21', '059_jh-001_PO4_D21',
     'medianInt_Condition5', 'cv_Condition5', 'log2ratio_Condition5', 'fracNAFeatures.Condition5',
     'pValue_Condition5', 'FTestPValue', 'FTestQValue']]

afb2_five_hr = phospho_time_course_afb2[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '006_jh-001_PO4r_D21', '024_jh-001_PO4_D21', '042_jh-001_PO4_D21', '060_jh-001_PO4_D21',
     'medianInt_Condition6', 'cv_Condition6', 'log2ratio_Condition6', 'fracNAFeatures.Condition6',
     'pValue_Condition6', 'FTestPValue', 'FTestQValue']]

# Avo3 time course data

phospho_time_course_avo = pd.read_csv(
    r".\TOR_AID_Mass_spec\Ctrl_Avo3_Kog1_Phosphoproteome_time_course_Raw\SQ_P317_MH_JH_jh001_timecoures_Avo3\SQ_P317_MH_JH_jh001_timecoures_Avo3_PEPTIDE.tsv",
    sep="\t")

# Separating each time point data from Avo3 data

avo_zero_hr = phospho_time_course_avo[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '013_jh-001_PO4_D21', '031_jh-001_PO4_D21', '049_jh-001_PO4_D21', '067_jh-001_PO4_D21',
     'medianInt_Condition1', 'cv_Condition1', 'fracNAFeatures.Condition1']]

avo_30_min = phospho_time_course_avo[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '014_jh-001_PO4_D21', '032_jh-001_PO4_D21', '050_jh-001_PO4_D21', '068_jh-001_PO4_D21',
     'medianInt_Condition2', 'cv_Condition2', 'log2ratio_Condition2', 'fracNAFeatures.Condition2',
     'pValue_Condition2', 'FTestPValue', 'FTestQValue']]

avo_one_hr = phospho_time_course_avo[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '015_jh-001_PO4_D21', '033_jh-001_PO4_D21', '051_jh-001_PO4_D21', '069_jh-001_PO4_D21',
     'medianInt_Condition3', 'cv_Condition3', 'log2ratio_Condition3', 'fracNAFeatures.Condition3',
     'pValue_Condition3', 'FTestPValue', 'FTestQValue']]

avo_two_hr = phospho_time_course_avo[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '016_jh-001_PO4_D21', '034_jh-001_PO4_D21', '052_jh-001_PO4_D21', '070_jh-001_PO4_D21',
     'medianInt_Condition4', 'cv_Condition4', 'log2ratio_Condition4', 'fracNAFeatures.Condition4',
     'pValue_Condition4', 'FTestPValue', 'FTestQValue']]

avo_three_hr = phospho_time_course_avo[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '017_jh-001_PO4_D21', '035_jh-001_PO4_D21', '053_jh-001_PO4_D21', '071_jh-001_PO4_D21',
     'medianInt_Condition5', 'cv_Condition5', 'log2ratio_Condition5', 'fracNAFeatures.Condition5',
     'pValue_Condition5', 'FTestPValue', 'FTestQValue']]

avo_five_hr = phospho_time_course_avo[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '018_jh-001_PO4_D21', '036_jh-001_PO4_D21', '054_jh-001_PO4_D21', '072_jh-001_PO4_D21',
     'medianInt_Condition6', 'cv_Condition6', 'log2ratio_Condition6', 'fracNAFeatures.Condition6',
     'pValue_Condition6', 'FTestPValue', 'FTestQValue']]

# Kog1 time course data

phospho_time_course_kog = pd.read_csv(
    r".\TOR_AID_Mass_spec\Ctrl_Avo3_Kog1_Phosphoproteome_time_course_Raw\SQ_P317_MH_JH_jh001_timecoures_kog1\SQ_P317_MH_JH_jh001_timecoures_Kog1_PEPTIDE.tsv",
    sep="\t")

# Separating each time point data from Kog1 data

kog_zero_hr = phospho_time_course_kog[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '007_jh-001_PO4r_D21', '025_jh-001_PO4_D21', '043_jh-001_PO4_D21', '061_jh-001_PO4_D21',
     'medianInt_Condition1', 'cv_Condition1', 'fracNAFeatures.Condition1']]

kog_30_min = phospho_time_course_kog[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '008_jh-001_PO4r_D21', '026_jh-001_PO4_D21', '044_jh-001_PO4_D21', '062_jh-001_PO4_D21',
     'medianInt_Condition2', 'cv_Condition2', 'log2ratio_Condition2', 'fracNAFeatures.Condition2',
     'pValue_Condition2', 'FTestPValue', 'FTestQValue']]

kog_one_hr = phospho_time_course_kog[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '009_jh-001_PO4r_D21', '027_jh-001_PO4_D21', '045_jh-001_PO4_D21', '063_jh-001_PO4_D21',
     'medianInt_Condition3', 'cv_Condition3', 'log2ratio_Condition3', 'fracNAFeatures.Condition3',
     'pValue_Condition3', 'FTestPValue', 'FTestQValue']]

kog_two_hr = phospho_time_course_kog[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '010_jh-001_PO4r_D21', '028_jh-001_PO4_D21', '046_jh-001_PO4_D21', '064_jh-001_PO4_D21',
     'medianInt_Condition4', 'cv_Condition4', 'log2ratio_Condition4', 'fracNAFeatures.Condition4',
     'pValue_Condition4', 'FTestPValue', 'FTestQValue']]

kog_three_hr = phospho_time_course_kog[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '011_jh-001_PO4r_D21', '029_jh-001_PO4_D21', '047_jh-001_PO4_D21', '065_jh-001_PO4_D21',
     'medianInt_Condition5', 'cv_Condition5', 'log2ratio_Condition5', 'fracNAFeatures.Condition5',
     'pValue_Condition5', 'FTestPValue', 'FTestQValue']]

kog_five_hr = phospho_time_course_kog[
    ['peptide', 'proteinName', 'ac', 'geneName', 'proteinDescription', 'idScore', 'idQValue',
     'retentionTime', 'ptm', 'nbPtmsPerPeptide', 'motifX', 'modifCoord', 'allAccessions',
     '012_jh-001_PO4r_D21', '030_jh-001_PO4_D21', '048_jh-001_PO4_D21', '066_jh-001_PO4_D21',
     'medianInt_Condition6', 'cv_Condition6', 'log2ratio_Condition6', 'fracNAFeatures.Condition6',
     'pValue_Condition6', 'FTestPValue', 'FTestQValue']]

fold_change = 0

"""
Function to remove the redundant phosphosites and consider the only with highest intensity.
Remove entries where no phosphosites obtained. 
consider the only significant entries (p value <=0.05) and considering differentially expressed phosphosites. 
If control dataset given then remove similarly altered phosphosite from test dataset
"""


def signi_data(df, time_stamp, ctr=None, var=str()):
    data_1 = df[['ac', 'geneName', 'modifCoord', 'motifX', 'log2ratio_Condition' + var, 'pValue_Condition' + var,
                 'medianInt_Condition' + var]]
    # sorting data based on median intensity of condition2
    data_1.sort_values('medianInt_Condition' + var, ascending=False, ignore_index=True, inplace=True)
    # Removing NA values from coordinates values
    data_1.dropna(subset=['modifCoord'], inplace=True)
    # Converting double entries into single
    data_1['modifCoord'] = data_1['modifCoord'].astype(str)
    # Creating new accession ID by combing Gene name and AA coordinates
    data_1['phospho_site'] = data_1['geneName'] + '_' + data_1['modifCoord']
    data_1['time_stamp'] = str(time_stamp)
    # Get only significant entries
    data_signi = data_1[data_1['pValue_Condition' + var] <= 0.05]
    # dropping duplicated phosphosites
    data_signi_1 = data_signi.drop_duplicates('phospho_site', keep='first')
    # Considering differentially expressed phosphosites
    up_down_prt = data_signi_1[(data_signi_1['log2ratio_Condition' + var] > fold_change) |
                               (data_signi_1['log2ratio_Condition' + var] < fold_change)]
    if ctr is not None:  # removing auxin effect
        # ctr_data = ctr.rename(columns={'log2ratio_Condition' + var: 'log2ratio_Condition2'})
        df_auxin_data_merge = up_down_prt.merge(ctr, how='inner', on='phospho_site', suffixes=("_df", "_aux"))
        df_auxin_data_merge['potential_auxin_effect'] = df_auxin_data_merge["log2ratio_Condition"+str(var)+"_df"] * \
                                                        df_auxin_data_merge['log2ratio_Condition'+str(var)+"_aux"]
        aux_effect_prot = df_auxin_data_merge[df_auxin_data_merge['potential_auxin_effect'] > 0]['phospho_site']
        aux_removed_data = up_down_prt[~(up_down_prt['phospho_site'].isin(aux_effect_prot))]
        return aux_removed_data
    else:
        return up_down_prt


# Afb2 significant data
path_to_store_signi_afb2_data = r"./re_analysis/phospho_proteome_time_course/afb2_time_course_phospho/"

afb2_30min_signi = signi_data(afb2_30_min, "30_min", ctr=None, var=str(2))
afb2_30min_signi.to_csv(path_to_store_signi_afb2_data + "afb2_30min_signi_time_course.tsv", header=True, index=False, sep="\t")

afb2_one_hr_signi = signi_data(afb2_one_hr, "1_hr", ctr=None, var=str(3))
afb2_one_hr_signi.to_csv(path_to_store_signi_afb2_data + "afb2_1hr_signi_time_course.tsv",  header=True, index=False, sep="\t")

afb2_two_hr_signi = signi_data(afb2_two_hr, "2_hr", ctr=None, var=str(4))
afb2_two_hr_signi.to_csv(path_to_store_signi_afb2_data + "afb2_2hr_signi_time_course.tsv",  header=True, index=False, sep="\t")

afb2_three_hr_signi = signi_data(afb2_three_hr, "3_hr", ctr=None, var=str(5))
afb2_three_hr_signi.to_csv(path_to_store_signi_afb2_data + "afb2_3hr_signi_time_course.tsv",  header=True, index=False, sep="\t")

afb2_five_hr_signi = signi_data(afb2_five_hr, "5_hr", ctr=None, var=str(6))
afb2_five_hr_signi.to_csv(path_to_store_signi_afb2_data + "afb2_5hr_signi_time_course.tsv",  header=True, index=False, sep="\t")

# Avo3 normalize significant data (from this auxin effect removed)

path_to_store_signi_avo_data = r"./re_analysis/phospho_proteome_time_course/torc2_time_course_phospho/"

avo_30_min_signi = signi_data(avo_30_min, "30_min", ctr=afb2_30min_signi, var=str(2))
avo_30_min_signi.to_csv(path_to_store_signi_avo_data + "TORC2_30_min_signi_time_course.tsv", header=True, index=False, sep="\t")

avo_1_hr_signi = signi_data(avo_one_hr, "1_hr", ctr=afb2_one_hr_signi, var=str(3))
avo_1_hr_signi.to_csv(path_to_store_signi_avo_data + "TORC2_1_hr_signi_time_course.tsv", header=True, index=False, sep="\t")

avo_2_hr_signi = signi_data(avo_two_hr, "2_hr", ctr=afb2_two_hr_signi, var=str(4))
avo_2_hr_signi.to_csv(path_to_store_signi_avo_data + "TORC2_2_hr_signi_time_course.tsv", header=True, index=False, sep="\t")

avo_3_hr_signi = signi_data(avo_three_hr, "3_hr", ctr=afb2_three_hr_signi, var=str(5))
avo_3_hr_signi.to_csv(path_to_store_signi_avo_data + "TORC2_3_hr_signi_time_course.tsv", header=True, index=False, sep="\t")

avo_5_hr_signi = signi_data(avo_five_hr, "5_hr", ctr=afb2_five_hr_signi, var=str(6))
avo_5_hr_signi.to_csv(path_to_store_signi_avo_data + "TORC2_5_hr_signi_time_course.tsv", header=True, index=False, sep="\t")

# Kog normalize significant data (from this auxin effect removed)

path_to_store_signi_kog_data = r"./re_analysis/phospho_proteome_time_course/torc1_time_course_phospho/"

kog_30_min_signi = signi_data(kog_30_min, "30_min", ctr=afb2_30min_signi, var=str(2))
kog_30_min_signi.to_csv(path_to_store_signi_kog_data + "TORC1_30_min_signi_time_course.tsv", header=True, index=False, sep="\t")

kog_1_hr_signi = signi_data(kog_one_hr, "1_hr", ctr=afb2_one_hr_signi, var=str(3))
kog_1_hr_signi.to_csv(path_to_store_signi_kog_data + "TORC1_1_hr_signi_time_course.tsv", header=True, index=False, sep="\t")

kog_2_hr_signi = signi_data(kog_two_hr, "2_hr", ctr=afb2_two_hr_signi, var=str(4))
kog_2_hr_signi.to_csv(path_to_store_signi_kog_data + "TORC1_2_hr_signi_time_course.tsv", header=True, index=False, sep="\t")

kog_3_hr_signi = signi_data(kog_three_hr, "3_hr", ctr=afb2_three_hr_signi, var=str(5))
kog_3_hr_signi.to_csv(path_to_store_signi_kog_data + "TORC1_3_hr_signi_time_course.tsv", header=True, index=False, sep="\t")

kog_5_hr_signi = signi_data(kog_five_hr, "5_hr", ctr=afb2_five_hr_signi, var=str(6))
kog_5_hr_signi.to_csv(path_to_store_signi_kog_data + "TORC1_5_hr_signi_time_course.tsv", header=True, index=False, sep="\t")

# Upset plots

# Upset plot for Avo3 data all deregulated
avo3_dergulated_phospho = from_contents({"0.5hr": avo_30_min_signi['phospho_site'],
                                        "1hr": avo_1_hr_signi['phospho_site'],
                                        "2hr": avo_2_hr_signi["phospho_site"],
                                        "3hr": avo_3_hr_signi["phospho_site"],
                                        "5hr": avo_5_hr_signi["phospho_site"]})
avo3_dergulated_phospho_table = avo3_dergulated_phospho.sort_values('id').reset_index()
avo3_dergulated_phospho_table.to_excel(path_to_store_signi_avo_data + "TORC2_deregulated_upset_table.xlsx", header=True,
                                       index=False)
avo3_deregulated_plot = UpSet(avo3_dergulated_phospho, subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
avo3_deregulated_plot.plot()
plt.suptitle("TORC2 time course deregulated phospho-peptides")
plt.savefig(path_to_store_signi_avo_data + "TORC2_time_course_deregulated_upset_plot.tif", dpi=600)
plt.close()

# Upset plot for Avo3 hypo phosphorylated

avo3_hypo_deregulated_phospho = from_contents({
    "0.5hr_hypo": avo_30_min_signi[avo_30_min_signi['log2ratio_Condition2'] < 0]['phospho_site'],
    "1hr_hypo": avo_1_hr_signi[avo_1_hr_signi['log2ratio_Condition3'] < 0]['phospho_site'],
    "2hr_hypo": avo_2_hr_signi[avo_2_hr_signi['log2ratio_Condition4'] < 0]['phospho_site'],
    "3hr_hypo": avo_3_hr_signi[avo_3_hr_signi['log2ratio_Condition5'] < 0]['phospho_site'],
    "5hr_hypo": avo_5_hr_signi[avo_5_hr_signi['log2ratio_Condition6'] < 0]['phospho_site']
})
avo3_hypo_deregulated_phospho_table = avo3_hypo_deregulated_phospho.sort_values('id').reset_index()
avo3_hypo_deregulated_phospho_table.to_excel(path_to_store_signi_avo_data + "TORC2_hypo_time_course_upset_table.xlsx",
                                             header=True, index=False)
avo3_hypo_deregulated_phospho_plot = UpSet(avo3_hypo_deregulated_phospho,  subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
avo3_hypo_deregulated_phospho_plot.plot()
plt.suptitle("TORC2 time course hypo phospho-peptides")
plt.savefig(path_to_store_signi_avo_data + "TORC2_time_course_hypo_upset_plot.tif", dpi=600)
plt.close()

# Upset plot for Avo3 hypo phosphorylated

avo3_hyper_deregulated_phospho = from_contents({
    "0.5hr_hyper": avo_30_min_signi[avo_30_min_signi['log2ratio_Condition2'] > 0]['phospho_site'],
    "1hr_hyper": avo_1_hr_signi[avo_1_hr_signi['log2ratio_Condition3'] > 0]['phospho_site'],
    "2hr_hyper": avo_2_hr_signi[avo_2_hr_signi['log2ratio_Condition4'] > 0]['phospho_site'],
    "3hr_hyper": avo_3_hr_signi[avo_3_hr_signi['log2ratio_Condition5'] > 0]['phospho_site'],
    "5hr_hyper": avo_5_hr_signi[avo_5_hr_signi['log2ratio_Condition6'] > 0]['phospho_site']
})
avo3_hyper_deregulated_phospho_table = avo3_hyper_deregulated_phospho.sort_values('id').reset_index()
avo3_hyper_deregulated_phospho_table.to_excel(path_to_store_signi_avo_data + "TORC2_hyper_time_course_upset_table.xlsx",
                                             header=True, index=False)
avo3_hyper_deregulated_phospho_plot = UpSet(avo3_hyper_deregulated_phospho,  subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
avo3_hyper_deregulated_phospho_plot.plot()
plt.suptitle("TORC2 time course hyper phospho-peptides")
plt.savefig(path_to_store_signi_avo_data + "TORC2_time_course_hyper_upset_plot.tif", dpi=600)
plt.close()

# Upset plot for Kog1 data all deregulated
kog1_deregulated_phospho = from_contents({"0.5hr": kog_30_min_signi['phospho_site'],
                                        "1hr": kog_1_hr_signi['phospho_site'],
                                        "2hr": kog_2_hr_signi["phospho_site"],
                                        "3hr": kog_3_hr_signi["phospho_site"],
                                        "5hr": kog_5_hr_signi["phospho_site"]})
kog1_deregulated_phospho_table = kog1_deregulated_phospho.sort_values('id').reset_index()
kog1_deregulated_phospho_table.to_excel(path_to_store_signi_kog_data + "TORC1_deregulated_upset_table.xlsx", header=True,
                                        index=False)
kog1_deregulated_plot = UpSet(kog1_deregulated_phospho, subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
kog1_deregulated_plot.plot()
plt.suptitle("TORC1 time course deregulated phospho-peptides")
plt.savefig(path_to_store_signi_kog_data + "TORC1_time_course_deregulated_upset_plot", dpi=600)
plt.close()

# Upset plot for Kog1 hypo phosphorylated

kog1_hypo_deregulated_phospho = from_contents({
    "0.5hr_hypo": kog_30_min_signi[kog_30_min_signi['log2ratio_Condition2'] < 0]['phospho_site'],
    "1hr_hypo": kog_1_hr_signi[kog_1_hr_signi['log2ratio_Condition3'] < 0]['phospho_site'],
    "2hr_hypo": kog_2_hr_signi[kog_2_hr_signi['log2ratio_Condition4'] < 0]['phospho_site'],
    "3hr_hypo": kog_3_hr_signi[kog_3_hr_signi['log2ratio_Condition5'] < 0]['phospho_site'],
    "5hr_hypo": kog_5_hr_signi[kog_5_hr_signi['log2ratio_Condition6'] < 0]['phospho_site']
})
kog1_hypo_deregulated_phospho_table = kog1_hypo_deregulated_phospho.sort_values('id').reset_index()
kog1_hypo_deregulated_phospho_table.to_excel(path_to_store_signi_kog_data + "TORC1_hypo_time_course_upset_table.xlsx",
                                             header=True, index=False)
kog1_hypo_deregulated_phospho_plot = UpSet(kog1_hypo_deregulated_phospho,  subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
kog1_hypo_deregulated_phospho_plot.plot()
plt.suptitle("TORC1 time course hypo phospho-peptides")
plt.savefig(path_to_store_signi_kog_data + "TORC1_time_course_hypo_upset_plot.tif", dpi=600)
plt.close()

# Upset plot for KOg1 hypo phosphorylated

kog1_hyper_deregulated_phospho = from_contents({
    "0.5hr_hyper": kog_30_min_signi[kog_30_min_signi['log2ratio_Condition2'] > 0]['phospho_site'],
    "1hr_hyper": kog_1_hr_signi[kog_1_hr_signi['log2ratio_Condition3'] > 0]['phospho_site'],
    "2hr_hyper": kog_2_hr_signi[kog_2_hr_signi['log2ratio_Condition4'] > 0]['phospho_site'],
    "3hr_hyper": kog_3_hr_signi[kog_3_hr_signi['log2ratio_Condition5'] > 0]['phospho_site'],
    "5hr_hyper": kog_5_hr_signi[kog_5_hr_signi['log2ratio_Condition6'] > 0]['phospho_site']
})
kog1_hyper_deregulated_phospho_table = kog1_hyper_deregulated_phospho.sort_values('id').reset_index()
kog1_hyper_deregulated_phospho_table.to_excel(path_to_store_signi_kog_data + "TORC1_hyper_time_course_upset_table.xlsx",
                                             header=True, index=False)
kog1_hyper_deregulated_phospho_plot = UpSet(kog1_hyper_deregulated_phospho,  subset_size="count", show_counts=True, sort_by='degree', sort_categories_by='-input')
kog1_hyper_deregulated_phospho_plot.plot()
plt.suptitle("TORC1 time course hyper phospho-peptides")
plt.savefig(path_to_store_signi_kog_data + "TORC1_time_course_hyper_upset_plot", dpi=600)
plt.close()

