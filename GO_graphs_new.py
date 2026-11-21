import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
import textwrap

# For TORC1 alltime points GO

# Loading alltime points GO data
torc1_go_file = r'.\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_GO_list_all_time1.xlsx'
torc1_go_data = pd.ExcelFile(torc1_go_file)


def get_go_data (data, sheet_name, time, hits):
    empty_df = pd.DataFrame(columns=['GO biological process complete', 'Fold Enrichment', '-log10(pValue)'])
    for i in range(len(sheet_name)):
        df = pd.read_excel(data, sheet_name=sheet_name[i], header=0)
        df1 = df[['GO biological process complete', 'Fold Enrichment', 'raw P value']]
        df1['-log10(pValue)'] = -np.log10(df1['raw P value'])
        df2 = df1.drop(columns=['raw P value'])
        df3 = df2.head(hits)
        df3.rename(columns={'Fold Enrichment': 'Fold Enrichment'+'_'+time[i], '-log10(pValue)': '-log10(pValue)'+'_'+time[i]}, inplace=True)
        empty_df = empty_df.merge(df3, on='GO biological process complete', how='outer')
    df4 = empty_df.drop(columns=['Fold Enrichment', '-log10(pValue)'])
    df4["NA_count"] = df4.isna().sum(axis=1)
    return df4.sort_values(ascending=True, by='NA_count')


torc1_go_common = get_go_data(torc1_go_file, torc1_go_data.sheet_names, time=['30min', '1hr', '2hr', '3hr', '5hr'], hits=15)

torc1_go_common_fold_enrich = torc1_go_common.melt(id_vars='GO biological process complete',
value_vars=['Fold Enrichment_30min', 'Fold Enrichment_1hr', 'Fold Enrichment_2hr', 'Fold Enrichment_3hr', 'Fold Enrichment_5hr'],
                                                   var_name='Time', value_name="Fold_Enrichment")

torc1_go_common_fold_enrich['Time'].replace(regex={'Fold Enrichment_30min':'30min', 'Fold Enrichment_1hr':'1hr', 'Fold Enrichment_2hr':'2hr',
                                                   'Fold Enrichment_3hr':'3hr', 'Fold Enrichment_5hr':'5hr'}, inplace=True)

torc1_go_common_p_value = torc1_go_common.melt(id_vars='GO biological process complete',
value_vars=['-log10(pValue)_30min', '-log10(pValue)_1hr', '-log10(pValue)_2hr', '-log10(pValue)_3hr', '-log10(pValue)_5hr'],
                                                   var_name='Time', value_name="-log10(pValue)")

torc1_go_common_p_value['Time'].replace(regex={r'-log10\(pValue\)_30min':'30min', r'-log10\(pValue\)_1hr': '1hr', r'-log10\(pValue\)_2hr': '2hr',
                                               r'-log10\(pValue\)_3hr': '3hr', r'-log10\(pValue\)_5hr': '5hr'}, inplace=True)

torc1_go_common1 = pd.merge(torc1_go_common_fold_enrich, torc1_go_common_p_value, on=['GO biological process complete','Time'])
torc1_go_common1.to_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_top15_GO_term.xlsx', header=True, index=False)

torc1_go_manually_curated = pd.read_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_top15_GO_term_mannually_curated.xlsx', header=0)
# plot
plt.figure(figsize=(12, 12))
ax = sns.scatterplot(data=torc1_go_manually_curated, x="Time", y="GO biological process complete", size="Fold_Enrichment", hue="-log10(pValue)",
                sizes=(50, 350), palette="flare", edgecolors="black", legend=True)
# wrap Y-axis tick labels
wrapped_labels = [
    textwrap.fill(label.get_text(), width=60)
    for label in ax.get_yticklabels()
]
ax.set_yticklabels(wrapped_labels)
ax.tick_params(axis='y', which='major', labelsize=8)
sns.move_legend(ax, loc='upper left', bbox_to_anchor=(1, 1))
plt.subplots_adjust(left=0.5)
plt.subplots_adjust(right=0.77)
plt.subplots_adjust(top=0.99)
plt.subplots_adjust(bottom=0.048)
plt.savefig(r'.\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_GO_manual_scattergraph.svg', dpi=600)
plt.close()

# For TORC2
torc2_go_file = r'.\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_GO_pairwise1.xlsx'
torc2_go_data = pd.ExcelFile(torc2_go_file)
torc2_go_common = get_go_data(torc2_go_file, torc2_go_data.sheet_names, time=['30min', '1hr', '2hr', '3hr', '5hr'], hits=15)

torc2_go_common_fold_enrich = torc2_go_common.melt(id_vars='GO biological process complete',
value_vars=['Fold Enrichment_30min', 'Fold Enrichment_1hr', 'Fold Enrichment_2hr', 'Fold Enrichment_3hr', 'Fold Enrichment_5hr'],
                                                   var_name='Time', value_name="Fold_Enrichment")

torc2_go_common_fold_enrich['Time'].replace(regex={'Fold Enrichment_30min':'30min', 'Fold Enrichment_1hr':'1hr', 'Fold Enrichment_2hr':'2hr',
                                                   'Fold Enrichment_3hr':'3hr', 'Fold Enrichment_5hr':'5hr'}, inplace=True)

torc2_go_common_p_value = torc2_go_common.melt(id_vars='GO biological process complete',
value_vars=['-log10(pValue)_30min', '-log10(pValue)_1hr', '-log10(pValue)_2hr', '-log10(pValue)_3hr', '-log10(pValue)_5hr'],
                                                   var_name='Time', value_name="-log10(pValue)")

torc2_go_common_p_value['Time'].replace(regex={r'-log10\(pValue\)_30min':'30min', r'-log10\(pValue\)_1hr': '1hr', r'-log10\(pValue\)_2hr': '2hr',
                                               r'-log10\(pValue\)_3hr': '3hr', r'-log10\(pValue\)_5hr': '5hr'}, inplace=True)

torc2_go_common2 = pd.merge(torc2_go_common_fold_enrich, torc2_go_common_p_value, on=['GO biological process complete','Time'])
torc2_go_common2.to_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_top15_GO_term.xlsx', index=False, header=True)

torc2_go_manually_curated = pd.read_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_top15_GO_term_manually_curated.xlsx', header=0)

plt.figure(figsize=(12,12))
ax = sns.scatterplot(data=torc2_go_manually_curated, x="Time", y="GO biological process complete", size="Fold_Enrichment", hue="-log10(pValue)",
                          sizes=(15, 330), palette="flare", edgecolors="black", legend=True)

# wrap Y-axis tick labels
wrapped_labels = [
    textwrap.fill(label.get_text(), width=60)
    for label in ax.get_yticklabels()
]
ax.set_yticklabels(wrapped_labels)
ax.tick_params(axis='y', which='major', labelsize=8)
sns.move_legend(ax, loc='upper left', bbox_to_anchor=(1, 1))
plt.subplots_adjust(left=0.42)
plt.subplots_adjust(right=0.7)
plt.subplots_adjust(top=0.99)
plt.subplots_adjust(bottom=0.048)
plt.savefig(r'.\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_GO_manual_scattergraph.svg', dpi=600)
plt.close()

# Plotting bar graph for torc1_30min_1hr common GO
torc1_30min_1hr_com_go = pd.read_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_30min_1hr_com_GO.xlsx', header=0)
torc1_30min_1hr_com_go['-log10(pValue)'] = -np.log10(torc1_30min_1hr_com_go['raw P value'])


def wrap_labels(labels, width=60):
    return ["<br>".join(textwrap.wrap(label, width)) for label in labels]


categories = torc1_30min_1hr_com_go['GO biological process complete'].to_list()
wrapped_categories = wrap_labels(categories, width=70)
fig = px.bar(torc1_30min_1hr_com_go, x='Fold Enrichment', y=wrapped_categories,
             color='-log10(pValue)', color_continuous_scale='brwnyl')
# Set background to white
fig.update_layout(
    plot_bgcolor='light grey',  # Background of plotting area
    # paper_bgcolor='grey', # Background of entire figure
    title='TORC1 30 min & 1 hr common GO',
    yaxis_title="GO biological process"
)
fig.update_layout(height=800, width=1200)
fig.write_image(r".\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_30m_1hr_com_GO.png", scale=6)

# TORC1 30 min - 1hr common hyper and hypo graphs
torc1_30min_1hr_com_go_hyper = pd.read_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_30min_1hr_com_GO.xlsx', header=0, sheet_name='hyper')
torc1_30min_1hr_com_go_hyper['-log10(pValue)'] = -np.log10(torc1_30min_1hr_com_go_hyper['raw P value'])
categories = torc1_30min_1hr_com_go_hyper['GO biological process complete'].to_list()
wrapped_categories = wrap_labels(categories, width=70)
fig = px.bar(torc1_30min_1hr_com_go_hyper, x='Fold Enrichment', y=wrapped_categories,
             color='-log10(pValue)', color_continuous_scale='oranges')
# Set background to white
fig.update_layout(
    plot_bgcolor='light grey',  # Background of plotting area
    title='TORC1 30 min & 1 hr common GO Hyper',
    yaxis_title="GO biological process"
)
fig.update_layout(height=800, width=1200)
fig.write_image(r".\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_30m_1hr_com_GO_hyper.png", scale=6)

# Hypo
torc1_30min_1hr_com_go_hypo = pd.read_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_30min_1hr_com_GO.xlsx', header=0, sheet_name='hypo')
torc1_30min_1hr_com_go_hypo['-log10(pValue)'] = -np.log10(torc1_30min_1hr_com_go_hypo['raw P value'])
categories = torc1_30min_1hr_com_go_hypo['GO biological process complete'].to_list()
wrapped_categories = wrap_labels(categories, width=70)
fig = px.bar(torc1_30min_1hr_com_go_hypo, x='Fold Enrichment', y=wrapped_categories,
             color='-log10(pValue)', color_continuous_scale='teal')
# Set background to white
fig.update_layout(
    plot_bgcolor='light grey',  # Background of plotting area
    title='TORC1 30 min & 1 hr common GO Hypo',
    yaxis_title="GO biological process"
)
fig.update_layout(height=800, width=1200)
fig.write_image(r".\re_analysis\phospho_proteome_pair_wise\torc1_pairwise_phospho\torc1_30m_1hr_com_GO_hypo.png", scale=6)



# Plotting bar graph for torc1_30min_1hr common GO

torc2_30min_1hr_com_go = pd.read_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_30min_1hr_com_GO.xlsx', header=0)
torc2_30min_1hr_com_go['-log10(pValue)'] = -np.log10(torc2_30min_1hr_com_go['raw P value'])

# categories = torc2_30min_1hr_com_go['GO biological process complete'].to_list()
# wrapped_categories = wrap_labels(categories, width=70)

fig = px.bar(torc2_30min_1hr_com_go, x='Fold Enrichment', y='GO biological process complete',
             color='-log10(pValue)', color_continuous_scale='brwnyl')
# Set background to white
fig.update_layout(
    plot_bgcolor='light grey',  # Background of plotting area
    # paper_bgcolor='grey', # Background of entire figure
    title='TORC2 30 min & 1 hr common GO',
    yaxis_title="GO biological process"
)
fig.update_layout(height=800, width=1200)
fig.write_image(r".\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_30m_1hr_com_GO.png", scale=6)

# For Hyper
torc2_30min_1hr_com_go_hyper = pd.read_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_30min_1hr_com_GO.xlsx', header=0, sheet_name='hyper')
torc2_30min_1hr_com_go_hyper['-log10(pValue)'] = -np.log10(torc2_30min_1hr_com_go_hyper['raw P value'])

fig = px.bar(torc2_30min_1hr_com_go_hyper, x='Fold Enrichment', y='GO biological process complete',
             color='-log10(pValue)', color_continuous_scale='oranges')
# Set background to white
fig.update_layout(
    plot_bgcolor='light grey',  # Background of plotting area
    # paper_bgcolor='grey', # Background of entire figure
    title='TORC2 30 min & 1 hr common GO Hyper',
    yaxis_title="GO biological process"
)
fig.update_layout(height=800, width=1200)
fig.write_image(r".\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_30m_1hr_com_GO_hyper.png", scale=6)

# For Hypo
torc2_30min_1hr_com_go_hypo = pd.read_excel(r'.\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_30min_1hr_com_GO.xlsx', header=0, sheet_name='hypo')
torc2_30min_1hr_com_go_hypo['-log10(pValue)'] = -np.log10(torc2_30min_1hr_com_go_hypo['raw P value'])

fig = px.bar(torc2_30min_1hr_com_go_hypo, x='Fold Enrichment', y='GO biological process complete',
             color='-log10(pValue)', color_continuous_scale='teal')
# Set background to white
fig.update_layout(
    plot_bgcolor='light grey',  # Background of plotting area
    # paper_bgcolor='grey', # Background of entire figure
    title='TORC2 30 min & 1 hr common GO Hypo',
    yaxis_title="GO biological process"
)
fig.update_layout( width=1200)
fig.write_image(r".\re_analysis\phospho_proteome_pair_wise\torc2_pairwise_phospho\torc2_30m_1hr_com_GO_hypo.png", scale=6)

# YPK2 independent potential TORC2 substrates GO graph

ypk2_independent_torc2_go = pd.read_excel(r'./re_analysis/ypk2_results/ypk2_independent_torc2_GO.xlsx', header=0)
ypk2_independent_torc2_go['-log10(pValue)'] = -np.log10(ypk2_independent_torc2_go['raw P value'])
fig = px.bar(ypk2_independent_torc2_go, x='Fold Enrichment', y='GO biological process complete',
             color='-log10(pValue)', color_continuous_scale='brwnyl')
# Set background to white
fig.update_layout(
    plot_bgcolor='light grey',  # Background of plotting area
    # paper_bgcolor='grey', # Background of entire figure
    title='YPK2 independent TORC2 substrates 30min & 1hr common',
    yaxis_title="GO biological process"
)
fig.update_layout( width=1200)
fig.write_image(r".\re_analysis\ypk2_results/ypk2_independent_torc2_GO.png", scale=6)

# TORC1-TORC2 30m-1hr common GO

torc1_torc2_common = pd.read_excel(r'./re_analysis/phospho_proteome_graphs/TORC1_TORC2_30m_1hr_common_GO_list.xlsx', header=0)
torc1_torc2_common['-log10(pValue)'] = -np.log10(torc1_torc2_common['raw P value'])
fig = px.bar(torc1_torc2_common, x='Fold Enrichment', y='GO biological process complete',
             color='-log10(pValue)', color_continuous_scale='brwnyl')
# Set background to white
fig.update_layout(
    plot_bgcolor='light grey',  # Background of plotting area
    # paper_bgcolor='grey', # Background of entire figure
    title='TORC1 and TORC2 early responders common',
    yaxis_title="GO biological process"
)
fig.update_layout( width=1200)
fig.write_image(r".\re_analysis\phospho_proteome_graphs/TORC1_TORC2_30m_1hr_common_GO.png", scale=6)













# Scatter plot having size and color
test1 = torc2_go_manually_curated.copy()
test1.fillna(0, inplace=True)
fig = px.scatter(test1, x="Time", y="GO biological process complete", size="Fold_Enrichment", color="-log10(pValue)", size_max=15)
fig.show()

torc1_go_common1['-log10(pValue)1'] = torc1_go_common1['-log10(pValue)']**2
torc1_go_common1['Fold_Enrichment S'] = torc1_go_common1['Fold_Enrichment']**2

fig, ax = plt.subplots()
scatter = ax.scatter(data=torc1_go_common1, x="Time", y="GO biological process complete", c="Fold_Enrichment", s="-log10(pValue)1",
                 cmap="flare")

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.2)  # adjust pad for spacing
cbar = plt.colorbar(scatter)
cbar.set_label("Fold Enrichment")


size_values = np.linspace(np.sqrt(torc1_go_common1['-log10(pValue)1'].min()),
                          np.sqrt(torc1_go_common1['-log10(pValue)1'].max()), 3)
size_handles = [
    plt.scatter([], [], s=val**2, color="gray", alpha=0.6)
    for val in size_values
]
size_labels = [f"{val:.0f}" for val in size_values]
# place legend manually outside plot
ax.legend(
    handles=size_handles,
    labels=size_labels,
    title="-log10(pValue)",
    loc='upper left',           # inside the axes, will reposition next
    bbox_to_anchor=(1.05, 0.2), # x=1.05 moves right outside plot, y=0.6 is vertical pos
    borderaxespad=0
)
plt.show()
