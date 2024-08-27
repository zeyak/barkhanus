import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read CSV
df = pd.read_csv("/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/diplo_virulence.blastp", sep="\t",
                 header=None)
out_file = "/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/virulence_heatmap_grouped.png"
out_file2 = "/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/virulence_heatmap_grouped.xlsx"

# Make a sp column with the prefix of the ids in the second column and add it next as a separate column
df['sp'] = df[1].str.split("_").str[0]
# Replace "contig" with "bark"
df["sp"] = df["sp"].str.replace("contig", "bark")

# Group by virulence gene and species, and count the number of hits
df3 = df.groupby([0, 'sp']).size().reset_index(name='count')

# Make the first column the index
df3.set_index(0, inplace=True)
# Make a pivot table
df3 = df3.pivot(columns='sp', values='count').fillna(0)

# Rename the columns
df3 = df3.rename(columns={
    'HIN': "H. inflata",
    'TPC1': "Trepomonas pc1",
    'bark': "S. barkhanus",
    'SS50377': "S. salmonicida",
    'GL50803': "G. intestinalis",
    'GMRT': "G. muris"
})

# Reorder the columns as needed
df3 = df3[["H. inflata", "Trepomonas pc1", "S. salmonicida", "S. barkhanus", "G. intestinalis", "G. muris"]]


# Define the gene categories with lists of genes
gene_categories = {
    'Surface': ['GL50803_113797', 'GL50803_115066', 'GL50803_10330', 'GL50803_103887', 'GL50803_2902'],
    'Alpha-giardin': ['GL50803_11654', 'GL50803_11683', 'GL50803_17153', 'GL50803_15097', 'GL50803_15101'],
    'Disc': ['GL50803_4812', 'GL50803_86676', 'GL50803_17230', 'GL50803_4410'],
    'Cysteine': ['GL50803_14019', 'GL50803_16160', 'GL50803_16779'],
    'Enzymes': ['GL50803_11118', 'GL50803_112103', 'GL50803_10311'],
    'Hydrogenosome': ['SS50377_23299', 'SS50377_25409', 'SS50377_27457', 'SS50377_28003'],
    'Encystation': ['GL50803_5638', 'GL50803_5435', 'GL50803_2421', 'GL50803_8245', 'GL50803_14259', 'GL50803_16069'],
    'Meiosis': ['GL50803_15279', 'GL50803_4084', 'GL50803_13104', 'GL50803_6626']
}

# Invert the dictionary to map each gene to its category
gene_to_category = {gene: category for category, genes in gene_categories.items() for gene in genes}

# Map the genes to the categories using the inverted dictionary
df3['Category'] = df3.index.map(gene_to_category)

# Sort by 'Category' to organize genes by their categories
df3 = df3.sort_values(by='Category')
# Save to Excel
df3.to_excel(out_file2)

# Drop the 'Category' column before plotting
df3_numeric = df3.drop(columns=['Category'])

# Replace 0 values with NaN to leave them blank in the heatmap
df3_numeric = df3_numeric.replace(0, np.nan)


# Plot the heatmap with each gene under its respective category
sns.set_theme(style="white")

f, ax = plt.subplots(figsize=(10, 15))  # Adjust the figsize to accommodate the number of genes

# Plot the heatmap and adjust cbar location
heatmap = sns.heatmap(df3_numeric,
                      cmap=sns.color_palette("gray_r", as_cmap=True),
                      square=False,  # Set square=False to handle uneven row lengths better
                      fmt='g',
                      linewidths=.4,
                      annot=True,
                      cbar_kws={"shrink": .5, "location": "right", "pad": 0.15}  # Adjust the pad to move the colorbar
                      )
heatmap.set_ylabel('')  # Removes the y-axis title
heatmap.set_xlabel('')

# Group the heatmap by categories visually by adding horizontal lines and labels
for category in df3['Category'].unique():
    first_gene_idx = df3.index.get_loc(df3[df3['Category'] == category].index[0])
    last_gene_idx = df3.index.get_loc(df3[df3['Category'] == category].index[-1])
    ax.hlines(last_gene_idx + 1, *ax.get_xlim(), color='black', linewidth=2)

    # Adjust the position of the category labels (closer to the heatmap)
    ax.text(len(df3.columns) - 0.7, (first_gene_idx + last_gene_idx) / 2, category,
            verticalalignment='center', horizontalalignment='center',
            rotation=90, fontsize=12, fontweight='bold')

    # Adjust the position of the color bar (further away from the heatmap)
    cbar_kws = {"shrink": .5, "location": "right", "pad": 0.2}

# Adjust the plot to make space for the labels
plt.subplots_adjust(right=0.85)  # Adjusting the right parameter to move the colorbar further away

# Save the plot
plt.savefig(out_file, format="png", bbox_inches='tight', dpi=1200)
plt.show()
