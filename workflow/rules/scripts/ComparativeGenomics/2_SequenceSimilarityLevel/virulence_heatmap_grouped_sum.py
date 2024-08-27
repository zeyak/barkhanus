import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#read csv
df = pd.read_csv("/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/diplo_virulence.blastp", sep="\t", header=None)
out_file = "/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/virulence_heatmap_grouped_sum.png"
out_file2 = "/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/virulence_heatmap_grouped_sum.xlsx"

#make a sp column with the prefix of the ids in teh secodn column and add it next as a sperate column
df['sp'] = df[1].str.split("_").str[0]
#replae contig with bark
df["sp"] = df["sp"].str.replace("contig", "bark")


# ref virluence genes
df1 = pd.DataFrame(df[0].unique())
#diplo virluence hits
df2 = pd.DataFrame(df[1].unique())


#virulence gene count table
df2['count'] = df[1].value_counts().sort_index().values

# groupby virulence gene and species and count the number of hits
df3 = df.groupby([0, 'sp']).size().reset_index(name='count')

#make firsr column the index
df3.set_index(0, inplace=True)
#make a pivot table
df3 = df3.pivot(columns='sp', values='count').fillna(0)
#renaem the columns
df3 = df3.rename(columns={ 'HIN': "H. inflata",
                'TPC1': "Trepomonas pc1",
                'bark': "S. barkhanus",
                'SS50377': "S. salmonicida",
                'GL50803': "G. intestinalis",
               'GMRT': "G. muris"})

df3 = df3[["H. inflata", "Trepomonas pc1" ,"S. salmonicida","S. barkhanus", "G. intestinalis",
             "G. muris"]]
#to excel file
df3.to_excel(out_file2)

# Define the gene categories with lists of genes
gene_categories = {
    'Surface': [
        'GL50803_113797', 'GL50803_115066', 'GL50803_10330',
        'GL50803_103887', 'GL50803_2902'
    ],
    'Alpha-giardins': [
        'GL50803_11654', 'GL50803_11683', 'GL50803_17153',
        'GL50803_15097', 'GL50803_15101'
    ],
    'Disc proteins': [
        'GL50803_4812', 'GL50803_86676', 'GL50803_17230',
        'GL50803_4410'
    ],
    'Cysteine proteases': [
        'GL50803_14019', 'GL50803_16160', 'GL50803_16779'
    ],
    'Enzymes': [
        'GL50803_11118', 'GL50803_112103', 'GL50803_10311'
    ],
    'Hydrogenosome': [
        'SS50377_23299', 'SS50377_25409', 'SS50377_27457', 'SS50377_28003'
    ],
    'Encystation': [
        'GL50803_5638', 'GL50803_5435', 'GL50803_2421',
        'GL50803_8245', 'GL50803_14259', 'GL50803_16069'
    ],
    'Meiosis': [
        'GL50803_15279', 'GL50803_4084', 'GL50803_13104',
        'GL50803_6626'
    ]
}

# Invert the dictionary to map each gene to its category
gene_to_category = {gene: category for category, genes in gene_categories.items() for gene in genes}

# Assuming df3 is already prepared from your previous steps
# Map the genes to the categories using the inverted dictionary
df3['Category'] = df3.index.map(gene_to_category)

# Group by the new categories and sum the counts
df3_grouped = df3.groupby('Category').sum()

# Replace 0 values with NaN (already in your code)
df3_grouped = df3_grouped.replace(0, np.nan)

# Plot the heatmap with the grouped data
sns.set_theme(style="white")

f, ax = plt.subplots(figsize=(10, 10))
heatmap = sns.heatmap(df3_grouped,
                     cmap=sns.color_palette("gray_r", as_cmap=True),
                     square=True,
                     fmt='g',
                     linewidths=.4,
                     annot=True,
                     cbar_kws={"shrink": .5}
                     )
heatmap.set_ylabel('')  # Removes the y-axis title
heatmap.set_xlabel('')

# Save the plot
out_file = "/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/diplo_virulence_grouped.blastp.png"
plt.savefig(out_file, format="png", bbox_inches='tight', dpi=1200)
plt.show()
