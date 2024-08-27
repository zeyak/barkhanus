import pandas as pd
#read csv
df = pd.read_csv("/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/diplo_virulence.blastp", sep="\t", header=None)
out_file = "/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/diplo_virulence.blastp.png"


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
#get group for each species
df4 = df.groupby('sp').size().reset_index(name='count')
#get group for bark and GL50803_10311
df5 = df.groupby([0, 'sp']).get_group(( 'GL50803_10311', 'bark'))


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
df3.to_excel("/Users/zeyku390/PycharmProjects/barkhanus/results/diamond/virulence/diplo_virulence.blastp.xlsx")

#makea heatmap of the pivot table
#leave zero values as white
import numpy as np
df3 = df3.replace(0, np.nan)
#plot
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme(style="white")

f, ax = plt.subplots(figsize=(10, 10))
heatmap = sns.heatmap(df3,
                     #norm=LogNorm(),
                     cmap=sns.color_palette("gray_r", as_cmap=True),
                     square=True,
                     fmt='g',
                     linewidths=.4,
                     annot=True,
                     cbar_kws={"shrink": .5}
                     )
heatmap.set_ylabel('')  # Removes the y-axis title
heatmap.set_xlabel('')
plt.savefig(out_file, format="png", bbox_inches='tight', dpi=1200)
plt.show()









