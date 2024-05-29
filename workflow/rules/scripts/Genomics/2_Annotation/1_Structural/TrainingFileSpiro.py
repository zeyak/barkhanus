import pandas as pd

df =pd.read_csv("resources/TrainGlimmerHMM/S_salmonicida/ncbi_dataset/data/GCA_000497125.2/genomic.gff",
                sep='\t', header=None, comment='#')

#get exons on teh 2nd column
#make df with 3 colummns: 0,2, 3, 4
df_exons = df[df[2] == "exon"]
df_exons = df_exons.iloc[:, [0, 3, 4]]

df_exons.to_csv("resources/TrainGlimmerHMM/spiro_exons.csv", sep='\t', header=False, index=False)
