import pandas as pd

# Read the GFF file, ignoring comment lines starting with '#'
df = pd.read_csv("resources/TrainGlimmerHMM/S_salmonicida/ncbi_dataset/data/GCA_000497125.2/genomic.gff",
                 sep='\t', header=None, comment='#')

# Filter to get only exon rows
df_exons = df[df[2] == "exon"]

# Select columns 0, 3, and 4
df_exons = df_exons.iloc[:, [0, 3, 4]]

df_exons.to_csv("resources/TrainGlimmerHMM/spiro_exons.csv", sep='\t', header=False, index=False)
