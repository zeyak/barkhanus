from snakemake.shell import shell

mfasta = snakemake.input.mfasta
genome = snakemake.input.genome
exon = snakemake.input.exon

gff = snakemake.output.gff
train_dir = snakemake.output.train_dir

shell(f"trainGlimmerHMM {mfasta} {exon} -n 150 -v 50 -d {train_dir}")
#shell(f"python glimmerhmm.py")
shell(f"glimmerhmm -g {genome} -o {gff} {train_dir}")
