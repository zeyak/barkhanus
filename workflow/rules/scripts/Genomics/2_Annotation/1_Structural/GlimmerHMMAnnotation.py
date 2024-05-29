from snakemake.shell import shell

mfasta_file = snakemake.input.mfasta_file
genome = snakemake.input.genome
exon_file = snakemake.input.exon_file

gff = snakemake.output.gff
train_dir = snakemake.output.train_dir

shell(f"trainGlimmerHMM {mfasta_file} {exon_file} -n 150 -v 50 -d {train_dir}")
#shell(f"python glimmerhmm.py")
shell(f"glimmerhmm -g {genome} -o {gff} {train_dir}")
