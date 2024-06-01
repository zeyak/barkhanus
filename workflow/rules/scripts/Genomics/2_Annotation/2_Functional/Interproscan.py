from snakemake.shell import shell

proteome = snakemake.input.proteome
threads = snakemake.params.threads
out = snakemake.output.out

#remove the stars representing the stop codons
shell(f"""sed 's/*//' {proteome} > {proteome}_.faa""")
shell(f"""interproscan.sh -i {proteome}_.faa -o {out} -f tsv -iprlookup -goterms --pathways -cpu {threads} """)

"""
/data/zeynep/interproscan-5.47-82.0/interproscan.sh -i results/Genomics/2_Annotation/1_Structural/glimmerhmm/flye/raw/genome.faa_.faa -o results/Genomics/2_Annotation/interproscan_results.tsv -f tsv -iprlookup -goterms --pathways -cpu 32

"""