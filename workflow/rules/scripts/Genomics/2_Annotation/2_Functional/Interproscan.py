from snakemake.shell import shell

proteome = snakemake.input.proteome
threads = snakemake.params.threads
out = snakemake.output.out

#remove the stars representing the stop codons
shell(f"""sed 's/*//' {proteome} > {proteome}_.faa""")
shell(f"""interproscan.sh -i {proteome}_.faa -o {out} -f gff3,tsv -iprlookup -goterms --pathways -cpu {threads} """)