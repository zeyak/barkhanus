from snakemake.shell import shell

genome = snakemake.input.genome
out = snakemake.output[0]

#seq_identity = snakemake.params.seq_identity
threads = snakemake.params.threads
identity = snakemake.wildcards.n


#shell(f"""cd-hit -i {genome} -o {out} -c {wildcards.n} -T {threads}""")
shell(f"""cd-hit -i {genome} -o {out} -c {identity} -T {threads} -g 1 -M 5000""")