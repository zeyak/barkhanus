from snakemake.shell import shell

genome= snakemake.input.genome
db = snakemake.input.db

outfmt = snakemake.params.outfmt
threads = snakemake.params.threads
max_target_seqs = snakemake.params.max_target_seqs
max_hsps = snakemake.params.max_hsps
more_sensitive = snakemake.params.more_sensitive
evalue = snakemake.params.evalue

output = snakemake.output

shell(f"""diamond blastp --query {genome} \
--db {db} \
--out {output} \
--outfmt {outfmt} \
--threads {threads} \
--evalue {evalue} \
--max-target-seqs {max_target_seqs} \
--max-hsps {max_hsps} \
--more-sensitive {more_sensitive} \
""")
