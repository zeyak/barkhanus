from snakemake.shell import shell

threads = snakemake.params.threads
db_name = snakemake.params.db_name
engine= snakemake.params.engine
out_dir = snakemake.params.out_dir

shell(f"cd {out_dir}")
shell(f"RepeatModeler -database {db_name} -engine {engine} -threads {threads}")

