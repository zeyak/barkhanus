import os
import subprocess
from Bio import SeqIO

# Function to create a directory if it doesn't exist
def create_directory(dir_path):
    try:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(f"Directory created: {dir_path}")
        else:
            print(f"Directory already exists: {dir_path}")
    except Exception as e:
        print(f"Error creating directory {dir_path}: {e}")
        raise

def trainGlimmerHMM(mfasta_file, exon_file, args):
    ''' Train GlimmerHMM with core genes '''
    train_dir = args['-d']
    create_directory(train_dir)  # Ensure train_dir exists

    command = ['trainGlimmerHMM', mfasta_file, exon_file] + sum(map(list, zip(args.keys(), args.values())), [])
    print(f"Running command: {' '.join(command)}")

    cline = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, stderr = cline.communicate()

    if cline.returncode != 0:
        print(f"Error running trainGlimmerHMM: {stderr.decode('utf-8')}")
        raise RuntimeError(f"trainGlimmerHMM failed with return code {cline.returncode}")
    print("trainGlimmerHMM completed successfully.")

def glimmerhmm(genome_file, train_dir, args):
    ''' Run glimmerHMM to predict genes/exons '''
    command = ['glimmerhmm', genome_file, train_dir] + sum(map(list, zip(args.keys(), args.values())), [])
    print(f"Running command: {' '.join(command)}")

    cline = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, stderr = cline.communicate()

    if cline.returncode != 0:
        print(f"Error running glimmerhmm: {stderr.decode('utf-8')}")
        raise RuntimeError(f"glimmerhmm failed with return code {cline.returncode}")
    print("glimmerhmm completed successfully.")

# Extracted from the provided script
genome_file = snakemake.input.genome_file
exon_file = snakemake.input.exon_file
mfasta_file = snakemake.input.mfasta_file

home_dir = snakemake.params.home_dir
train_dir = snakemake.params.train_dir

conf_file = snakemake.output.conf_file
gff_output = snakemake.output.gff_output

# Run trainGlimmerHMM
trainGlimmerHMM(mfasta_file, exon_file, {'-d': train_dir, '-n': '150', '-v': '50'})

# Modify the train config file to add onlytga=1
with open(conf_file, 'a') as conf_h:
    conf_h.write("BoostSgl 10\n")
    conf_h.write("onlytga 1\n")

# Loop through the scaffolds, and predict separately, then concatenate
with open(gff_output, 'w') as gff_h:
    gff_h.write("##gff-version 3\n")

    with open(genome_file, 'r') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seqid = record.id
            seqfile = f"{seqid}.fa"
            with open(seqfile, 'w') as seqh:
                SeqIO.write(record, seqh, 'fasta')

            tmp_output = os.path.join(home_dir, f"{seqid}.gff")
            glimmerhmm(seqfile, train_dir, {'-o': tmp_output, '-g': '', '-v': ''})

            with open(tmp_output, 'r') as tmp_h:
                for i, line in enumerate(tmp_h):
                    if i != 0:
                        gff_h.write(line)

            subprocess.check_call(f'rm {seqfile} {tmp_output}', shell=True)
