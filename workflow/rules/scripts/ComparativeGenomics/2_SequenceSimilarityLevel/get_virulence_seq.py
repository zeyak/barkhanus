# List of virulence gene IDs
giardia_db_ids = [
    "GL50803_113797",
    "GL50803_115066",
    "GL50803_10330",
    "GL50803_103887",
    "GL50803_2902",
    "GL50803_11654",
    "GL50803_11683",
    "GL50803_17153",
    "GL50803_15097",
    "GL50803_15101",
    "GL50803_4812",
    "GL50803_86676",
    "GL50803_17230",
    "GL50803_4410",
    "GL50803_14019",
    "GL50803_16160",
    "GL50803_16779",
    "GL50803_11118",
    "GL50803_112103",
    "GL50803_10311",
    "SS50377_23299",
    "SS50377_25409",
    "SS50377_27457",
    "SS50377_28003",
    "GL50803_5638",
    "GL50803_5435",
    "GL50803_2421",
    "GL50803_8245",
    "GL50803_14259",
    "GL50803_16069",
    "GL50803_15279",
    "GL50803_4084",
    "GL50803_13104",
    "GL50803_6626"
]

# File paths
input_fasta = "/Users/zeyku390/PycharmProjects/barkhanus/resources/DiploProteoms/diplo_prot.fa"
output_fasta = "/Users/zeyku390/PycharmProjects/barkhanus/results/ComparativeGenomics/virulence_sequences.fa"


# Read the fasta file and write matching sequences to output
with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
    write_sequence = False
    for line in infile:
        if line.startswith(">"):
            # Extract the ID and remove any whitespace or newlines
            sequence_id = line.split(" ")[0][1:].strip()
            if sequence_id in giardia_db_ids:
                write_sequence = True
                outfile.write(line)  # Write header
                print(f"Matched ID: {sequence_id}")  # Debugging print statement
            else:
                write_sequence = False
        elif write_sequence:
            outfile.write(line)  # Write sequence



