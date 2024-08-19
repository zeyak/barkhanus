# List of virulence gene IDs
giardia_db_ids = [
    "GL50803_00113797",
    "GL50803_00115066",
    "GL50803_0010330",
    "GL50803_00103887",
    "GL50803_002902",
    "GL50803_0011654",
    "GL50803_0011683",
    "GL50803_17153",
    "GL50803_0015097",
    "GL50803_15101",
    "GL50803_004812",
    "GL50803_0086676",
    "GL50803_0017230",
    "GL50803_004410",
    "GL50803_0014019",
    "GL50803_0016160",
    "GL50803_0016779",
    "GL50803_0011118",
    "GL50803_00112103",
    "GL50803_0010311",
    "SS50377_23299",
    "SS50377_25409",
    "SS50377_27457",
    "SS50377_28003",
    "GL50803_005638",
    "GL50803_005435",
    "GL50803_002421",
    "GL50803_008245",
    "GL50803_0014259",
    "GL50803_0016069",
    "GL50803_0015279",
    "GL50803_004084",
    "GL50803_0013104",
    "GL50803_006626"
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



