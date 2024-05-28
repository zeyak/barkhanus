rule fastqc_before_trimming:
    input:
         input_dir="/data/zeynep/barkhanus_data/DNA/raw",
    params:
          threads=32,
    output:
          out_dir=directory("results/Genomics/1_Assembly/1_Preprocessing/fastqc_before_trimming/"),
    conda:
         "envs/genomics.yaml",
    script:
          "scripts/Genomics/1_Assembly/1_Preprocessing/ReadQualityCheck.py"

#Assembly
rule flye:
    input:
         reads="/data/zeynep/barkhanus_data/DNA/{process}/nanopore.fastq.gz",
    params:
          genome_size="114m",
          threads=32,
    output:
          out_dir=directory("results/Genomics/1_Assembly/2_Assemblers/flye/{process}/")
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/2_Assemblers/FlyeAssembler.py"

rule setup_nr_db:  #FIX how to actuvste this before running blastn
    output:
        outdir = protected(directory("/data/zeynep/databases"))
    conda:
        "envs/genomics.yaml"
    script:
        "scripts/Genomics/1_Assembly/3_Evaluation/setup_nr_db.py"

rule blastn:
    input:
        query="results/Genomics/1_Assembly/2_Assemblers/{assembler}/assembly.fasta",
        db="/data/zeynep/databases"
    output:
        "results/Genomics/1_Assembly/3_Evaluation/blastn/{assembler}/{db}/assembly.blastn"
    params:
        outfmt= "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle",
        threads=32,
        evalue=1e-10,
        db_prefix="/data/zeynep/databases/{db}"
    conda:
        "envs/genomics.yaml"
    script:
        "scripts/Genomics/1_Assembly/3_Evaluation/blastn.py"

rule meryl:
    input:
         genome="results/Genomics/1_Assembly/2_Assemblers/{assembler}/assembly.fasta",
    output:
          merylDB=directory("results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/merlyDB"),
          repetitive_k15="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/repetitive_k15.txt",
    params:
          threads=30,
          nanopore=True
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/3_Evaluation/CalculateKmerLongReads.py"

rule winnowmap:
    input:
         genome="results/Genomics/1_Assembly/2_Assemblers/{assembler}/assembly.fasta",
         long_read="/data/zeynep/barkhanus_data/DNA/clean/{long_read}.fastq.gz",
         merylDB="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/merlyDB",
         repetitive_k15="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/repetitive_k15.txt",
    output:
          sorted_bam="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/{long_read}.bam",
    #bai="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/{long_read}.bai"
    params:
          threads=32,
          nanopore=True
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/3_Evaluation/MapLongReadsToAssembly.py"

#Evaluation
rule quast:
    input:
         assembly="results/Genomics/1_Assembly/2_Assemblers/{assembler}/{process}/assembly.fasta",
    params:
          threads=32
    output:
          report_dir=directory("results/Genomics/1_Assembly/3_Evaluation/quast/{assembler}/{process}/")
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/3_Evaluation/AssemblyQualityCheck.py"

rule multiqc:
    input:
         input_dir="results/Genomics/1_Assembly/3_Evaluation/",
    params:
          threads=32
    output:
          out_dir=directory("results/Genomics/1_Assembly/3_Evaluation/multiqc/{assembler}")
    conda:
         "envs/genomics.yaml"
    shell:
         'multiqc {input.input_dir} -o {output.out_dir}'

"""rule plot_coverage_cont:
    input:
         #coverage on assembley
         run1="results/Genomics/1_Assembly/3_Evaluation/bowtie2/paired/{assembler}/illumina_run1.bam",
         run2="results/Genomics/1_Assembly/3_Evaluation/bowtie2/paired/{assembler}/illumina_run2.bam",
         run3="results/Genomics/1_Assembly/3_Evaluation/bowtie2/paired/{assembler}/illumina_run3.bam",
         pac="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/pacbio.bam",
         nano="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/nanopore.bam"
    output:
          out="results/Genomics/1_Assembly/3_Evaluation/deeptools/{assembler}.png",
          outraw="results/Genomics/1_Assembly/3_Evaluation/deeptools/{assembler}/outRawCounts.txt"
    params:
          threads=32,
    # ill1_P = "ill1_P",
    # ill2_P = "ill2_P",
    # ill3_P = "ill3_P",
    #pac = "pac",
    #nano = "nano",
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/3_Evaluation/PlotCoverage.py""""
