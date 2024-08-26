rule orthofinder:
    input:
         proteome="resources/Comparison"
    output:
        directory('results/ComparativeGenomics/2_SequenceSimilarityLevel/')
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/ComparativeGenomics/2_SequenceSimilarityLevel/orthofinder.py"

rule make_diamond_db:
    input:
        "resources/{db}/{db}.faa"
    output:
        "resources/{db}/{db}.db.dmnd"
    conda:
        "envs/genomics.yaml"
    shell:
        "diamond makedb --in {input} --db {output}"

rule diamond_blastp_virulence:
    input:
        virulence="resources/virulence/virulence_sequences.faa",
        db = "resources/{db}/{db}.db.dmnd"
    output:
        'results/ComparativeGenomics/2_SequenceSimilarityLevel/diamond_blastp/'
    conda:
        "envs/genomics.yaml"
    script:
        "scripts/Genomics/2_Annotation/2_Functional/Diamond.py"

rule cdhit_faa:
    input:
        genome = "results/Genomics/2_Annotation/1_Structural/{annotation}/{assembler}/genome.faa"
    params:
        threads=32
    output: "results/ComparativeGenomics/2_SequenceSimilarityLevel/cdhit/{annotation}/{assembler}/genome_{n}.cdhit"
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/ComparativeGenomics/1_GenomeStructureLevel/cdhit.py"