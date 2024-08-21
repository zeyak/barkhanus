rule orthofinder:
    input:
         proteome="resources/comparison/"
    output:
        directory('results/ComparativeGenomics/2_SequenceSimilarityLevel/')
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/ComparativeGenomics/2_SequenceSimilarityLevel/Orthofinder.py"