rule orthofinder:
    input:
         proteome="/opt/zeynep/barkhanus/resources/Comparison"
    output:
        directory('results/ComparativeGenomics/2_SequenceSimilarityLevel/')
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/ComparativeGenomics/2_SequenceSimilarityLevel/orthofinder.py"