# GenoDiploBasic Pipeline

## Overview

The implementation of the Snakemake workflow management system revealed the need for custom adjustments and the integration of additional software for each genome project. Managing command-line tools sequentially introduced significant repetitive tasks, making it challenging to maintain comprehensive records and reproduce results. The project structure frequently evolved based on final outcomes, emphasizing the importance of reproducibility for efficient time management and keeping the project current for future advancements.

## GenoDiplo Pipeline

The GenoDiplo pipeline was initially tailored for the large diplomonad genome of *H. inflata*, characterized by sparse introns. However, the pipeline's structure is adaptable for other genome projects with the incorporation of additional software tools. In this study, a simplified version of the GenoDiplo pipeline was applied to *S. barkhanus*.

## Genome Assembly

The GenoDiploBasic pipeline focused on genome assembly using only Nanopore long reads, excluding Illumina polishing and contamination processes due to the bacteria-free cultivation of *S. barkhanus*. Nanopore reads, being long and accurate, facilitated a compact assembly given the genome size compared to *S. salmonicida*. The FLYE assembler processed the high-quality Nanopore sequencing data, yielding an assembly of XXX Mbp and XXX contigs for the commensal *S. barkhanus*, the second assembled *Spironucleus* genome in the Hexamitinae branch.

## Gene Prediction and Functional Annotation

### Gene Prediction

The genome annotation pipeline employed custom workflows previously used for other diplomonad genomes. The pipeline utilized GlimmerHMM, trained on *S. salmonicida*, alongside Prodigal for prokaryotic gene calling to compare single-exon genes.

### Functional Annotation

The functional annotation pipeline adopted a hierarchical approach similar to the GenoDiplo pipeline for *H. inflata*. The process began with sequence similarity searches across diplomonad genomes and transcriptomes (e.g., *Trepomonas*), followed by InterProScan to identify functional genes not found in previous diplomonad annotations. This resulted in the identification of XXX functional genes and XXX hypothetical genes in *S. barkhanus*. Additionally, tRNA, rRNA, and RepeatMasker software were employed to annotate non-coding regions of the genome.

## Getting Started

### Prerequisites

- Snakemake

### Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/genodiplo.git
cd genodiplo
```

Install the required dependencies:

```bash
# Example for installing dependencies
conda env create -f environment.yml
conda activate genodiplo
```

### Usage

Run the pipeline:

```bash
snakemake --cores <number_of_cores>
```

Customize the configuration file (`config.yaml`) to match your project requirements.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to discuss your ideas.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
