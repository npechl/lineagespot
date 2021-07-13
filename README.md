# Lineagespot

Lineagespot is a framework written in [R](https://www.r-project.org/), and aims to identify SARS-CoV-2 related mutations based on a single (or a list) of variant(s) file(s) (i.e., [variant calling format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format)). 

## Installation

### Prerequisites

The following [CRAN](https://cran.r-project.org/) packages need to be installed in order to run the project:

```R
install.packages(c("data.table", "vcfR", "stringr", "readr", "stringdist", "seqinr", "openxlsx"))
```

### Installing

The project can be downloaded as a ZIP file or by using git:

```
git clone https://github.com/BiodataAnalysisGroup/lineagespot.git
```

## Usage

The project consists of four main scripts:

- `01_merge_vcf.R`: A script for merging all identified mutations into an integrated table.
- `02a-compare_with_pangolin.R`: A script for comparing mutations with Pangolin tool.
- `02b-collapse.R`: A script for collapsing the table produced by the previous step (02a).
- `02c-compare-files.R`: A script for comparing different output files from previous steps (02a - 02b).

### Inputs

- `01_merge_vcf.R`:
    - `path_to_vcfFolder`: path to the folder containing VCF files
    - `path_to_metaData`: path to a tab-delimited meta data file
    
- `02a-compare_with_pangolin.R`: 
    - `reference.path`: Path to the referece SARS-CoV-2 genome file (in FASTA format).
    - `vcf.path`: Path to the inputs VCF file.
    - `decision.rules.path`: Path of the Pangolin decision rules file (by default this file is downloaded to the given path).
    - `nreads`: Number of total reads.

**It should be noted that the first script applys on VCF files which contain FORMAT or sample-specific information. Examples of this kind of format can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format) or [here](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/02_variant-calling.html).**

- `02b-collapse.R`:
    - `analysis_output_filename`: Input file path.

- `02c-compare-files.R` (optional):
    - `input.folder`: Path to folder with multiple Lineagespot output files, that need to be compared

*before running the ```03-compare-files.R``` script, make sure you have run the ```01-find-lineages-v3.R``` script for all the different vcf files.

### Outputs

Lineagespot provides an output tab-delimited table with all the identified mutations for all of the provided samples.

## Raw data analysis

The processing steps of the raw fastq files can be found [here](docs/raw-data-analysis.md).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use the tool, please cite the following work:

Nikolaos Pechlivanis, Maria Tsagiopoulou, Maria Christina Maniou, Anastasis Togkousidis, Evangelia Mouchtaropoulou, Taxiarchis Chassalevris, Serafeim Chaintoutis, Chrysostomos Dovas, Maria Petala, Margaritis Kostoglou, Thodoris Karapantsios, Stamatia Laidou, Elisavet Vlachonikola, Anastasia Chatzidimitriou, Agis Papadopoulos, Nikolaos Papaioannou, Anagnostis Argiriou, Fotis Psomopoulos, "_Detecting SARS-CoV-2 lineages and mutational load in municipal wastewater; a use-case in the metropolitan area of Thessaloniki, Greece_", medRxiv 2021.03.17.21252673; doi: [https://doi.org/10.1101/2021.03.17.21252673](https://doi.org/10.1101/2021.03.17.21252673)




