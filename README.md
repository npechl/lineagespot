# Lineagespot

Lineagespot is a framework written in [R](https://www.r-project.org/), and aims to identify and assign different SARS-CoV-2 lineages based on a single variant file (i.e., variant calling format). 

## Installation

### Prerequisites

The following [CRAN](https://cran.r-project.org/) packages need to be installed in order to run the project:

```
install.packages(c("data.table", 
                   "stringr", 
                   "readr, 
                   "stringdist, 
                   "seqinr", 
                   "vcfR", 
                   "openxlsx"))
```

### Installing

The project can be downloaded as a ZIP file or by using git:

```
git clone https://github.com/BiodataAnalysisGroup/lineagespot.git
```

## Usage

The project consists of three main scripts:

- `01-find-lineages-v3.R`
- `02-map.R`
- `03-compare-files.R`

### Inputs

- `01-find-lineages-v3.R`: 
    - `reference.path`: Path to the referece SARS-CoV-2 genome file (in FASTA format).
    - `vcf.path`: Path to the inputs VCF file.
    - `decision.rules.path`: Path of the Pangolin decision rules file (by default this file is downloaded to the given path).
    - `nreads`: Number of total reads.

- `02-map.R`:
    - `total_ratio_threshold`: Total ratio minimum threshold.
    - `tree_ratio_threshold`: Tree ratio minimum threshold.
    - `analysis_output_folder`: Path to output folder.
    - `analysis_output_filename`: Input file path.
    - `samples_detected_filepath`: Path to XLSX file with relevant samples to be compared.

- `03-compare-files.R` (optional):
    - `input.folder`: Path to folder with multiple Lineagespot output files, that need to be compared

### Outputs

Two output files are produced by the Lineagespot:

#### A tab-delimited file (TSV) containing the most probable lineages that have been found:

|Lineage |Rules |Total |Tree Overlap |Total Overlap |Tree Ratio |Total Ratio |Tree Av. DP |Total Av. AD |Av. DP |Total Run Reads |
| ------ | ---- | ---- | ----------- | ------------ | --------- | ---------- | ---------- | ----------- | ----- | -------------- |
|B.1.177.17 |28931!='A',... |4 |4 |4 |1 |1 |7 |7 |21.3780 |69706 |
|B.1.177 |28931!='A',... |5 |5 |5 |1 |1 |1 |1 |21.3780 |69706 |
|B.1.177 |28931!='A',... |6 |4 |4 |0.66 |0.66 |1 |1 |21.3780 |69706 |



#### A collapsed table:

|Lineage |Max Tree Ratio |Max Total Ratio |Max Tree Av. AD |Max Total Av. AD |Average of Tree Av. AD |Average of Total Av. AD |
| ------ | ------------- | -------------- | -------------- | --------------- | --------------------- | ---------------------- | 
|B.1.177.17 |1 |1 |7 |7 |7 |7 |
|B.1.177 |1 |1 |7 |10 |1.01 |5.56 |

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.