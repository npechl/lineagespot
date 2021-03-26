# Details

## Details regarding the execution of the scripts
As menttioned on the main Readme of the project, the workflow consists of three R scripts:
- `01-find-lineages-v3.R`: This is probably the script that you're mostly interested in, responsible for taking taking as input a single VCF file and, based on the rules derived from Pangoline and the reference SARS-CoV-2 genome, constructs the [first output table](https://github.com/BiodataAnalysisGroup/lineagespot#a-tab-delimited-file-tsv-containing-the-most-probable-lineages-that-have-been-found).
  - **input**: A single TSV file
  - **Parameters to specify**:
   1. `reference.path`: Path to the referece SARS-CoV-2 genome file (in FASTA format). (Line 16)
   2. `vcf.path`: Path to the VCF input file. (Line 18)
   3. `decision.rules.path`: Path of the Pangolin decision rules file (by default this file is downloaded to the given path). (Line 20 - No need to change this)
   4. `nreads`: Number of total reads. (Line 22)
  - **Output**: A single TSV file, containing the [first output table](https://github.com/BiodataAnalysisGroup/lineagespot#a-tab-delimited-file-tsv-containing-the-most-probable-lineages-that-have-been-found).
   1. **Output folder**: Same as input folder, in our case `vcf-files` directory.
   2. **Output file**: A single `.tsv` file having the same main name as input file, with a suffix `-lineagespot.tsv` at the end. For example, if the input file is `L1_S22_L001_freebayes.ann.vcf`, the output file will be `L1_S22_L001_freebayes.ann-lineagespot.tsv`.

- `02-collapse.R`: This script is responsible for merging the rows of the [first output table](https://github.com/BiodataAnalysisGroup/lineagespot#a-tab-delimited-file-tsv-containing-the-most-probable-lineages-that-have-been-found) that refer to the same Sars-CoV-2 lineage into a single row. It constructs the [second output table](https://github.com/BiodataAnalysisGroup/lineagespot#a-collapsed-table).
  - **input**: A single TSV file, occured as output from the `01-find-lineages-v3.R` script.
  - **Parameters to specify**:
   1. `analysis_output_file`: Path to the first TSV file, occured from `01-find-lineages-v3.R` script. (Line9)
  - **Output**: A single TSV file, containing the [second output table](https://github.com/BiodataAnalysisGroup/lineagespot#a-collapsed-table).
   1. **Output folder**: Same as input folder, in our case `vcf-files` directory.
   2. **Output file**: A single `.tsv` file having the same main name as input file, with a suffix `-collapsed.tsv` at the end. For example, if the input file is `L1_S22_L001_freebayes.ann-lineagespot.tsv`, the output file will be `L1_S22_L001_freebayes.ann-lineagespot-collapsed.tsv`.

- `03-compare-files.R`: This script is responsible for comparing the results from the `01-find-lineages-v3.R` script, when the latter is given multiple inputs. Practically, it compares the TSV tables and detects differences. It's useful when dealing with VCF files occured from different analysis methodologies, or VCF files from the same analysis method but from different samples (in that way we can compare samples within a time period and observe how specific mutations variate over time).
  - **input**: A set of TSV files occured as outputs from the `01-find-lineages-v3.R` script.
  - **Parameters to specify**:
   1. `input.folder`: The directory where the TSV tables are placed (here in `vcf-files` directory. (Line 6)
   2. `outputs`: The outputs folder, which is going to be created and where all comparison files are placed (Here we place the outputs inside `vcf-files/Outputs_compared` directory). (Line 12)
  - **Output**: A directory, containing all CSV comparison files, occured from pairwise comparisons of input files (in our case, the output directory is `vcf-files/Outputs_compared`). Each CSV output file corresponds to a specific pair of TSV tables.
 

## Details regarding differences between output table 1 and table 2
1. Each row in table 1 corresponds to a single rule in Pangoline's Decision Tree, meaning that it contains the majority of lineages more than once.
2. In contrast, table 2 contains all lineages once. Table 2 occurs from merging rules that end up to the same lineage in a single row.
3. Most significant columns in table 1 are:
 - `Lineage`: The lineage derived from this specific rule.
 - `Rules`: The corresponding rules
 - `Total`: Number of sub-rules in a rule (Basically, every rule corresponds to a single path in Pangoline's Decision Tree, so it definitely consists of sub-rules at every decision node).
 - `Tree Overlap`: Number of sub-rules that are met in our data, if we examine them serially.
 - `Total Overlap`: Number of sub-rules that are met in our data, if we don't examine them serially, but as a general set.
 - `Tree Ratio`, `Total Ratio`: The corresponding ratios, that lie within [0,1]. For example, if the number of sub-rules is 5 and four of them are verified (non-serially), the total ratio is `0.8`.
 - `Tree Av. DP`: The average value of DP in the positions of the sub-rules that are verified in our data, in a serial way (Tree-way).
 - `Total Av. DP`: The average value of DP in the positions of the sub-rules that are verified in our data, in a non-serial way (Total sub-rules).
 - `Av. DP`: The average value of DP parameter in our sample.
 - `Number of reads`: The number of reads in our sample.
 - `Sig`: Significance, indicates at which rate we find the lineage.
4. In contrast, the columns in table 2 are:
 - `Lineage`: The lineage
 - `Mean Tree Ratio`: The mean value of the Tree ratio column from all rules that end up to a specific lineage in table 1.
 - `Mean Total Ratio`: The mean value of the Total ratio column from all rules that end up to a specific lineage in table 1.
 - `Mean Total Ratio Var`: The variance of Total Ratio values from all rules that end up to a specific lineage in table 1.
 - `Mean Tree Av DP`: The mean value of Tree AV. DP values from all rules that end up to a specific lineage in table 1. In our case, it is a column of zeros, but it strictly depends on data.
 - `Mean Total Av DP`: The mean value of Total AV. DP values from all rules that end up to a specific lineage in table 1.

