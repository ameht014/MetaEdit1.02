# MetaEdit - Computational Identification of RNA editing in
Microbiomes

MetaEdit is a computational pipeline designed to detect and characterize RNA editing events within bacterial microbiomes using metagenomic and metatranscriptomic data. By identifying discrepancies between DNA and RNA sequencing reads, MetaEdit provides insights into post-transcriptional modifications that influence microbial gene regulation, protein function, and adaptation within complex communities like the human gut microbiome.

## Methodology

MetaEdit analyzes sequencing reads in six steps:

* Read Filtering: High-quality reads are retained after trimming adapters and filtering host-derived sequences.
* Alignment: Reads are mapped to microbial and host reference genomes.
* Base Distribution Calculation: Base frequencies are computed for each genomic position.
* Coverage Analysis: Ensures sufficient sequencing depth for reliable editing detection.
* RNA Editing Prediction: Consistent discrepancies between DNA and RNA are aggregated across samples to identify robust editing sites.
* Functional Annotation: RNA editing sites are linked to corresponding genomic and protein changes.


## Source Code
[src/](src/) contains the source code for the pipeline. It includes [aimap_modified](src/aimap_modified), which is a modified version of the [AIMAP pipeline](https://github.com/bioone/aimap) used to predict RNA editing of a single bacteria in a single sample. The [src/metaEdit.py](src/metaEdit.py) contains the main mMtaEdit function to merge the results from the single-sample analysis and filter the editing sites based on criteria such as sequencing depth, consistency across samples, and the absence of corresponding variants in the DNA reads. 

## Usage

The [singularity environment file](./env/metaEdit.def) contains all dependencies required to run the pipeline.

The [experiments](experiments/) folder provides an example of running MetaEdit along with the results reported in the paper.

[experiments/1-select_samples.ipynb](experiments/1-select_samples.ipynb) shows the cohort analysis used to select samples from iHMP study.

[experiments/2-download_dna.sh](experiments/2-download_dna.sh) provides the shell script used to download the DNA sequencing data from iHMP study.

[experiments/3-download_rna.sh](experiments/3-download_rna.sh) provides the shell script used to download the RNA sequencing data from iHMP study.

[experiments/4-submit_DNA.sh](experiments/4-submit_DNA.sh) provides the shell script used to submit the jobs for analysis of DNA sequencing data, including trimming, alignment, base distribution calculation, coverage analysis, and SNPs calling.

[experiments/5-submit_RNA.sh](experiments/5-submit_RNA.sh) provides the shell script used to submit the jobs for analysis of RNA sequencing data, including trimming, alignment, base distribution calculation, coverage analysis, and RNA editing prediction.

[experiments/6-submit_filter_RNA.sh](experiments/6-submit_filter_RNA.sh) provides the shell script used to submit the jobs for filtering the RNA editing sites and output final RNA editing results.