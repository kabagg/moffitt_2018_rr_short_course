README
================
Keith Baggerly
2018-05-23

-   [Overview](#overview)
-   [Brief Results](#brief-results)
-   [Running the Analysis](#running-the-analysis)
    -   [Required Libraries](#required-libraries)

Overview
========

We want to illustrate assembly of a reproducible analysis using a dataset we care about. Our workflow closely follows that of Jenny Bryan's [packages-report-EXAMPLE](https://github.com/jennybc/packages-report-EXAMPLE) on GitHub.

Several years ago, [Potti et al](https://www.nature.com/articles/nm1491) claimed to have found a way to use microarray profiles of a specific panel of cell lines (the NCI60) to predict cancer patient response to chemotherapeutics from a similar profile of the patient's tumor. Using different subsets of cell lines, they made predictions for several different drugs. We wanted to apply their method, so we asked them to send us lists of which cell lines were used to make predictions for which drugs. The method doesn't work; we describe our full analyses [here](https://projecteuclid.org/euclid.aoas/1267453942).

The first dataset we received from Potti et al didn't have the cell lines labeled. We want to see if we can identify where the numbers came from and see if there were any oddities that should have raised red flags early on.

Brief Results
=============

-   [01\_gather\_data](results/01_gather_data.md) downloads the raw datasets used from the web.
-   [02\_describing\_raw\_data](results/02_describing_raw_data.md) both summarizes background information about the datasets in question and describes what we can see by looking at the structure of the data with minimal processing.
-   [03\_preprocessing\_data\_base\_r\_version](results/03_preprocessing_data_base_r_version.md) reorganizes the datasets into matrices of expression values and data frames of annotation using consistent names.
-   [04\_check\_sample\_matches\_and\_corrs](results/04_check_sample_matches_and_corrs.md) establishes the mappings desired.
-   [05\_report\_matches\_to\_potti\_columns](results/05_report_matches_to_potti_columns.md) summarizes the results of our analyses in more formal "report" format.

Running the Analysis
====================

Roughly, our analyses involve running the R and Rmd files in [R](R) in the order they appear.

Run [R/95\_make\_clean.R](R/95_make_clean.R) to clear out any downstream products.

Run [R/99\_make\_all.R](R/99_make_all.R) to re-run the analysis from beginning to end, including generating this README.

Raw data from the web is stored in \[data\]\[data\].

Reports and interim results are stored in \[results\]\[results\].

Required Libraries
------------------

These analyses were performed in RStudio 1.1.447 using R version 3.5.0 (2018-04-23), and use (in alphabetical order):

-   downloader 0.4
-   GEOquery 2.47.18
-   here 0.1
-   lattice 0.20.35
-   magrittr 1.5
-   readr 1.1.1
-   rmarkdown 1.9
-   tidyr 0.8.1

All of the above except `GEOquery` are available from [CRAN](https://cran.r-project.org).

`GEOquery` is available from [BioConductor](https://www.bioconductor.org), and can be installed from the R command line with

``` r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
```
