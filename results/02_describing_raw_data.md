Describing the Raw Data
================
Keith Baggerly
2018-05-23

-   [Overview](#overview)
-   [Libraries](#libraries)
-   [The Potti Data](#the-potti-data)
    -   [Text Description](#text-description)
    -   [File Examination](#file-examination)
-   [The NCI60 Data](#the-nci60-data)
    -   [Text Description](#text-description-1)
    -   [File Examination](#file-examination-1)
-   [The GEO Data](#the-geo-data)
    -   [Text Description](#text-description-2)
    -   [File Examination](#file-examination-2)

Overview
========

The three raw datasets we're working with (potti, nci60, and geo) all contain gene expression data from Affymetrix U95 arrays, but they also have associated metadata. For review, we briefly describe each of the three below, and look at matrix arrangements of the data files.

Libraries
=========

``` r
library(here)
library(readr)
library(GEOquery)
```

The Potti Data
==============

Text Description
----------------

The figshare page,

<https://figshare.com/s/66603862d770b4c73146>

has some description of the dataset which isn't included in the downloaded version. Some relevant excerpts:

> In late 2006, [Potti et al](https://www.nature.com/articles/nm1491) published an article in Nature Medicine in which they claimed to have found a way to use drug sensitivity information and genomic profiles of cell lines to infer likely patient response to treatment from a patient's genomic profile.

> The cell lines in question were those from the NCI60 panel, which has been maintained by the National Cancer Institute (NCI) for many decades.

> Since investigators at MD Anderson wanted to use the approach to improve care, they asked us for help evaluating the approach. We asked the authors if they could be more specific about precisely which genomic profiles were being used, and which cell lines were being treated as sensitive and resistant for each of the 7 drugs they examined.

> The first file we got back in response was the "chemo.zip" file included here.

File Examination
----------------

``` r
potti_tibble <- 
  read_delim(here("data", "chemo.zip"), delim = "\t")
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   probe_set = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
potti <- as.matrix(potti_tibble[, 2:ncol(potti_tibble)])
rownames(potti) <- potti_tibble[[1]]

rm(potti_tibble)

potti[1:5, 1:5]
```

    ##             Adria0         0       0_1       0_2        0_3
    ## 36460_at  41.67195  21.82034 125.79484  93.45925  79.063210
    ## 36461_at 171.39024 122.09759 218.35896 179.27033  77.453125
    ## 36462_at 147.49791 203.84113 211.10683 208.76929 216.542786
    ## 36463_at 151.96898 122.31647 128.46779 186.08681 111.650917
    ## 36464_at  20.99779  35.43748  12.16202  19.52815   8.648886

The uncompressed data file, "Chemo predictors (U-95) All - FINAL.txt", is a tab-delimited file with a 12559 rows: 1 header row and 12558 data rows, where each row corresponds to a probeset id (roughly a gene) from the Affymetrix U-95A gene chip (sequences used to query the transcriptome were built using the 95th assembly of the Unigene consensus definitions). There are actually 12625 probesets on the U-95A, but 67 of these are "control probes" used to assess general functioning of the assay. The control probe names begin with "AFFX", and these do not appear to be present here.

There are 135 columns in the dataset. The first column lists the probeset ids (rownames), and the rest give (we presume) expression measurements from one of the cell lines in the NCI-60 panel. The column names suggest the columns are grouped by drug, as there are 7 "name blocks" which start with

"(drug name)0 0 0..." and end with "...1 1 1 (drug name)1"

For a given drug, we suspect 0/1 labels indicate sensitive/resistant status, with 0 indicating one group and 1 indicating the other. We do not see anything in the dataset itself specifying which group is which.

Total numbers of cell lines vary by drug, as do relative numbers of 0's and 1's. It's not immediately clear what rule was used to select these just these lines.

Since the total number of data columns (134) exceeds the number of cell lines in the NCI60 panel (59), we suspect some columns should be repeated, indicating some cell lines supplied information about more than one drug.

Looking at the data matrix with View, we see the expression values are recorded to 6 decimal places, so matching just one probeset may serve to match the entire column. Looking at the first row (probeset 36460\_at) and first column (Adria0), we see a value of 41.671947. Skimming across columns, there are 4 with precisely this value.

``` r
c(1:ncol(potti))[which(potti[1,] == potti[1,1])]
```

    ## [1]   1  43  95 119

We don't (yet) know the names of the specific cell lines involved, but we think we have some idea what the values represent.

The NCI60 Data
==============

Text Description
----------------

Again, some of this comes from the figshare page above.

> There are 59 distinct cell lines in the NCI60 panel, and both drug sensitivity and genomic profile information for these cell lines are publicly available from the NCI's Developmental Therapeutics Program (DTP).

A list of the cell lines in the panel is available here

<https://dtp.cancer.gov/discovery_development/nci-60/cell_list.htm>

This list also includes links for a few of the cell lines which have been discovered to be something other that what they were thought to be in many earlier studies. There were 60 cell lines when the NCI60 panel was first developed, but cell line MDA-N was soon found to have been derived from MDA-MB-435, so MDA-N has rarely been separately profiled since.

> The raw NCI60 microarray data is available from the DTP:

> <https://wiki.nci.nih.gov/display/ncidtpdata/molecular+target+data>

> This page has a link to data from Affymetrix U95A arrays run on the panel in triplicate by Novartis, "WEB\_DATA\_NOVARTIS\_ALL.ZIP"

> The current data on growth inhibition is here:

> <https://wiki.nci.nih.gov/display/NCIDTPdata/NCI-60+Growth+Inhibition+Data>

The molecular target data page above also contains some brief descriptions of the types of information reported in the various table columns

> When uncompressed the file is comma delimited in the following format:

> File Format: Probe Set Name, ID (composite of the moltid derived from this measurement, and a letter to distinguish individual arrays), GENE, cellname, pname, PANELNBR, CELLNBR, Signal, Detection, P value Gene assignments are based on Unigene Build \#U225 (August 2010)

File Examination
----------------

``` r
nci60_tibble <- 
  read_csv(here("data", "WEB_DATA_NOVARTIS_ALL.zip"))
```

    ## Parsed with column specification:
    ## cols(
    ##   `Probe Set Name` = col_character(),
    ##   ID = col_character(),
    ##   GENE = col_character(),
    ##   panelnbr = col_integer(),
    ##   cellnbr = col_integer(),
    ##   pname = col_character(),
    ##   cellname = col_character(),
    ##   Signal = col_double(),
    ##   Detection = col_character(),
    ##   `P Value` = col_double()
    ## )

``` r
nci60_tibble[1:13, ]
```

    ## # A tibble: 13 x 10
    ##    `Probe Set Name` ID      GENE   panelnbr cellnbr pname  cellname Signal
    ##    <chr>            <chr>   <chr>     <int>   <int> <chr>  <chr>     <dbl>
    ##  1 36460_at         GC2685… POLR1C        7       5 Leuke… K-562     120. 
    ##  2 36460_at         GC2685… POLR1C        7       6 Leuke… MOLT-4    114. 
    ##  3 36460_at         GC2685… POLR1C        7       3 Leuke… CCRF-CEM   67.4
    ##  4 36460_at         GC2685… POLR1C        7      10 Leuke… RPMI-82…  148. 
    ##  5 36460_at         GC2685… POLR1C        7       8 Leuke… HL-60(T…  105. 
    ##  6 36460_at         GC2685… POLR1C        7      19 Leuke… SR         41.8
    ##  7 36460_at         GC2685… POLR1C       12      14 CNS    SF-268     67.6
    ##  8 36460_at         GC2685… POLR1C       12      15 CNS    SF-295     70.8
    ##  9 36460_at         GC2685… POLR1C       12      16 CNS    SF-539     59.9
    ## 10 36460_at         GC2685… POLR1C       12       2 CNS    SNB-19     62.4
    ## 11 36460_at         GC2685… POLR1C       12       5 CNS    SNB-75     28.8
    ## 12 36460_at         GC2685… POLR1C       12       9 CNS    U251       61.7
    ## 13 36460_at         GC2685… POLR1C        7       5 Leuke… K-562      95.8
    ## # ... with 2 more variables: Detection <chr>, `P Value` <dbl>

``` r
length(unique(nci60_tibble$`Probe Set Name`))
```

    ## [1] 12625

``` r
nrow(nci60_tibble) / length(unique(nci60_tibble$`Probe Set Name`))
```

    ## [1] 180.0143

``` r
length(unique(nci60_tibble$cellname))
```

    ## [1] 59

Each cell line was intended to be profiled in triplicate, so we expected to see data on 3 \* 59 = 177 arrays. As it happens, we appear to have data on a *fraction* over 180, which suggests (a) some cell lines were run more that 3 times and (b) we may have some duplicate entries. The cell line names we actually want are in the "cellname" column. The replicate index for the array is given by the one-letter suffix in the "ID" column - the first 12 rows are reporting values from the "B" set of array replicates for the cell lines named. Row 13, by contrast, is reporting the value of probeset 36460\_at for the "A" set replicate of cell line K-562, which is a leukemia cell line. Many different tumor types are represented in the NCI60 panel, so cell lines of a given type are often grouped into subpanels. Looking at the first row, we see that K-562 is a Leukemia cell line, and the Leukemia panel is number 7. K-562 is the 5th cell line of a larger leukemia panel; this larger panel includes some cell lines which are not included in the NCI60. Looking at the Signal values with View suggests these are reported to 6 decimal places, as were the values in the table we got from Potti et al. That, together with the fact that the same probeset (36460\_at) is at the top of both lists, makes us optimistic that the Potti et al values may come from here. As a quick check, let's see if we can match the first entry of the Potti et al data matrix.

``` r
which(nci60_tibble$Signal == potti[1,1])
```

    ## [1] 21

``` r
nci60_tibble[which(nci60_tibble$Signal == potti[1,1]), ]
```

    ## # A tibble: 1 x 10
    ##   `Probe Set Name` ID        GENE   panelnbr cellnbr pname cellname Signal
    ##   <chr>            <chr>     <chr>     <int>   <int> <chr> <chr>     <dbl>
    ## 1 36460_at         GC26855_A POLR1C       12      16 CNS   SF-539     41.7
    ## # ... with 2 more variables: Detection <chr>, `P Value` <dbl>

We get one (and only one!) hit, suggesting the decimal place accuracy may be detailed enough to make matching pretty easy. Here, this match suggests that one of the cell lines in Potti et al's "0 group" for Adriamycin is SF-539, which is derived from a tumor of the central nervous system (CNS). Googling suggests this is a glioblastoma (brain tumor) cell line. The value used for SF-539 is for the "A" replicate assay of the cell line, not from some aggregate across replicates.

The GEO Data
============

Text Description
----------------

Potti et al describe using the NCI60 data to predict response for various patient cohorts for which corresponding microarray data were publicly available. One of these cohorts was a group of 24 women with breast cancer who were treated with single agent docetaxel; the first report on this cohort was by [Chang et al](https://www.ncbi.nlm.nih.gov/pubmed/12907009) in the Lancet in 2003. Chang et al dichotomized the patients into "Sensitive" and "Resistant" subcohorts of sizes 11 and 13, respectively. These patient tumors were interrogated with the same type of microarray (Affymetrix U-95Av2) as the NCI60 panel, so matching gene measurements will be easier for this dataset than for others where genes were interrogated using different strands of cDNA. The gene expression data were posted to the Gene Expression Omnibus (GEO) as GSE349 (the resistant cohort) and GSE350 (the sensitive cohort). As [supplementary information](http://www.thelancet.com/pb-assets/Lancet/extras/01art11086webtable.pdf), the authors provide a table (in pdf) of the expression values by sample for each of the 92 genes they found to be particularly important for distinguishing sensitive from resistant.

Interestingly, GSE349 has data for 14 samples, not 13, and GSE350 has data for 10 samples, not 11. Personal communication with one of the Chang et al authors (Sue Hilsenbeck) confirmed that one of sensitive samples had been uploaded with the resistant group by mistake.

> Sample \#377 is mislabeled in the GEO DB. It is listed there as resistant, but in reality it was sensitive.

Checking the GEO annotation, sample 377 corresponds to GSM4913.

File Examination
----------------

Data are available from GEO in a variety of formats; in general, if all we really want are the expression values, we the "Series Matrix" files are streamlined for this and faster to download. Here we've downloaded the SOFT format files for the two Gene Set Series (GSEs) in order to capture more of the annotation and structure of the data as they are presented on GEO itself. In this case, both GSEs are comprised of samples (GSMs) from the same type of assay, and the values should be in the same order. We'll check this for GSE349 here.

``` r
gse349_gq <-
  getGEO(filename = here("data", "GSE349_family.soft.gz"))

names(gse349_gq@gsms)
```

    ##  [1] "GSM4901" "GSM4902" "GSM4904" "GSM4905" "GSM4906" "GSM4909" "GSM4910"
    ##  [8] "GSM4911" "GSM4912" "GSM4913" "GSM4916" "GSM4918" "GSM4922" "GSM4924"

We have data on one array platform (GPL8300, the U-95Av2 array) and 14 individual samples (with GSM ids). Now we can check what's recorded for each GSM, and preallocate space for a matrix to hold the expression values assuming the ordering is the same.

``` r
Columns(gse349_gq@gsms[[1]])
```

    ##   Column Description
    ## 1 ID_REF            
    ## 2  VALUE

``` r
Table(gse349_gq@gsms[[1]])[1:5,]
```

    ##            ID_REF      VALUE
    ## 1  AFFX-MurIL2_at   54.11840
    ## 2 AFFX-MurIL10_at   81.42730
    ## 3  AFFX-MurIL4_at    5.89822
    ## 4  AFFX-MurFAS_at   68.90200
    ## 5  AFFX-BioB-5_at 1385.53000

``` r
dim(Table(gse349_gq@gsms[[1]]))
```

    ## [1] 12625     2

``` r
n_genes <- nrow(Table(gse349_gq@gsms[[1]]))
n_samples <- length(gse349_gq@gsms)

gse349 <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(gse349) <- Table(gse349_gq@gsms[[1]])[, "ID_REF"]
colnames(gse349) <- names(gse349_gq@gsms)
```

Now we can load in the data one GSM at a time, only loading values if the probeset ordering is the same as what we've seen before.

``` r
for(i in 1:n_samples){
  if(all(rownames(gse349) == Table(gse349_gq@gsms[[i]])[, "ID_REF"])){
    gse349[, i] <- Table(gse349_gq@gsms[[i]])[, "VALUE"]
  }
}

gse349[1,]
```

    ## GSM4901 GSM4902 GSM4904 GSM4905 GSM4906 GSM4909 GSM4910 GSM4911 GSM4912 
    ## 54.1184 49.3722 49.1230 75.3204 54.2316 52.3357 41.1346 28.7636 43.5179 
    ## GSM4913 GSM4916 GSM4918 GSM4922 GSM4924 
    ## 46.4989 37.4464 33.4237 35.5482 32.2541

All columns have nonzero entries, so probeset ordering is maintained and everything loaded successfully. Let's take a look at the first few values.

``` r
gse349[1:4, 1:4]
```

    ##                  GSM4901  GSM4902  GSM4904  GSM4905
    ## AFFX-MurIL2_at  54.11840 49.37220 49.12300 75.32040
    ## AFFX-MurIL10_at 81.42730 65.21570 65.17600 96.05920
    ## AFFX-MurIL4_at   5.89822  6.07494  5.89822  5.89822
    ## AFFX-MurFAS_at  68.90200 47.21600 37.82200 73.67270

``` r
gse349[65:70, 1:4]
```

    ##                       GSM4901   GSM4902  GSM4904  GSM4905
    ## AFFX-YEL018w/_at      11.9622   8.83069  11.3278  18.9825
    ## AFFX-YEL024w/RIP1_at  72.0352  73.25030  57.5064 121.7420
    ## AFFX-YEL021w/URA3_at 141.2240 119.19100 124.5780 146.8590
    ## 31307_at              88.9291 113.68900  80.4923  93.0181
    ## 31308_at             199.5690 162.03200 153.5350 224.4540
    ## 31309_r_at            15.5754  19.49600  22.6177  23.4324

``` r
grep("^AFFX", rownames(gse349))
```

    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
    ## [24] 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46
    ## [47] 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67

The first few probeset ids don't match what we saw for the other datasets. These are control probes (their names begin with AFFX). There are 67 control probes on the U-95A arrays, so it seems likely these are in the first 67 positions; a quick check with grep confirms this. Even after the control probes, however, the probeset ordering differs from that used in the other datasets, as the first non-control probeset is 31307\_at as opposed to 36460\_at.

Going back to the control probes, let's take a closer look at that third row.

``` r
gse349[3, ]
```

    ##  GSM4901  GSM4902  GSM4904  GSM4905  GSM4906  GSM4909  GSM4910  GSM4911 
    ##  5.89822  6.07494  5.89822  5.89822  8.19463  8.33974 16.07990 21.64230 
    ##  GSM4912  GSM4913  GSM4916  GSM4918  GSM4922  GSM4924 
    ## 30.90970 25.21470 20.16330 10.87110 27.79290 22.98790

``` r
apply(gse349, 2, min)
```

    ## GSM4901 GSM4902 GSM4904 GSM4905 GSM4906 GSM4909 GSM4910 GSM4911 GSM4912 
    ## 5.89822 5.89822 5.89822 5.89822 5.89822 5.89822 5.89822 5.89822 5.89822 
    ## GSM4913 GSM4916 GSM4918 GSM4922 GSM4924 
    ## 5.89822 5.89822 5.89822 5.89822 5.89822

There are tied values in this dataset - we see several occurrences of 5.89822. This is because the quantification method used for this dataset is different than the one used for the other data, and in particular it allows for ties near the extremes. Checking with apply shows the value we see repeated is the minimum value for each sample, corresponding to effectively no expression.
