Preprocessing Raw Data, Base R Version
================
Keith Baggerly
2018-05-23

-   [Outline](#outline)
-   [Libraries](#libraries)
-   [The Potti Data](#the-potti-data)
-   [The NCI60 Data](#the-nci60-data)
-   [The GEO Data](#the-geo-data)

Outline
=======

We want to organize the three raw datasets (Potti, NCI60, and GEO) in roughly consistent fashion to allow for analysis and cross-checking. Here, we define this to mean producing a gene by sample expression matrix for each dataset (with rownames and colnames), and data frame of sample information for each. Some challenges for each dataset include

-   parsing the header row of the Potti dataset to extract information about drug (e.g., Adria, Doce, etc), contrast group (0 or 1), and column index within the overall matrix for reference.
-   arranging the table of information from the NCI into matrix form with columns using the cell line name from cell name and the replicate array suffix from ID, while also looking for and removing duplicate row entries (since there are more rows than the 180 arrays we suspect were used would generate)
-   combining the results from the sensitive and resistant cohorts, and correcting the mislabeling of one sample in the annotation.

We save the data and annotation matrices in results, one RData file per dataset.

Libraries
=========

``` r
library(here)
library(readr)
library(GEOquery)
library(tidyr)
library(magrittr)
```

The Potti Data
==============

First, we load the expression data into a matrix.

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

Next, we parse the column names to extract sample information.

``` r
## Strip the suffixes R added to duplicate names on loading

raw_column_names <- 
  gsub("_[[:digit:]]+$", "", colnames(potti))

## The last character is either 0 or 1

contrast_group <- 
  substr(raw_column_names, 
         nchar(raw_column_names), nchar(raw_column_names))

## Fill in the gaps, since names are given just
## for the first and last entries for a drug

drug_name <- rep(NA, ncol(potti))
drug_name[nchar(raw_column_names) > 1] <- 
  substr(raw_column_names[nchar(raw_column_names) > 1],
         1, nchar(raw_column_names)[nchar(raw_column_names) > 1] - 1)
for(i in 2:length(drug_name)){
  if(is.na(drug_name[i])){
    drug_name[i] <- drug_name[i - 1]
  }
}

## Recast the column index as a 0-padded string of width 3

raw_index <- sprintf("%03d", c(1:ncol(potti)))

## Join all components together for a more useful name

full_names <- paste(drug_name, 
                    contrast_group, 
                    raw_index,
                    sep = "_")

## Combine into a data frame

potti_info <- 
  data.frame(drug_name = drug_name,
             contrast_group = contrast_group,
             raw_index = raw_index,
             row.names = full_names,
             stringsAsFactors = FALSE)

## Replace the initial (less helpful) column names

colnames(potti) <- full_names

## Clean up

rm(i, drug_name, contrast_group, raw_index, full_names)

## Save the results

save(potti, potti_info,
     file = here("results", "potti_data.RData"))
```

The NCI60 Data
==============

First, we load the data.

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
id_width <- nchar(nci60_tibble$ID[[1]])
all(nchar(nci60_tibble$ID) == nchar(nci60_tibble$ID[[1]]))
```

    ## [1] TRUE

Next, we look for duplicate entries. From our previous explorations, we know there are results for *about* 180 arrays, assuming 12625 probesets per array.

``` r
length(unique(nci60_tibble$`Probe Set Name`))
```

    ## [1] 12625

``` r
nrow(nci60_tibble) / 180
```

    ## [1] 12626

``` r
temp_tab <- table(nci60_tibble$`Probe Set Name`)
temp_tab[1:5]
```

    ## 
    ##  100_g_at   1000_at   1001_at 1002_f_at 1003_s_at 
    ##       360       180       180       180       180

``` r
which(temp_tab != 180)
```

    ## 100_g_at 
    ##        1

Looking more closely, we see that the number of rows is 12626 \* 180, so there are precisely 180 more rows than we expect. Checking the frequencies of the probesets, we see there's just one which is present more than 180 times: "100\_g\_at", which may well be the first probeset in alphabetical order. This is present exactly 360 times, so we suspect the values were repeated for some reason.

``` r
temp_tib <- nci60_tibble[nci60_tibble$`Probe Set Name` == "100_g_at", ]
temp_tib[1:3,]
```

    ## # A tibble: 3 x 10
    ##   `Probe Set Name` ID       GENE   panelnbr cellnbr pname  cellname Signal
    ##   <chr>            <chr>    <chr>     <int>   <int> <chr>  <chr>     <dbl>
    ## 1 100_g_at         GC32976… RABGG…        7       5 Leuke… K-562     192. 
    ## 2 100_g_at         GC32976… RABGG…        7       6 Leuke… MOLT-4     46.9
    ## 3 100_g_at         GC32976… RABGG…        7       3 Leuke… CCRF-CEM  145. 
    ## # ... with 2 more variables: Detection <chr>, `P Value` <dbl>

``` r
temp_tib[181:183,]
```

    ## # A tibble: 3 x 10
    ##   `Probe Set Name` ID       GENE   panelnbr cellnbr pname  cellname Signal
    ##   <chr>            <chr>    <chr>     <int>   <int> <chr>  <chr>     <dbl>
    ## 1 100_g_at         GC32977… RABGG…        7       5 Leuke… K-562     192. 
    ## 2 100_g_at         GC32977… RABGG…        7       6 Leuke… MOLT-4     46.9
    ## 3 100_g_at         GC32977… RABGG…        7       3 Leuke… CCRF-CEM  145. 
    ## # ... with 2 more variables: Detection <chr>, `P Value` <dbl>

``` r
all.equal(temp_tib[1:180, c(1, 3:10)], temp_tib[181:360, c(1, 3:10)])
```

    ## [1] TRUE

For some reason, this probeset was matched to two gene cluster (GC) ids, so the values were repeated twice. Since we don't care that much about the gene cluster, we can simply drop the rows associated with the (arbitrarily chosen) latter of the two, GC32977.

First, we split the ID column, which is a composite of gene associated information (the GC) and sample associated information (the replicate) into separate columns.

``` r
nci60_tibble <- 
  nci60_tibble %>%
  separate(ID, into = c("gene_id", "replicate"), sep = "_")

nci60_tibble[1:3, ]
```

    ## # A tibble: 3 x 11
    ##   `Probe Set Name` gene_id replicate GENE  panelnbr cellnbr pname cellname
    ##   <chr>            <chr>   <chr>     <chr>    <int>   <int> <chr> <chr>   
    ## 1 36460_at         GC26855 B         POLR…        7       5 Leuk… K-562   
    ## 2 36460_at         GC26855 B         POLR…        7       6 Leuk… MOLT-4  
    ## 3 36460_at         GC26855 B         POLR…        7       3 Leuk… CCRF-CEM
    ## # ... with 3 more variables: Signal <dbl>, Detection <chr>, `P
    ## #   Value` <dbl>

Then we trim off the unwanted rows.

``` r
nci60_tibble <- nci60_tibble[nci60_tibble$gene_id != "GC32977", ]
```

The resulting table is now properly dimensioned to be reorganized in matrix form.

Of course, before performing such a reorganization we'd like to check whether the sample results are always presented in the same order within a probeset, because if they are (which we suspect) we can organize things faster. The Signal, Detection, and P Value column entries may differ, but we expect the rest of the values (annotation) to be consistent.

Now we check whether all of the sample info is consistently ordered.

``` r
n_samples <- 180
n_genes <- 12625

sample_cols <- c("replicate", "panelnbr", "cellnbr", "pname", "cellname")
for(i in 1:length(sample_cols)){
  temp <- nci60_tibble[[sample_cols[i]]]
  if( all( temp == rep(temp[1:n_samples], n_genes) ) ){
    cat(sample_cols[i], "matches\n")
  }else{
    cat(sample_cols[i], "doesn't match\n")
  }
}
```

    ## replicate matches
    ## panelnbr matches
    ## cellnbr matches
    ## pname matches
    ## cellname matches

Now we check whether all of the gene info is consistently ordered.

``` r
n_samples <- 180
n_genes <- 12625

gene_cols <- c("Probe Set Name", "gene_id", "GENE")
for(i in 1:length(gene_cols)){
  temp <- nci60_tibble[[gene_cols[i]]]
  temp[is.na(temp)] <- ""
  if( all( temp == 
           rep(temp[seq(from = 1, to = n_genes * n_samples, 
                        by = n_samples)], each = n_samples) ) ){
    cat(gene_cols[i], "matches\n")
  }else{
    cat(gene_cols[i], "doesn't match\n")
  }
}
```

    ## Probe Set Name matches
    ## gene_id matches
    ## GENE matches

Ok, everything matches. Let's reorganize what remains into gene information, sample information, and expression values, discarding the Detection and P Value columns.

``` r
key_gene_rows <- seq(1, n_genes * n_samples, n_samples)

nci60_gene_info <- 
  as.data.frame(nci60_tibble[key_gene_rows, c("gene_id", "GENE")])
names(nci60_gene_info) <- c("gene_cluster", "gene_name")
rownames(nci60_gene_info) <- nci60_tibble[["Probe Set Name"]][key_gene_rows]

nci60_sample_info <- 
  as.data.frame(
    nci60_tibble[1:n_samples, 
                 c("cellname", "replicate", "pname",
                   "panelnbr", "cellnbr")])
rownames(nci60_sample_info) <- 
  paste(nci60_sample_info[["cellname"]],
        nci60_sample_info[["replicate"]], sep = "_")

nci60 <- 
  matrix(nci60_tibble[["Signal"]], 
         nrow = n_genes, ncol = n_samples, byrow = TRUE)
rownames(nci60) <- rownames(nci60_gene_info)
colnames(nci60) <- rownames(nci60_sample_info)
```

Now we clean up and save the results.

``` r
rm(nci60_tibble, temp_tib, temp_tab,
   key_gene_rows, n_genes, n_samples, i, 
   temp, gene_cols, sample_cols)

save(nci60, nci60_gene_info, nci60_sample_info,
     file = here("results", "nci60_data.RData"))
```

The GEO Data
============

Here we simply load in the matrices of expression values from the SOFT files of the two patient cohorts, combine them if their row names match, and record which samples are resistant and which are sensitive.

We load GSE349 (the resistant samples and one mislabeled sensitive sample) first.

``` r
gse349_gq <-
  getGEO(filename = here("data", "GSE349_family.soft.gz"))

n_genes <- nrow(Table(gse349_gq@gsms[[1]]))
n_samples <- length(gse349_gq@gsms)

gse349 <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(gse349) <- Table(gse349_gq@gsms[[1]])[, "ID_REF"]
colnames(gse349) <- names(gse349_gq@gsms)

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

Now we repeat the process for GSE350 (the sensitive samples).

``` r
gse350_gq <-
  getGEO(filename = here("data", "GSE350_family.soft.gz"))

n_genes <- nrow(Table(gse350_gq@gsms[[1]]))
n_samples <- length(gse350_gq@gsms)

gse350 <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(gse350) <- Table(gse350_gq@gsms[[1]])[, "ID_REF"]
colnames(gse350) <- names(gse350_gq@gsms)

for(i in 1:n_samples){
  if(all(rownames(gse350) == Table(gse350_gq@gsms[[i]])[, "ID_REF"])){
    gse350[, i] <- Table(gse350_gq@gsms[[i]])[, "VALUE"]
  }
}

gse350[1,]
```

    ## GSM4903 GSM4907 GSM4908 GSM4914 GSM4915 GSM4917 GSM4919 GSM4920 GSM4921 
    ## 38.0593 44.0123 41.8354 35.4666 43.0591 34.2387 28.6799 30.1072 45.7519 
    ## GSM4923 
    ## 47.4269

Check whether the tables can be safely combined...

``` r
all(rownames(gse349) == rownames(gse350))
```

    ## [1] TRUE

Now we bundle the matrices and record the response status.

``` r
geo <- cbind(gse349, gse350)

geo_sample_info <- 
  data.frame(sample_status = 
               c(rep("Resistant", ncol(gse349)),
                 rep("Sensitive", ncol(gse350))),
             geo_dataset = 
                c(rep("GSE349", ncol(gse349)),
                  rep("GSE350", ncol(gse350))),
             row.names = colnames(geo),
             stringsAsFactors = FALSE)
geo_sample_info["GSM4913", "sample_status"] <- "Sensitive"
```

Now we reorder the columns matrix columns to put all of the resistant samples first and all of the sensitive samples next, sorting by GSM id within a cohort.

``` r
sorted_geo_names <- 
  c(sort(rownames(geo_sample_info)[
    geo_sample_info$sample_status == "Resistant"]),
    sort(rownames(geo_sample_info)[
    geo_sample_info$sample_status == "Sensitive"]))

geo_sample_info <- 
  geo_sample_info[sorted_geo_names, ]

geo <- geo[, sorted_geo_names]
```

Now we clean up and save the results.

``` r
rm(gse349_gq, gse349, gse350_gq, gse350, i, 
   n_genes, n_samples, sorted_geo_names)

save(geo, geo_sample_info,
     file = here("results", "geo_data.RData"))
```
