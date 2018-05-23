Check Sample Matches and Correlations
================
Keith Baggerly
2018-05-23

-   [Overview](#overview)
-   [Librares and Data](#librares-and-data)
    -   [Libraries](#libraries)
    -   [Data](#data)
-   [Matching the Cell Lines Used](#matching-the-cell-lines-used)
-   [Check Sample Correlations in the Potti et al Data](#check-sample-correlations-in-the-potti-et-al-data)
-   [Check Other Summaries for Potti et al Data](#check-other-summaries-for-potti-et-al-data)
-   [Matching Potti to GEO](#matching-potti-to-geo)
-   [Combine all Match Info Into a Data Frame](#combine-all-match-info-into-a-data-frame)
-   [Save the Results](#save-the-results)

Overview
========

Given what we've seen of the Potti et al and NCI60 datasets, we want to see whether we can nail down the precise mappings used.

We load the data matrices for Potti et al, the NCI60, and the docetaxel patient cohort posted to GEO. We check for exact matches between expression values in the Potti et al and NCI60 datasets, and between the Potti et al and GEO datasets if the former fails. We also examine the pairwise sample correlations between columns to see whether that clarifies any of the structure present.

Precise matching on data values lets us match all of the Potti et al columns for 5 of the 7 drugs to columns in the NCI60 dataset; we know which cell lines are being used in each of the contrast groups. We find no matches for columns associated with Docetaxel (Doce) and Cytoxan (Cytox); correlation plots show the contents of these columns are wholly uncorrelated with those in the rest of the data matrix. Examining minimums shows the minimum values in these columns are (a) all the same, and (b) the same as the minimum values we saw in earlier examination of the GEO data. Again looking for exact matches shows the Doce and Cytox columns match columns from the GEO datasets after allowing for an offset of 67 rows to exclude the Affymetrix control probes reported in GEO. These matches are exact for the first 12535 of 12558 probesets reported; the last 23 use a slightly different order. Importantly, the values in the Potti table match the GEO values *in the order they are presented at GEO*, but these probeset ids do not match those used in the Potti et al table, effectively producing a random scrambling of the data values.

Librares and Data
=================

Libraries
---------

``` r
library(here)
library(lattice)
```

Data
----

``` r
load(here("results", "potti_data.RData"))
load(here("results", "nci60_data.RData"))
load(here("results", "geo_data.RData"))
```

Matching the Cell Lines Used
============================

Our first question is whether we can match the data values supplied by Potti et al to the NCI60 data, and thus infer the identities of the cell lines used to characterize response for each drug. Based on our earlier description of the raw data tables, we think checking for matches of the first one or two probeset values might suffice.

``` r
potti_nci60_names <- rep(NA, ncol(potti))
for(i1 in 1:ncol(potti)){
  temp <- which(nci60["36460_at", ] == potti["36460_at", i1])
  if(length(temp) == 1){
    potti_nci60_names[i1] <- colnames(nci60)[temp]
  }
}

drug_names <- unique(potti_info$drug_name)
contrast_levels <- c(0, 1)
for(i1 in 1:length(drug_names)){
  for(i2 in 1:length(contrast_levels)){
    cat("drug: ", drug_names[i1], " level ", contrast_levels[i2],
        "\n", 
        potti_nci60_names[ 
          (potti_info$drug_name == drug_names[i1]) & 
          (potti_info$contrast_group == contrast_levels[i2])], "\n",
        fill = TRUE)
  }
}
```

    ## drug:  Adria  level  0 
    ##  SF-539_A SNB-75_A MDA-MB-435_A NCI-H23_A M14_A 
    ## MALME-3M_A SK-MEL-2_A SK-MEL-28_A SK-MEL-5_A UACC-62_A 
    ## 
    ## drug:  Adria  level  1 
    ##  NCI/ADR-RES_A HCT-15_A HT29_A EKVX_A NCI-H322M_A 
    ## IGROV1_A OVCAR-3_A OVCAR-4_A OVCAR-5_A OVCAR-8_A SK-OV-3_A CAKI-1_A 
    ## 
    ## drug:  Doce  level  0 
    ##  NA NA NA NA NA NA NA NA NA NA 
    ## 
    ## drug:  Doce  level  1 
    ##  NA NA NA NA NA NA NA NA NA NA 
    ## 
    ## drug:  Etopo  level  0 
    ##  SF-539_A BT-549_A MDA-MB-231/ATCC_A NCI/ADR-RES_A 
    ## HOP-62_A NCI-H226_A SK-MEL-28_A UACC-257_A 786-0_A 
    ## 
    ## drug:  Etopo  level  1 
    ##  MCF7_A HCC-2998_A HCT-15_A SW-620_A NCI-H322M_A 
    ## PC-3_A OVCAR-4_A OVCAR-5_A 
    ## 
    ## drug:  5-FU  level  0 
    ##  MCF7_A COLO 205_A HCT-116_A NCI-H460_A LOX IMVI_A 
    ## SK-MEL-5_A A498_A UO-31_A 
    ## 
    ## drug:  5-FU  level  1 
    ##  NCI/ADR-RES_A MDA-MB-435_A SW-620_A EKVX_A M14_A 
    ## OVCAR-8_A SN12C_A 
    ## 
    ## drug:  Cytox  level  0 
    ##  NA NA NA NA NA NA NA NA NA NA 
    ## 
    ## drug:  Cytox  level  1 
    ##  NA NA NA NA NA NA NA NA NA NA 
    ## 
    ## drug:  Topo  level  0 
    ##  SF-539_A SNB-75_A U251_A HS 578T_A HOP-62_A 
    ## NCI-H226_A NCI-H23_A LOX IMVI_A OVCAR-8_A A498_A ACHN_A CAKI-1_A UO-31_A 
    ## 
    ## drug:  Topo  level  1 
    ##  K-562_A RPMI-8226_A MDA-MB-435_A SK-MEL-5_A 
    ## HCC-2998_A HCT-116_A HCT-15_A NCI-H322M_A SK-MEL-28_A COLO 205_A 
    ## 
    ## drug:  Taxol  level  0 
    ##  SF-295_A SF-539_A HS 578T_A MDA-MB-435_A 
    ## COLO 205_A HCC-2998_A HT29_A OVCAR-3_A DU-145_A 
    ## 
    ## drug:  Taxol  level  1 
    ##  CCRF-CEM_A SW-620_A A549/ATCC_A EKVX_A MALME-3M_A 
    ## SK-MEL-28_A OVCAR-8_A 786-0_A

We get perfect matches for all cell lines used for 5 of the 7 drugs. We get no matches at all for two of drugs: Doce (Docetaxel) and Cytox (Cytoxan). Failing for all cell lines for these drugs suggests something more involved might be going on there.

All of the matches we do see are for the "A" set of cell line replicates, suggesting the way they dealt with having triplicate measurements of each line was to simply take the first in each instance.

Check Sample Correlations in the Potti et al Data
=================================================

``` r
potti_cors <- cor(potti)
```

Let's try plotting the pairwise correlations for all columns in the Potti et al dataset.

``` r
levelplot(potti_cors,
          xlab = "", ylab = "",
          scales = list(alternating=0, tck=0),
          main = "Pairwise Sample Correlations\nin Potti et al Data")
```

![](/Users/kabaggerly/Professional/Talks/2018/2018_05_31_Moffitt/Git/moffitt_2018_rr_short_course/results/04_check_sample_matches_and_corrs_files/figure-markdown_github/plot_raw_potti_cors-1.png)

The columns for the two drugs (Docetaxel and Cytoxan) stand out clearly. They're effectively uncorrelated with the columns of the NCI60 data matrix.

We can also look for sample reuse by checking where the pairwise correlations are quite high (say above 0.99).

``` r
levelplot(potti_cors > 0.99,
          xlab = "", ylab = "",
          scales = list(alternating=0, tck=0),
          main = "Pairwise Sample Correlations > 0.99\nin Potti et al Data")
```

![](/Users/kabaggerly/Professional/Talks/2018/2018_05_31_Moffitt/Git/moffitt_2018_rr_short_course/results/04_check_sample_matches_and_corrs_files/figure-markdown_github/plot_potti_high_cors-1.png)

As it happens, there are a few off-diagonal contiguous stretches of "hits", and these are coming from the Docetaxel and Cytoxan columns. We can see this by looking at the values for the first probeset.

``` r
## Doce 0, Cytox 1

potti["36460_at", (potti_info$drug_name == "Doce") &
        (potti_info$contrast_group == 0)]
```

    ## Doce_0_023 Doce_0_024 Doce_0_025 Doce_0_026 Doce_0_027 Doce_0_028 
    ##    79.9480    81.3032   100.7450    73.0558    24.7488    57.2766 
    ## Doce_0_029 Doce_0_030 Doce_0_031 Doce_0_032 
    ##    36.8400    42.7512    24.8366    80.3417

``` r
potti["36460_at", (potti_info$drug_name == "Cytox") &
        (potti_info$contrast_group == 1)]
```

    ## Cytox_1_085 Cytox_1_086 Cytox_1_087 Cytox_1_088 Cytox_1_089 Cytox_1_090 
    ##     79.9480     81.3032    100.7450     73.0558     24.7488     57.2766 
    ## Cytox_1_091 Cytox_1_092 Cytox_1_093 Cytox_1_094 
    ##     36.8400     42.7512     24.8366     80.3417

``` r
## Doce 1, Cytox 0

potti["36460_at", (potti_info$drug_name == "Doce") &
        (potti_info$contrast_group == 0)]
```

    ## Doce_0_023 Doce_0_024 Doce_0_025 Doce_0_026 Doce_0_027 Doce_0_028 
    ##    79.9480    81.3032   100.7450    73.0558    24.7488    57.2766 
    ## Doce_0_029 Doce_0_030 Doce_0_031 Doce_0_032 
    ##    36.8400    42.7512    24.8366    80.3417

``` r
potti["36460_at", (potti_info$drug_name == "Cytox") &
        (potti_info$contrast_group == 1)]
```

    ## Cytox_1_085 Cytox_1_086 Cytox_1_087 Cytox_1_088 Cytox_1_089 Cytox_1_090 
    ##     79.9480     81.3032    100.7450     73.0558     24.7488     57.2766 
    ## Cytox_1_091 Cytox_1_092 Cytox_1_093 Cytox_1_094 
    ##     36.8400     42.7512     24.8366     80.3417

The Docetaxel and Cytoxan columns are identical but for the reversal of the contrast group labels. This makes little sense biologically.

Check Other Summaries for Potti et al Data
==========================================

Given the discrepancy, we also checked whether other summary statistics applied to the Potti et al data might clarify the situation. We looked at values typically returned by `summary` (most not shown), of which the minimum was most enlightening.

``` r
plot(apply(potti, 2, min),
     ylab = "Minimum Expression Value",
     xlab = "Raw Column Index",
     main = "Minimum Expression for Potti et al Data")

drug_boundaries <- 
  which(potti_info$drug_name[1:(ncol(potti) - 1)] != 
          potti_info$drug_name[2:ncol(potti)])
abline(v = drug_boundaries + 0.5)

drug_centers <- 
  (c(1, drug_boundaries+1) + c(drug_boundaries, ncol(potti))) / 2
drug_center_names <- 
  potti_info$drug_name[floor(drug_centers)]
mtext(drug_center_names, side = 3, 
      at = drug_centers)
```

![](/Users/kabaggerly/Professional/Talks/2018/2018_05_31_Moffitt/Git/moffitt_2018_rr_short_course/results/04_check_sample_matches_and_corrs_files/figure-markdown_github/plot_potti_minimums-1.png)

The minimums for the Doce and Cytox samples all look to be about the same. Let's look at the actual values.

``` r
apply(potti[, potti_info$drug_name == "Doce"], 2, min)
```

    ## Doce_0_023 Doce_0_024 Doce_0_025 Doce_0_026 Doce_0_027 Doce_0_028 
    ##    5.89822    5.89822    5.89822    5.89822    5.89822    5.89822 
    ## Doce_0_029 Doce_0_030 Doce_0_031 Doce_0_032 Doce_1_033 Doce_1_034 
    ##    5.89822    5.89822    5.89822    5.89822    5.89822    5.89822 
    ## Doce_1_035 Doce_1_036 Doce_1_037 Doce_1_038 Doce_1_039 Doce_1_040 
    ##    5.89822    5.89822    5.89822    5.89822    5.89822    5.89822 
    ## Doce_1_041 Doce_1_042 
    ##    5.89822    5.89822

``` r
apply(potti[, potti_info$drug_name == "Cytox"], 2, min)
```

    ## Cytox_0_075 Cytox_0_076 Cytox_0_077 Cytox_0_078 Cytox_0_079 Cytox_0_080 
    ##     5.89822     5.89822     5.89822     5.89822     5.89822     5.89822 
    ## Cytox_0_081 Cytox_0_082 Cytox_0_083 Cytox_0_084 Cytox_1_085 Cytox_1_086 
    ##     5.89822     5.89822     5.89822     5.89822     5.89822     5.89822 
    ## Cytox_1_087 Cytox_1_088 Cytox_1_089 Cytox_1_090 Cytox_1_091 Cytox_1_092 
    ##     5.89822     5.89822     5.89822     5.89822     5.89822     5.89822 
    ## Cytox_1_093 Cytox_1_094 
    ##     5.89822     5.89822

All of the minimum values are 5.89822. We've seen this value before. This is the minimum value for all of the samples in the datasets from GEO, which are not cell lines, but rather patient samples. Since this was one of the cohorts used to test their predictor for docetaxel, this may be the test data as opposed to the training data. It's not at all clear why these should be present for cytoxan, nor is it clear why the correlations beteen these samples and the cell lines. would all be near zero.

Matching Potti to GEO
=====================

We can explore the hypothesis that some of the Potti data might be coming from the GEO data by looking for matching values for a few entries. We'll start with the first probeset in the Potti data, 36460\_at.

``` r
potti["36460_at", "Doce_0_023"]
```

    ## [1] 79.948

``` r
which(geo == potti["36460_at", "Doce_0_023"], arr.ind = TRUE)
```

    ##          row col
    ## 31307_at  68  14

``` r
which(geo == potti["36460_at", "Doce_0_024"], arr.ind = TRUE)
```

    ##          row col
    ## 31307_at  68  15

``` r
which(geo == potti["36460_at", "Doce_0_025"], arr.ind = TRUE)
```

    ##             row col
    ## 39191_at   9479   1
    ## 31307_at     68  16
    ## 34027_f_at 2074  21

It appears row 1 of the Potti et al columns is mapping to row 68 of the GEO data. While row 68 of the GEO data is the first non-control probe row, the probeset id is *not* the same; the probeset ordering is different between the NCI60 data and the GEO data.

We can check this further by looking at the next few rows for some of the samples that match.

``` r
colnames(geo)[14]
```

    ## [1] "GSM4903"

``` r
potti[1:5, "Doce_0_023"]
```

    ## 36460_at 36461_at 36462_at 36463_at 36464_at 
    ##  79.9480 105.1320  12.8672  66.4968 329.8380

``` r
geo[67 + c(1:5), "GSM4903"]
```

    ##   31307_at   31308_at 31309_r_at   31310_at   31311_at 
    ##    79.9480   105.1320    12.8672    66.4968   329.8380

The values match exactly, but since the probesets are in a different order, the values have essentially been randomly scrambled with respect to each other. Random scrambling would explain why we correlations near zero. We can check how extensive the offset is for this sample (for which there are 12558 data rows in the Potti et al table):

``` r
sum(potti[, "Doce_0_023"] == geo[67 + c(1:nrow(potti)), "GSM4903"])
```

    ## [1] 12535

``` r
which(potti[, "Doce_0_023"] != geo[67 + c(1:nrow(potti)), "GSM4903"])
```

    ##   36435_at   36436_at 36437_s_at   36438_at   36439_at   36440_at 
    ##      12534      12535      12536      12537      12538      12539 
    ##   36441_at 36442_g_at   36443_at 36444_s_at   36445_at 36446_s_at 
    ##      12540      12541      12542      12543      12544      12545 
    ##   36448_at 36449_s_at   36450_at   36451_at   36452_at   36453_at 
    ##      12547      12548      12549      12550      12551      12552 
    ##   36454_at   36455_at   36456_at   36457_at   36459_at 
    ##      12553      12554      12555      12556      12558

``` r
all( sort(potti[12534:nrow(potti), "Doce_0_023"]) == 
       sort(geo[67 + c(12534:nrow(potti)), "GSM4903"]) )
```

    ## [1] TRUE

We match all but the last 23 of the 12558 entries using the offset, and the remaining entries involve the same values, just scrambled in a different way. This is where this column in the Potti data comes from.

Now, there are just 20 columns in the Potti data labeled as belonging to Docetaxel (and Cytoxan), but there are 24 columns of data in the GEO dataset.

We're simply going to try brute force matching here.

``` r
n_doce <- 20
n_geo <- 24

match_count_matrix <- 
  matrix(0, nrow = n_doce, ncol = n_geo)
rownames(match_count_matrix) <- 
  colnames(potti)[potti_info$drug_name == "Doce"]
colnames(match_count_matrix) <- 
  colnames(geo)

for(i1 in 1:n_doce){
  for(i2 in 1:n_geo){
    match_count_matrix[i1, i2] <- 
      sum(potti[, rownames(match_count_matrix)[i1]] ==
            geo[67 + c(1:nrow(potti)), colnames(match_count_matrix)[i2]])
  }
}

table(match_count_matrix[match_count_matrix > 1000])
```

    ## 
    ## 12535 
    ##    20

We have matches of the same extent for all 20 of the Doce columns.
Let's eyeball roughly where these are.

``` r
levelplot(match_count_matrix,
          main = "N Matches of Potti Doce Columns\nand GEO Data Columns")
```

![](/Users/kabaggerly/Professional/Talks/2018/2018_05_31_Moffitt/Git/moffitt_2018_rr_short_course/results/04_check_sample_matches_and_corrs_files/figure-markdown_github/plot_potti_doce_geo_matches-1.png)

These are essentially in linear order, so the sensitive and resistant patient groups map pretty well to the 0/1 contrast group labels.

Let's collect all of the name matches into a vector. We can add these for Cytoxan too, given that the Doce and Cytox values are the same except for the contrast group label.

``` r
potti_geo_matches <- rep(NA, ncol(potti))
names(potti_geo_matches) <- colnames(potti)
for(i1 in 1:nrow(match_count_matrix)){
  potti_geo_matches[rownames(match_count_matrix)[i1]] <-
    colnames(match_count_matrix)[which(match_count_matrix[i1, ] > 1000)]
}

potti_geo_matches[(potti_info$drug_name == "Cytox") &
                    (potti_info$contrast_group == 0)] <-
  potti_geo_matches[(potti_info$drug_name == "Doce") &
                      (potti_info$contrast_group == 1)]

potti_geo_matches[(potti_info$drug_name == "Cytox") &
                    (potti_info$contrast_group == 1)] <-
  potti_geo_matches[(potti_info$drug_name == "Doce") &
                      (potti_info$contrast_group == 0)]
```

Similarly, let's grab a vector defining the sample status for the GEO samples used.

``` r
geo_status <- rep(NA, ncol(potti))
geo_status[!is.na(potti_geo_matches)] <- 
  geo_sample_info[potti_geo_matches[!is.na(potti_geo_matches)],
                  "sample_status"]
```

Combine all Match Info Into a Data Frame
========================================

We now know where all of the columns in the original data matrix came from. Here, we combine this into a single data frame for later use.

``` r
potti_matches <- potti_info
potti_matches$nci60_match <- potti_nci60_names
potti_matches$geo_match <- potti_geo_matches
potti_matches$geo_status <- geo_status
```

Save the Results
================

Now we save the results for later use, both as an RData object and as a csv file.

``` r
save(potti_matches,
     file = here("results", "potti_matches.RData"))
write.csv(potti_matches,
          file = here("results", "potti_matches.csv"))
```
