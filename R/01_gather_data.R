#' ---
#' title: Gathering Raw Data
#' author: Keith Baggerly
#' date: "`r Sys.Date()`"
#' output: github_document
#' ---

#' Here we're assembling the raw data from the three
#' datasets we want to use (potti, nci60, and geo) 
#' from various repositories on the web, and collecting
#' them in data/.

library(here)
library(downloader)

#' Set up the data directory if it doesn't exist

if (!dir.exists(here("data"))) {
  cat("creating data/ directory\n")
  dir.create(here("data"))
}

#'
#' Location of the expression matrix from 
#' Potti et al on figshare
#'

potti_url <- 
  "https://ndownloader.figshare.com/files/10615624?private_link=66603862d770b4c73146"
potti_data_file <-
  "chemo.zip"

##download(potti_url, 
##         destfile = here("data", potti_data_file),
##         mode = "wb")

#' 
#' Location of the NCI60 Gene Expression Data from the NCI's 
#' Developmental Therapeutics Program (DTP) page
#' 

#' In general, the file is linked to from 
#' https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data.
#' We've used the direct link below, but I worry this may 
#' change over time. 

nci60_url <- 
  paste0("https://wiki.nci.nih.gov/download/attachments/",
         "155845004/WEB_DATA_NOVARTIS_ALL.zip?version=1&",
         "modificationDate=1378406329000&api=v2&download=true")
nci60_data_file <- "WEB_DATA_NOVARTIS_ALL.zip"

##if(!file.exists(here("data", nci60_data_file))){
##  cat("downloading ", nci60_data_file, "\n")
##  download(nci60_url, 
##           destfile = here("data", nci60_data_file),
##           mode = "wb")
##}

#' 
#' Locations of the soft files for GSE349 and GSE350
#' from the Gene Expression Omnibus (GEO)
#' 

gse349_soft_url <- 
  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE349/soft/"
gse349_soft_file <- "GSE349_family.soft.gz" 

gse350_soft_url <- 
  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE350/soft/"
gse350_soft_file <- "GSE350_family.soft.gz"

##download(paste0(gse349_soft_url, gse349_soft_file), 
##         destfile = here("data", gse349_soft_file),
##         mode = "wb")
##
##download(paste0(gse350_soft_url, gse350_soft_file), 
##         destfile = here("data", gse350_soft_file),
##         mode = "wb")

#' 
#' Loop through all of the above, and acquire
#' files we don't have yet
#' 

url_list <- 
  c(potti_url,
    nci60_url,
    paste0(gse349_soft_url, gse349_soft_file),
    paste0(gse350_soft_url, gse350_soft_file))

data_file_list <-
  c(potti_data_file,
    nci60_data_file,
    gse349_soft_file,
    gse350_soft_file)

for(i1 in 1:length(data_file_list)){
  
  if(!file.exists(here("data", data_file_list[i1]))){
    cat("downloading ", data_file_list[i1], "\n")

    download(url_list[i1], 
             destfile = here("data", data_file_list[i1]),
             mode = "wb")
  }
}

