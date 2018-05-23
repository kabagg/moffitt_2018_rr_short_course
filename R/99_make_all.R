#' ---
#' title: Make All
#' author: Keith Baggerly
#' date: "`r Sys.Date()`"
#' output: github_document
#' ---

library(here)
library(rmarkdown)

if(!dir.exists(here("results"))){
  dir.create(here("results"))
}

files_in_r_to_run <- 
  c("01_gather_data.R",
    "02_describing_raw_data.Rmd",
    "03_preprocessing_data_base_r_version.Rmd",
    "04_check_sample_matches_and_corrs.Rmd",
    "05_report_matches_to_potti_columns.Rmd")

for(i1 in 1:length(files_in_r_to_run)){
  
  rmarkdown::render(here("R", files_in_r_to_run[i1]),
                    output_format = 
                      github_document(html_preview = TRUE, toc = TRUE),
                    output_dir = here("results"))
  
}

rmarkdown::render(here("README.Rmd"),
                  output_format = 
                    github_document(html_preview = TRUE, toc = TRUE),
                  output_dir = here())
