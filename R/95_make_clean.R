#' ---
#' title: Make Clean
#' author: Keith Baggerly
#' date: "`r Sys.Date()`"
#' output: github_document
#' ---

library(here)

dirs_to_clean <- c("results", "data")

for(i1 in 1:length(dirs_to_clean)){
  temp_file_list <- 
    dir(here(dirs_to_clean[i1]), recursive = TRUE)
  file.remove(here(dirs_to_clean[i1], temp_file_list))
}
