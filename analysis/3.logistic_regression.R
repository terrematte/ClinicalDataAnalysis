#' ---
#' title: "A logistic regression analysis of TCGA-KIRC"
#' output: 
#'   html_document: 
#'     default
#'   github_document: 
#'     df_print: paged
#'     html_preview: FALSE
#'     keep_html: TRUE
#'   pdf_document:
#'     latex_engine: xelatex
#' knit: (function(inputFile, encoding) {
#'   rmarkdown::render(inputFile, encoding = encoding, output_format = "all") })    
#' ---
#' 
#' This project contains a pipeline of clinical analysis of the Cancer Genome Atlas Kidney Renal Clear Cell Carcinoma (TCGA-KIRC) data of patients from [Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-KIRC%22%5D%7D%7D%5D%7D) and [cBioPortal](https://www.cbioportal.org/study/summary?id=kirp_tcga).
#' 
#' Previously, we presented a [An exploratory preprocessing analysis](1.preprocessing.md), and [a correlation analysis](2.correlation.md).
#' 
#' In this final section, we present a logistic regression analysis of each clinical variable weight for TCGA-KIRC.
#' 

#' 
#' 
## ----message=FALSE, warning=FALSE, paged.print=FALSE--------------------------
# Set the packages of interest
packages = c("tidyverse","skimr","finalfit")

# if a package is installed, it will be loaded
# otherwise, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

suppressMessages(library("tidyverse"))
setwd(".")

#' 
#' 
#' ## 1. Importing data
#' 
## ----message=FALSE, warning=FALSE, paged.print=FALSE, echo = FALSE------------

kirc_clinic <- read_csv("data/kirc_clinic.csv")


#' 
#' ## TO DO
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
