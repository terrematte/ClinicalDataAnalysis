#' ---
#' title: "A Preprocessing analysis of clinical data of TCGA-KIRC patients"
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
#' This project contains a pipeline for analysis of The Cancer Genome Atlas Kidney - Renal Clear Cell Carcinoma (TCGA-KIRC) clinical data, from [Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-KIRC%22%5D%7D%7D%5D%7D) and [cBioPortal](https://www.cbioportal.org/study/summary?id=kirp_tcga).
#' 
#' In this section, the initial preprocessing is applied to clean the data and arrange following the Tidyverse philosophy. Exploratory Data Analysis summarizes their main characteristics. 
#' 
#' 

#' 
#' 
## ----message=FALSE, warning=FALSE, echo = FALSE-------------------------------
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
## ----message=FALSE, warning=FALSE, paged.print=TRUE---------------------------
kirc_clin_raw <- read_delim("data/kirc_tcga_clinical_data.tsv", "\t", 
                            escape_double = FALSE, 
                            trim_ws = TRUE)

#' 
#' 
## ----echo=FALSE, message=FALSE, results='hide', paged.print=TRUE--------------
class(kirc_clin_raw) 
dim(kirc_clin_raw) 
names(kirc_clin_raw) 
glimpse(kirc_clin_raw)
skim(kirc_clin_raw) 
#View(kirc_clin_raw)

#' 
#' 
#' ## 2. Cleaning data
#' 
#' Select variables based on NA count (> 50% complete is a good choice!).
#' 
#' <!-- # TO DO @PATRICK: simplify code NA_sum? -->
#' <!-- # kirc_clean <- kirc_clin_raw %>% -->
#' <!-- #     summarise_all(~ sum(is.na(.)))  -->
#' 
## -----------------------------------------------------------------------------
NA_fifty <- dim(kirc_clin_raw)[1]/2

NA_sum <- colSums(is.na(kirc_clin_raw))
NA_sum <- as.data.frame(NA_sum)
NA_sum <- tibble::rownames_to_column(NA_sum, "variables")
NA_sum <- NA_sum %>%
     filter(NA_sum < NA_fifty)

kirc_clean <- kirc_clin_raw %>%
     select(one_of(NA_sum$variables))

#' 
#' Remove duplicate observations:
#' 
## -----------------------------------------------------------------------------
kirc_clean0 <- kirc_clean %>%
     distinct_at('Patient ID', .keep_all = TRUE)

#' 
#' Remove nuneric variables with unique observations:  
#' <!-- # TO DO @PATRICK: function to select variables with unique observations? -->
#' <!-- # kirc_cleanX <- kirc_clean1 %>% -->
#' <!-- #     summarise_if(is.numeric, ~ n=unique(.)) -->
#' 
## ----message=FALSE, warning=FALSE, paged.print=TRUE---------------------------
kirc_clean0 %>%
     select_if(is.numeric) %>%
     skim()

kirc_clean1 <-  kirc_clean0  %>%
     select(!c('Last Alive Less Initial Pathologic Diagnosis Date Calculated Day Value', 
               'Number of Samples Per Patient', 
               'Sample type id'))

#' 
#' Remove character variables with unique observations:
#' 
## ----message=FALSE, warning=FALSE, paged.print=TRUE---------------------------
kirc_clean1 %>%
     select_if(is.character) %>%
     skim()

#' 
## ----message=FALSE, warning=FALSE, paged.print=TRUE---------------------------
kirc_clean2 <- kirc_clean1  %>%
     select(!c('Study ID', 'Cancer Type', 'Cancer Type Detailed', 
               'Neoplasm Histologic Type Name', 'ICD-10 Classification', 
               'International Classification of Diseases for Oncology, Third Edition ICD-O-3 Site Code', 
               'Informed consent verified', 'Is FFPE', 'Oncotree Code', 'Sample Type', 'Tumor Tissue Site'))

#' 
#' Remove character variables with similar information - check each one!
#' 
## -----------------------------------------------------------------------------
kirc_clean2 %>%
     select_if(is.character) %>%
     skim()

#' 
## -----------------------------------------------------------------------------
table(kirc_clean2$`Overall Survival Status`, exclude = NULL)
table(kirc_clean2$`Patient's Vital Status`, exclude = NULL)

#' 
## -----------------------------------------------------------------------------
kirc_clean3 <- kirc_clean2  %>%
     select(!c('Sample ID', 'Other Patient ID', 'Other Sample ID', 'Pathology Report File Name', 'Pathology report uuid', "Patient's Vital Status"))

#' 
#' Remove other variables not directly related to patient - check each one!
#' 
## -----------------------------------------------------------------------------
kirc_clean2 %>%
     select_if(is.character) %>%
     skim()

kirc_clean4 <- kirc_clean3  %>%
     select(!c('Form completion date','International Classification of Diseases for Oncology, Third Edition ICD-O-3 Histology Code','Vial number'))

#' 
#' ## 3. Changing variables names
#' 
#' Using snake_style 
#' 
## -----------------------------------------------------------------------------
kirc_clean4 <- kirc_clean4 %>%
     rename(patient_id = 'Patient ID',
            age = 'Diagnosis Age',
            metastasis_stg = 'American Joint Committee on Cancer Metastasis Stage Code',
            neoplasm_ln_stg = 'Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code',
            neoplasm_stg = 'Neoplasm Disease Stage American Joint Committee on Cancer Code',
            tumor_stg = 'American Joint Committee on Cancer Tumor Stage Code',
            disease_free_mth = 'Disease Free (Months)',
            disease_free_stt = 'Disease Free Status',
            ethnicity = 'Ethnicity Category', 
            frac_genome_alter = 'Fraction Genome Altered',
            histology_grd = 'Neoplasm Histologic Grade',
            hemoglobin = 'Hemoglobin level',
            neoadj_therapy = 'Neoadjuvant Therapy Type Administered Prior To Resection Text',
            prior_cancer = 'Prior Cancer Diagnosis Occurence',
            year_diagnose = 'Year Cancer Initial Diagnosis',
            tumor_lateral = 'Primary Tumor Laterality',
            long_dim = 'Longest Dimension',
            primer_ln_ind3 = 'Primary Lymph Node Presentation Assessment Ind-3',
            mutation_cnt = 'Mutation Count',
            over_surv_mth = 'Overall Survival (Months)',
            over_surv_stt = 'Overall Survival Status',
            platelet = 'Platelet count',
            tissue_prospect = 'Tissue Prospective Collection Indicator',
            race = 'Race Category',
            tissue_retrospect = 'Tissue Retrospective Collection Indicator',
            serum_ca = 'Serum calcium level',
            gender = 'Sex',
            short_dim = 'Shortest Dimension',
            second_long_dim = 'Specimen Second Longest Dimension',
            tissue_site = 'Tissue Source Site',
            person_neoplasm_stt = 'Person Neoplasm Status',
            wbc = 'WBC')

#' 
#' ## 4. Taming data
#' 
#' Use lubridate for dates
#' 
## -----------------------------------------------------------------------------
kirc_clean4 <- kirc_clean4 %>%
     mutate_if(is.character, as.factor) %>%
     mutate(patient_id = as.character(patient_id))

#' 
#' ## 5. Checking NA patterns 
#' 
#' Check distincts types of NAs: MCAR, MAR, MNAR
#' 
## -----------------------------------------------------------------------------
kirc_clean4  %>%
     missing_plot()

missing_glimpse(kirc_clean4)

#' 
#' ## 6. Checking numeric variables
#' 
#' Check data distribution, plausible ranges, outliers;
#' Thinking about deleting outliers from dataset? Need to evaluate carefully each one!
#' 
#' <!-- # TO DO @PATRICK: codigo para analizar todas as variaveis numericas? -->
#' <!-- # kirc_clean6 <-  kirc_clean4 %>% -->
#' <!-- #      select_if(is.numeric) %>% -->
#' <!-- #      ggplot(aes(funs(.)) + -->
#' <!-- #      geom_boxplot(width = .5) + -->
#' <!-- #      geom_jitter(width = 0.05, alpha = 0.2, color = "orange")om_boxplot(width = .5) + -->
#' <!-- #      geom_jitter(width = 0.05, alpha = 0.2, color = "orange") -->
#' 
## -----------------------------------------------------------------------------
kirc_clean4 %>%
     select_if(is.numeric) %>%
     summary()

#' 
## -----------------------------------------------------------------------------
ggplot(kirc_clean4, aes(age)) +
     geom_histogram(bins = 20, alpha = 0.8, color = "red")

#' 
## -----------------------------------------------------------------------------
ggplot(kirc_clean4, aes(year_diagnose)) +
     geom_histogram(bins = 20, alpha = 0.8, color = "red")

#' 
## -----------------------------------------------------------------------------
ggplot(kirc_clean4, aes(x ='', y=disease_free_mth)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$disease_free_mth)
# filter(disease_free_mth >= 0) 

#' 
## -----------------------------------------------------------------------------
ggplot(kirc_clean4, aes(x ='', y=frac_genome_alter)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$frac_genome_alter)

#' 
## -----------------------------------------------------------------------------
ggplot(kirc_clean4, aes(x ='', y=long_dim)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$long_dim)

#' 
## -----------------------------------------------------------------------------
ggplot(kirc_clean4, aes(x ='', y=mutation_cnt)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$mutation_cnt)

#' 
## -----------------------------------------------------------------------------
ggplot(kirc_clean4, aes(x ='', y=over_surv_mth)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$over_surv_mth)

#' 
## -----------------------------------------------------------------------------
ggplot(kirc_clean4, aes(x ='', y=short_dim)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$short_dim)

#' 
## -----------------------------------------------------------------------------
ggplot(kirc_clean4, aes(x ='', y=second_long_dim)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$second_long_dim)

#' 
#' ## 7. Checking categorical variables
#' 
#' Check frequency, lables and levels 
#' 
## -----------------------------------------------------------------------------
kirc_clean4 %>%
     select_if(is.factor) %>%
     summary() 

# agregating levels
kirc_clinic <- kirc_clean4 %>%
     mutate(tumor_stg = fct_collapse(tumor_stg,
                             T1 = c('T1', 'T1a', 'T1b'),
                             T2 = c('T2', 'T2a', 'T2b'),
                             T3 = c('T3', 'T3a', 'T3b', 'T3c')))

kirc_clinic <- kirc_clinic %>%
     mutate(prior_cancer = fct_collapse(prior_cancer, 
               Yes = c('Yes', 'Yes, History of Prior Malignancy', 'Yes, History of Synchronous/Bilateral Malignancy')))

kirc_clinic <- kirc_clinic %>%
     mutate(gender = fct_collapse(gender, Male = c('MALE', 'Male')))
                                        
kirc_clinic <- kirc_clinic %>%
     mutate(tissue_site = fct_collapse(tissue_site,
                         A = c('A3', 'AK', 'AS'),
                         B = c('B0', 'B2', 'B4', 'B8', 'BP'),
                         C = c('CJ', 'CW', 'CZ'),
                         OTHERS = c('G6', 'GK', 'MM', 'MW',
                                    '3Z', '6D', 'DV', 'EU', 'T7')))

# droping levels ??? What about others?? 
# check chunk bellow!
kirc_clinic <- kirc_clinic %>%
     mutate(race = fct_recode(race, NULL = 'ASIAN'))

# kirc_clinic <- kirc_clinic %>%
#     mutate(race = fct_drop(race, only = 'ASIAN'))

# recoding levels ??
# 
# kirc_clinic <- kirc_clinic %>%
#      mutate(gender = fct_recode(gender, '1'='Male', '2'='Female'))
# 
# kirc_clinic <- kirc_clinic %>%
#      mutate(gender = if_else(gender %in% c('Male', 'Female'), 1, 0))


#' 
## -----------------------------------------------------------------------------
kirc_clinic %>%
     select_if(is.factor) %>%
     summary()

#' 
#' ## 8. Saving dataset
#' 
## -----------------------------------------------------------------------------
write_csv(kirc_clinic, path = "data/kirc_clinic.csv")

rm(kirc_clean4, kirc_clean3, kirc_clean2, kirc_clean1, kirc_clean0, kirc_clean, NA_sum, NA_fifty)

#' 
#' ## Further analysis
#' 
#' - [Comparison and Hyphotesis test](2.correlation.md) 
#' - [Logistic Regression Model](3.logistic_regression.md)
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
