#' ---
#' title: "A correlation analysis of clinical variables of TCGA-KIRC patients"
#' output: 
#'   html_document: 
#'     default
#'   github_document: 
#'     df_print: paged
#'   pdf_document:
#'     latex_engine: xelatex
#' knit: (function(inputFile, encoding) {
#'   rmarkdown::render(inputFile, encoding = encoding, output_format = "all") })     
#' ---
#' 
#' This project contains a pipeline of clinical analysis of the Cancer Genome Atlas Kidney Renal Clear Cell Carcinoma (TCGA-KIRC) data of patients from [Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-KIRC%22%5D%7D%7D%5D%7D) and [cBioPortal](https://www.cbioportal.org/study/summary?id=kirp_tcga).
#' 
#' Previously, we presented [an exploratory preprocessing analysis](1.preprocessing.md). In this section, we present a correlation analysis with t-test and ANOVA test to investigate significative distinctions between clinical variables according to their vital status. 
#' 
#' 

#' 
#' 
## ----message=FALSE, warning=FALSE, paged.print=FALSE, echo = FALSE------------
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
#' ## 2. Taming data 
#' 
## -----------------------------------------------------------------------------
kirc_clinic <- kirc_clinic %>%
  mutate_if(is.character, as.factor) %>%
  mutate(patient_id = as.character(patient_id))

#' 
#' ## 3. Checking categorical variables
#' 
#' check frequency, lables and levels 
#' 
## -----------------------------------------------------------------------------
kirc_clinic %>%
  select_if(is.factor) %>%
  summary() 

# agregating levels
kirc_clinic <- kirc_clinic %>%
  mutate(tumor_stg = fct_collapse(tumor_stg,
                                  T1 = c('T1', 'T1a', 'T1b'),
                                  T2 = c('T2', 'T2a', 'T2b'),
                                  T3 = c('T3', 'T3a', 'T3b', 'T3c')))

kirc_clinic <- kirc_clinic %>%
  mutate(prior_cancer = fct_collapse(prior_cancer, 
                                     Yes = c('Yes', 'Yes, History of Prior Malignancy', 'Yes, History of Synchronous/Bilateral Malignancy')))

kirc_clinic <- kirc_clinic %>%
  mutate(sex = fct_collapse(sex, Male = c('MALE', 'Male')))

kirc_clinic <- kirc_clinic %>%
  mutate(tissue_site = fct_collapse(tissue_site,
                                    A = c('A3', 'AK', 'AS'),
                                    B = c('B0', 'B2', 'B4', 'B8', 'BP'),
                                    C = c('CJ', 'CW', 'CZ'),
                                    G = c('G6', 'GK'),
                                    M = c('MM', 'MW')))

# droping levels
kirc_clinic <- kirc_clinic %>%
  mutate(race = fct_recode(race, NULL = 'ASIAN'))

kirc_clinic <- kirc_clinic %>%
  mutate(tissue_site = fct_recode(tissue_site, NULL = '3Z', NULL='6D', NULL='DV', NULL='EU', NULL='G', NULL='M', NULL='T7'))

#' 
#' ## 4. Checking variables
#' 
## -----------------------------------------------------------------------------
glimpse(kirc_clinic)
skim(kirc_clinic) 
#View(kirc_clinic)

#' 
#' ## 5. Numeric variables vs. over_surv_stt
#' graphic visualization and t-test
#' 
## -----------------------------------------------------------------------------
# PATRICK: codigo para analizar todas as variaveis numericas?
kirc_clinic %>%
  select_if(is.numeric) %>%
  summary()

ggplot(kirc_clinic, aes(age, fill= over_surv_stt)) +
  geom_histogram(bins = 15, position = "dodge")
t.test(kirc_clinic$age ~ kirc_clinic$over_surv_stt) 

ggplot(kirc_clinic, aes(x=over_surv_stt, y=disease_free_mth)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
t.test(kirc_clinic$disease_free_mth ~ kirc_clinic$over_surv_stt) 

ggplot(kirc_clinic, aes(x=over_surv_stt, y=frac_genome_alter)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
t.test(kirc_clinic$frac_genome_alter ~ kirc_clinic$over_surv_stt)

ggplot(kirc_clinic, aes(year_diagnose, fill= over_surv_stt)) +
  geom_histogram(bins = 15, position = "dodge")
t.test(kirc_clinic$year_diagnose ~ kirc_clinic$over_surv_stt) 

ggplot(kirc_clinic, aes(x=over_surv_stt, y=long_dim)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
t.test(kirc_clinic$long_dim ~ kirc_clinic$over_surv_stt)

ggplot(kirc_clinic, aes(x=over_surv_stt, y=mutation_cnt)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
t.test(kirc_clinic$mutation_cnt ~ kirc_clinic$over_surv_stt)

ggplot(kirc_clinic, aes(x=over_surv_stt, y=over_surv_mth)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
t.test(kirc_clinic$over_surv_mth ~ kirc_clinic$over_surv_stt)

ggplot(kirc_clinic, aes(x=over_surv_stt, y=short_dim)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
t.test(kirc_clinic$short_dim ~ kirc_clinic$over_surv_stt)

ggplot(kirc_clinic, aes(x=over_surv_stt, y=second_long_dim)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
t.test(kirc_clinic$second_long_dim ~ kirc_clinic$over_surv_stt)

# fazer uma table com as variaveis dependentes, indpendente e p-valores

#' 
#' ## 4. Categorical variables vs. over_surv_stt
#' 
#' Tabulation and chi-square test
#' 
## -----------------------------------------------------------------------------
# talvez isso possa sair uma vez que ja tem a mesma analise com tablefit
kirc_clinic %>%
  select_if(is.factor) %>%
  summary() 

t_metas_stg <- table(kirc_clinic$metastasis_stg, kirc_clinic$over_surv_stt, exclude = NULL)
t_metas_stg <- addmargins(round(100*prop.table(t_metas_stg)))
t_metas_stg
chisq.test(x = kirc_clinic$metastasis_stg, y = kirc_clinic$over_surv_stt) 

t_lynph <- table(kirc_clinic$neoplasm_ln_stg, kirc_clinic$over_surv_stt, exclude = NULL)
t_lynph <- addmargins(round(100*prop.table(t_lynph)))
t_lynph
chisq.test(x = kirc_clinic$neoplasm_ln_stg, y = kirc_clinic$over_surv_stt) 

t_neop <- table(kirc_clinic$neoplasm_stg, kirc_clinic$over_surv_stt, exclude = NULL)
t_neop <- addmargins(round(100*prop.table(t_neop)))
t_neop
chisq.test(x = kirc_clinic$neoplasm_stg, y = kirc_clinic$over_surv_stt) 

t_tumor <- table(kirc_clinic$tumor_stg, kirc_clinic$over_surv_stt, exclude = NULL)
t_tumor <- addmargins(round(100*prop.table(t_tumor)))
t_tumor
chisq.test(x = kirc_clinic$tumor_stg, y = kirc_clinic$over_surv_stt) 

t_free <- table(kirc_clinic$disease_free_stt, kirc_clinic$over_surv_stt, exclude = NULL)
t_free <- addmargins(round(100*prop.table(t_free)))
t_free
chisq.test(x = kirc_clinic$disease_free_stt, y = kirc_clinic$over_surv_stt) 

t_prior <- table(kirc_clinic$prior_cancer, kirc_clinic$over_surv_stt, exclude = NULL)
t_prior <- addmargins(round(100*prop.table(t_prior)))
t_prior
chisq.test(x = kirc_clinic$prior_cancer, y = kirc_clinic$over_surv_stt) 

t_neo <- table(kirc_clinic$neoadj_therapy, kirc_clinic$over_surv_stt, exclude = NULL)
t_neo <- addmargins(round(100*prop.table(t_neo)))
t_neo
chisq.test(x = kirc_clinic$neoadj_therapy, y = kirc_clinic$over_surv_stt) 

t_platelet <- table(kirc_clinic$platelet, kirc_clinic$over_surv_stt, exclude = NULL)
t_platelet <- addmargins(round(100*prop.table(t_platelet)))
t_platelet
chisq.test(x = kirc_clinic$platelet, y = kirc_clinic$over_surv_stt)

t_prospect <- table(kirc_clinic$tissue_prospect, kirc_clinic$over_surv_stt, exclude = NULL)
t_prospect <- addmargins(round(100*prop.table(t_prospect)))
t_prospect
chisq.test(x = kirc_clinic$tissue_prospect, y = kirc_clinic$over_surv_stt)

t_race <- table(kirc_clinic$race, kirc_clinic$over_surv_stt, exclude = NULL)
t_race <- addmargins(round(100*prop.table(t_race)))
t_race
chisq.test(x = kirc_clinic$race, y = kirc_clinic$over_surv_stt) 

t_retros <- table(kirc_clinic$tissue_retrospect, kirc_clinic$over_surv_stt, exclude = NULL)
t_retros <- addmargins(round(100*prop.table(t_retros)))
t_retros
chisq.test(x = kirc_clinic$tissue_retrospect, y = kirc_clinic$over_surv_stt)

t_ca <- table(kirc_clinic$serum_ca, kirc_clinic$over_surv_stt, exclude = NULL)
t_ca <- addmargins(round(100*prop.table(t_ca)))
t_ca
chisq.test(x = kirc_clinic$serum_ca, y = kirc_clinic$over_surv_stt) 

t_sex <- table(kirc_clinic$sex, kirc_clinic$over_surv_stt, exclude = NULL)
t_sex <- addmargins(round(100*prop.table(t_sex)))
t_sex
chisq.test(x = kirc_clinic$sex, y = kirc_clinic$over_surv_stt) 

t_site <- table(kirc_clinic$tissue_site, kirc_clinic$over_surv_stt, exclude = NULL)
t_site <- addmargins(round(100*prop.table(t_site)))
t_site
chisq.test(x = kirc_clinic$tissue_site, y = kirc_clinic$over_surv_stt) 

t_neop_st <- table(kirc_clinic$person_neoplasm_stt, kirc_clinic$over_surv_stt, exclude = NULL)
t_neop_st <- addmargins(round(100*prop.table(t_neop_st)))
t_neop_st
chisq.test(x = kirc_clinic$person_neoplasm_stt, y = kirc_clinic$over_surv_stt) 

t_wbc <- table(kirc_clinic$wbc, kirc_clinic$over_surv_stt, exclude = NULL)
t_wbc <- addmargins(round(100*prop.table(t_wbc)))
t_wbc
chisq.test(x = kirc_clinic$wbc, y = kirc_clinic$over_surv_stt) 

#' 
#' ## 7. FinalFit
#' 
#' summarise variables/factors by a categorical variable
#' 
## -----------------------------------------------------------------------------
explanatory <- names(kirc_clinic %>%
              select(-over_surv_stt) %>%
              select_if(is.factor))
dependent <-  'over_surv_stt'

table_fit <- kirc_clinic %>%
  summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)
table_fit
warnings()

#knitr::kable(table_fit, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

#' 
#' ## Further analysis
#' 
#' - [A logistic regression analysis](3.logistic_regression.md) of each clinical variable weight.
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 