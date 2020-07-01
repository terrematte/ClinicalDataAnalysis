# Clinical Data Analysis

This project contains a pipeline for analysis of The Cancer Genome Atlas Kidney - Renal Clear Cell Carcinoma (TCGA-KIRC) clincal data, from [Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-KIRC%22%5D%7D%7D%5D%7D) and [cBioPortal](https://www.cbioportal.org/study/summary?id=kirp_tcga).

- [ ] [Preprocessing and Exploratory Analysis](analysis/1.preprocessing.md). Initial preprocessing is applied to clean the data and arrange following the Tdyverse philosophy. Than,  Exploratory Data Analysis summarize their main characteristics.
- [ ] [Comparisons and Hypothesis test](analysis/2.correlation.md). Chi-squared test is applied to compare two or more proportions of categorical variables and T-student test to compare the means of numeric ones regardind the levels of 'Overall_Survival_Status'. The Hypoteis test is performed and p-value indicates the strength of evidence in supportting the null hypothesis.
- [ ] [Logistic Regression Model](analysis/3.logistic_regression.md). The logistic regression is apllied to investigate the relation between a set of variables and a binary outcome, i.e, the 'Overall_Survival_Status'.


#### Authors :busts_in_silhouette:

 :bust_in_silhouette: [Beatriz Stransky](https://github.com/bia-stransky) [(Associate Professor, Phd, UFRN)](http://lattes.cnpq.br/3142264445097872)
 
 :bust_in_silhouette: [Patrick Terrematte](https://github.com/terrematte) [(Assistant Professor, UFERSA; PhD Student, UFRN)](http://lattes.cnpq.br/4283045850342312)
 
 :bust_in_silhouette: [Josivan Justino](https://github.com/Josivan-br) [(Assistant Professor, UNIR; PhD Student, UFRN)](http://lattes.cnpq.br/6470296449367089)







