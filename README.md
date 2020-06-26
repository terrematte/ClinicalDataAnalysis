# Clinical Data Analysis

This project contains a pipeline of clinical analysis of the Cancer Genome Atlas Kidney Renal Clear Cell Carcinoma (TCGA-KIRC) data of patients from [Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-KIRC%22%5D%7D%7D%5D%7D) and [cBioPortal](https://www.cbioportal.org/study/summary?id=kirp_tcga).

- [x] [An exploratory preprocessing analysis](analysis/1.preprocessing.md) of clinical variables in R with [tidyverse](https://www.tidyverse.org/), [skimr](https://github.com/ropensci/skimr) and [finalfit](https://github.com/ewenharrison/finalfit) packages.
- [ ] [A correlation analysis](analysis/2.correlation.md) with t-test and ANOVA checking significant distinction  between variables acording their vital status.
- [ ] [A logistic regression analysis](analysis/3.logistic_regression.md) of each clinical variable weight.


#### Authors :busts_in_silhouette:

 :bust_in_silhouette: [Beatriz Stransky](https://github.com/bia-stransky) [(Associate Professor, Phd, UFRN)](http://lattes.cnpq.br/3142264445097872)
 
 :bust_in_silhouette: [Patrick Terrematte](https://github.com/terrematte) [(Assistant Professor, UFERSA; PhD Student, UFRN)](http://lattes.cnpq.br/4283045850342312)
 
 :bust_in_silhouette: [Josivan Justino](https://github.com/Josivan-br) [(Assistant Professor, UNIR; PhD Student, UFRN)](http://lattes.cnpq.br/6470296449367089)







