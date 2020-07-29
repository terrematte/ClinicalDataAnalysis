########## TCGA Clinical Data Analysis ##########
# logistic regression model for "over_surv_stt" 
# using Tidyverse whenever possible


install.packages("tidyverse", dependencies = TRUE)
library(tidyverse)

# setwd('.')


## 1. Importing data ---------------------

kirc_glm <- read_csv("kirc_glm.csv")


## 2. Taming data -------------------------
# use lubridate for dates
kirc_glm <- kirc_glm %>%
     mutate_if(is.character, as.factor) %>%
     mutate(patient_id = as.character(patient_id),
            age = as.integer(age),
            year_diagnose = as.integer(year_diagnose))

glimpse(kirc_glm)
View(kirc_glm)

# Option: make one variable per column, rather than refer to the data set every time.


## 3. Simple Logistic Regression -------------

# R organise values (called levels) of categorical variables alphabetically, by default.
contrasts(kirc_glm$over_surv_stt)

# Changing the order of levels 
kirc_glm$over_surv_stt <- relevel(kirc_glm$over_surv_stt, ref = "LIVING") 

# The simplest model is one with no predictors whatsoever in it. 
# This is called the empty or null model and is rarely useful. 
# It has the assumption that everyone has the same odds of having the outcome - to die. 

m_nul <- glm(kirc_glm$over_surv_stt ~ 1, family=binomial(link=logit))
summary(m_nul)

# The intercept coefficient means the the log odds of dying
# To know the odds of die, just exponentiate it
exp(-0.7099) # = 0.4916934

# To interpret odds as probabilities, just divide the odds by 1 plus the odds.
0.49 / 1.49  # = 0.3288591, or 33%.

# Check how R has interpreted the binary outcome 'over_surv_stt'
table(m$y)
table(kirc_glm$over_surv_stt)


# 4. Fitting GLM to a continuous variable -------------

# GLM assumes a linear relation of continuous variables with the outcome, i.e. the log odds of dying. 
# Check this premisse plotting one against the other.

# create a cross tabulation of age and sruvival status  
surv_by_age <- table(kirc_glm$age, kirc_glm$over_surv_stt) 

# output the frequencies of sruvival status by age 
freq_table <- prop.table(surv_by_age, margin = 1) 

# calculate the odds of dying 
odds <- freq_table[, "DECEASED"]/freq_table[, "LIVING"] 

# calculate the log odds 
logodds_die <- log(odds) 

# plot the ages found in the sample against the log odds of dying
plot(rownames(freq_table), logodds_die, xlab = 'age')

m1 <- glm(kirc_glm$over_surv_stt ~ kirc_glm$age, family=binomial(link=logit))
summary(m1)

# Log odds of dying = intercept + (coefficient for gender) * age
#                   = -3.04 + 0.04 * age
# when you are born (0 age) you have 5% of probablility to dye. 
# when you are 50 years old, you have 25%.

## odds_die vs. frac_genome_alter ( => OUT)
surv_by_frac <- table(kirc_glm$frac_genome_alter, kirc_glm$over_surv_stt) 
freq_table <- prop.table(surv_by_frac, margin = 1) 
odds <- freq_table[, "DECEASED"]/freq_table[, "LIVING"] 
logodds_die <- log(odds) 
plot(rownames(freq_table), logodds_die, xlab = 'frac_genome_alter')

# odds_die vs year_diagnose  ( => Ok)
surv_by_year <- table(kirc_glm$year_diagnose, kirc_glm$over_surv_stt) 
freq_table <- prop.table(surv_by_year, margin = 1) 
odds <- freq_table[, "DECEASED"]/freq_table[, "LIVING"] 
logodds_die <- log(odds) 
plot(rownames(freq_table), logodds_die, xlab = 'year_diagnose')

# odds_die vs. long_dim  ( => Ok)
surv_by_long <- table(kirc_glm$long_dim, kirc_glm$over_surv_stt) 
freq_table <- prop.table(surv_by_long, margin = 1) 
odds <- freq_table[, "DECEASED"]/freq_table[, "LIVING"] 
logodds_die <- log(odds) 
plot(rownames(freq_table), logodds_die, xlab = 'long_dim')

# odds_die vs. mutation_cnt ( => OUT)
surv_by_mut <- table(log(kirc_glm$mutation_cnt), kirc_glm$over_surv_stt) 
freq_table <- prop.table(surv_by_mut, margin = 1) 
odds <- freq_table[, "DECEASED"]/freq_table[, "LIVING"] 
logodds_die <- log(odds) 
plot(rownames(freq_table), logodds_die, xlab = 'mutation_cnt')

# odds_die vs. short_dim ( => OUT, colinear with long_dim, fewer observations)
surv_by_short <- table(kirc_glm$short_dim, kirc_glm$over_surv_stt) 
freq_table <- prop.table(surv_by_short, margin = 1) 
odds <- freq_table[, "DECEASED"]/freq_table[, "LIVING"] 
logodds_die <- log(odds) 
plot(rownames(freq_table), logodds_die, xlab = 'short_dim')

# odds_die vs. second_long_dim ( => OUT, colinear with long_dim, fewer observations)
surv_by_second <- table(kirc_glm$second_long_dim, kirc_glm$over_surv_stt) 
freq_table <- prop.table(surv_by_second , margin = 1) 
odds <- freq_table[, "DECEASED"]/freq_table[, "LIVING"] 
logodds_die <- log(odds) 
plot(rownames(freq_table), logodds_die, xlab = 'second_long_dim')


# MULTIPLE LOGISTIC (numeric)
# Removed: over_surv_mth, disease_free_mth
multi_nun <- glm(kirc_glm$over_surv_stt ~ kirc_glm$age + kirc_glm$year_diagnose + kirc_glm$long_dim, family=binomial(link=logit))
summary(multi_nun)

# Coefficients:
#      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            635.126452  98.948477   6.419 1.37e-10 ***
#      kirc_glm$age             0.039580   0.008889   4.453 8.48e-06 ***
#      kirc_glm$year_diagnose  -0.318564   0.049338  -6.457 1.07e-10 ***
#      kirc_glm$long_dim        0.371952   0.153008   2.431   0.0151 *  
#      ---
#     
# AIC: 565.05


## 5. Fitting GLM to a categorical variable -------------

# Fit a logistic regression with “gender” as the predictor variable.
m2 <- glm(kirc_glm$over_surv_stt ~ kirc_glm$gender, family=binomial(link=logit))
summary(m2)

# Log odds of dying = intercept + (coefficient for gender) * gendertype
#                   = -0.662 + 0.075 * be a male
contrasts(kirc_glm$gender)

# The odds of surviving being female compared male (to two decimal places)
m2$coefficients
round(exp(m2$coefficients), digits = 2)

# MULTIPLE LOGISTIC (character)
# Neoplasm stage agregates information from metastasis, linphnode and tumor stages (TMN).
# Risc factors: calcium serum high* (MSKCC and IMDC), hemoglobin low (MSKCC and IMDC), white blod cells (neutrophils) high (IMDC), platelet high (IMDC).

predict_char <- kirc_glm %>% 
     select_if(is.factor) %>%
     select(-c(over_surv_stt, metastasis_stg, lymph_stg, tumor_stg, disease_free_stt, person_neoplasm_stt, ethnicity)) %>%
     names()

multi_char <-glm(over_surv_stt ~.,
                 data = kirc_glm[ , c(predict_char, "over_surv_stt")],
                 family = "binomial")
summary(multi_char)

# neoplasm_stgStage III-IV  -1.1144     0.2922  -3.814 0.000137 ***
# histology_grdG3-G4        -0.5071     0.2935  -1.728 0.083965 .  
# hemoglobinLow             -0.8747     0.3157  -2.771 0.005589 ** 
# neoadj_therapyYes         -1.5878     0.7165  -2.216 0.026692 *  
# plateletNormal             1.1110     0.5204   2.135 0.032782 * 
#
# 7 variables not significant
# AIC: 381.42


# 6. Full model ----------
# Excluding all numerical variables presenting non-linear relations with outcome, or colinearities
# Removing one variable at time, selecting the ones with higher p-value

predict <- kirc_glm %>% 
     select(-c(over_surv_stt, metastasis_stg, lymph_stg, tumor_stg, disease_free_stt, person_neoplasm_stt, ethnicity, tissue_site, prior_cancer, gender, race, platelet, histology_grd, tumor_lateral, serum_ca,wbc )) %>%
     select(-c(over_surv_mth, disease_free_mth, frac_genome_alter, mutation_cnt, short_dim, second_long_dim, long_dim)) %>%
     select(-patient_id)  %>%
     names()

m_full <-glm(over_surv_stt ~.,
            data = kirc_glm[ , c(predict, "over_surv_stt")],
            family = "binomial")
summary(m_full)

# Even the model has converged, inspect the standart error of odds radtios and the sizes of odds ratios thenselves.  
# Std. Error > 1 (>10?) is not aceptable.
# Std. Error intercept ???


# Null Deviance and Residual Deviance --------------------
# analyse table of deviance 
anova(m_full, test = "Chisq") 

# Resid. Dev shows the deviances of the models compared with the saturated model.
# Adding each new variable to the model explains the data better, reducing the deviance.
# The p-value indicates that the corresponding added variable causes a significant change in deviance, and thus is a better fitting model.


# 7. Testing ------------------
# The deviance is given by default, other tests have to be requested in R.

# McFadden’s r-squared -------------------------
# design your logistic regression: mfull 
# run a null model: m0 
R2 <- 1 - logLik(m_full)/logLik(m_null) 
R2 


# c-statistic ------------------------
install.packages("DescTools") 
require(DescTools) 

# design your logistic regression: mfull 
# generate the c-statistic 
Cstat(m_full) 
