# Estimating-Mean-Total-Costs-under-Moderate-Censoring


#### 📝 Read the paper [here](https://github.com/JonasWeinert/Estimating-Mean-Total-Costs-under-Moderate-Censoring/blob/main/Weinert_MeanCostUnderModerateCensoring.docx)
---

***In this MSC project, I benchmarked different parametric and non-parametric approaches to medical cost data survival analysis by***
- Implementing several models from the current academic discourse in stata (some self operationalised, w/o packages)
- Writing my own bootstrapping program to obtain comparable standard errors for model comparison.


## Abstract
###Background
To inform health technology reimbursement policy on which technology should be adapted, the costs associated with a technology in relation to its benefits to recipients must be estimated. However, the analysis of costs that are associated with a given treatment is more complicated as it faces several issues with the nature of cost data in clinical studies. The asymmetry of the distribution of costs, the presence of outliers, and censoring are challenges that need to be overcome in the statistical analysis of health care costs [11].
This study will test and compare three estimation approaches with a total of seven specifications on cost data that presents moderate levels of censoring, strong skewness and does not provide information on cost histories.

###Methods
Seven possible approaches (one nonparametric, three parametric complete case estimators and three two stage parametric estimators) to model the mean cost to time of death under censoring. Strengths and weaknesses of each approach will be illustrated based on the assumptions made by the respective models. Afterwards, measures to evaluate and compare the approaches will be introduced.

###Findings
Ignoring censored observations leads to an expected bias in the total cost estimate. The nonparametric and adjusted parametric models produce the lowest uncertainty around the estimate. The lack of cost histories and covariates poses the biggest obstacle to more robust estimations.



