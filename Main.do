**# Cleaning
label define deathlabel 0 "censored" 1 "uncensored", replace
label values death deathlabel
encode patientid, generate(idnum)
generate lnCost=ln(cost)

**#Declare the data as survival data by stset, 
stset time, failure(death)
//stci, emean 


**# Inspecting Censoring and Skewness
summarize
summarize cost, detail
summarize cost if death==1, detail
summarize cost if death==0
summarize time if death==1
histogram cost if death==1 & cost>0


browse

*plot  survival curve and  mean survival time.
sts graph
stci,rmean 

**# Data inspection
* plot expenditure by survival times likke huang
// censored
scatter cost time if death==1, ytitle(Total Cost) xtitle(Total Time until death) saving("scatter cost time censored")
// uncensored
scatter cost time if death==0, color(red) ytitle(Total Cost) xtitle(Total Time until censoring) saving("scatter cost time uncensored")
// in comparison
twoway (scatter lnCost time), by(death)
* inspect times
hist time if death==1
hist time if death==0
twoway (hist time), by(death)
** spikeplot censoring over studytime
xtile pctTime = time, n(100)
gen cens = 1 if death == 0
replace cens = 0 if death == 0
graph bar (mean) cens, over(pctTime)


**# 1a B&T Non parametrics estimation
hcost idnum cost, l(8360) method(0) 
eststo
** 8360 = max survival time

**# 2a Naive OLS //only complete cases
reg cost time if death==1
linktest
capture drop naiveOLS
predict naiveOLS if death==1, xb
capture drop naiveOLSresid
predict naiveOLSresid if death==1, resid
sum naiveOLS // Mean = mean cost on time
eststo



**# 2b OLS with Logged Costs + smear

regress lnCost time if death==1
linktest
estat hettest //Beusch Pagan Test for heteroskedasticity
estimates store ols_ident_gauss
predict LnOLSlog if death==1, xb
predict residlnCost if death==1, residuals 
pnorm residsmeared
//swilk residlnCost // determine whether smearing is necessary

**keep the betas
matrix betaln=e(b)
svmat betaln

egen b0ln=min(betaln2)
replace b0ln=. if death==0

egen b1ln=min(betaln1)
replace b1ln=. if death==0

**compute the smearing estimator
generate eresidlnCost=exp(residlnCost)
egen seresidlnCost=sum(eresidlnCost)
generate smearlnCost=(1/18576)*seresidlnCost if death==1

gen residsmeared  = residlnCost * smearlnCost

** retransform fitted values with smearing factor and obtain mean
generate OLSlog=(exp(b0ln+b1ln*LnOLSlog))*smearlnCost
summarize OLSlog if death==1

**Hosmer–Lemeshow test
xtile pctOLSLOG = OLSlog, n(10)
tab pctOLSLOG, gen(il)
global il il*
reg eresidlnCost $il
test $il

**# 2c GLM

** generate multiple GLMs with different link functions
* Log + Gamma
quietly glm cost time if death==1, link(log) family(gamma)
estimates store glm_log_gam
//predict glm_log_gam2, deviance
*pnorm glm_log_gam2
*graph export "/Users/jonas/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Courses/HP425/Summative/Graphical output/glm_log_gam.png", replace
*sqrt + gamma
quietly glm cost time if death==1, link(power .5) family(gamma)
estimates store glm_sqrt_gam
*predict glm_sqrt_gam, deviance
*pnorm glm_sqrt_gam
*graph export "/Users/jonas/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Courses/HP425/Summative/Graphical output/glm_sqrt_gam.png", replace
* log + gaussian
quietly glm cost time if death==1, link(log) family(gaussian)
estimates store glm_log_gau
*predict glm_log_gau, deviance
*pnorm glm_log_gau
*graph export "/Users/jonas/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Courses/HP425/Summative/Graphical output/glm_log_gau.png", replace
* log + inverted gaussian
quietly glm cost time if death==1, link(log) family(igaussian)
estimates store glm_log_igau
*predict glm_log_igau, deviance
*pnorm glm_log_igau
*graph export "/Users/jonas/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Courses/HP425/Summative/Graphical output/glm_log_igau.png", replace
* sqrt + gaussian
quietly glm cost time if death==1, link(power .5) family(gaussian)
estimates store glm_sqrt_gau
*predict glm_sqrt_gau, deviance
*pnorm glm_sqrt_gau
*graph export "/Users/jonas/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Courses/HP425/Summative/Graphical output/glm_sqrt_gau.png", replace
* log + poisson
quietly glm cost time if death==1, link(log) family(poisson) scale(x2) // x2 added to correct stadard errors
estimates store glm_log_poi
*predict glm_log_poi, deviance
*pnorm glm_log_poi
*graph export "/Users/jonas/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Courses/HP425/Summative/Graphical output/glm_log_poi.png", replace
* sqrt + poisson
quietly glm cost time if death==1, link(power .5) family(poisson) scale(x2)
estimates store glm_sqrt_poi
*predict glm_sqrt_poi, deviance
*pnorm glm_sqrt_poi
*graph export "/Users/jonas/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Courses/HP425/Summative/Graphical output/glm_sqrt_poi.png", replace

* OLS: Gaiss+identity
quietly glm lnCost time if death==1, link(identity) family(gaussian)
estimates store glm_gau_ident_ln



** compare models based on AIC & BIC
estimates stats *

** Running Log-Gamma GLM again to compute mean estimate//Based on AIC and BIC Log + gamma performs best

glm cost time if death==1, link(log) family(gamma)
linktest, link(log) family(gamma)
predict GLM if death==1
sum GLM
predict residGLM, deviance
//qnorm residGLM
//scatter GLM residGLM
//hist residGLM
// export graph

**Hosmer–Lemeshow test
xtile pctGLM = GLM, n(10)
tab pctGLM, gen(igam)
global igam igam*
reg residGLM $igam
test $igam



// alternative to get mean: predict double resid

**#2d Carides et al. // accounting for censoring through inverse probability weighting
*Carides with Naive OLS
**Stage 1: linear model - betas estimated though OLS on only uncensored costs observations
regress cost time if death==1


*keep the coefficients
matrix beta=e(b)
svmat beta

egen b0=min(beta2)
replace b0=. if death==0

egen b1=min(beta1)
replace b1=. if death==0


**2n stage: compute the mean costs from the betas and the mean survival //times mean survival time from KM estimator
*mean survival with KM
summarize time
stci, emean 
generate meanClinear=b0+b1*3369.838
summarize meanClinear

**#2eCarides - Log transformation. 



**1st stage: predictions, Betas and smearing estimator for the first stage are used from the previous Log Transformed Naive OLS

**2n stage - mean costs
stci,rmean 
generate meanClnCostsm=(exp(b0ln+b1ln*3369.838))*smearlnCost
summarize meanClnCostsm

**#2f Carides et al with Log Gamma GLMs
**Stage 1: GLM betas estimated through GLM with log link on gamma distribution
glm cost time if death==1, link(log) family(gamma)
matrix betacglm=e(b)
svmat betacglm

egen b01=min(betacglm2)
replace b01=. if death==0

egen b11=min(betacglm1)
replace b11=. if death==0

**Stage 2: ompute the mean costs from the betas and the mean survival //times mean survival time from KM estimator

*mean survival with KM
summarize time
stci, emean 
generate meanCGLM4=exp(b01+b11*3369.838) // since betas are obtained from Gamma regression, exponentiation is necessary to obtain additive effect of change in time


summarize meanCGLM4
