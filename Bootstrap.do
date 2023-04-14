1.	Bootstrapping
**#Bootstrap BT
matrix wiu = (r(mean))

capture program drop BTboot
program define BTboot, rclass
preserve
bsample 
stset time, failure(death)
	hcost idnum cost, l(`r(max)') method(0) 
matrix bt = e(b)
capture drop bt1 bt2
svmat bt
capture drop btest
		egen btest=min(bt1)
		sum btest	
	return scalar wiu_1 = r(mean)
	restore
end

capture drop BootBT
bootstrap BootBT = r(wiu_1), rep(1000) saving(BootAll, replace): BTboot


**#Bootstrap Simple OLS

matrix NOLSMmean = (r(mean))

capture program drop NaiveOLSboot
program define NaiveOLSboot, rclass
preserve
		bsample 
		reg cost time if death==1
		capture drop naiveOLS
		predict naiveOLS, xb
		sum naiveOLS if death==1 
		return scalar NOLSMmean_1 = r(mean)
restore
end

capture drop BootNaiveOLS
bootstrap BootNaiveOLS = r(NOLSMmean_1), rep(1000) saving(BootNOLS, replace): NaiveOLSboot
estat bootstrap, all
eststo 


**#Bootstrap Simple Log

matrix LogGLMmean = (r(mean))

capture program drop NaiveLogboot
program define NaiveLogboot, rclass
preserve
	bsample 
			count if death ==1
			matrix observations = (r(N))
			capture drop observations
			svmat observations
			capture drop obsN
			egen obsN=min(observations)
		regress lnCost time if death==1
		capture drop LnOLSlog
		predict LnOLSlog, xb 
		capture drop residlnCost
		predict residlnCost if death==1, residuals 
		
		matrix betaln=e(b)
		svmat betaln

		capture drop b0ln
		egen b0ln=min(betaln2)
		replace b0ln=. if death==0

		capture drop b1ln
		egen b1ln=min(betaln1)
		replace b1ln=. if death==0

		capture drop eresidlnCost
		generate eresidlnCost=exp(residlnCost) if death==1
		capture drop seresidlnCost
		egen seresidlnCost=sum(eresidlnCost)
		capture drop smearlnCost
		generate smearlnCost=(1/obsN)*seresidlnCost if death==1

		capture drop OLSlog
		generate OLSlog=(exp(b0ln+b1ln*LnOLSlog))*smearlnCost
		summarize OLSlog if death==1
	
		return scalar LogGLMmean_1 = r(mean)
restore		
end

capture drop BootNaivelog
bootstrap BootNaivelog = r(LogGLMmean_1), rep(1000) saving(BootOLSLog, replace): NaiveLogboot
estat bootstrap, all
eststo

**# Bootstrap Simple GLM 

matrix GLMmean = (r(mean))

capture program drop GLMboot
program define GLMboot, rclass
preserve
	bsample 
		glm cost time if death==1, link(log) family(gamma)		
		capture drop GLMbootfit2
		predict GLMbootfit2 if death==1
		sum GLMbootfit2 if death==1
		return scalar GLMmean_1 = r(mean)
restore		
end

bootstrap BootGLMSimple = r(GLMmean_1), rep(1000) saving(BootGLM, replace): GLMboot
estat bootstrap, all
eststo


**#Bootstrap Carides Simple OLS

matrix CarNaive = (r(mean))

capture program drop CarNaiveOLS
program define CarNaiveOLS, rclass
preserve
	bsample 
		regress cost time if death==1
		*keep the coefficients
		capture drop beta
		matrix betaa=e(b)
		capture drop betaa
		svmat betaa
		capture drop b0
		egen b0=min(betaa2)
		replace b0=. if death==0
		capture drop b1
		egen b1=min(betaa1)
		replace b1=. if death==0
		**2n stage: compute the mean costs from the betas and the mean survival //times mean survival time from KM estimator
		*mean survival with KM
		stset time, failure(death)
		stci, emean
			matrix Kaplan = (r(emean))
			capture drop Kaplan
			svmat Kaplan
			capture drop KMmean
			egen KMmean=min(Kaplan)
		summarize time
		stci, emean 
		capture drop meanClinear
		generate meanClinear=b0+b1*KMmean
		summarize meanClinear
		return scalar CarNaive_1 = r(mean)
restore		
end
capture drop BootCarNaiveOLSS
bootstrap BootCarNaiveOLSS = r(CarNaive_1), rep(1000) saving(BootGLM, replace): CarNaiveOLS
estat bootstrap, all
eststo

**#Bootstrap Carides Log

matrix LogCar = (r(mean))

capture program drop CarLogboot
program define CarLogboot, rclass
	preserve
		bsample 
		count if death ==1
		matrix observations = (r(N))
			capture drop observations
			svmat observations
			capture drop obsN
			egen obsN=min(observations)
			regress lnCost time if death==1
			capture drop LnOLSlog
			predict LnOLSlog, xb 
			capture drop residlnCost
			predict residlnCost if death==1, residuals 	
			matrix betaln=e(b)
			svmat betaln
			capture drop b0ln
			egen b0ln=min(betaln2)
			replace b0ln=. if death==0
			capture drop b1ln
			egen b1ln=min(betaln1)
			replace b1ln=. if death==0
			capture drop eresidlnCost
			generate eresidlnCost=exp(residlnCost) if death==1
			capture drop seresidlnCost
			egen seresidlnCost=sum(eresidlnCost)		
			capture drop smearlnCost
			generate smearlnCost=(1/obsN)*seresidlnCost if death==1	
			stset time, failure(death)
			stci, emean
			matrix Kaplan = (r(emean))
			capture drop Kaplan
			svmat Kaplan
			capture drop KMmean
			egen KMmean=min(Kaplan)
			capture drop meanClnCostsm
			generate meanClnCostsm=(exp(b0ln+b1ln*KMmean))*smearlnCost
			summarize meanClnCostsm 
			return scalar LogCar_1 = r(mean)
	restore		
end

capture drop BootCarlog
bootstrap BootCarlog = r(LogCar_1), rep(1000) saving(BootCarLog, replace): CarLogboot
estat bootstrap, all
eststo


**# Bootstrap Carides with GLM

matrix CGLMmean = (r(mean))

capture program drop CGLMboot
program define CGLMboot, rclass
preserve
	bsample 
		glm cost time if death==1, link(log) family(gamma)
		matrix betacglm=e(b)
		svmat betacglm
		capture drop b01
		egen b01=min(betacglm2)
		replace b01=. if death==0
		capture drop b11
		egen b11=min(betacglm1)
		replace b11=. if death==0	
		stset time, failure(death)
		stci, emean
			matrix Kaplan = (r(emean))
			capture drop Kaplan
			svmat Kaplan
			capture drop KMmean
			egen KMmean=min(Kaplan)
		capture drop meanCGLM4
		generate meanCGLM4=exp(b01+b11*KMmean) 
		summarize meanCGLM4
		return scalar CGLMmean_1 = r(mean)
restore		
end

capture drop BootCaridesGLM
bootstrap BootCaridesGLM = r(CGLMmean_1), rep(100) saving(BootCarGLM, replace): CGLMboot
estat bootstrap, all
eststo


**# Compare bootstrapped estimators 

use BootAll.dta
capture merge 1:1 _n using BootNOLS , update
capture merge 1:1 _n using BootOLSLog , update
capture merge 1:1 _n using BootGLM , update
capture merge 1:1 _n using CarNaiveOLS , update
capture merge 1:1 _n using BootCarLog , update
capture merge 1:1 _n using BootCarGLM , update

summarize
