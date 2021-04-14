//install packages
net install merlin, from(https://www.mjcrowther.co.uk/software/merlin)
net install stmixed, from(https://www.mjcrowther.co.uk/software/stmixed)

//load data
use /Users/michael/Desktop/tramme-ipd-rep-MJC/dat,clear

gen died = alive1_death2==2
drop if sui==0
stset sui, f(died)

//reference test model, fixed effects only, same baseline
merlin (sui ages mmrcs fev1pps, family(rp, df(3) fail(died)))

//model 1:
//random intercept
//proportional hazards
//random slopes
	
	timer clear
	timer on 1
	merlin (sui ages mmrcs fev1pps 						///
				ages#M1[cohort]@1 						///
				mmrcs#M2[cohort]@1 						///
				fev1pps#M3[cohort]@1 					///
				M4[cohort]@1							///
				, 										///
				family(rp, df(3) fail(died))			///
			)											///
			, covariance(unstructured)					///
			restartvalues(M1 1 M2 0.01 M3 1 M4 0.01) 	///
			adaptopts(log)
	timer off 1
	timer list
	//or use wrapper
	stmixed ages mmrcs fev1pps || cohort: ages mmrcs fev1pps, distribution(rp) df(3) covariance(unstructuted)
		
		
//model 2:
//stratified baseline
//proportional hazards
//random slopes

	//generate dummies for cohort membership
	tab cohort, gen(cstub)
	forvalues i = 1/22 {
		local base `base' cstub`i'#rcs(sui, df(3) orthog event)
	}
	merlin (sui ages mmrcs fev1pps 				///
				`base'							///
				ages#M1[cohort]@1 				///
				mmrcs#M2[cohort]@1 				///
				fev1pps#M3[cohort]@1 			///
				, 								///
				family(rp, df(3) fail(died))	///
			)									///
			, covariance(unstructured) 			///
			restartvalues(M1 1 M2 0.01 M3 1)
	
	//or use wrapper
	stmixed ages mmrcs fev1pps `base' || cohort: ages mmrcs fev1pps			///
			, noconstant distribution(rp) df(3) covariance(unstructuted)
			
			
//model 3:
//random intercept
//non-proportional hazards
//random slopes

	merlin (sui ages mmrcs fev1pps 						///
				ages#rcs(sui, df(3) orthog event)		///
				mmrcs#rcs(sui, df(3) orthog event)		///
				fev1pps#rcs(sui, df(3) orthog event)	///
				ages#M1[cohort]@1 						///
				mmrcs#M2[cohort]@1 						///
				fev1pps#M3[cohort]@1 					///
				M4[cohort]@1							///
				, 										///
				family(rp, df(3) fail(died))			///
			)											///
			, intpoints(3) cov(unstr) 					///
			restartvalues(M1 1 M2 0.01 M3 1 M4 0.1)
	
	//or use wrapper
	stmixed ages mmrcs fev1pps || cohort: ages mmrcs fev1pps			///
			, distribution(rp) df(3) tvc(ages mmrcs fev1pps) dftvc(3)	///
			 covariance(unstructuted)

			 
//model 4:
//stratified risks
//non-proportional hazards
//random slopes

	forvalues i = 1/22 {
		local base `base' cstub`i'#rcs(sui, df(3) orthog event)
	}
	merlin (sui ages mmrcs fev1pps 						///
				ages#rcs(sui, df(3) orthog event)		///
				mmrcs#rcs(sui, df(3) orthog event)		///
				fev1pps#rcs(sui, df(3) orthog event)	///
				`base'									///
				ages#M1[cohort]@1 						///
				mmrcs#M2[cohort]@1 						///
				fev1pps#M3[cohort]@1 					///
				, 										///
				family(rp, df(3) fail(died))			///
			)											///
			, intpoints(3) cov(unstr) 					///
			restartvalues(M1 1 M2 0.01 M3 1)
			
	//or use wrapper
	stmixed ages mmrcs fev1pps `base' || cohort: ages mmrcs fev1pps		///
			, noconstant distribution(rp) df(3) 						///
			tvc(ages mmrcs fev1pps) dftvc(3)							///
			covariance(unstructuted)
	
	
//model 5:
//stratified baseline
//proportional hazards
//fixed effects only

	forvalues i = 1/22 {
		local base `base' cstub`i'#rcs(sui, df(3) orthog event)
	}
	merlin (sui ages mmrcs fev1pps `base', family(rp, df(3) fail(died)))
	
	//or use wrapper
	stmerlin ages mmrcs fev1pps, dist(rp) df(3)

	
//model 6:
//stratified baseline risks
//non-proportional hazards
//fixed effects only

	forvalues i = 1/22 {
		local base `base' cstub`i'#rcs(sui, df(3) orthog event)
	}
	merlin (sui ages mmrcs fev1pps 						///
				ages#rcs(sui, df(3) orthog event)		///
				mmrcs#rcs(sui, df(3) orthog event)		///
				fev1pps#rcs(sui, df(3) orthog event)	///
				`base'									///
				, 										///
				family(rp, df(3) fail(died))			///
			)
	//or use wrapper
	forvalues i = 1/22 {
		local basevars `basevars' cstub`i'
	}
	stmerlin ages mmrcs fev1pps , dist(rp) df(3) 			///
			tvc(ages mmrcs fev1pps `basevars') dftvc(3)
			
