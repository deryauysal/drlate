

/*Ado-file for the LATE estimators*/
/* DR-LATE project!*/
/*Version: 09082022*/
/*Author: DU
ThisVersion: IPWRA LATE & LATT
RA LATE & LATT
NIPW LATE & LATT
AIPW -LATE
NAIPW - LATE
*Comments from 22.08.2022 Tymon
(3) Standardization done!  standardize only if not binary  
(5) corrected

AIPW - LATT 
NAIPW -LATT done 09/03/2023


To do:
Testing
add cbps? 
*/

capture program drop late

program late, eclass 

version 17 

syntax  namelist(min=4) [if] [, xvarsy(varlist) xvarst(varlist) xvarsips(varlist) omodel(string) tmodel(string) zmodel(string) vce(string) samplew(varname) stat(string)]
	

	
	local orig `0'
	tokenize `namelist'
	local method `1'
	macro shift
	local yvar `1'  //outcome
	macro shift
	local tvar `1'  //treatment
	macro shift 
	local zvar `1'	 //instrument	
	macro shift
	local preds `*'		//predetermined variables (if models use same set of covariates, otherwise empty)

	
	

		// Give values to options that were not set
	if "`xvarsy'" == "" {
		local xvarsy `preds'
	}
	if "`xvarst'" == "" {
		local xvarst `preds'
	}
	if "`xvarsips'" == "" {
		local xvarsips `preds'
	}
	if "`omodel'" == "" {
		local omodel regress
	}
	if "`omodel'" == "linear" {
		local omodel regress
	}
	if "`tmodel'" == "" {
		local tmodel logit
	}
	if "`tmodel'" == "linear" {
		local tmodel regress
	}
	if "`zmodel'" == "" {
		local zmodel logit
	}
	if "`vce'" == "" {
		local vce robust
	}
	if "`samplew'" == "" {
		tempvar samplew 
		gen `samplew' = 1
	}
	if "`vce'" == "" {
		local vce robust
	}
		
	if "`stat'" == "" {
		local stat late
	}	
		marksample touse			
	quietly {
/*
		tempvar troot 
		tempname valuest
		qui tab `tvar' if `touse', gen(`troot') matrow(`valuest')
		if rowsof(`valuest') != 2 {
			noi di as error "The treatment variable must only take 2 values in the sample."
			error 450
			exit
		}
*/
	//standardize the variables only if they are taking more than two values
	//Outcome Model
		foreach x in `xvarsy' {
		tempvar  `x'_y_st
		tempvar troot 
		tempname valuest
		qui tab `x' if `touse', gen(`troot') matrow(`valuest')
		if rowsof(`valuest') != 2 {
		egen  ``x'_y_st' = std(`x')
		}
		else {
		gen ``x'_y_st' = `x'	
		}
	} 

	local xvarsy_STD
	foreach item in `xvarsy' {
	local xvarsy_STD `xvarsy_STD' ``item'_y_st'

	}
	//	display "`xvarsy_STD'"
		local xvarsy "`xvarsy_STD'"
	//	display "`xvarsy'"
	
	//Treatment Model	
		foreach x in `xvarst' {
		tempvar  `x'_t_st
		tempvar troot 
		tempname valuest
		qui tab `x' if `touse', gen(`troot') matrow(`valuest')
		if rowsof(`valuest') != 2 {
		egen  ``x'_t_st' = std(`x')
		}
		else {
			gen  ``x'_t_st' = `x'
		}
	} 

	local xvarst_STD
	foreach item in `xvarst' {
	local xvarst_STD `xvarst_STD' ``item'_t_st'

	}
	//	display "`xvarst_STD'"
		local xvarst "`xvarst_STD'"
	//	display "`xvarst'"		
		
	//IPS Model	
		foreach x in `xvarsips' {
		tempvar  `x'_ips_st
				tempvar troot 
		tempname valuest
		qui tab `x' if `touse', gen(`troot') matrow(`valuest')
		if rowsof(`valuest') != 2 {
		egen  ``x'_ips_st' = std(`x')
		}
		else {
		gen  ``x'_ips_st' = `x'			
		}
	} 

	local xvarsips_STD
	foreach item in `xvarsips' {
	local xvarsips_STD `xvarsips_STD' ``item'_ips_st'

	}
		local xvarsips "`xvarsips_STD'"
		
				
		
		
		
		tempvar troot 
		tempname valuest
		qui tab `tvar' if `touse', gen(`troot') matrow(`valuest')
		if rowsof(`valuest') != 2 {
			noi di as error "The treatment variable must only take 2 values in the sample."
			error 450
			exit
		}
	    quietly summarize `tvar' if `touse'
	    if r(min)!=0 | r(max)!=1 {
		noi display as error "Treatment must only take on values zero or one"
		error 450
		exit
	    }
		
		tempvar zroot 
		tempname valuesz
		qui tab `zvar' if `touse', gen(`zroot') matrow(`valuesz')
		if rowsof(`valuesz') != 2 {
			noi di as error "The instrumental variable must only take 2 values in the sample."
			error 450
			exit
		}			
	    quietly summarize `zvar' if `touse'
	    if r(min)!=0 | r(max)!=1 {
		noi display as error "Instrument must only take on values zero or one"
		error 450
		exit
	    }	
		
if "`omodel'" == "logit" {
		tempvar yroot 
		tempname valuesy
		qui tab `yvar' if `touse', gen(`yroot') matrow(`valuesy')
		if rowsof(`valuesy') != 2 {
			noi di as error "The outcome variable must only take 2 values in the sample."
			error 450
			exit
		}
	    quietly summarize `yvar' if `touse'
	    if r(min)!=0 | r(max)!=1 {
		noi display as error "Outcome must only take on values zero or one"
		error 450
		exit
	    }
	
}
if "`omodel'" == "poisson" {

	    quietly summarize `yvar' if `touse'
	    if r(min) < 0 {
		noi display as error "Outcome must be greater than or equal to zero"
		error 450
		exit
	    }
	
}		
tempname dmeanz1 dmeanz0
//Check the compliance:
	    quietly summarize `tvar' if `zvar' == 1 & `touse'==1
		scalar `dmeanz1' = r(mean)
		quietly summarize `tvar' if `zvar'  == 0 & `touse'==1
		scalar `dmeanz0' = r(mean)
//Check this! 
	tempvar y1hat y0hat d1hat d0hat ips ipsxb  wr1 wr0 w0s w1 w0 w1s invw0 invw1 omega1 omega0 term1 term2 term3 term4 
	tempname bips by1 by0 bd1 bd0 denom1s denom0s num1s num0s nums  denoms lates

	if "`method'" == "ipwra" {
	//IPS
set maxiter 300
if "`zmodel'" == "logit" {
  `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix `bips' = e(b)
 qui predict double `ips'
 }
else if "`zmodel'" == "cbps" {
	if "`stat'" == "late" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew')
  }
  	else if "`stat'" == "latt" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew') stat(att)
  }
   qui matrix `bips' = e(b)
   qui predict double `ipsxb'
   qui ge `ips' = logistic(`ipsxb')
  
}	


  


//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar' `xvarsy' if `zvar'==1 & `touse'==1 [pw = `samplew'/`ips']
qui matrix `by1' = e(b)
qui predict double `y1hat'

 
 
 qui `omodel' `yvar' `xvarsy' if `zvar'==0 & `touse'==1 [pw = `samplew'/(1-`ips')]
qui matrix `by0' = e(b)
qui predict double `y0hat'

//Conditional mean of d1 & d0#


	    if `dmeanz1'==1 {
        ge double `d1hat'= 1
	    }
		else {

qui `tmodel' `tvar' `xvarst' if `zvar'==1  & `touse'==1  [pw =  `samplew'/`ips']
matrix `bd1' = e(b)
qui predict double `d1hat'
}
	    
	    if `dmeanz0'==0 {
        ge double `d0hat' = 0
	    }
		else {
qui `tmodel' `tvar' `xvarst' if `zvar'==0  & `touse'==1 [pw =  `samplew'/(1-`ips')]
matrix `bd0' = e(b)
qui predict double `d0hat'
		}
		
		
qui sum `d1hat' if `touse'==1 [iw = `samplew']
matrix `denom1s' = r(mean)
qui sum `d0hat'  if `touse'==1 [iw = `samplew']
matrix `denom0s' = r(mean)
qui sum `y1hat' if `touse'==1  [iw = `samplew']
matrix `num1s' = r(mean)
qui sum `y0hat' if `touse'==1 [iw = `samplew']
matrix `num0s' = r(mean)
matrix `nums' = `num1s'-`num0s'
matrix `denoms' = `denom1s'-`denom0s'
matrix `lates' =`nums'*invsym(`denoms')

//drop y1hat y0hat d1hat d0hat ips
	
//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((exp({zhat:})/(1+exp({zhat:})))*(1/(1+exp({zhat:})))))

local eqips_inst instruments(eqips: `xvarsips')
}


if "`tmodel'" == "logit" {
local eqd1 (eqd1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`tvar' - exp({d1: `xvarst' _cons})/(1+exp({d1:})))))
           
local eqd1_inst instruments(eqd1: `xvarst')

local eqd0 (eqd0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`tvar'-exp({d0: `xvarst' _cons})/(1+exp({d0:}))))) 

local eqd0_inst instruments(eqd0: `xvarst' )
}


if "`tmodel'" =="regress" {
local eqd1 (eqd1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`tvar' - ({d1: `xvarst' _cons}))))
           
local eqd1_inst instruments(eqd1: `xvarst')

local eqd0 (eqd0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`tvar'-({d0: `xvarst' _cons})))) 

local eqd0_inst instruments(eqd0: `xvarst' )
}

if "`tmodel'" == "logit" {
	if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}- ((exp({d1:})/(1+exp({d1:})))))
	}
	else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-((exp({d0:})/(1+exp({d0:}))))))
	}
	else {
local denom (denom: {denom}- (((exp({d1:})/(1+exp({d1:}))))-((exp({d0:})/(1+exp({d0:}))))))
	}
}
else if "`tmodel'" =="regress"  {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}- ({d1:}))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-({d0:})))
		}
		else {
local denom (denom: {denom}- (({d1:})-({d0:})))
		}
}



if "`omodel'" == "regress" {
local eqy0 (eqy0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`yvar'-({y0: `xvarsy'  _cons})))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`yvar'-({y1: `xvarsy'  _cons})))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num (num: {num}-(({y1:})-({y0:})))
}
else if "`omodel'"  == "logit" {
local eqy0 (eqy0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`yvar'-exp({y0: `xvarsy' _cons})/(1+exp({y0:}))))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`yvar'-exp({y1: `xvarsy' _cons})/(1+exp({y1:}))))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num (num: {num}-(((exp({y1: })/(1+exp({y1:}))))-((exp({y0: })/(1+exp({y0:}))))))

}
else if "`omodel'"  == "poisson"  {
local eqy0 (eqy0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`yvar'-exp({y0: `xvarsy'  _cons})))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`yvar'-exp({y1: `xvarsy'  _cons})))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num (num: {num}-((exp({y1: }))-(exp({y0: }))))

}

local late (late: ({late} - {num}/{denom}))
} 

	

	/*Estimation of LATT*/
if "`stat'" == "latt" {
	
	if `dmeanz1'==1 {
		ge double `d1hat'= 1
	}
	else {
		
		qui regress  `tvar'  if `zvar' == 1 & `touse'==1 [pw = `samplew']
		matrix `bd1' = e(b)
		qui predict double `d1hat'
	}
	
	if `dmeanz0'==0 {
		ge double `d0hat' = 0
	}
	else {
		qui `tmodel' `tvar' `xvarst'  if `zvar' == 0 & `touse'==1 [pw = `samplew'*(`ips'/(1-`ips'))]
		matrix `bd0' = e(b)
		qui predict double `d0hat'
	}
	
	qui regress `yvar' if `zvar' == 1 & `touse'==1 [pw = `samplew']
	matrix `by1' = e(b)
	qui predict `y1hat'
	
	qui `omodel' `yvar' `xvarsy'  if `zvar' == 0 & `touse'==1 [pw = `samplew'*(`ips'/(1-`ips'))]
	matrix `by0' = e(b)
	qui predict double `y0hat'
	
	
	qui sum `d1hat'  if `zvar' == 1 [iw = `samplew']
	matrix `denom1s' = r(mean)
	qui sum `d0hat'  if `zvar' == 1 [iw = `samplew']
	matrix `denom0s' = r(mean)
	qui sum `y1hat'  if `zvar' == 1 [iw = `samplew']
	matrix `num1s' = r(mean)
	qui sum `y0hat'  if `zvar' == 1 [iw = `samplew']
	matrix `num0s' = r(mean)
	matrix `nums' = `num1s'-`num0s'
	matrix `denoms' = `denom1s'-`denom0s'
	
	
	matrix `lates' =`nums'*invsym(`denoms')
	
	//drop y1hat y0hat d1hat d0hat ips
	
	
	
	//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - ((exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*(1-`zvar')))

local eqips_inst instruments(eqips: `xvarsips')
}

	


	local eqd1 (eqd1: ((1.`zvar'*(`tvar' - ({d1: _cons})))))
	
	local eqd1_inst instruments(eqd1: )
	
	
	if "`tmodel'" == "logit" {
		
		local eqd0 (eqd0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`tvar'-exp({d0: `xvarst' _cons})/(1+exp({d0:}))))) 
		
		local eqd0_inst instruments(eqd0: `xvarst' )
	}
	
	
	if "`tmodel'" =="regress" {
		
		local eqd0 (eqd0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`tvar'-({d0: `xvarst' _cons})))) 
		
		local eqd0_inst instruments(eqd0: `xvarst' )
	}
	

	
	
	if "`tmodel'" == "logit" {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-(exp({d0: })/(1+exp({d0:})))))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}-(({d1: })-(exp({d0: })/(1+exp({d0:})))))))
		}
	}
	else if "`tmodel'" =="regress"  {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-({d0: })))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}-(({d1: })-({d0: })))))
		}
	}
	

	
	local eqy1 (eqy1: (1.`zvar'*(`yvar'-({y1:  _cons})))) 
	
	local eqy1_inst instruments(eqy1:  )
	
	
	
	if "`omodel'" == "regress" {

		local eqy0 (eqy0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`yvar'-({y0: `xvarsy'  _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )
		
		
		local num (num: (1.`zvar'*({num}-(({y1: })-({y0: })))))
		
	}
	else if "`omodel'"  == "logit" {
		
		
		local eqy0 (eqy0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`yvar'-exp({y0: `xvarsy' _cons})/(1+exp({y0:}))))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )
		
		local num (num: (1.`zvar'*({num}-(({y1: })-(exp({y0: })/(1+exp({y0:})))))))
		
		
	}
	else if "`omodel'"  == "poisson"  {
		local eqy0 (eqy0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`yvar'-exp({y0: `xvarsy' _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )
		
		local num (num: (1.`zvar'*({num}-(({y1: })-(exp({y0: }))))))
		
		
	}
	
	local late (late: ({late} - {num}/{denom}))
}


if "`stat'" == "late" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	 
if "`stat'" == "latt" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm   `eqips'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom'  if `touse' [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}
	}
	else if "`method'" == "ipw" {
			//IPS
set maxiter 300

if "`zmodel'" == "logit" {
  `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix `bips' = e(b)
 qui predict double `ips'
 }
else if "`zmodel'" == "cbps" {
		if "`stat'" == "late" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew')
  }
  	else if "`stat'" == "latt" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew') stat(att)
  }
   qui matrix `bips' = e(b)
   qui predict double `ipsxb'
   qui ge `ips' = logistic(`ipsxb')
  
}	
  

		
//Conditional mean of y1 & y0#
if "`stat'" == "late" {

qui gen `y1hat' = ((`zvar')/`ips')*`yvar' if `touse'==1

qui gen `y0hat' = (((1-`zvar'))/(1-`ips'))*`yvar' if `touse'==1



//Conditional mean of d1 & d0#


	    if `dmeanz1'==1 {
        ge double `d1hat'= 1
	    }
		else {
qui gen `d1hat' = ((`zvar')/`ips')*`tvar' if `touse'==1

}
	    
	    if `dmeanz0'==0 {
        ge double `d0hat' = 0
	    }
		else {
		
qui gen `d0hat' = (((1-`zvar'))/(1-`ips'))*`tvar' if `touse'==1


		}
		
		
qui sum `d1hat' if `touse'==1 [iw = `samplew']
matrix `denom1s' = r(mean)
qui sum `d0hat'  if `touse'==1 [iw = `samplew']
matrix `denom0s' = r(mean)
qui sum `y1hat' if `touse'==1  [iw = `samplew']
matrix `num1s' = r(mean)
qui sum `y0hat' if `touse'==1 [iw = `samplew']
matrix `num0s' = r(mean)
matrix `nums' = `num1s'-`num0s'
matrix `denoms' = `denom1s'-`denom0s'
matrix `lates' =`nums'*invsym(`denoms')

//drop y1hat y0hat d1hat d0hat ips
	
//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((exp({zhat:})/(1+exp({zhat:})))*(1/(1+exp({zhat:})))))

local eqips_inst instruments(eqips: `xvarsips')
}



local eqd1 (eqd1: ((`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*`tvar') - ({d1:  _cons}))
local eqd0 (eqd0: (((1-`zvar')/(1-(exp({zhat:})/(1+exp({zhat:})))))*`tvar')-({d0:  _cons})) 



		if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}- ({d1:}))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-({d0:})))
		}
		else {
local denom (denom: {denom}- (({d1:})-({d0:})))
		}




local eqy0 (eqy0: (((1-`zvar')/(1-(exp({zhat:})/(1+exp({zhat:})))))*`yvar')-({y0:   _cons})) 

local eqy0_inst instruments(eqy0:  )

local eqy1 (eqy1: ((`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*`yvar')-({y1:   _cons})) 

local eqy1_inst instruments(eqy1:  )

local num (num: {num}-(({y1:})-({y0:})))



local late (late: ({late} - {num}/{denom}))

} 
	 
if "`stat'" == "latt" {
	

	 sum `zvar' if `touse'==1
	matrix `w1s' = r(mean)
	 sum `zvar' if `touse'==1
	scalar `w1' = r(mean)
	
	 ge double `omega1'= (`zvar')/`w1'
	 ge double `omega0'= (((1-`zvar')*`ips')/(1-`ips'))/`w1'
	

			
	

	
	if `dmeanz1'==1 {
		 ge double `d1hat'= 1
	}
	else {
		 gen double `d1hat' = `omega1'*`tvar' if `touse'==1

	}
	
	if `dmeanz0'==0 {
		 ge double `d0hat'= 0
	}
	else {
		
		 gen double `d0hat' = `omega0'*`tvar' if `touse'==1
		
	}
	
	
	 gen double `y1hat' = `omega1'*`yvar' if `touse'==1
	
	 gen double `y0hat' = `omega0'*`yvar' if `touse'==1	

	
	 sum `d1hat'  if `touse'==1 [iw = `samplew']
	matrix `denom1s' = r(mean)
	 sum `d0hat'  if `touse'==1 [iw = `samplew']
	matrix `denom0s' = r(mean)
	 sum `y1hat'  if `touse'==1 [iw = `samplew']
	matrix `num1s' = r(mean)
	 sum `y0hat'  if `touse'==1 [iw = `samplew']
	matrix `num0s' = r(mean)
	matrix `nums' = `num1s'-`num0s'
	matrix `denoms' = `denom1s'-`denom0s'
	
	
	matrix `lates' =`nums'*invsym(`denoms')

	
	//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - ((exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*(1-`zvar')))

local eqips_inst instruments(eqips: `xvarsips')
}

	
	local eqw1 (eqw1: (`zvar' - {w1: _cons}))
    local eqw1_inst instruments(eqw1: )

	


	local eqd1 (eqd1: ((`zvar'*`tvar')/{w1:} - ({d1: _cons})))
	
	local eqd1_inst instruments(eqd1: )
	
		local eqd0 (eqd0: (((1-`zvar')*(exp({zhat:})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*`tvar')/{w1:}-({d0:  _cons}))		
		local eqd0_inst instruments(eqd0:  )
	

		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: (({denom}-(({d1: })-(0)))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (({denom}-(1-({d0: })))))
		}
		else {
			local denom (denom: (({denom}-(({d1: })-({d0: })))))
		}
	
	

	
		local eqy1 (eqy1: ((`zvar'*`yvar')/{w1:} - ({y1: _cons})))
	
	    local eqy1_inst instruments(eqy1: )
	
		local eqy0 (eqy0: (((1-`zvar')*(exp({zhat:})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*`yvar')/{w1:}-({y0:  _cons}))		
		local eqy0_inst instruments(eqy0:  )
	

		
		local num (num: (({num}-(({y1: })-({y0: })))))
		
	
	
	
	local late (late: ({late} - {num}/{denom}))
}	 
	 
	 
	 if "`stat'" == "late" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `denom0s', `denom1s', `num0s', `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `denom1s', `num0s', `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `denom0s', `num0s', `num1s', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom'  if `touse' [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	 
if "`stat'" == "latt" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial =  (`bips', `w1s', `denom0s', `denom1s', `num0s', `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqw1' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `w1s', `denom1s', `num0s', `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm    `eqips' `eqw1' `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `w1s', `denom0s', `num0s', `num1s', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm  `eqips' `eqw1' `eqd0' `eqy0' `eqy1'  `num' `denom'   if `touse' [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	
	}
	else if "`method'" == "nipw" {
			//IPS
set maxiter 300

if "`zmodel'" == "logit" {
  `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix `bips' = e(b)
 qui predict double `ips'
 }
else if "`zmodel'" == "cbps" {
	if "`stat'" == "late" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew')
  }
  	else if "`stat'" == "latt" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew') stat(att)
  }
  qui matrix `bips' = e(b)
   qui predict double `ipsxb'
   qui ge `ips' = logistic(`ipsxb')
  
}	

  

		
//Conditional mean of y1 & y0#
if "`stat'" == "late" {
qui regress `yvar'  if `zvar'==1 & `touse'==1 [pw = `samplew'/`ips']
qui matrix `by1' = e(b)
qui predict double `y1hat'

 
 
 qui regress `yvar'  if `zvar'==0 & `touse'==1 [pw = `samplew'/(1-`ips')]
qui matrix `by0' = e(b)
qui predict double `y0hat'

//Conditional mean of d1 & d0#


	    if `dmeanz1'==1 {
        ge double `d1hat'= 1
	    }
		else {

qui regress `tvar'  if `zvar'==1  & `touse'==1  [pw =  `samplew'/`ips']
matrix `bd1' = e(b)
qui predict double `d1hat'
}
	    
	    if `dmeanz0'==0 {
        ge double `d0hat' = 0
	    }
		else {
qui regress `tvar'  if `zvar'==0  & `touse'==1 [pw =  `samplew'/(1-`ips')]
matrix `bd0' = e(b)
qui predict double `d0hat'
		}
		
		
qui sum `d1hat' if `touse'==1 [iw = `samplew']
matrix `denom1s' = r(mean)
qui sum `d0hat'  if `touse'==1 [iw = `samplew']
matrix `denom0s' = r(mean)
qui sum `y1hat' if `touse'==1  [iw = `samplew']
matrix `num1s' = r(mean)
qui sum `y0hat' if `touse'==1 [iw = `samplew']
matrix `num0s' = r(mean)
matrix `nums' = `num1s'-`num0s'
matrix `denoms' = `denom1s'-`denom0s'
matrix `lates' =`nums'*invsym(`denoms')

//drop y1hat y0hat d1hat d0hat ips
	
//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((exp({zhat:})/(1+exp({zhat:})))*(1/(1+exp({zhat:})))))

local eqips_inst instruments(eqips: `xvarsips')
}



local eqd1 (eqd1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`tvar' - ({d1:  _cons}))))
           
local eqd1_inst instruments(eqd1: )

local eqd0 (eqd0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`tvar'-({d0:  _cons})))) 

local eqd0_inst instruments(eqd0:  )



		if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}- ({d1:}))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-({d0:})))
		}
		else {
local denom (denom: {denom}- (({d1:})-({d0:})))
		}




local eqy0 (eqy0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`yvar'-({y0:   _cons})))) 

local eqy0_inst instruments(eqy0:  )

local eqy1 (eqy1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`yvar'-({y1:   _cons})))) 

local eqy1_inst instruments(eqy1:  )

local num (num: {num}-(({y1:})-({y0:})))


local late (late: ({late} - {num}/{denom}))

} 
	 
if "`stat'" == "latt" {
	
	

	
	if `dmeanz1'==1 {
		ge double `d1hat'= 1
	}
	else {
		
		qui regress  `tvar'  if `zvar' == 1 & `touse'==1 [pw = `samplew']
		matrix `bd1' = e(b)
		qui predict double `d1hat'
	}
	
	if `dmeanz0'==0 {
		ge double `d0hat'= 0
	}
	else {
		qui regress `tvar'   if `zvar' == 0 & `touse'==1 [pw = `samplew'*(`ips'/(1-`ips'))]
		matrix `bd0' = e(b)
		qui predict double `d0hat'
	}
	
	qui regress `yvar' if `zvar' == 1 & `touse'==1 [pw = `samplew']
	matrix `by1' = e(b)
	qui predict double `y1hat'
	
	qui regress `yvar'  if `zvar' == 0 & `touse'==1 [pw = `samplew'*(`ips'/(1-`ips'))]
	matrix `by0' = e(b)
	qui predict double `y0hat'
	
	
	qui sum `d1hat'  if `zvar' == 1 [iw = `samplew']
	matrix `denom1s' = r(mean)
	qui sum `d0hat'  if `zvar' == 1 [iw = `samplew']
	matrix `denom0s' = r(mean)
	qui sum `y1hat'  if `zvar' == 1 [iw = `samplew']
	matrix `num1s' = r(mean)
	qui sum `y0hat'  if `zvar' == 1 [iw = `samplew']
	matrix `num0s' = r(mean)
	matrix `nums' = `num1s'-`num0s'
	matrix `denoms' = `denom1s'-`denom0s'
	
	
	matrix `lates' =`nums'*invsym(`denoms')
	
	//drop y1hat y0hat d1hat d0hat ips
	
	
	
	//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - ((exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*(1-`zvar')))

local eqips_inst instruments(eqips: `xvarsips')
}

	


	local eqd1 (eqd1: ((1.`zvar'*(`tvar' - ({d1: _cons})))))
	
	local eqd1_inst instruments(eqd1: )
	

		
		local eqd0 (eqd0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`tvar'-({d0:  _cons})))) 
		
		local eqd0_inst instruments(eqd0:  )

	

	

		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-({d0: })))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}-(({d1: })-({d0: })))))
		}
	
	
	local eqy1 (eqy1: (1.`zvar'*(`yvar'-({y1:  _cons})))) 
	
	local eqy1_inst instruments(eqy1:  )
	
	
	


		local eqy0 (eqy0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`yvar'-({y0:   _cons})))) 
		
		local eqy0_inst instruments(eqy0:  )
		
		
		local num (num: (1.`zvar'*({num}-(({y1: })-({y0: })))))
		
	

	
	local late (late: ({late} - {num}/{denom}))
}	 
	 
	 
	 if "`stat'" == "late" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	 
if "`stat'" == "latt" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm   `eqips'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	
	}
	else if "`method'" == "ra" {
	//IPS
set maxiter 300


//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar' `xvarsy' if `zvar'==1 & `touse'==1 [pw = `samplew']
qui matrix `by1' = e(b)
qui predict double `y1hat'

 
 
 qui `omodel' `yvar' `xvarsy' if `zvar'==0 & `touse'==1 [pw = `samplew']
qui matrix `by0' = e(b)
qui predict double `y0hat'

//Conditional mean of d1 & d0#


	    if `dmeanz1'==1 {
        ge double `d1hat'= 1
	    }
		else {

qui `tmodel' `tvar' `xvarst' if `zvar'==1  & `touse'==1  [pw =  `samplew']
matrix `bd1' = e(b)
qui predict double `d1hat'
}
	    
	    if `dmeanz0'==0 {
        ge double `d0hat' = 0
	    }
		else {
qui `tmodel' `tvar' `xvarst' if `zvar'==0  & `touse'==1 [pw =  `samplew']
matrix `bd0' = e(b)
qui predict double `d0hat'
		}
		
		
qui sum `d1hat' if `touse'==1 [iw = `samplew']
matrix `denom1s' = r(mean)
qui sum `d0hat'  if `touse'==1 [iw = `samplew']
matrix `denom0s' = r(mean)
qui sum `y1hat' if `touse'==1  [iw = `samplew']
matrix `num1s' = r(mean)
qui sum `y0hat' if `touse'==1 [iw = `samplew']
matrix `num0s' = r(mean)
matrix `nums' = `num1s'-`num0s'
matrix `denoms' = `denom1s'-`denom0s'
matrix `lates' =`nums'*invsym(`denoms')

//drop y1hat y0hat d1hat d0hat 
	
//Moment Functions

if "`tmodel'" == "logit" {
local eqd1 (eqd1: ((1.`zvar')*(`tvar' - exp({d1: `xvarst' _cons})/(1+exp({d1:})))))
           
local eqd1_inst instruments(eqd1: `xvarst')

local eqd0 (eqd0: ((0.`zvar')*(`tvar'-exp({d0: `xvarst' _cons})/(1+exp({d0:}))))) 

local eqd0_inst instruments(eqd0: `xvarst' )
}


if "`tmodel'" =="regress" {
local eqd1 (eqd1: ((1.`zvar')*(`tvar' - ({d1: `xvarst' _cons}))))
           
local eqd1_inst instruments(eqd1: `xvarst')

local eqd0 (eqd0: ((0.`zvar')*(`tvar'-({d0: `xvarst' _cons})))) 

local eqd0_inst instruments(eqd0: `xvarst' )
}

if "`tmodel'" == "logit" {
	if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}- ((exp({d1:})/(1+exp({d1:})))))
	}
	else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-((exp({d0:})/(1+exp({d0:}))))))
	}
	else {
local denom (denom: {denom}- (((exp({d1:})/(1+exp({d1:}))))-((exp({d0:})/(1+exp({d0:}))))))
	}
}
else if "`tmodel'" =="regress"  {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}- ({d1:}))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-({d0:})))
		}
		else {
local denom (denom: {denom}- (({d1:})-({d0:})))
		}
}



if "`omodel'" == "regress" {
local eqy0 (eqy0: ((0.`zvar')*(`yvar'-({y0: `xvarsy'  _cons})))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar')*(`yvar'-({y1: `xvarsy'  _cons})))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num (num: {num}-(({y1:})-({y0:})))
}
else if "`omodel'"  == "logit" {
local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})/(1+exp({y0:}))))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar')*(`yvar'-exp({y1: `xvarsy' _cons})/(1+exp({y1:}))))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num (num: {num}-(((exp({y1: })/(1+exp({y1:}))))-((exp({y0: })/(1+exp({y0:}))))))

}
else if "`omodel'"  == "poisson"  {
local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy'  _cons})))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar')*(`yvar'-exp({y1: `xvarsy'  _cons})))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num (num: {num}-((exp({y1: }))-(exp({y0: }))))

}

local late (late: ({late} - {num}/{denom}))
} 

	

	/*Estimation of LATT*/
if "`stat'" == "latt" {
	
	if `dmeanz1'==1 {
		ge double `d1hat'= 1
	}
	else {
		
		qui regress  `tvar'  if `zvar' == 1 & `touse'==1 [pw = `samplew']
		matrix `bd1' = e(b)
		qui predict double `d1hat'
	}
	
	if `dmeanz0'==0 {
		ge double `d0hat' = 0
	}
	else {
		qui `tmodel' `tvar' `xvarst'  if `zvar' == 0 & `touse'==1 [pw = `samplew']
		matrix `bd0' = e(b)
		qui predict double `d0hat'
	}
	
	qui regress `yvar' if `zvar' == 1 & `touse'==1 [pw = `samplew']
	matrix `by1' = e(b)
	qui predict double `y1hat'
	
	qui `omodel' `yvar' `xvarsy'  if `zvar' == 0 & `touse'==1 [pw = `samplew']
	matrix `by0' = e(b)
	qui predict double `y0hat'
	
	
	qui sum `d1hat'  if `zvar' == 1 [iw = `samplew']
	matrix `denom1s' = r(mean)
	qui sum `d0hat'  if `zvar' == 1 [iw = `samplew']
	matrix `denom0s' = r(mean)
	qui sum `y1hat'  if `zvar' == 1 [iw = `samplew']
	matrix `num1s' = r(mean)
	qui sum `y0hat'  if `zvar' == 1 [iw = `samplew']
	matrix `num0s' = r(mean)
	matrix `nums' = `num1s'-`num0s'
	matrix `denoms' = `denom1s'-`denom0s'
	
	
	matrix `lates' =`nums'*invsym(`denoms')
	
	//drop y1hat y0hat d1hat d0hat 
	
	
	
	//Moment Functions
	local eqd1 (eqd1: ((1.`zvar'*(`tvar' - ({d1: _cons})))))
	
	local eqd1_inst instruments(eqd1: )
	
	
	if "`tmodel'" == "logit" {
		
		local eqd0 (eqd0: ((0.`zvar')*(`tvar'-exp({d0: `xvarst' _cons})/(1+exp({d0:}))))) 
		
		local eqd0_inst instruments(eqd0: `xvarst' )
	}
	
	
	if "`tmodel'" =="regress" {
		
		local eqd0 (eqd0: ((0.`zvar')*(`tvar'-({d0: `xvarst' _cons})))) 
		
		local eqd0_inst instruments(eqd0: `xvarst' )
	}
	

	
	
	if "`tmodel'" == "logit" {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-(exp({d0: })/(1+exp({d0:})))))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}-(({d1: })-(exp({d0: })/(1+exp({d0:})))))))
		}
	}
	else if "`tmodel'" =="regress"  {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-({d0: })))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}-(({d1: })-({d0: })))))
		}
	}
	

	
	local eqy1 (eqy1: (1.`zvar'*(`yvar'-({y1:  _cons})))) 
	
	local eqy1_inst instruments(eqy1:  )
	
	
	
	if "`omodel'" == "regress" {

		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-({y0: `xvarsy'  _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )
		
		
		local num (num: (1.`zvar'*({num}-(({y1: })-({y0: })))))
		
	}
	else if "`omodel'"  == "logit" {
		
		
		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})/(1+exp({y0:}))))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )
		
		local num (num: (1.`zvar'*({num}-(({y1: })-(exp({y0: })/(1+exp({y0:})))))))
		
		
	}
	else if "`omodel'"  == "poisson"  {
		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )
		
		local num (num: (1.`zvar'*({num}-(({y1: })-(exp({y0: }))))))
		
		
	}
	
	local late (late: ({late} - {num}/{denom}))
}

if "`stat'" == "late" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bd0', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm  `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  if `touse' [pw = `samplew'],   `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],    `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqd0' `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],   `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	 
if "`stat'" == "latt" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bd0', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],   `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = ( `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm    `eqy0' `eqy1'  `num' `denom' if `touse' [pw = `samplew'],   `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bd0', `by0', `by1', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqd0' `eqy0' `eqy1'  `num' `denom'  if `touse' [pw = `samplew'],   `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}		
		
		
	}	
	else if "`method'" == "aipw" {
		
			//IPS
set maxiter 300
if "`zmodel'" == "logit" {
  `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix `bips' = e(b)
 qui predict double `ips'
 }
else if "`zmodel'" == "cbps" {
	if "`stat'" == "late" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew')
  }
  	else if "`stat'" == "latt" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew') stat(att)
  }
   qui matrix `bips' = e(b)
   qui predict double `ipsxb'
   qui ge `ips' = logistic(`ipsxb')
  
}	

  

		
//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar' `xvarsy' if `zvar'==1 & `touse'==1 [pw = `samplew']
qui matrix `by1' = e(b)
qui predict double `y1hat'

 
 
 qui `omodel' `yvar' `xvarsy' if `zvar'==0 & `touse'==1 [pw = `samplew']
qui matrix `by0' = e(b)
qui predict `y0hat'

//Conditional mean of d1 & d0#


	    if `dmeanz1'==1 {
        ge double `d1hat'= 1
	    }
		else {

qui `tmodel' `tvar' `xvarst' if `zvar'==1  & `touse'==1  [pw =  `samplew']
matrix `bd1' = e(b)
qui predict double `d1hat'
}
	    
	    if `dmeanz0'==0 {
        ge double `d0hat'= 0
	    }
		else {
qui `tmodel' `tvar' `xvarst' if `zvar'==0  & `touse'==1 [pw =  `samplew']
matrix `bd0' = e(b)
qui predict double `d0hat'
		}
		
	qui gen double `term1' = ((`zvar'*`yvar')-((`zvar'-`ips')*`y1hat'))/`ips'
	qui gen double `term2'  = (((1-`zvar')*`yvar')+((`zvar'-`ips')*`y0hat'))/(1-`ips')
	qui gen double `term3'  = ((`zvar'*`tvar')-((`zvar'-`ips')*`d1hat'))/`ips'
	qui gen double `term4'  = (((1-`zvar')*`tvar')+((`zvar'-`ips')*`d0hat'))/(1-`ips')
	
	

qui sum `term1'  if  `touse'==1  [iw = `samplew']
matrix `num1s' = r(mean)
qui sum `term2'  if  `touse'==1 [iw = `samplew']
matrix `num0s' = r(mean)
matrix `nums' = `num1s'-`num0s'
qui sum `term3'  if  `touse'==1 [iw = `samplew']
matrix `denom1s' = r(mean)
qui sum `term4'  if  `touse'==1 [iw = `samplew']
matrix `denom0s' = r(mean)
matrix `denoms' = `denom1s'-`denom0s'
matrix `lates' =`nums'*invsym(`denoms')

//drop y1hat y0hat d1hat d0hat term* ips

//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((exp({zhat:})/(1+exp({zhat:})))*(1/(1+exp({zhat:})))))

local eqips_inst instruments(eqips: `xvarsips')
}



if "`tmodel'" == "logit" {
local eqd1 (eqd1: ((1.`zvar')*(`tvar' - exp({d1: `xvarst' _cons})/(1+exp({d1:})))))
           
local eqd1_inst instruments(eqd1: `xvarst')

local eqd0 (eqd0: ((0.`zvar')*(`tvar'-exp({d0: `xvarst' _cons})/(1+exp({d0:}))))) 

local eqd0_inst instruments(eqd0: `xvarst' )
}


if "`tmodel'" =="regress" {
local eqd1 (eqd1: ((1.`zvar')*(`tvar' - ({d1: `xvarst' _cons}))))
           
local eqd1_inst instruments(eqd1: `xvarst')

local eqd0 (eqd0: ((0.`zvar')*(`tvar'-({d0: `xvarst' _cons})))) 

local eqd0_inst instruments(eqd0: `xvarst' )
}

if "`tmodel'" == "logit" {
	if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}- ((((`zvar'*`tvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({d1: })/(1+exp({d1:}))))))/(exp({zhat:})/(1+exp({zhat:})))) )
	}
	else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-((((1-`zvar')*`tvar')+((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({d0: })/(1+exp({d0:}))))) / (1-(exp({zhat:})/(1+exp({zhat:})))) ) ))
	}
	else {
local denom1 ( denom1: {denom1}- ((((`zvar'*`tvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({d1: })/(1+exp({d1:}))))))/(exp({zhat:})/(1+exp({zhat:})))) )
local denom0 ( denom0: {denom0}- ((((1-`zvar')*`tvar')+((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({d0: })/(1+exp({d0:}))))) / (1-(exp({zhat:})/(1+exp({zhat:})))) ) )
local denom (denom: {denom}-({denom1}-{denom0}))
	}
}
else if "`tmodel'" =="regress"  {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}- (((`zvar'*`tvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*{d1:}))/(exp({zhat:})/(1+exp({zhat:})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-((((1-`zvar')*`tvar')+((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*{d0:}))/(1-(exp({zhat:})/(1+exp({zhat:}))))) ))
		}
		else {
local denom1 (denom1: {denom1}- (((`zvar'*`tvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*{d1:}))/(exp({zhat:})/(1+exp({zhat:})))))
local denom0 (denom0: {denom0}- ((((1-`zvar')*`tvar')+((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*{d0:}))/(1-(exp({zhat:})/(1+exp({zhat:}))))) )
local denom (denom: {denom}-({denom1}-{denom0}))
		}
}



if "`omodel'" == "regress" {
local eqy0 (eqy0: ((0.`zvar')*(`yvar'-({y0: `xvarsy'  _cons})))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar')*(`yvar'-({y1: `xvarsy'  _cons})))) 

local eqy1_inst instruments(eqy1: `xvarsy' )


local num1 (num1: {num1}- (((`zvar'*`yvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*{y1:}))/(exp({zhat:})/(1+exp({zhat:})))))
local num0 (num0: {num0}-  ((((1-`zvar')*`yvar')+((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*{y0:}))/(1-(exp({zhat:})/(1+exp({zhat:}))))) )

}
else if "`omodel'"  == "logit" {
local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})/(1+exp({y0:}))))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar')*(`yvar'-exp({y1: `xvarsy' _cons})/(1+exp({y1:}))))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num1 ( num1: {num1}- ((((`zvar'*`yvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({y1: })/(1+exp({y1:}))))))/(exp({zhat:})/(1+exp({zhat:})))) )
local num0 (num0: {num0}-  ((((1-`zvar')*`yvar')+((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({y0: })/(1+exp({y0:}))))) / (1-(exp({zhat:})/(1+exp({zhat:})))) ) )

}
else if "`omodel'"  == "poisson"  {
local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy'  _cons})))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar')*(`yvar'-exp({y1: `xvarsy'  _cons})))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num1 ( num1: {num1}- (((`zvar'*`yvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({y1: }))))/(exp({zhat:})/(1+exp({zhat:})))) )
local num0 (num0: {num0}-  ((((1-`zvar')*`yvar')+((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({y0: }))))/(1-(exp({zhat:})/(1+exp({zhat:}))))) )
}
local num (num: {num}-({num1}-{num0}) )


local late (late: ({late} - {num}/{denom}))
}
else if "`stat'" == "latt" {
	
	
	 sum `zvar' if `touse'==1
	matrix `w1s' = r(mean)
	 sum `zvar' if `touse'==1
	scalar `w1' = r(mean)
	
	 ge double `omega1'= (`zvar')/`w1'
	 ge double `omega0'= (((1-`zvar')*`ips')/(1-`ips'))/`w1'
	
	
	if `dmeanz1'==1 {
		 ge double `d1hat'= 1
	}
	else {
		 gen double `d1hat' = `omega1'*`tvar' if `touse'==1

	}
	  	if `dmeanz0'==0 {
		ge double `d0hat' = 0
	}
	else {
		 `tmodel' `tvar' `xvarst'  if `zvar' == 0 & `touse'==1 [pw = `samplew']
		matrix `bd0' = e(b)
		qui predict double `d0hat'
	}
	
	gen double `y1hat' = `omega1'*`yvar' if `touse'==1
	

	 `omodel' `yvar' `xvarsy'  if `zvar' == 0 & `touse'==1 [pw = `samplew']
	matrix `by0' = e(b)
	qui predict double `y0hat'
	
	
	


	

	
	//qui gen double `term1' = ((`zvar'*`yvar')-((`zvar'-`ips')*`y1hat'))/`ips'
	 gen double `term2'  = (`omega0'*(`yvar'-`y0hat'))+`zvar'*`y0hat'/`w1'
	//qui gen double `term3'  = ((`zvar'*`tvar')-((`zvar'-`ips')*`d1hat'))/`ips'
	 gen double `term4'  = (`omega0'*(`tvar'-`d0hat'))+`zvar'*`d0hat'/`w1'

	 sum `d1hat'  if `touse'==1 [iw = `samplew']
	matrix `denom1s' = r(mean)
	 sum `term4'  if `touse'==1 [iw = `samplew']
	matrix `denom0s' = r(mean)
	 sum `y1hat'  if `touse'==1 [iw = `samplew']
	matrix `num1s' = r(mean)
	 sum `term2'  if `touse'==1 [iw = `samplew']
	matrix `num0s' = r(mean)
	matrix `nums' = `num1s'-`num0s'
	matrix `denoms' = `denom1s'-`denom0s'
	
	
	matrix `lates' =`nums'*invsym(`denoms')
	


//Moment Functions

	//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - ((exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*(1-`zvar')))

local eqips_inst instruments(eqips: `xvarsips')
}

	
	local eqw1 (eqw1: (`zvar' - {w1: _cons}))
    local eqw1_inst instruments(eqw1: )

	local eqd1 (eqd1: ((`zvar'*`tvar')/{w1:} - ({d1: _cons})))
	
	local eqd1_inst instruments(eqd1: )
	
	local eqy1 (eqy1: ((`zvar'*`yvar')/{w1:} - ({y1: _cons})))
	
	local eqy1_inst instruments(eqy1: )
	
	
	if "`tmodel'" == "logit" {
		
		local eqd0 (eqd0: ((0.`zvar')*(`tvar'-exp({d0: `xvarst' _cons})/(1+exp({d0:}))))) 
		
		local eqd0_inst instruments(eqd0: `xvarst' )
	}
	
	
	if "`tmodel'" =="regress" {
		
		local eqd0 (eqd0: ((0.`zvar')*(`tvar'-({d0: `xvarst' _cons})))) 
		
		local eqd0_inst instruments(eqd0: `xvarst' )
	}
	

	
	
	if "`tmodel'" == "logit" {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((({d1:} - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (({denom}-(1-(-`zvar'*((exp({d0: })/(1+exp({d0:}))))+ ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-(exp({d0: })/(1+exp({d0:})))))/(1/(1+exp({zhat:}))))/{w1:}))))
		}
		else {
			local denom (denom: (({denom}  -(`zvar'*(`tvar'-(exp({d0: })/(1+exp({d0:}))))- ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-(exp({d0: })/(1+exp({d0:})))))/(1/(1+exp({zhat:}))))/{w1:})))
		}
	}
	else if "`tmodel'" =="regress"  {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((({d1:} - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (({denom}-(1-(-`zvar'*({d0:})+ ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-{d0:}))/(1/(1+exp({zhat:}))))/{w1:}))))
		}
		else {
			local denom (denom: (({denom}  -(`zvar'*(`tvar'-{d0:})- ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-{d0:}))/(1/(1+exp({zhat:}))))/{w1:})))
	
	}
	}
	

	
	
	if "`omodel'" == "regress" {

		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-({y0: `xvarsy'  _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )
		
		local num (num: (({num}  -(`zvar'*(`yvar'-{y0:})- ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`yvar'-{y0:}))/(1/(1+exp({zhat:}))))/{w1:})))		
	}
	else if "`omodel'"  == "logit" {
		
		
		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})/(1+exp({y0:}))))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )

		local num (num: (({num}  -(`zvar'*(`yvar'-(exp({y0: })/(1+exp({y0:}))))- ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`yvar'-(exp({y0: })/(1+exp({y0:})))))/(1/(1+exp({zhat:}))))/{w1:})))
		
	}
	else if "`omodel'"  == "poisson"  {
		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )

		local num (num: (({num}  -(`zvar'*(`yvar'-(exp({y0: })))- ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`yvar'-(exp({y0: }))))/(1/(1+exp({zhat:}))))/{w1:})))		
		
	}
	
	local late (late: ({late} - {num}/{denom}))
	
} 
	if "`stat'" == "late" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', `num1s', `num0s', `denom1s', `denom0s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm  `eqips' `eqd0' `eqd1'  `eqy0' `eqy1' `num1'  `num0' `denom1' `denom0'  `num' `denom' if `touse'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial =  (`bips', `bd1', `by0', `by1',  `num1s', `num0s', `nums', `denoms')
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1' `num1'  `num0'   `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1', `num1s', `num0s', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqips' `eqd0'  `eqy0' `eqy1' `num1'  `num0'  `num' `denom'  if `touse' [pw = `samplew'],   `eqips_inst' `eqd0_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}
else if "`stat'" == "latt" {  

	if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `w1s', `bd0', `denom1s', `by0', `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqw1' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  if `touse' [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `w1s', `denom1s', `by0', `num1s', `nums', `denoms')
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm    `eqips' `eqw1' `eqd1'  `eqy0' `eqy1'  `num' `denom'  if `touse' [pw = `samplew'],  `eqips_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `w1s', `bd0', `by0', `num1s', `nums', `denoms')
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqips' `eqw1' `eqd0' `eqy0' `eqy1'  `num' `denom'  if `touse' [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	
	}
	
		else if "`method'" == "naipw" {
			
			//IPS
set maxiter 300

if "`zmodel'" == "logit" {
  `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix `bips' = e(b)
 qui predict double `ips'
 }
else if "`zmodel'" == "cbps" {
	if "`stat'" == "late" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew')
  }
  	else if "`stat'" == "latt" {
  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew') stat(att)
  }
  
   qui matrix `bips' = e(b)
   qui predict double `ipsxb'
   qui ge `ips' = logistic(`ipsxb')
  
}	

  

		
//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar' `xvarsy' if `zvar'==1 & `touse'==1 [pw = `samplew']
qui matrix `by1' = e(b)
qui predict double `y1hat'

 
 
 qui `omodel' `yvar' `xvarsy' if `zvar'==0 & `touse'==1 [pw = `samplew']
qui matrix `by0' = e(b)
qui predict double `y0hat'

//Conditional mean of d1 & d0#


	    if `dmeanz1'==1 {
        ge double `d1hat'= 1
	    }
		else {

qui `tmodel' `tvar' `xvarst' if `zvar'==1  & `touse'==1  [pw =  `samplew']
matrix `bd1' = e(b)
qui predict double `d1hat'
}
	    
	    if `dmeanz0'==0 {
        ge double `d0hat'= 0
	    }
		else {
qui `tmodel' `tvar' `xvarst' if `zvar'==0  & `touse'==1 [pw =  `samplew']
matrix `bd0' = e(b)
qui predict double `d0hat'
		}
		
	qui gen double `invw1' = `zvar'/`ips'
	qui gen double `invw0' = (1-`zvar')/(1-`ips')
	
	qui sum `invw1'  if  `touse'==1  
	scalar `wr1' = r(mean)
	qui sum `invw0'  if  `touse'==1 
	scalar `wr0' = r(mean)
	qui sum `invw1'  if  `touse'==1 
	matrix `w1s' = r(mean)
	qui sum `invw0'  if  `touse'==1   
	matrix `w0s' = r(mean)
	
	
	
	qui ge double `omega1'= (`zvar'/`ips')/`wr1'
	qui ge double `omega0'= ((1-`zvar')/(1-`ips'))/`wr0'
	
	qui ge double `term1'= `omega1'*( `yvar' -`y1hat')+`y1hat'
	qui sum `term1'  if  `touse'==1  [iw = `samplew']
	matrix `num1s' =  r(mean)
	qui ge double `term2'= `omega0'*( `yvar' -`y0hat')+`y0hat'
	qui sum `term2'  if  `touse'==1  [iw = `samplew']
	matrix `num0s' =  r(mean)
	qui ge double `term3' = `omega1'*( `tvar' -`d1hat')+`d1hat'
	qui sum `term3'  if  `touse'==1  [iw = `samplew']
	matrix `denom1s' =  r(mean)	
	qui ge double `term4' = `omega0'*( `tvar' -`d0hat')+`d0hat'
	qui sum `term4'  if  `touse'==1  [iw = `samplew']
	matrix `denom0s' =  r(mean)
	
	matrix `nums' = `num1s'-`num0s'
    matrix `denoms' = `denom1s'-`denom0s'
    matrix `lates' =`nums'*invsym(`denoms')
		

		
//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((exp({zhat:})/(1+exp({zhat:})))*(1/(1+exp({zhat:})))))

local eqips_inst instruments(eqips: `xvarsips')
}


local eqw1 (eqw1: (((`zvar'/(exp({zhat:})/(1+exp({zhat:})))) - {w1: _cons})))
local eqw1_inst instruments(eqw1: )
local eqw0 (eqw0: ((((1-`zvar')/(1-(exp({zhat:})/(1+exp({zhat:})))) - {w0: _cons}))))
local eqw0_inst instruments(eqw0: )

if "`tmodel'" == "logit" {
local eqd1 (eqd1: ((1.`zvar')*(`tvar' - exp({d1: `xvarst' _cons})/(1+exp({d1:})))))
           
local eqd1_inst instruments(eqd1: `xvarst')

local eqd0 (eqd0: ((0.`zvar')*(`tvar'-exp({d0: `xvarst' _cons})/(1+exp({d0:}))))) 

local eqd0_inst instruments(eqd0: `xvarst' )
}


if "`tmodel'" =="regress" {
local eqd1 (eqd1: ((1.`zvar')*(`tvar' - ({d1: `xvarst' _cons}))))
           
local eqd1_inst instruments(eqd1: `xvarst')

local eqd0 (eqd0: ((0.`zvar')*(`tvar'-({d0: `xvarst' _cons})))) 

local eqd0_inst instruments(eqd0: `xvarst' )
}

if "`tmodel'" == "logit" {
	if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}- ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`tvar'-((exp({d1: })/(1+exp({d1:})))))-((exp({d1: })/(1+exp({d1:})))))
	}
	else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-((0.`zvar'/(1-exp({zhat:})/(1+exp({zhat:}))))/{w0:})*(`tvar'-((exp({d0: })/(1+exp({d0:})))))-((exp({d0: })/(1+exp({d0:}))))))
	}
	else {
local denom1 (denom1: {denom1}- ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`tvar'-((exp({d1: })/(1+exp({d1:})))))-((exp({d1: })/(1+exp({d1:})))))
local denom0 (denom0: {denom0}- ((0.`zvar'/(1-exp({zhat:})/(1+exp({zhat:}))))/{w0:})*(`tvar'-((exp({d0: })/(1+exp({d0:})))))-((exp({d0: })/(1+exp({d0:})))))
local denom (denom: {denom}-({denom1}-{denom0}))
	}
}
else if "`tmodel'" =="regress"  {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
local denom (denom: {denom}-  ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`tvar'-({d1:}))-({d1:}))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
local denom (denom: {denom}- (1-((0.`zvar'/(1-exp({zhat:})/(1+exp({zhat:}))))/{w0:})*(`tvar'-({d0:}))-({d0:})))
		}
		else {
local denom1 (denom1: {denom1}- ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`tvar'-({d1:}))-({d1:}))

local denom0 (denom0: {denom0}- ((0.`zvar'/(1-exp({zhat:})/(1+exp({zhat:}))))/{w0:})*(`tvar'-({d0:}))-({d0:}))

local denom (denom: {denom}-({denom1}-{denom0}))
		}
}



if "`omodel'" == "regress" {
local eqy0 (eqy0: ((0.`zvar')*(`yvar'-({y0: `xvarsy'  _cons})))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar')*(`yvar'-({y1: `xvarsy'  _cons})))) 

local eqy1_inst instruments(eqy1: `xvarsy' )


local num1 (num1: {num1}- ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`yvar'-({y1:}))-({y1:}))

local num0 (num0: {num0}- ((0.`zvar'/(1-exp({zhat:})/(1+exp({zhat:}))))/{w0:})*(`yvar'-({y0:}))-({y0:}))

}
else if "`omodel'"  == "logit" {
local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})/(1+exp({y0:}))))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar')*(`yvar'-exp({y1: `xvarsy' _cons})/(1+exp({y1:}))))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num1 (num1: {num1}- ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`yvar'-((exp({y1: })/(1+exp({y1:})))))-((exp({y1: })/(1+exp({y1:})))))

local num0 (num0: {num0}- ((0.`zvar'/(1-exp({zhat:})/(1+exp({zhat:}))))/{w0:})*(`yvar'-((exp({y0: })/(1+exp({y0:})))))-((exp({y0: })/(1+exp({y0:})))))

}
else if "`omodel'"  == "poisson"  {
local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy'  _cons})))) 

local eqy0_inst instruments(eqy0: `xvarsy' )

local eqy1 (eqy1: ((1.`zvar')*(`yvar'-exp({y1: `xvarsy'  _cons})))) 

local eqy1_inst instruments(eqy1: `xvarsy' )

local num1 (num1: {num1}- ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`yvar'-(exp({y1: })))-(exp({y1: })))

local num0 (num0: {num0}- ((0.`zvar'/(1-exp({zhat:})/(1+exp({zhat:}))))/{w0:})*(`yvar'-(exp({y0: })))-(exp({y0: })))

}
local num (num: {num}-({num1}-{num0}) )


local late (late: ({late} - {num}/{denom}))
} 
else if "`stat'" == "latt" {
		
	 sum `zvar' if `touse'==1
	matrix `w1s' = r(mean)
	 sum `zvar' if `touse'==1
	scalar `w1' = r(mean)
	
	 ge double `omega1'= (`zvar')/`w1'
	 ge double `omega0'= (((1-`zvar')*`ips')/(1-`ips'))
	sum `omega0' if `touse'==1
	scalar `w0'= r(mean)    
	matrix `w0s' = r(mean) 
	replace `omega0' =  `omega0'/`w0'


	
	if `dmeanz1'==1 {
		 ge double `d1hat'= 1
	}
	else {
		 gen double `d1hat' = `omega1'*`tvar' if `touse'==1

	}
	  	if `dmeanz0'==0 {
		ge double `d0hat' = 0
	}
	else {
		 `tmodel' `tvar' `xvarst'  if `zvar' == 0 & `touse'==1 [pw = `samplew']
		matrix `bd0' = e(b)
		qui predict double `d0hat'
	}
	
	gen double `y1hat' = `omega1'*`yvar' if `touse'==1
	

	 `omodel' `yvar' `xvarsy'  if `zvar' == 0 & `touse'==1 [pw = `samplew']
	matrix `by0' = e(b)
	qui predict double `y0hat'
	
	
	//qui gen double `term1' = ((`zvar'*`yvar')-((`zvar'-`ips')*`y1hat'))/`ips'
	 gen double `term2'  = (`omega0'*(`yvar'-`y0hat'))+`zvar'*`y0hat'/`w1'
	//qui gen double `term3'  = ((`zvar'*`tvar')-((`zvar'-`ips')*`d1hat'))/`ips'
	 gen double `term4'  = (`omega0'*(`tvar'-`d0hat'))+`zvar'*`d0hat'/`w1'

	 sum `d1hat'  if `touse'==1 [iw = `samplew']
	matrix `denom1s' = r(mean)
	 sum `term4'  if `touse'==1 [iw = `samplew']
	matrix `denom0s' = r(mean)
	 sum `y1hat'  if `touse'==1 [iw = `samplew']
	matrix `num1s' = r(mean)
	 sum `term2'  if `touse'==1 [iw = `samplew']
	matrix `num0s' = r(mean)
	matrix `nums' = `num1s'-`num0s'
	matrix `denoms' = `denom1s'-`denom0s'
	
	
	matrix `lates' =`nums'*invsym(`denoms')
	
	

//Moment Functions

	//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - ((exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*(1-`zvar')))

local eqips_inst instruments(eqips: `xvarsips')
}

	
	local eqw1 (eqw1: (`zvar' - {w1: _cons}))
    local eqw1_inst instruments(eqw1: )

	local eqd1 (eqd1: ((`zvar'*`tvar')/{w1:} - ({d1: _cons})))
	
	local eqd1_inst instruments(eqd1: )
	
	local eqy1 (eqy1: ((`zvar'*`yvar')/{w1:} - ({y1: _cons})))
	
	local eqy1_inst instruments(eqy1: )
	
	
	local eqw0 (eqw0: (((1-`zvar')/(1-(exp({zhat:})/(1+exp({zhat:})))))*(exp({zhat:})/(1+exp({zhat:}))) - {w0: _cons}))
    local eqw0_inst instruments(eqw0: )
	
	if "`tmodel'" == "logit" {
		
		local eqd0 (eqd0: ((0.`zvar')*(`tvar'-exp({d0: `xvarst' _cons})/(1+exp({d0:}))))) 
		
		local eqd0_inst instruments(eqd0: `xvarst' )
	}
	
	
	if "`tmodel'" =="regress" {
		
		local eqd0 (eqd0: ((0.`zvar')*(`tvar'-({d0: `xvarst' _cons})))) 
		
		local eqd0_inst instruments(eqd0: `xvarst' )
	}
	

	
	
	if "`tmodel'" == "logit" {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((({d1:} - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (({denom}-(1-(-(`zvar'*((exp({d0: })/(1+exp({d0:})))))/{w1:}+ (((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-(exp({d0: })/(1+exp({d0:})))))/(1/(1+exp({zhat:}))))/{w0:})))))
		}
		else {
		local denom (denom: (({denom}  -((`zvar'*(`tvar'-(exp({d0: })/(1+exp({d0:})))))/{w1:} - (((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-(exp({d0: })/(1+exp({d0:})))))/(1/(1+exp({zhat:}))))/{w0:}))))	


		}
	}
	else if "`tmodel'" =="regress"  {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((({d1:} - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (({denom}-(1-(-`zvar'*({d0:})/{w1:}+ (((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-{d0:}))/(1/(1+exp({zhat:}))))/{w1:})))))
		}
		else {
			local denom (denom: (({denom}  -((`zvar'*(`tvar'-{d0:}))/{w1:} - (((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-{d0:}))/(1/(1+exp({zhat:}))))/{w0:}))))	

	}
	}
	

	
	
	if "`omodel'" == "regress" {

		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-({y0: `xvarsy'  _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )
		
		local num (num: (({num}  -((`zvar'*(`yvar'-{y0:}))/{w1:} - (((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`yvar'-{y0:}))/(1/(1+exp({zhat:}))))/{w0:}))))	
	}
	else if "`omodel'"  == "logit" {
		
		
		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})/(1+exp({y0:}))))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )

		local num (num: (({num}  -((`zvar'*(`yvar'-(exp({y0: })/(1+exp({y0:})))))/{w1:} - (((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`yvar'-(exp({y0: })/(1+exp({y0:})))))/(1/(1+exp({zhat:}))))/{w0:}))))	

	}
	else if "`omodel'"  == "poisson"  {
		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )

		local num (num: (({num}  -((`zvar'*(`yvar'-(exp({y0: }))))/{w1:} - (((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`yvar'-(exp({y0: }))))/(1/(1+exp({zhat:}))))/{w0:}))))	
		
	}
	
	local late (late: ({late} - {num}/{denom}))


	}		
	if "`stat'" == "late" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', `w1s', `w0s', `num1s', `num0s', `denom1s', `denom0s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm  `eqips' `eqd0' `eqd1'  `eqy0' `eqy1' `eqw1' `eqw0' `num1'  `num0' `denom1' `denom0'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial =  (`bips', `bd1', `by0', `by1', `w1s', `w0s',  `num1s', `num0s', `nums', `denoms')
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1' `eqw1' `eqw0' `num1'  `num0'   `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1', `w1s', `w0s', `num1s', `num0s', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqips' `eqd0'  `eqy0' `eqy1' `eqw1' `eqw0' `num1'  `num0'  `num' `denom' if `touse'  [pw = `samplew'],   `eqips_inst' `eqd0_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}
else if "`stat'" == "latt" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `w1s', `w0s', `bd0', `by0',`denom1s',  `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm  `eqips' `eqw1' `eqw0' `eqd0'  `eqy0' `eqd1'  `eqy1'  `num' `denom' if `touse' [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqy0_inst' `eqd1_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial =  (`bips', `w1s', `w0s', `by0',`denom1s',  `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm  `eqips' `eqw1' `eqw0' `eqy0' `eqd1'  `eqy1'  `num' `denom' if `touse' [pw = `samplew'], `eqips_inst' `eqy0_inst' `eqd1_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial =  (`bips', `w1s', `w0s', `bd0', `by0', `num1s', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips' `eqw1' `eqw0' `eqd0'  `eqy0' `eqy1'  `num' `denom'  if `touse'  [pw = `samplew'],   `eqips_inst' `eqd0_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}
	}
set maxiter 300
scalar nob =  e(N)
scalar ncl = e(N_clust)
matrix vc= e(V)
local ry = rownumb(vc, "num:_cons")
local rd = rownumb(vc, "denom:_cons")
ereturn scalar ated = _b[denom:_cons]
ereturn scalar atey = _b[num:_cons]
ereturn scalar vard = vc[`rd', `rd'] 
ereturn scalar vary = vc[`ry', `ry'] 
scalar covdy= vc[`ry', `rd']
ereturn scalar seatey = sqrt(vc[`ry', `ry'] )
ereturn scalar seated = sqrt( vc[`rd', `rd'] )
scalar ated = _b[denom:_cons]
scalar atey = _b[num:_cons]
scalar vard = vc[`rd', `rd'] 
scalar vary = vc[`ry', `ry'] 
ereturn scalar ipwra_late = atey/ated
ereturn scalar ipwra_late_se = sqrt((1/ated)^2*vary+(atey/ated^2)^2*vard-2*(atey/ated^3)*covdy)
ereturn scalar ipwra_late_var =((1/ated)^2*vary+(atey/ated^2)^2*vard-2*(atey/ated^3)*covdy)
 scalar ipwra_late = atey/ated
 scalar ipwra_late_se = sqrt((1/ated)^2*vary+(atey/ated^2)^2*vard-2*(atey/ated^3)*covdy)
 scalar ipwra_late_var =((1/ated)^2*vary+(atey/ated^2)^2*vard-2*(atey/ated^3)*covdy)
	}  
ereturn clear
	// Display results
if "`stat'" == "late" {
	tempname b V
	matrix `b' = (atey, ated, ipwra_late)
	matrix `V' = (vary, covdy, 0 \ covdy, vard, 0  \ 0, 0,ipwra_late_var)
	/*matrix `b' = J(1,1,ipwra_late)
	matrix `V' = J(1,1,ipwra_late_var)*/
	matrix rownames `b' = " "
	matrix colnames `b' = "Numerator" "Denominator" "LATE"
	matrix rownames `V' = "Numerator" "Denominator" "LATE"
	matrix colnames `V' = "Numerator" "Denominator" "LATE"
	
	di as text "LATE Estimation"
	di as text "Outcome        : "  "`yvar'"
	di as text "Treatment      : "  "`tvar'"
	di as text "Instrument     : "  "`zvar'"
	if "`method'" == "ipwra" {
	di as text "Estimator      : DR LATE"
		    if "`omodel'" == "regress" {
		di as text "Outcome model  : linear"
	}
	else if "`omodel'" == "logit" {
		di as text "Outcome model  : logit"
	}
	else if "`omodel'" == "poisson" {
		di as text "Outcome model  : poisson"
	}
	    if "`tmodel'" == "regress" {
		di as text "Treatment model: linear"
	}
	else if "`tmodel'" == "logit" {
		di as text "Treatment model: logit"
	}
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
		if "`method'" == "ipw" {
	di as text "Estimator      : inverse-probability weights"
	di as text "Outcome model  : weighted mean"
	di as text "Treatment model: weighted mean"
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
	if "`method'" == "nipw" {
	di as text "Estimator      : normalized inverse-probability weights"
	di as text "Outcome model  : weighted mean"
	di as text "Treatment model: weighted mean"
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
		if "`method'" == "ra" {
	di as text "Estimator      : regression adjustment"
		    if "`omodel'" == "regress" {
		di as text "Outcome model  : linear"
	}
	else if "`omodel'" == "logit" {
		di as text "Outcome model  : logit"
	}
	else if "`omodel'" == "poisson" {
		di as text "Outcome model  : poisson"
	}
	    if "`tmodel'" == "regress" {
		di as text "Treatment model: linear"
	}
	else if "`tmodel'" == "logit" {
		di as text "Treatment model: logit"
	}
	 di as text "IPS Model      : none"

	}
	if "`method'" == "aipw" {
	di as text "Estimator      : augmented IPW"
		    if "`omodel'" == "regress" {
		di as text "Outcome model  : linear"
	}
	else if "`omodel'" == "logit" {
		di as text "Outcome model  : logit"
	}
	else if "`omodel'" == "poisson" {
		di as text "Outcome model  : poisson"
	}
	    if "`tmodel'" == "regress" {
		di as text "Treatment model: linear"
	}
	else if "`tmodel'" == "logit" {
		di as text "Treatment model: logit"
	}
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
		if "`method'" == "naipw" {
	di as text "Estimator      : Normalized augmented IPW"
		    if "`omodel'" == "regress" {
		di as text "Outcome model  : linear"
	}
	else if "`omodel'" == "logit" {
		di as text "Outcome model  : logit"
	}
	else if "`omodel'" == "poisson" {
		di as text "Outcome model  : poisson"
	}
	    if "`tmodel'" == "regress" {
		di as text "Treatment model: linear"
	}
	else if "`tmodel'" == "logit" {
		di as text "Treatment model: logit"
	}
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
    di as text "Number of obs = " nob
	
	
	di
	ereturn clear
	ereturn post `b' `V', esample(`touse')
	ereturn display
	
}
if "`stat'" == "latt" {  
	tempname b V
	matrix `b' = (atey, ated, ipwra_late)

	matrix `V' = (vary, covdy, 0 \ covdy, vard, 0  \ 0, 0,ipwra_late_var)
	matrix rownames `b' = " "
	matrix colnames `b' = "Numerator" "Denominator" "LATT"
	matrix rownames `V' = "Numerator" "Denominator" "LATT"
	matrix colnames `V' = "Numerator" "Denominator" "LATT"

	di as text "LATT Estimation"
	di as text "Outcome        : "  "`yvar'"
	di as text "Treatment      : "  "`tvar'"
	di as text "Instrument     : "  "`zvar'"
	if "`method'" == "ipwra" {
	di as text "Estimator      : DR LATT"
		    if "`omodel'" == "regress" {
		di as text "Outcome model  : linear"
	}
	else if "`omodel'" == "logit" {
		di as text "Outcome model  : logit"
	}
	else if "`omodel'" == "poisson" {
		di as text "Outcome model  : poisson"
	}
	    if "`tmodel'" == "regress" {
		di as text "Treatment model: linear"
	}
	else if "`tmodel'" == "logit" {
		di as text "Treatment model: logit"
	}
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
		if "`method'" == "ipw" {
	di as text "Estimator      : inverse-probability weights"
	di as text "Outcome model  : weighted mean"
	di as text "Treatment model: weighted mean"
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
	if "`method'" == "nipw" {
	di as text "Estimator      : normalized inverse-probability weights"
	di as text "Outcome model  : weighted mean"
	di as text "Treatment model: weighted mean"
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
		if "`method'" == "ra" {
	di as text "Estimator      : regression adjustment"
		    if "`omodel'" == "regress" {
		di as text "Outcome model  : linear"
	}
	else if "`omodel'" == "logit" {
		di as text "Outcome model  : logit"
	}
	else if "`omodel'" == "poisson" {
		di as text "Outcome model  : poisson"
	}
	    if "`tmodel'" == "regress" {
		di as text "Treatment model: linear"
	}
	else if "`tmodel'" == "logit" {
		di as text "Treatment model: logit"
	}
	 di as text "IPS Model      : none"

	}
	if "`method'" == "aipw" {
	di as text "Estimator      : augmented IPW"
		    if "`omodel'" == "regress" {
		di as text "Outcome model  : linear"
	}
	else if "`omodel'" == "logit" {
		di as text "Outcome model  : logit"
	}
	else if "`omodel'" == "poisson" {
		di as text "Outcome model  : poisson"
	}
	    if "`tmodel'" == "regress" {
		di as text "Treatment model: linear"
	}
	else if "`tmodel'" == "logit" {
		di as text "Treatment model: logit"
	}
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
		if "`method'" == "naipw" {
	di as text "Estimator      : Normalized augmented IPW"
		    if "`omodel'" == "regress" {
		di as text "Outcome model  : linear"
	}
	else if "`omodel'" == "logit" {
		di as text "Outcome model  : logit"
	}
	else if "`omodel'" == "poisson" {
		di as text "Outcome model  : poisson"
	}
	    if "`tmodel'" == "regress" {
		di as text "Treatment model: linear"
	}
	else if "`tmodel'" == "logit" {
		di as text "Treatment model: logit"
	}
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	}
    di as text "Number of obs = " nob
	
	 
	
	di
	ereturn clear
	ereturn post `b' `V', esample(`touse')
	ereturn display
}
	

end 
