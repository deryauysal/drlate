

/*Ado-file for the LATE estimators*/
/* DR-LATE project!*/
/*Version: 10082022*/
/*Author: DU
ThisVersion: IPWRA LATE & LATT
RA LATE & LATT
IPW LATE & LATT
AIPW -LATE
NAIPW - LATE
To do:
AIPW - LATT
NAIPW -LATT
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
		local samplew 1
	}
	if "`vce'" == "" {
		local vce robust
	}
		
	if "`stat'" == "" {
		local stat late
	}	
		marksample touse			
	quietly {
		
		tempvar troot 
		tempname valuest
		tab `tvar' if `touse', gen(`troot') matrow(`valuest')
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
		tab `zvar' if `touse', gen(`zroot') matrow(`valuesz')
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
		tab `yvar' if `touse', gen(`yroot') matrow(`valuesy')
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
//Check the compliance:
	    quietly summarize `tvar' if `zvar' == 1 & `touse'==1
		scalar dmeanz1 = r(mean)
		quietly summarize `tvar' if `zvar'  == 0 & `touse'==1
		scalar dmeanz0 = r(mean)
		

	if "`method'" == "ipwra" {
	//IPS
set maxiter 300

 qui `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix bips = e(b)
 qui predict ips

  


//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar' `xvarsy' if `zvar'==1 & `touse'==1 [pw = `samplew'/ips]
qui matrix by1 = e(b)
qui predict y1hat

 
 
 qui `omodel' `yvar' `xvarsy' if `zvar'==0 & `touse'==1 [pw = `samplew'/(1-ips)]
qui matrix by0 = e(b)
qui predict y0hat

//Conditional mean of d1 & d0#


	    if dmeanz1==1 {
        ge d1hat= 1
	    }
		else {

qui `tmodel' `tvar' `xvarst' if `zvar'==1  & `touse'==1  [pw =  `samplew'/ips]
matrix bd1 = e(b)
qui predict d1hat
}
	    
	    if dmeanz0==0 {
        ge d0hat= 0
	    }
		else {
qui `tmodel' `tvar' `xvarst' if `zvar'==0  & `touse'==1 [pw =  `samplew'/(1-ips)]
matrix bd0 = e(b)
qui predict d0hat
		}
		
		
qui sum d1hat if `touse'==1 [iw = `samplew']
matrix denom1s = r(mean)
qui sum d0hat  if `touse'==1 [iw = `samplew']
matrix denom0s = r(mean)
qui sum y1hat if `touse'==1  [iw = `samplew']
matrix num1s = r(mean)
qui sum y0hat if `touse'==1 [iw = `samplew']
matrix num0s = r(mean)
matrix nums = num1s-num0s
matrix denoms = denom1s-denom0s
matrix lates =nums*invsym(denoms)

drop y1hat y0hat d1hat d0hat ips
	
//Moment Functions
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')




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
	if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}- ((exp({d1:})/(1+exp({d1:})))))
	}
	else if dmeanz0 != 0 & dmeanz1 ==1 {
local denom (denom: {denom}- (1-((exp({d0:})/(1+exp({d0:}))))))
	}
	else {
local denom (denom: {denom}- (((exp({d1:})/(1+exp({d1:}))))-((exp({d0:})/(1+exp({d0:}))))))
	}
}
else if "`tmodel'" =="regress"  {
		if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}- ({d1:}))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
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
	
	if dmeanz1==1 {
		ge d1hat= 1
	}
	else {
		
		qui regress  `tvar'  if `zvar' == 1 & `touse'==1 [pw = `samplew']
		matrix bd1 = e(b)
		qui predict d1hat
	}
	
	if dmeanz0==0 {
		ge d0hat= 0
	}
	else {
		qui `tmodel' `tvar' `xvarst'  if `zvar' == 0 & `touse'==1 [pw = `samplew'*(ips/(1-ips))]
		matrix bd0 = e(b)
		qui predict d0hat
	}
	
	qui regress `yvar' if `zvar' == 1 & `touse'==1 [pw = `samplew']
	matrix by1 = e(b)
	qui predict y1hat
	
	qui `omodel' `yvar' `xvarsy'  if `zvar' == 0 & `touse'==1 [pw = `samplew'*(ips/(1-ips))]
	matrix by0 = e(b)
	qui predict y0hat
	
	
	qui sum d1hat  if `zvar' == 1 [iw = `samplew']
	matrix denom1s = r(mean)
	qui sum d0hat  if `zvar' == 1 [iw = `samplew']
	matrix denom0s = r(mean)
	qui sum y1hat  if `zvar' == 1 [iw = `samplew']
	matrix num1s = r(mean)
	qui sum y0hat  if `zvar' == 1 [iw = `samplew']
	matrix num0s = r(mean)
	matrix nums = num1s-num0s
	matrix denoms = denom1s-denom0s
	
	
	matrix lates =nums*invsym(denoms)
	
	drop y1hat y0hat d1hat d0hat ips
	
	
	
	//Moment Functions
	local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))
	
	local eqips_inst instruments(eqips: `xvarsips')
	
	


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
		if dmeanz0 == 0 & dmeanz1 !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-(exp({d0: })/(1+exp({d0:})))))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}-(({d1: })-(exp({d0: })/(1+exp({d0:})))))))
		}
	}
	else if "`tmodel'" =="regress"  {
		if dmeanz0 == 0 & dmeanz1 !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
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
if dmeanz0 !=0 & dmeanz1 !=1  {
matrix initial = (bips, bd0, bd1, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
qui gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if dmeanz0 == 0 & dmeanz1 !=1 {
matrix initial = (bips, bd1, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
qui gmm   `eqips'  `eqd1'  `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if dmeanz0 !=0 & dmeanz1 ==1 {
matrix initial = (bips, bd0, by0, by1, nums, denoms)	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
qui gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	 
if "`stat'" == "latt" {
if dmeanz0 !=0 & dmeanz1 !=1  {
matrix initial = (bips, bd0, bd1, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
qui gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if dmeanz0 == 0 & dmeanz1 !=1 {
matrix initial = (bips, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm   `eqips'  `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],  `eqips_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if dmeanz0 !=0 & dmeanz1 ==1 {
matrix initial = (bips, bd0, by0, by1, nums, denoms)	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
qui gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}
	}

	else if "`method'" == "ipw" {
			//IPS
set maxiter 300

 qui `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix bips = e(b)
 qui predict ips

		
//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar'  if `zvar'==1 & `touse'==1 [pw = `samplew'/ips]
qui matrix by1 = e(b)
qui predict y1hat

 
 
 qui `omodel' `yvar'  if `zvar'==0 & `touse'==1 [pw = `samplew'/(1-ips)]
qui matrix by0 = e(b)
qui predict y0hat

//Conditional mean of d1 & d0#


	    if dmeanz1==1 {
        ge d1hat= 1
	    }
		else {

qui `tmodel' `tvar'  if `zvar'==1  & `touse'==1  [pw =  `samplew'/ips]
matrix bd1 = e(b)
qui predict d1hat
}
	    
	    if dmeanz0==0 {
        ge d0hat= 0
	    }
		else {
qui `tmodel' `tvar'  if `zvar'==0  & `touse'==1 [pw =  `samplew'/(1-ips)]
matrix bd0 = e(b)
qui predict d0hat
		}
		
		
qui sum d1hat if `touse'==1 [iw = `samplew']
matrix denom1s = r(mean)
qui sum d0hat  if `touse'==1 [iw = `samplew']
matrix denom0s = r(mean)
qui sum y1hat if `touse'==1  [iw = `samplew']
matrix num1s = r(mean)
qui sum y0hat if `touse'==1 [iw = `samplew']
matrix num0s = r(mean)
matrix nums = num1s-num0s
matrix denoms = denom1s-denom0s
matrix lates =nums*invsym(denoms)

drop y1hat y0hat d1hat d0hat ips
	
//Moment Functions
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')




if "`tmodel'" == "logit" {
local eqd1 (eqd1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`tvar' - exp({d1:  _cons})/(1+exp({d1:})))))
           
local eqd1_inst instruments(eqd1: )

local eqd0 (eqd0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`tvar'-exp({d0:  _cons})/(1+exp({d0:}))))) 

local eqd0_inst instruments(eqd0:  )
}


if "`tmodel'" =="regress" {
local eqd1 (eqd1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`tvar' - ({d1:  _cons}))))
           
local eqd1_inst instruments(eqd1: )

local eqd0 (eqd0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`tvar'-({d0:  _cons})))) 

local eqd0_inst instruments(eqd0:  )
}

if "`tmodel'" == "logit" {
	if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}- ((exp({d1:})/(1+exp({d1:})))))
	}
	else if dmeanz0 != 0 & dmeanz1 ==1 {
local denom (denom: {denom}- (1-((exp({d0:})/(1+exp({d0:}))))))
	}
	else {
local denom (denom: {denom}- (((exp({d1:})/(1+exp({d1:}))))-((exp({d0:})/(1+exp({d0:}))))))
	}
}
else if "`tmodel'" =="regress"  {
		if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}- ({d1:}))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
local denom (denom: {denom}- (1-({d0:})))
		}
		else {
local denom (denom: {denom}- (({d1:})-({d0:})))
		}
}



if "`omodel'" == "regress" {
local eqy0 (eqy0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`yvar'-({y0:   _cons})))) 

local eqy0_inst instruments(eqy0:  )

local eqy1 (eqy1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`yvar'-({y1:   _cons})))) 

local eqy1_inst instruments(eqy1:  )

local num (num: {num}-(({y1:})-({y0:})))
}
else if "`omodel'"  == "logit" {
local eqy0 (eqy0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`yvar'-exp({y0:  _cons})/(1+exp({y0:}))))) 

local eqy0_inst instruments(eqy0:  )

local eqy1 (eqy1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`yvar'-exp({y1:  _cons})/(1+exp({y1:}))))) 

local eqy1_inst instruments(eqy1:  )

local num (num: {num}-(((exp({y1: })/(1+exp({y1:}))))-((exp({y0: })/(1+exp({y0:}))))))

}
else if "`omodel'"  == "poisson"  {
local eqy0 (eqy0: ((0.`zvar'/(1-(exp({zhat:})/(1+exp({zhat:})))))*(`yvar'-exp({y0:   _cons})))) 

local eqy0_inst instruments(eqy0:  )

local eqy1 (eqy1: ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))*(`yvar'-exp({y1:   _cons})))) 

local eqy1_inst instruments(eqy1:  )

local num (num: {num}-((exp({y1: }))-(exp({y0: }))))

}

local late (late: ({late} - {num}/{denom}))

} 
	 
if "`stat'" == "latt" {
	
	

	
	if dmeanz1==1 {
		ge d1hat= 1
	}
	else {
		
		qui regress  `tvar'  if `zvar' == 1 & `touse'==1 [pw = `samplew']
		matrix bd1 = e(b)
		qui predict d1hat
	}
	
	if dmeanz0==0 {
		ge d0hat= 0
	}
	else {
		qui `tmodel' `tvar'   if `zvar' == 0 & `touse'==1 [pw = `samplew'*(ips/(1-ips))]
		matrix bd0 = e(b)
		qui predict d0hat
	}
	
	qui regress `yvar' if `zvar' == 1 & `touse'==1 [pw = `samplew']
	matrix by1 = e(b)
	qui predict y1hat
	
	qui `omodel' `yvar'  if `zvar' == 0 & `touse'==1 [pw = `samplew'*(ips/(1-ips))]
	matrix by0 = e(b)
	qui predict y0hat
	
	
	qui sum d1hat  if `zvar' == 1 [iw = `samplew']
	matrix denom1s = r(mean)
	qui sum d0hat  if `zvar' == 1 [iw = `samplew']
	matrix denom0s = r(mean)
	qui sum y1hat  if `zvar' == 1 [iw = `samplew']
	matrix num1s = r(mean)
	qui sum y0hat  if `zvar' == 1 [iw = `samplew']
	matrix num0s = r(mean)
	matrix nums = num1s-num0s
	matrix denoms = denom1s-denom0s
	
	
	matrix lates =nums*invsym(denoms)
	
	drop y1hat y0hat d1hat d0hat ips
	
	
	
	//Moment Functions
	local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))
	
	local eqips_inst instruments(eqips: `xvarsips')
	
	


	local eqd1 (eqd1: ((1.`zvar'*(`tvar' - ({d1: _cons})))))
	
	local eqd1_inst instruments(eqd1: )
	
	
	if "`tmodel'" == "logit" {
		
		local eqd0 (eqd0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`tvar'-exp({d0:  _cons})/(1+exp({d0:}))))) 
		
		local eqd0_inst instruments(eqd0:  )
	}
	
	
	if "`tmodel'" =="regress" {
		
		local eqd0 (eqd0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`tvar'-({d0:  _cons})))) 
		
		local eqd0_inst instruments(eqd0:  )
	}
	

	
	
	if "`tmodel'" == "logit" {
		if dmeanz0 == 0 & dmeanz1 !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-(exp({d0: })/(1+exp({d0:})))))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}-(({d1: })-(exp({d0: })/(1+exp({d0:})))))))
		}
	}
	else if "`tmodel'" =="regress"  {
		if dmeanz0 == 0 & dmeanz1 !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-({d0: })))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}-(({d1: })-({d0: })))))
		}
	}
	

	
	local eqy1 (eqy1: (1.`zvar'*(`yvar'-({y1:  _cons})))) 
	
	local eqy1_inst instruments(eqy1:  )
	
	
	
	if "`omodel'" == "regress" {

		local eqy0 (eqy0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`yvar'-({y0:   _cons})))) 
		
		local eqy0_inst instruments(eqy0:  )
		
		
		local num (num: (1.`zvar'*({num}-(({y1: })-({y0: })))))
		
	}
	else if "`omodel'"  == "logit" {
		
		
		local eqy0 (eqy0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`yvar'-exp({y0:  _cons})/(1+exp({y0:}))))) 
		
		local eqy0_inst instruments(eqy0:  )
		
		local num (num: (1.`zvar'*({num}-(({y1: })-(exp({y0: })/(1+exp({y0:})))))))
		
		
	}
	else if "`omodel'"  == "poisson"  {
		local eqy0 (eqy0: ((0.`zvar'*((exp({zhat:})/(1+exp({zhat:})))/(1-(exp({zhat:})/(1+exp({zhat:}))))))*(`yvar'-exp({y0:  _cons})))) 
		
		local eqy0_inst instruments(eqy0:  )
		
		local num (num: (1.`zvar'*({num}-(({y1: })-(exp({y0: }))))))
		
		
	}
	
	local late (late: ({late} - {num}/{denom}))
}	 
	 
	 
	 if "`stat'" == "late" {
if dmeanz0 !=0 & dmeanz1 !=1  {
matrix initial = (bips, bd0, bd1, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
qui gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if dmeanz0 == 0 & dmeanz1 !=1 {
matrix initial = (bips, bd1, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
qui gmm   `eqips'  `eqd1'  `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if dmeanz0 !=0 & dmeanz1 ==1 {
matrix initial = (bips, bd0, by0, by1, nums, denoms)	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
qui gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	 
if "`stat'" == "latt" {
if dmeanz0 !=0 & dmeanz1 !=1  {
matrix initial = (bips, bd0, bd1, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
qui gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if dmeanz0 == 0 & dmeanz1 !=1 {
matrix initial = (bips, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm   `eqips'  `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],  `eqips_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if dmeanz0 !=0 & dmeanz1 ==1 {
matrix initial = (bips, bd0, by0, by1, nums, denoms)	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
qui gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	
	}
	else if "`method'" == "ra" {
	//IPS
set maxiter 300


//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar' `xvarsy' if `zvar'==1 & `touse'==1 [pw = `samplew']
qui matrix by1 = e(b)
qui predict y1hat

 
 
 qui `omodel' `yvar' `xvarsy' if `zvar'==0 & `touse'==1 [pw = `samplew']
qui matrix by0 = e(b)
qui predict y0hat

//Conditional mean of d1 & d0#


	    if dmeanz1==1 {
        ge d1hat= 1
	    }
		else {

qui `tmodel' `tvar' `xvarst' if `zvar'==1  & `touse'==1  [pw =  `samplew']
matrix bd1 = e(b)
qui predict d1hat
}
	    
	    if dmeanz0==0 {
        ge d0hat= 0
	    }
		else {
qui `tmodel' `tvar' `xvarst' if `zvar'==0  & `touse'==1 [pw =  `samplew']
matrix bd0 = e(b)
qui predict d0hat
		}
		
		
qui sum d1hat if `touse'==1 [iw = `samplew']
matrix denom1s = r(mean)
qui sum d0hat  if `touse'==1 [iw = `samplew']
matrix denom0s = r(mean)
qui sum y1hat if `touse'==1  [iw = `samplew']
matrix num1s = r(mean)
qui sum y0hat if `touse'==1 [iw = `samplew']
matrix num0s = r(mean)
matrix nums = num1s-num0s
matrix denoms = denom1s-denom0s
matrix lates =nums*invsym(denoms)

drop y1hat y0hat d1hat d0hat 
	
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
	if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}- ((exp({d1:})/(1+exp({d1:})))))
	}
	else if dmeanz0 != 0 & dmeanz1 ==1 {
local denom (denom: {denom}- (1-((exp({d0:})/(1+exp({d0:}))))))
	}
	else {
local denom (denom: {denom}- (((exp({d1:})/(1+exp({d1:}))))-((exp({d0:})/(1+exp({d0:}))))))
	}
}
else if "`tmodel'" =="regress"  {
		if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}- ({d1:}))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
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
	
	if dmeanz1==1 {
		ge d1hat= 1
	}
	else {
		
		qui regress  `tvar'  if `zvar' == 1 & `touse'==1 [pw = `samplew']
		matrix bd1 = e(b)
		qui predict d1hat
	}
	
	if dmeanz0==0 {
		ge d0hat= 0
	}
	else {
		qui `tmodel' `tvar' `xvarst'  if `zvar' == 0 & `touse'==1 [pw = `samplew']
		matrix bd0 = e(b)
		qui predict d0hat
	}
	
	qui regress `yvar' if `zvar' == 1 & `touse'==1 [pw = `samplew']
	matrix by1 = e(b)
	qui predict y1hat
	
	qui `omodel' `yvar' `xvarsy'  if `zvar' == 0 & `touse'==1 [pw = `samplew']
	matrix by0 = e(b)
	qui predict y0hat
	
	
	qui sum d1hat  if `zvar' == 1 [iw = `samplew']
	matrix denom1s = r(mean)
	qui sum d0hat  if `zvar' == 1 [iw = `samplew']
	matrix denom0s = r(mean)
	qui sum y1hat  if `zvar' == 1 [iw = `samplew']
	matrix num1s = r(mean)
	qui sum y0hat  if `zvar' == 1 [iw = `samplew']
	matrix num0s = r(mean)
	matrix nums = num1s-num0s
	matrix denoms = denom1s-denom0s
	
	
	matrix lates =nums*invsym(denoms)
	
	drop y1hat y0hat d1hat d0hat 
	
	
	
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
		if dmeanz0 == 0 & dmeanz1 !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-(exp({d0: })/(1+exp({d0:})))))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}-(({d1: })-(exp({d0: })/(1+exp({d0:})))))))
		}
	}
	else if "`tmodel'" =="regress"  {
		if dmeanz0 == 0 & dmeanz1 !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
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
if dmeanz0 !=0 & dmeanz1 !=1  {
matrix initial = (bd0, bd1, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm  `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],   `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if dmeanz0 == 0 & dmeanz1 !=1 {
matrix initial = (bd1, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqd1'  `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],    `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if dmeanz0 !=0 & dmeanz1 ==1 {
matrix initial = (bips, bd0, by0, by1, nums, denoms)	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqd0' `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],   `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	 
if "`stat'" == "latt" {
if dmeanz0 !=0 & dmeanz1 !=1  {
matrix initial = (bd0, bd1, by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],   `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if dmeanz0 == 0 & dmeanz1 !=1 {
matrix initial = ( by0, by1, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm    `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],   `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if dmeanz0 !=0 & dmeanz1 ==1 {
matrix initial = (bd0, by0, by1, nums, denoms)	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqd0' `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],   `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}		
		
		
	}	
	else if "`method'" == "aipw" {
		
			//IPS
set maxiter 300

 qui `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix bips = e(b)
 qui predict ips
		
//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar' `xvarsy' if `zvar'==1 & `touse'==1 [pw = `samplew']
qui matrix by1 = e(b)
qui predict y1hat

 
 
 qui `omodel' `yvar' `xvarsy' if `zvar'==0 & `touse'==1 [pw = `samplew']
qui matrix by0 = e(b)
qui predict y0hat

//Conditional mean of d1 & d0#


	    if dmeanz1==1 {
        ge d1hat= 1
	    }
		else {

qui `tmodel' `tvar' `xvarst' if `zvar'==1  & `touse'==1  [pw =  `samplew']
matrix bd1 = e(b)
qui predict d1hat
}
	    
	    if dmeanz0==0 {
        ge d0hat= 0
	    }
		else {
qui `tmodel' `tvar' `xvarst' if `zvar'==0  & `touse'==1 [pw =  `samplew']
matrix bd0 = e(b)
qui predict d0hat
		}
		
	qui gen term1 = ((`zvar'*`yvar')-((`zvar'-ips)*y1hat))/ips
	qui gen term2 = (((1-`zvar')*`yvar')+((`zvar'-ips)*y0hat))/(1-ips)
	qui gen term3 = ((`zvar'*`tvar')-((`zvar'-ips)*d1hat))/ips
	qui gen term4 = (((1-`zvar')*`tvar')+((`zvar'-ips)*d0hat))/(1-ips)
	
	
		qui sum term1
	scalar t1 = r(mean)
	qui sum term2
	scalar t2 = r(mean)
	qui sum term3
	scalar t3 = r(mean)
	qui sum term4
	scalar t4 = r(mean)

qui sum term1
matrix num1s = r(mean)
qui sum term2
matrix num0s = r(mean)
matrix nums = num1s-num0s
qui sum term3
matrix denoms1 = r(mean)
qui sum term4
matrix denoms0 = r(mean)
matrix denoms = denoms1-denoms0
matrix lates =nums*invsym(denoms)

drop y1hat y0hat d1hat d0hat term* ips

//Moment Functions
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')


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
	if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}- ((((`zvar'*`tvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({d1: })/(1+exp({d1:}))))))/(exp({zhat:})/(1+exp({zhat:})))) )
	}
	else if dmeanz0 != 0 & dmeanz1 ==1 {
local denom (denom: {denom}- (1-((((1-`zvar')*`tvar')+((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({d0: })/(1+exp({d0:}))))) / (1-(exp({zhat:})/(1+exp({zhat:})))) ) ))
	}
	else {
local denom1 ( denom1: {denom1}- ((((`zvar'*`tvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({d1: })/(1+exp({d1:}))))))/(exp({zhat:})/(1+exp({zhat:})))) )
local denom0 ( denom0: {denom0}- ((((1-`zvar')*`tvar')+((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*(exp({d0: })/(1+exp({d0:}))))) / (1-(exp({zhat:})/(1+exp({zhat:})))) ) )
local denom (denom: {denom}-({denom1}-{denom0}))
	}
}
else if "`tmodel'" =="regress"  {
		if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}- (((`zvar'*`tvar')-((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))*{d1:}))/(exp({zhat:})/(1+exp({zhat:})))))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
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
	if "`stat'" == "late" {
if dmeanz0 !=0 & dmeanz1 !=1  {
matrix initial = (bips, bd0, bd1, by0, by1, num1s, num0s, denom1s, denom0s, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm  `eqips' `eqd0' `eqd1'  `eqy0' `eqy1' `num1'  `num0' `denom1' `denom0'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if dmeanz0 == 0 & dmeanz1 !=1 {
matrix initial =  (bips, bd1, by0, by1,  num1s, num0s, nums, denoms)
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1' `num1'  `num0'   `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if dmeanz0 !=0 & dmeanz1 ==1 {
matrix initial = (bips, bd0, by0, by1, num1s, num0s, nums, denoms)	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqips' `eqd0'  `eqy0' `eqy1' `num1'  `num0'  `num' `denom'   [pw = `samplew'],   `eqips_inst' `eqd0_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
if "`stat'" == "latt" {
	drop ips
	error 198	
}	
	}
		else if "`method'" == "naipw" {
			
			//IPS
set maxiter 300

 qui `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix bips = e(b)
 qui predict ips
		
//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar' `xvarsy' if `zvar'==1 & `touse'==1 [pw = `samplew']
qui matrix by1 = e(b)
qui predict y1hat

 
 
 qui `omodel' `yvar' `xvarsy' if `zvar'==0 & `touse'==1 [pw = `samplew']
qui matrix by0 = e(b)
qui predict y0hat

//Conditional mean of d1 & d0#


	    if dmeanz1==1 {
        ge d1hat= 1
	    }
		else {

qui `tmodel' `tvar' `xvarst' if `zvar'==1  & `touse'==1  [pw =  `samplew']
matrix bd1 = e(b)
qui predict d1hat
}
	    
	    if dmeanz0==0 {
        ge d0hat= 0
	    }
		else {
qui `tmodel' `tvar' `xvarst' if `zvar'==0  & `touse'==1 [pw =  `samplew']
matrix bd0 = e(b)
qui predict d0hat
		}
		
	qui gen invw1 = `zvar'/ips
	qui gen invw0 = (1-`zvar')/(1-ips)
	
	qui sum invw1
	scalar w1s = r(mean)
	qui sum invw0
	scalar w0s = r(mean)
	
	qui ge omega1= (`zvar'/ips)/w1s
	qui ge omega0= ((1-`zvar')/(1-ips))/w0s
	
	qui ge term1= omega1*( `yvar' -y1hat)+y1hat
	qui sum term1 
	matrix num1s =  r(mean)
	qui ge term2= omega0*( `yvar' -y0hat)+y0hat
	qui sum term2
	matrix num0s =  r(mean)
	qui ge term3= omega1*( `tvar' -d1hat)+d1hat
	qui sum term3 
	matrix denom1s =  r(mean)	
	qui ge term4= omega0*( `tvar' -d0hat)+d0hat
	qui sum term4
	matrix denom0s =  r(mean)
	
	matrix nums = num1s-num0s
    matrix denoms = denom1s-denom0s
matrix lates =nums*invsym(denoms)
		
drop ips invw0 invw1 y1hat y0hat omega1 omega0 term* d1hat
		
		
//Moment Functions
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

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
	if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}- ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`tvar'-((exp({d1: })/(1+exp({d1:})))))-((exp({d1: })/(1+exp({d1:})))))
	}
	else if dmeanz0 != 0 & dmeanz1 ==1 {
local denom (denom: {denom}- (1-((0.`zvar'/(1-exp({zhat:})/(1+exp({zhat:}))))/{w0:})*(`tvar'-((exp({d0: })/(1+exp({d0:})))))-((exp({d0: })/(1+exp({d0:}))))))
	}
	else {
local denom1 (denom1: {denom1}- ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`tvar'-((exp({d1: })/(1+exp({d1:})))))-((exp({d1: })/(1+exp({d1:})))))
local denom0 (denom0: {denom0}- ((0.`zvar'/(1-exp({zhat:})/(1+exp({zhat:}))))/{w0:})*(`tvar'-((exp({d0: })/(1+exp({d0:})))))-((exp({d0: })/(1+exp({d0:})))))
local denom (denom: {denom}-({denom1}-{denom0}))
	}
}
else if "`tmodel'" =="regress"  {
		if dmeanz0 == 0 & dmeanz1 !=1 {
local denom (denom: {denom}-  ((1.`zvar'/(exp({zhat:})/(1+exp({zhat:}))))/{w1:})*(`tvar'-({d1:}))-({d1:}))
		}
		else if dmeanz0 != 0 & dmeanz1 ==1 {
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
	if "`stat'" == "late" {
if dmeanz0 !=0 & dmeanz1 !=1  {
matrix initial = (bips, bd0, bd1, by0, by1,w1s, w0s, num1s, num0s, denom1s, denom0s, nums, denoms)	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm  `eqips' `eqd0' `eqd1'  `eqy0' `eqy1' `eqw1' `eqw0' `num1'  `num0' `denom1' `denom0'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if dmeanz0 == 0 & dmeanz1 !=1 {
matrix initial =  (bips, bd1, by0, by1,w1s, w0s,  num1s, num0s, nums, denoms)
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1' `eqw1' `eqw0' `num1'  `num0'   `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if dmeanz0 !=0 & dmeanz1 ==1 {
matrix initial = (bips, bd0, by0, by1,w1s, w0s, num1s, num0s, nums, denoms)	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqips' `eqd0'  `eqy0' `eqy1' `eqw1' `eqw0' `num1'  `num0'  `num' `denom'   [pw = `samplew'],   `eqips_inst' `eqd0_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}
if "`stat'" == "latt" {
	drop ips
	error 198	
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
	di as text "IPS Model      : logit"
	}
	if "`method'" == "ipw" {
	di as text "Estimator      : inverse-probability weights"
	di as text "Outcome model  : weighted mean"
	di as text "Treatment model: weighted mean"
	di as text "IPS Model      : logit"
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
	di as text "IPS Model      : logit"
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
	di as text "IPS Model      : logit"
	}
	if "`method'" == "ipw" {
	di as text "Estimator      : inverse-probability weights"
	di as text "Outcome model  : weighted mean"
	di as text "Treatment model: weighted mean"
	di as text "IPS Model      : logit"
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
    di as text "Number of obs = " nob
	 
	
	di
	ereturn clear
	ereturn post `b' `V', esample(`touse')
	ereturn display
}
	

end 