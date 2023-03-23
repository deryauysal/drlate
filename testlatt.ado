

/* Under one sided non compliance test latt vs att using dr methods

 */
capture program drop testlatt

program testlatt, eclass 

version 17 

syntax  namelist(min=3) [if] [, xvarsy(varlist) xvarst(varlist) xvarsips(varlist) omodel(string) tmodel(string) zmodel(string) vce(string) samplew(varname) stat(string)]
	

	local orig `0'
	tokenize `namelist'
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
		{
		noi display as error "Linear regression not allowed"
		error 450
		exit
	    }
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
		local stat latt
	}	
		marksample touse			
	//quietly {
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
		if `dmeanz0' != 0 {
		noi display as error "No one-sided noncpmpliance"
		error 450
		exit
	    }
//Check this! 
	tempvar y1hat y1hatd y0hat y0hatd d1hat d0hat ips ps ipsxb psxb wr1 wr0 w0s w1 w0 w1s invw0 invw1 omega1 omega0 term1 term2 term3 term4 
	tempname bips bps by1 by1d by0 by0d bd1 bd0 denom1s denom0s num1s num0s nums  att1s att0s atts denoms lates diffs

	
	
set maxiter 300
if "`zmodel'" == "logit" {
  `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix `bips' = e(b)
 qui predict double `ips'
 }
else if "`zmodel'" == "cbps" {

  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew') stat(att)
   
   qui matrix `bips' = e(b)
   qui predict double `ipsxb'
   qui ge `ips' = logistic(`ipsxb')
  
}


if "`tmodel'" == "logit" {
  `tmodel' `tvar' `xvarst' if `touse'==1 [pw = `samplew']
 qui matrix `bps' = e(b)
 qui predict double `ps'
 }
else if "`zmodel'" == "cbps" {

  `tmodel' `tvar' `xvarst' if `touse'==1, samplew(`samplew') stat(att)
   
   qui matrix `bps' = e(b)
   qui predict double `psxb'
   qui ge `ps' = logistic(`psxb')
  
}



	
	if `dmeanz1'==1 {
		ge double `d1hat'= 1

	}
	else {
		
		qui regress  `tvar'  if `zvar' == 1 & `touse'==1 [pw = `samplew']
		matrix `bd1' = e(b)
		qui predict double `d1hat'
	}
	

		ge double `d0hat' = 0
	
	
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
	
	//ATT:
	
	qui regress `yvar' if `tvar' == 1 & `touse'==1 [pw = `samplew']
	matrix `by1d' = e(b)
	qui predict `y1hatd'
	
	qui `omodel' `yvar' `xvarsy'  if `tvar' == 0 & `touse'==1 [pw = `samplew'*(`ps'/(1-`ps'))]
	matrix `by0d' = e(b)
	qui predict double `y0hatd'
	
	

	qui sum `y1hatd'  if `tvar' == 1 & `touse'==1 [iw = `samplew']
	matrix `att1s' = r(mean)
	qui sum `y0hatd'  if `tvar' == 1 & `touse'==1 [iw = `samplew']
	matrix `att0s' = r(mean)
	matrix `atts' = `att1s'-`att0s'
    
	matrix `diffs' = `lates'-`atts'

	
	
	//Moment Functions
if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - ((exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*(1-`zvar')))

local eqips_inst instruments(eqips: `xvarsips')
}

	
if "`tmodel'" == "logit" {
local eqps (eqps:`tvar' - exp({that: `xvarst' _cons})/(1+exp({that:})))

local eqps_inst instruments(eqps: `xvarst')

}
else if "`tmodel'" == "cbps" {
local eqps (eqps: (`tvar' - ((exp({that: `xvarst' _cons})/(1+exp({that:})))/((1/(1+exp({that:})))))*(1-`tvar')))

local eqps_inst instruments(eqps: `xvarst')
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
//Moments for ATT

	
	local eqy1d (eqy1d: (1.`tvar'*(`yvar'-({y1d:  _cons})))) 
	
	local eqy1d_inst instruments(eqy1d:  )
	
	
	if "`omodel'" == "regress" {

		local eqy0d (eqy0d: ((0.`tvar'*((exp({that:})/(1+exp({that:})))/(1-(exp({that:})/(1+exp({that:}))))))*(`yvar'-({y0d: `xvarsy'  _cons})))) 
		
		local eqy0d_inst instruments(eqy0d: `xvarsy' )
		
		
		local att (att: (1.`tvar'*({att}-(({y1d: })-({y0d: })))))
		
	}
	else if "`omodel'"  == "logit" {
		
		
		local eqy0d (eqy0d: ((0.`tvar'*((exp({that:})/(1+exp({that:})))/(1-(exp({that:})/(1+exp({that:}))))))*(`yvar'-exp({y0d: `xvarsy' _cons})/(1+exp({y0d:}))))) 
		
		local eqy0d_inst instruments(eqy0d: `xvarsy' )
		
		local att (att: (1.`tvar'*({att}-(({y1d: })-(exp({y0d: })/(1+exp({y0d:})))))))
		
		
	}
	else if "`omodel'"  == "poisson"  {
		local eqy0d (eqy0d: ((0.`tvar'*((exp({that:})/(1+exp({that:})))/(1-(exp({that:})/(1+exp({that:}))))))*(`yvar'-exp({y0d: `xvarsy' _cons})))) 
		
		local eqy0d_inst instruments(eqy0d: `xvarsy' )
		
		local att (att: (1.`tvar'*({att}-(({y1d: })-(exp({y0d: }))))))
		
		
	}

 	local diff (diff: ({diff}-({late}-{att})))	
	
	

matrix initial = (`bips', `bps', `by0', `by1', `nums', `denoms', `lates', `by1d', `by0d', `atts', `diffs')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm   `eqips' `eqps' `eqy0' `eqy1'  `num' `denom' `late' `eqy1d'  `eqy0d' `att' `diff' if `touse' [pw = `samplew'],  `eqips_inst'  `eqps_inst' `eqy0_inst' `eqy1_inst' `eqy1d_inst' `eqy0d_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	

//	}

end
