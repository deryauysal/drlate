	else if "`method'" == "aipw" {
		
			//IPS
set maxiter 300

 qui `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix `bips' = e(b)
 qui predict double `ips'
		
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
	
	
	qui sum `term1'
	scalar t1 = r(mean)
	qui sum `term2'
	scalar t2 = r(mean)
	qui sum `term3'
	scalar t3 = r(mean)
	qui sum `term4'
	scalar t4 = r(mean)

qui sum `term1'
matrix `num1s' = r(mean)
qui sum `term2'
matrix `num0s' = r(mean)
matrix `nums' = `num1s'-`num0s'
qui sum `term3'
matrix `denom1s' = r(mean)
qui sum `term4'
matrix `denom0s' = r(mean)
matrix `denoms' = `denom1s'-`denom0s'
matrix `lates' =`nums'*invsym(`denoms')

//drop y1hat y0hat d1hat d0hat term* ips

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
	
	//qui gen double `term1' = ((`zvar'*`yvar')-((`zvar'-`ips')*`y1hat'))/`ips'
	qui gen double `term2'  = (`ips'*(1-`zvar')*(`yvar'-`y0hat'))/(1-`ips')+`zvar'*`y0hat'
	//qui gen double `term3'  = ((`zvar'*`tvar')-((`zvar'-`ips')*`d1hat'))/`ips'
	qui gen double `term4'  = (`ips'*(1-`zvar')*(`tvar'-`d0hat'))/(1-`ips')+`zvar'*`d0hat'
	
	qui sum `d1hat'  if `zvar' == 1 [iw = `samplew']
	matrix `denom1s' = r(mean)
	qui sum `term4'  if `zvar' == 1 [iw = `samplew']
	matrix `denom0s' = r(mean)
	qui sum `y1hat'  if `zvar' == 1 [iw = `samplew']
	matrix `num1s' = r(mean)
	qui sum `term2'  if `zvar' == 1 [iw = `samplew']
	matrix `num0s' = r(mean)
	matrix `nums' = `num1s'-`num0s'
	matrix `denoms' = `denom1s'-`denom0s'
	
	
	matrix `lates' =`nums'*invsym(`denoms')
	
	

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
			local denom (denom: (1.`zvar'*({denom}-(1-(`zvar'*((exp({d0: })/(1+exp({d0:}))))+ ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-(exp({d0: })/(1+exp({d0:})))))/(1/(1+exp({zhat:}))))))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}  -(`zvar'*(`tvar'-(exp({d0: })/(1+exp({d0:}))))+ ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-(exp({d0: })/(1+exp({d0:})))))/(1/(1+exp({zhat:})))))))
		}
	}
	else if "`tmodel'" =="regress"  {
		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: ((1.`zvar'*(`tvar' - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (1.`zvar'*({denom}-(1-((`zvar'*({d0:})+ ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-{d0:}))/(1/(1+exp({zhat:})))))))))
		}
		else {
			local denom (denom: (1.`zvar'*({denom}  -(`zvar'*(`tvar'-{d0:})+ ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`tvar'-{d0:}))/(1/(1+exp({zhat:})))))))
		}
	}
	

	
	local eqy1 (eqy1: (1.`zvar'*(`yvar'-({y1:  _cons})))) 
	
	local eqy1_inst instruments(eqy1:  )
	
	
	
	if "`omodel'" == "regress" {

		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-({y0: `xvarsy'  _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )
		
		local num (num: (1.`zvar'*({num}  -(`zvar'*(`yvar'-{y0:})+ ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`yvar'-{y0:}))/(1/(1+exp({zhat:})))))))		
	}
	else if "`omodel'"  == "logit" {
		
		
		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})/(1+exp({y0:}))))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )

		local num (num: (1.`zvar'*({num}  -(`zvar'*(`yvar'-(exp({y0: })/(1+exp({y0:}))))+ ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`yvar'-(exp({y0: })/(1+exp({y0:})))))/(1/(1+exp({zhat:})))))))
		
	}
	else if "`omodel'"  == "poisson"  {
		local eqy0 (eqy0: ((0.`zvar')*(`yvar'-exp({y0: `xvarsy' _cons})))) 
		
		local eqy0_inst instruments(eqy0: `xvarsy' )

		local num (num: (1.`zvar'*({num}  -(`zvar'*(`yvar'-(exp({y0: })))+ ((exp({zhat:})/(1+exp({zhat:})))*(1-`zvar')*(`yvar'-(exp({y0: }))))/(1/(1+exp({zhat:})))))))		
		
	}
	
	local late (late: ({late} - {num}/{denom}))
	
} 
	if "`stat'" == "late" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', `num1s', `num0s', `denom1s', `denom0s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm  `eqips' `eqd0' `eqd1'  `eqy0' `eqy1' `num1'  `num0' `denom1' `denom0'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial =  (`bips', `bd1', `by0', `by1',  `num1s', `num0s', `nums', `denoms')
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1' `num1'  `num0'   `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1', `num1s', `num0s', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqips' `eqd0'  `eqy0' `eqy1' `num1'  `num0'  `num' `denom'   [pw = `samplew'],   `eqips_inst' `eqd0_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}
else if "`stat'" == "latt" {  //bence bd1, by1, ie eqd1 ve eqy1 are reluctant
	if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm   `eqips'  `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],  `eqips_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	
	}