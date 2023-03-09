	else if "`method'" == "ipw" {
			//IPS
set maxiter 300

 qui `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix `bips' = e(b)
 qui predict  double `ips'
		
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
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')


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
	
	
	qui gen double `invw1' = `zvar'
	
	qui sum `invw1' if `touse'==1
	matrix `w1s' = r(mean)
	
	
	qui ge double `omega1'= (`zvar')/`w1s'
	qui ge double `omega0'= (((1-`zvar')*`ips')/(1-`ips'))/`w1s'
	
	
	
	
	
	if `dmeanz1'==1 {
		qui ge double `d1hat'= 1
	}
	else {
		qui gen double `d1hat' = `omega1'*`tvar' if `touse'==1

	}
	
	if `dmeanz0'==0 {
		qui ge double `d0hat'= 0
	}
	else {
		
		qui gen double `d0hat' = `omega0'*`tvar' if `touse'==1
		
	}
	
	
	qui gen double `y1hat' = `omega1'*`yvar' if `touse'==1
	
	qui gen double `y0hat' = `omega0'*`yvar' if `touse'==1	

	
	qui sum `d1hat'  if `touse'==1 [iw = `samplew']
	matrix `denom1s' = r(mean)
	qui sum `d0hat'  if `touse'==1 [iw = `samplew']
	matrix `denom0s' = r(mean)
	qui sum `y1hat'  if `touse'==1 [iw = `samplew']
	matrix `num1s' = r(mean)
	qui sum `y0hat'  if `touse'==1 [iw = `samplew']
	matrix `num0s' = r(mean)
	matrix `nums' = `num1s'-`num0s'
	matrix `denoms' = `denom1s'-`denom0s'
	
	
	matrix `lates' =`nums'*invsym(`denoms')
	
	
	
	
	//Moment Functions
	local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))
	
	local eqips_inst instruments(eqips: `xvarsips')
	
	local eqw1 (eqw1: (`zvar' - {w1: _cons}))
    local eqw1_inst instruments(eqw1: )

	


	local eqd1 (eqd1: ((`zvar'*`tvar')/{w1:} - ({d1: _cons})))
	
	local eqd1_inst instruments(eqd1: )
	
		local eqd0 (eqd0: (((1-`zvar')*(exp({zhat:})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*`tvar')/{w1:}-({d0:  _cons}))		
		local eqd0_inst instruments(eqd0:  )
	

		if `dmeanz0' == 0 & `dmeanz1' !=1 {
			local denom (denom: (((`tvar' - ({denom: _cons})))))
		}
		else if `dmeanz0' != 0 & `dmeanz1' ==1 {
			local denom (denom: (({denom}-(1-({d0: })))))
		}
		else {
			local denom (denom: (({denom}-(({d1: })-({d0: })))))
		}
	
	

	
	local eqy1 (eqy1: (1.`zvar'*(`yvar'-({y1:  _cons})))) 
	
	local eqy1_inst instruments(eqy1:  )
	
		local eqy1 (eqy1: ((`zvar'*`yvar')/{w1:} - ({y1: _cons})))
	
	local eqy1_inst instruments(eqy1: )
	
		local eqy0 (eqy0: (((1-`yvar')*(exp({zhat:})/(1+exp({zhat:})))/((1/(1+exp({zhat:})))))*`tvar')/{w1:}-({y0:  _cons}))		
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
 gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `denom1s', `num0s', `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `denom0s', `num0s', `num1s', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd0' `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	 
if "`stat'" == "latt" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial =  (`bips', `w1s', `denom0s', `denom1s', `num0s', `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqw1' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `w1s', `denom1s', `num0s', `num1s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 1
 gmm    `eqips' `eqw1' `eqd1'  `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `w1s', `denom0s', `num0s', `num1s', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm  `eqips' `eqw1' `eqd0' `eqy0' `eqy1'  `num' `denom'    [pw = `samplew'],  `eqips_inst'  `eqd0_inst'  `eqy0_inst' `eqy1_inst'  onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}	
	
	}
