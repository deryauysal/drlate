	else if "`method'" == "ipw" {
			//IPS
set maxiter 300

 qui `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
 qui matrix `bips' = e(b)
 qui predict  double `ips'
		
//Conditional mean of y1 & y0#
if "`stat'" == "late" {
 qui `omodel' `yvar'  if `zvar'==1 & `touse'==1 [pw = `samplew'/`ips']
qui matrix `by1' = e(b)
qui predict double `y1hat'

 
 
 qui `omodel' `yvar'  if `zvar'==0 & `touse'==1 [pw = `samplew'/(1-`ips')]
qui matrix `by0' = e(b)
qui predict double `y0hat'

//Conditional mean of d1 & d0#


	    if `dmeanz1'==1 {
        ge double `d1hat'= 1
	    }
		else {

qui `tmodel' `tvar'  if `zvar'==1  & `touse'==1  [pw =  `samplew'/`ips']
matrix `bd1' = e(b)
qui predict double `d1hat'
}
	    
	    if `dmeanz0'==0 {
        ge double `d0hat' = 0
	    }
		else {
qui `tmodel' `tvar'  if `zvar'==0  & `touse'==1 [pw =  `samplew'/(1-`ips')]
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
		qui `tmodel' `tvar'   if `zvar' == 0 & `touse'==1 [pw = `samplew'*(`ips'/(1-`ips'))]
		matrix `bd0' = e(b)
		qui predict double `d0hat'
	}
	
	qui regress `yvar' if `zvar' == 1 & `touse'==1 [pw = `samplew']
	matrix `by1' = e(b)
	qui predict double `y1hat'
	
	qui `omodel' `yvar'  if `zvar' == 0 & `touse'==1 [pw = `samplew'*(`ips'/(1-`ips'))]
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
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm   `eqips' `eqd0' `eqd1'  `eqy0' `eqy1'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `bd1', `by0', `by1', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1'  `num' `denom' [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
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
	 
if "`stat'" == "latt" {
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