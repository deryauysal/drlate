		else if "`method'" == "naipw" {
			
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
	
	qui sum `invw1'
	scalar w1s = r(mean)
	qui sum `invw0'
	scalar w0s = r(mean)
	
	qui ge double `omega1'= (`zvar'/`ips')/w1s
	qui ge double `omega0'= ((1-`zvar')/(1-`ips'))/w0s
	
	qui ge double `term1'= `omega1'*( `yvar' -`y1hat')+`y1hat'
	qui sum `term1' 
	matrix `num1s' =  r(mean)
	qui ge double `term2'= `omega0'*( `yvar' -`y0hat')+`y0hat'
	qui sum `term2'
	matrix `num0s' =  r(mean)
	qui ge double `term3' = `omega1'*( `tvar' -`d1hat')+`d1hat'
	qui sum `term3' 
	matrix `denom1s' =  r(mean)	
	qui ge double `term4' = `omega0'*( `tvar' -`d0hat')+`d0hat'
	qui sum `term4'
	matrix `denom0s' =  r(mean)
	
	matrix `nums' = `num1s'-`num0s'
    matrix `denoms' = `denom1s'-`denom0s'
matrix `lates' =`nums'*invsym(`denoms')
		
//drop ips invw0 invw1 y1hat y0hat omega1 omega0 term* d1hat
		
		
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
	if "`stat'" == "late" {
if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `bd0', `bd1', `by0', `by1', w1s, w0s, `num1s', `num0s', `denom1s', `denom0s', `nums', `denoms')	
local k = colsof(initial )
matrix I = I(`k')  
set maxiter 5
 gmm  `eqips' `eqd0' `eqd1'  `eqy0' `eqy1' `eqw1' `eqw0' `num1'  `num0' `denom1' `denom0'  `num' `denom'  [pw = `samplew'],  `eqips_inst' `eqd0_inst' `eqd1_inst' `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')
		}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial =  (`bips', `bd1', `by0', `by1', w1s, w0s,  `num1s', `num0s', `nums', `denoms')
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm   `eqips'  `eqd1'  `eqy0' `eqy1' `eqw1' `eqw0' `num1'  `num0'   `num' `denom'  [pw = `samplew'],  `eqips_inst'  `eqd1_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}	
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `bd0', `by0', `by1',w1s, w0s, `num1s', `num0s', `nums', `denoms')	
matrix list initial
local k = colsof(initial )
matrix I = I(`k')  
qui set maxiter 5
 gmm    `eqips' `eqd0'  `eqy0' `eqy1' `eqw1' `eqw0' `num1'  `num0'  `num' `denom'   [pw = `samplew'],   `eqips_inst' `eqd0_inst'  `eqy0_inst' `eqy1_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce')	
}
}
if "`stat'" == "latt" {
	//drop ips
	error 198	
}
	}