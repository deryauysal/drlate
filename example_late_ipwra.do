
cd "C:\Users\derya.uysal\Dropbox\Derya\research\LATE\stata\gmm stata\late_ipwra_v20220808"


do "late.ado"

cls 
clear matrix
clear mata
set maxvar 120000

//use "Abadie2003_data.dta", clear
use "/Users/Derya/Dropbox/Derya/research/LATE/DR LATE/05-simulations/DR LATE (Abadie)/Abadie2003_data.dta"

egen incst= std(inc)
egen age25st= std(age25)
egen age25sqst= std(age25sq)

cls
/*LATE estimation and teffects estimation of the numerator */
late ipwra net_tfa p401 e401 incst age25st age25sqst marr fsize
late ipwra net_tfa p401 e401 incst age25st age25sqst marr fsize, zmodel(cbps)
//teffects ipwra (net_tfa incst age25st age25sqst marr fsize) (e401 incst age25st age25sqst marr fsize) , iterate(2)

late ipw net_tfa p401 e401 incst age25st age25sqst marr fsize
late ipw net_tfa p401 e401 incst age25st age25sqst marr fsize, zmodel(cbps)
late nipw net_tfa p401 e401 incst age25st age25sqst marr fsize
late nipw net_tfa p401 e401 incst age25st age25sqst marr fsize, zmodel(cbps)

//teffects ipw (net_tfa) (e401 incst age25st age25sqst marr fsize) , iterate(2)

late ra net_tfa p401 e401 incst age25st age25sqst marr fsize
late ra net_tfa p401 e401 incst age25st age25sqst marr fsize, zmodel(cbps)
//teffects ra (net_tfa incst age25st age25sqst marr fsize) (e401) , iterate(2)

late aipw net_tfa p401 e401 incst age25st age25sqst marr fsize
late aipw net_tfa p401 e401 incst age25st age25sqst marr fsize, zmodel(cbps)

//teffects aipw (net_tfa incst age25st age25sqst marr fsize) (e401 incst age25st age25sqst marr fsize) , iterate(2)

late naipw net_tfa p401 e401 incst age25st age25sqst marr fsize
late naipw net_tfa p401 e401 incst age25st age25sqst marr fsize, zmodel(cbps)


/*LATT estimation and teffects estimation of the numerator */
cls
late ipwra net_tfa p401 e401 incst age25st age25sqst marr fsize, stat(latt)

qui teffects ipwra (net_tfa incst age25st age25sqst marr fsize) (e401 incst age25st age25sqst marr fsize) , iterate(2) atet
display _b[ATET:r1vs0.e401]
display _se[ATET:r1vs0.e401]

late ipw net_tfa p401 e401 incst age25st age25sqst marr fsize, stat(latt)

late nipw net_tfa p401 e401 incst age25st age25sqst marr fsize, stat(latt)

qui teffects ipw (net_tfa) (e401 incst age25st age25sqst marr fsize) , iterate(2) atet
display _b[ATET:r1vs0.e401]
display _se[ATET:r1vs0.e401]

late ra net_tfa p401 e401 incst age25st age25sqst marr fsize, stat(latt)

qui teffects ra (net_tfa incst age25st age25sqst marr fsize) (e401) , iterate(2) atet
display _b[ATET:r1vs0.e401]
display _se[ATET:r1vs0.e401]

late aipw net_tfa p401 e401 incst age25st age25sqst marr fsize, stat(latt)
late naipw net_tfa p401 e401 incst age25st age25sqst marr fsize, stat(latt)


use "Finkelstein2012_er_num_mod_12m.dta", clear

late ipwra er_num_mod_12m ohp_all_ever_survey treatment ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2,  omodel(poisson) vce(cluster household_id) samplew(weight_12m)

teffects ipwra (er_num_mod_12m ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2, poisson) (treatment ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2)  [pweight=weight_12m] ,  vce(cluster household_id) iterate(2)

teffects ipwra (ohp_all_ever_survey ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2, poisson) (treatment ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2) [pweight=weight_12m] ,  vce(cluster household_id)  iterate(2)


use "Finkelstein2012_er_any_12m.dta", clear

late ipwra er_any_12m ohp_all_ever_survey treatment ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2,  omodel(logit) vce(cluster household_id) samplew(weight_12m)

teffects ipwra (er_any_12m ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2, logit) (treatment ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2) [pweight=weight_12m] ,  vce(cluster household_id) 

teffects ipwra (ohp_all_ever_survey ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2, logit) (treatment ddddraw_sur_2 ddddraw_sur_3 ddddraw_sur_4 ddddraw_sur_5 ddddraw_sur_6 ddddraw_sur_7 dddnumhh_li_2 dddnumhh_li_3 ddddraXnum_2_2 ddddraXnum_2_3 ddddraXnum_3_2 ddddraXnum_3_3 ddddraXnum_4_2 ddddraXnum_5_2 ddddraXnum_6_2 ddddraXnum_7_2) [pweight=weight_12m],  vce(cluster household_id)
