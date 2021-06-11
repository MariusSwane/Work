***********************************************************
*Creates Table 5: Population Density and GDP pc in Former
*European Colonies
***********************************************************

clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using maketable5, replace

/*Data Files Used
	maketable5
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
************************************************************************
***---Panel A: Log population density in 1500 as independent variable	
************************************************************************
	
use maketable5, clear
keep if ex2col==1

*---Column 1: Base Sample	
reg logpgp95 lpd1500s 

*---Column 2: No Africa
reg logpgp95 lpd1500s if africa!=1
	
*---Column 3: No Americas
reg logpgp95 lpd1500s if america!=1
	
*---Column 4: Just Americas
reg logpgp95 lpd1500s if america==1

*---Column 5: With continent dummies
reg logpgp95 lpd1500s america africa asia

*---Column 6: Without neo-Europes
replace neoeur=0 if shortnam=="ARG"
reg logpgp95 lpd1500s if neoeur!=1

*---Column 7: Controlling for latitude 
reg logpgp95 lpd1500s  lat_abst

*---Column 8: Controlling for climate 
reg logpgp95 lpd1500s  humid* temp* steplow  deslow stepmid desmid  drystep hiland drywint
testparm temp1-temp5
testparm humid1-humid4
testparm steplow-desmid
test (steplow=0) (deslow=0) (stepmid=0) (desmid=0) (drystep=0) (hiland=0) (drywint=0)

***---Column 9: Controlling for resources
regress logpgp95 lpd1500s coal landlock island goldm iron silv zinc oilres 
test (goldm=0) (iron=0) (silv=0) (zinc=0) (oilres=0)

***---Column 10: Controlling for colonial origin
regress logpgp95 lpd1500s f_french f_spain f_pothco f_dutch f_belg f_italy f_germ

***---Column 11: Controlling for religion
regress logpgp95 lpd1500s catho80 muslim80 notmcp80
test (catho80=0) (muslim80=0) (notmcp80=0)



***--------------Panel B: Log population and log land in 1500 as separate independent variables	
	
use maketable5, clear
keep if ex2col==1

*NOTES

*lpd1500s=ln(pop1500*100/(land95*adjland)) AND lland15=ln(land95*adjland)
*the lpop1500 or pop1500 vars in the dataset don't have enough obs to generate the # of observations in the table.

g lnpop1500=lpd1500s+lland15

*---Column 1: Base Sample	
reg logpgp95 lland15 lnpop1500 


*---Column 2: No Africa
reg logpgp95 lland15 lnpop1500 if africa!=1
	
*---Column 3: No Americas
reg logpgp95 lland15 lnpop1500 if america!=1
	
*---Column 4: Just Americas
reg logpgp95 lland15 lnpop1500 if america==1

*---Column 5: With continent dummies
reg logpgp95 lland15 lnpop1500 america africa asia

*---Column 6: Without neo-Europes
replace neoeur=0 if shortnam=="ARG"
reg logpgp95 lland15 lnpop1500 if neoeur!=1

*---Column 7: Controlling for latitude 
reg logpgp95 lland15 lnpop1500  lat_abst

*---Column 8: Controlling for climate 
reg logpgp95 lland15 lnpop1500  humid* temp* steplow  deslow stepmid desmid  drystep hiland drywint
testparm temp1-temp5
testparm humid1-humid4
testparm steplow-desmid
test (steplow=0) (deslow=0) (stepmid=0) (desmid=0) (drystep=0) (hiland=0) (drywint=0)

***---Column 9: Controlling for resources
regress logpgp95 lland15 lnpop1500 coal landlock island goldm iron silv zinc oilres 
test (goldm=0) (iron=0) (silv=0) (zinc=0) (oilres=0)

***---Column 10: Controlling for colonial origin
regress logpgp95 lland15 lnpop1500 f_french f_spain f_pothco f_dutch f_belg f_italy f_germ

***---Column 11: Controlling for religion
regress logpgp95 lland15 lnpop1500 catho80 muslim80 notmcp80 
test (catho80=0) (muslim80=0) (notmcp80=0)
	
	
***----Panel C: Using population density in 1000 A.D. as an instrument for population density in 1500 A.D.


use maketable5, clear
keep if ex2col==1

*---Column 1: Base Sample	
ivreg logpgp95 (lpd1500s=lpd1000s) 

*---Column 2: No Africa
ivreg logpgp95 (lpd1500s=lpd1000s) if africa!=1
	
*---Column 3: No Americas
ivreg logpgp95 (lpd1500s=lpd1000s) if america!=1
	
*---Column 4: Just Americas
ivreg logpgp95 (lpd1500s=lpd1000s) if america==1

*---Column 5: With continent dummies
ivreg logpgp95 america africa asia (lpd1500s=lpd1000s)

*---Column 6: Without neo-Europes
replace neoeur=0 if shortnam=="ARG"
ivreg logpgp95 (lpd1500s=lpd1000s) if neoeur!=1

*---Column 7: Controlling for latitude 
ivreg logpgp95  lat_abst (lpd1500s=lpd1000s)

*---Column 8: Controlling for climate 
ivreg logpgp95  humid* temp* steplow  deslow stepmid desmid  drystep hiland drywint (lpd1500s=lpd1000s)
testparm temp1-temp5
testparm humid1-humid4
testparm steplow-desmid
test (steplow=0) (deslow=0) (stepmid=0) (desmid=0) (drystep=0) (drywint=0)

***---Column 9: Controlling for resources
ivreg logpgp95 coal landlock island goldm iron silv zinc oilres  (lpd1500s=lpd1000s)
test (goldm=0) (iron=0) (silv=0) (zinc=0) (oilres=0) 

***---Column 10: Controlling for colonial origin
ivreg logpgp95 f_french f_spain f_pothco f_dutch f_belg f_italy f_germ (lpd1500s=lpd1000s)

***---Column 11: Controlling for religion
ivreg logpgp95 catho80 muslim80 notmcp80 (lpd1500s=lpd1000s) 
test (catho80=0) (muslim80=0) (notmcp80=0) 
