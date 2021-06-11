***********************************************************
*Creates Table 6: Robustness Checks for Urbanization and Log population Density
***********************************************************

clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using maketable6, replace

/*Data Files Used
	maketable6
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
use maketable6, clear
keep if baserf==1

***----------Panel A: Unweighted regressions

*Column 1: Base sample

reg logpgp95 sjb1500 

*Column 2: Assuming lower urbanization in the Americas

reg logpgp95 sj1b1500 

*Column 3: Assuming lower urbanization in North Africa
 
reg logpgp95 sj2b1500 

*Column 4: Assuming lower urbanization in Indian subcontinent

reg logpgp95 sj3b1500 

*Column 5: Using least favorable combination of assumptions

reg logpgp95 sj4b1500 

*Column 6: Using augmented Toynbee definition of civilization


*Column 7: Using land area in 1995 for population density

use maketable6, clear
keep if ex2col==1

reg logpgp95 lpd1500 

*Column 8: Alternative assumptions for log population density

reg logpgp95 lpd1500j

*Column 9: All countries never colonized by Europe

use maketable6, clear
regress logpgp95 sjb1500 if nevcol==1 

*Column 10: Europe (including Eastern Europe)

regress logpgp95 sjb1500 if europe==1 



****-----Panel B: Regressions weighted using log population in 1500
	
	
use maketable6, clear
keep if baserf==1

*Column 1: Base sample

reg logpgp95 sjb1500 [aw=lpa1500]

*Column 2: Assuming lower urbanization in the Americas

reg logpgp95 sj1b1500  [aw=lpa1500]

*Column 3: Assuming lower urbanization in North Africa
 
reg logpgp95 sj2b1500 [aw=lpa1500] 

*Column 4: Assuming lower urbanization in Indian subcontinent

reg logpgp95 sj3b1500 [aw=lpa1500] 

*Column 5: Using least favorable combination of assumptions

reg logpgp95 sj4b1500 [aw=lpa1500] 

*Column 6: Using augmented Toynbee definition of civilization


*Column 7: Using land area in 1995 for population density

use maketable6, clear
keep if ex2col==1

reg logpgp95 lpd1500 [aw=lpa1500] 

*Column 8: Alternative assumptions for log population density

reg logpgp95 lpd1500j [aw=lpa1500]

*Column 9: All countries never colonized by Europe

use maketable6, clear
regress logpgp95 sjb1500 if nevcol==1  [aw=lpa1500]

*Column 10: Europe (including Eastern Europe)

regress logpgp95 sjb1500 if europe==1 [aw=lpa1500] 
	
	
****--------------Panel C: Including both urbanization and log population density as independent variables


use maketable6, clear
keep if baserf==1

*Column 1: Base sample

reg logpgp95 sjb1500 lpd1500s

*Column 2: Assuming lower urbanization in the Americas

reg logpgp95 sj1b1500 lpd1500s

*Column 3: Assuming lower urbanization in North Africa
 
reg logpgp95 sj2b1500 lpd1500s

*Column 4: Assuming lower urbanization in Indian subcontinent

reg logpgp95 sj3b1500 lpd1500s

*Column 5: Using least favorable combination of assumptions

reg logpgp95 sj4b1500 lpd1500s

*Column 6: Using augmented Toynbee definition of civilization


*Column 7: Using land area in 1995 for population density

use maketable6, clear
keep if ex2col==1

reg logpgp95  sjb1500 lpd1500 

*Column 8: Alternative assumptions for log population density

reg logpgp95  sjb1500 lpd1500j

*Column 9: All countries never colonized by Europe

use maketable6, clear
regress logpgp95 sjb1500 lpd1500s if nevcol==1 

*Column 10: Europe (including Eastern Europe)

regress logpgp95 sjb1500 lpd1500s if europe==1 


***----Panel D: Instrumenting for urbanization in 1500 using log population density in 1500

use maketable6, clear
keep if baserf==1

*Column 1: Base sample

ivreg logpgp95 (sjb1500 =lpd1500s)

*Column 2: Assuming lower urbanization in the Americas

ivreg logpgp95 (sj1b1500 =lpd1500s)

*Column 3: Assuming lower urbanization in North Africa
 
ivreg logpgp95 (sj2b1500 =lpd1500s)

*Column 4: Assuming lower urbanization in Indian subcontinent

ivreg logpgp95 (sj3b1500 =lpd1500s)

*Column 5: Using least favorable combination of assumptions

ivreg logpgp95 (sj4b1500 =lpd1500s)

*Column 6: Using augmented Toynbee definition of civilization


*Column 7: Using land area in 1995 for population density

use maketable6, clear
keep if ex2col==1

ivreg logpgp95 (sjb1500=lpd1500)

*Column 8: Alternative assumptions for log population density

ivreg logpgp95 (sjb1500=lpd1500j)

*Column 9: All countries never colonized by Europe

use maketable6, clear
ivreg logpgp95 (sjb1500=lpd1500s) if nevcol==1 

*Column 10: Europe (including Eastern Europe)

ivreg logpgp95 (sjb1500=lpd1500s) if europe==1 


	