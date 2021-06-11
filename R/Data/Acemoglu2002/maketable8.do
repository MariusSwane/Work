***********************************************************
*Creates Table 8: GDP Per Capita and Institutions
***********************************************************

clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using maketable8, replace

/*Data Files Used
	maketable8
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
use maketable8, clear
keep if ex2col==1

***Panel A: Second-stage regressions and Panel B: First-stage regressions

*Column 1: Average protection against expropriation risk, 1985-1995

ivreg logpgp95 sjb1500 (avexpr=logem4), first

*Column 2: Average protection against expropriation risk, 1985-1995

ivreg logpgp95 lpd1500s (avexpr=logem4), first

*Column 3: Constraint on executive in 1990

ivreg logpgp95 sjb1500 (cons90=logem4), first

*Column 4: Constraint on executive in 1990

ivreg logpgp95 lpd1500s (cons90=logem4), first

*Column 5: Constraint on executive in first year of independence

ivreg logpgp95 sjb1500 indtime (cons1=logem4), first

*Column 6: Constraint on executive in first year of independence

ivreg logpgp95 lpd1500s indtime (cons1=logem4), first




***Panel C: Coefficient on institutions without urbanization or population density in 1500


*Column 1: Average protection against expropriation risk, 1985-1995

ivreg logpgp95 (avexpr=logem4) if  sjb1500!=.

*Column 2: Average protection against expropriation risk, 1985-1995

ivreg logpgp95 (avexpr=logem4) if  lpd1500s!=.

*Column 3: Constraint on executive in 1990

ivreg logpgp95 (cons90=logem4) if  sjb1500!=.

*Column 4: Constraint on executive in 1990

ivreg logpgp95 (cons90=logem4) if  lpd1500s!=.

*Column 5: Constraint on executive in first year of independence

ivreg logpgp95 indtime (cons1=logem4) if  sjb1500!=.

*Column 6: Constraint on executive in first year of independence

ivreg logpgp95 indtime (cons1=logem4) if  lpd1500s!=.
