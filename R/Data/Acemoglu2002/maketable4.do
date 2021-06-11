***********************************************************
*Creates Table 4: Alternative Measures of Urbanization
***********************************************************

clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using maketable4, replace

/*Data Files Used
	reversal_fortune
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
***----Panel A: Using AJR base measure of urbanization
	
use maketable4, clear
keep if baserf==1


***---Column 1: Base sample
regress logpgp95 sjb1500

***---Column 2: With continent dummies
regress logpgp95 sjb1500 america africa asia

***---Column 3: W/o neo-Europes

replace neoeur=0 if shortnam=="ARG"
regress logpgp95 sjb1500 if neoeur!=1

***---Column 4: Congrolling for latitude
regress logpgp95 sjb1500 lat_abst

***---Column 5: Controlling for resources
regress logpgp95 sjb1500 coal landlock island goldm iron silv zinc oilres 


***----Panel B: Using only Bairoch's estimates
	
use maketable4, clear
keep if baserf==1

***---Column 1
regress logpgp95 bonly

***---Column 2
regress logpgp95 bonly america africa asia

***---Column 3

replace neoeur=0 if shortnam=="ARG"
regress logpgp95 bonly if neoeur!=1

***---Column 4
regress logpgp95 bonly lat_abst

***---Column 5
regress logpgp95 bonly coal landlock island goldm iron silv zinc oilres 


***----Panel C: Using only Eggimann’s estimates
	
use maketable4, clear
keep if baserf==1

***---Column 1
regress logpgp95 eonly

***---Column 2
regress logpgp95 eonly america africa asia

***---Column 3

replace neoeur=0 if shortnam=="ARG"
regress logpgp95 eonly if neoeur!=1

***---Column 4
regress logpgp95 eonly lat_abst

***---Column 5
regress logpgp95 eonly coal landlock island goldm iron silv zinc oilres 




***----Panel D: Using only Chandler’s estimates
	
use maketable4, clear
keep if baserf==1

***---Column 1
regress logpgp95 cu1500

***---Column 2
regress logpgp95 cu1500 america africa asia

***---Column 3

replace neoeur=0 if shortnam=="ARG"
regress logpgp95 cu1500 if neoeur!=1

***---Column 4
regress logpgp95 cu1500 lat_abst

***---Column 5
regress logpgp95 cu1500 coal landlock island goldm iron silv zinc oilres 

	
	
***-------Panel E: Using Davis-Zipf Adjustment for Eggimann’s series

use maketable4, clear
keep if baserf==1

***---Column 1
regress logpgp95 eonlyal2

***---Column 2
regress logpgp95 eonlyal2 america africa asia

***---Column 3

replace neoeur=0 if shortnam=="ARG"
regress logpgp95 eonlyal2 if neoeur!=1

***---Column 4
regress logpgp95 eonlyal2 lat_abst

***---Column 5
regress logpgp95 eonlyal2 coal landlock island goldm iron silv zinc oilres 

	
	