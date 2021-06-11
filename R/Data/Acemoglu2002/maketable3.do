***********************************************************
*Creates Table 3: Urbanization in 1500 and GDP per Capita
*in 1995 for Former European Colonies
***********************************************************

clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using maketable3, replace

/*Data Files Used
	maketable3
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
use maketable3, clear
keep if baserf==1

***---Column 1: Base sample
regress logpgp95 sjb1500

***---Column 2: W/o North Africa
regress logpgp95 sjb1500 if nafrica!=1

***---Column 3: W/o the Americas
regress logpgp95 sjb1500 if america!=1

***---Column 4: Just the Americas
regress logpgp95 sjb1500 if america==1

***---Column 5: With continent dummies
regress logpgp95 sjb1500 america africa asia

***---Column 6: W/o neo-Europes

replace neoeur=0 if shortnam=="ARG"
regress logpgp95 sjb1500 if neoeur!=1

***---Column 7: Controlling for latitude
regress logpgp95 sjb1500 lat_abst

***---Column 8: Controlling for climate
regress logpgp95 sjb1500 humid* temp* steplow  deslow stepmid desmid 

testparm temp1-temp5
testparm humid1-humid4
testparm steplow-desmid

***---Column 9: Controlling for resources
regress logpgp95 sjb1500 coal landlock island goldm iron silv zinc oilres 
test (goldm=0) (iron=0) (silv=0) (zinc=0) (oilres=0)

***---Column 10: Controlling for colonial origin
regress logpgp95 sjb1500 f_french f_spain f_pothco f_dutch f_belg f_italy f_germ

***---Column 11: Controlling for religion
regress logpgp95 sjb1500 catho80 muslim80 notmcp80
test (catho80=0) (muslim80=0) (notmcp80=0)

	