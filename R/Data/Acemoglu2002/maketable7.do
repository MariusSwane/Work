***********************************************************
*Creates Table 7: Urbanization, Population Density, and Institutions
***********************************************************

clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using maketable7, replace

/*Data Files Used
	maketable7
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
use maketable7, clear
keep if ex2col==1

***---Column 1: Average protection against expropriation risk, 1985-1995

reg avexpr sjb1500 

***---Column 2: Average protection against expropriation risk, 1985-1995

reg avexpr lpd1500s

***---Column 3: Average protection against expropriation risk, 1985-1995

reg avexpr sjb1500 lpd1500s


***---Column 4: Constraint on executive in 1990

reg cons90 sjb1500 

***---Column 5: Constraint on executive in 1990

reg cons90 lpd1500s

***---Column 6: Constraint on executive in 1990

reg cons90 sjb1500 lpd1500s


***---Column 7: Constraint on executive in first year of independence

reg cons1 sjb1500 indtime 

***---Column 8:Constraint on executive in first year of independence

reg cons1 lpd1500s indtime 

***---Column 9: Constraint on executive in first year of independence

reg cons1 sjb1500 lpd1500s indtime 


****--------Panel B: Controlling for latitude


***---Column 1: Average protection against expropriation risk, 1985-1995

reg avexpr sjb1500 lat_abst

***---Column 2: Average protection against expropriation risk, 1985-1995

reg avexpr lpd1500s lat_abst

***---Column 3: Average protection against expropriation risk, 1985-1995

reg avexpr sjb1500 lpd1500s lat_abst

***---Column 4: Constraint on executive in 1990

reg cons90 sjb1500 lat_abst

***---Column 5: Constraint on executive in 1990

reg cons90 lpd1500s lat_abst

***---Column 6: Constraint on executive in 1990

reg cons90 sjb1500 lpd1500s lat_abst

***---Column 7: Constraint on executive in first year of independence

reg cons1 sjb1500 indtime lat_abst

***---Column 8:Constraint on executive in first year of independence

reg cons1 lpd1500s indtime lat_abst

***---Column 9: Constraint on executive in first year of independence

reg cons1 sjb1500 lpd1500s indtime lat_abst

