***********************************************************
*Creates Figure 4a: Urbanization Rate in India, the United States, and New World Countries with Low and High Urbanization, 800–1920
***********************************************************
clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using makefigure4, replace

/*Data Files Used
	makefigure4
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
use makefigure4, clear
keep if  shortnam=="USA"
reshape long cu, i(shortnam) j(year)
label var cu "urbanization (Chandler/Mitchell)"
label var year "year of urb. obs."


	

	
	