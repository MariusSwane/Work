**************************************************************************
*Creates Table 9: The Interaction of UK Industrialization and Institutions
**************************************************************************

clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using maketable9, replace

/*Data Files Used
	maketable9
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
***---Panel A: Dependent variable is industrial production per capita

use maketable9, clear

*Column 1: Former colonies, using only pre-1950 data

*Column 2: Former colonies, using data through 1980 (all data)

*Column 3: Former colonies, using only pre-1950 data

*Column 4: Former colonies, using only data pre-1950 and for independent countries

*Column 5: Former colonies, with average institutions for each country, using only pre-1950 data

*Column 6: Former colonies, with average institutions for each country, using only pre-1590 data

*Column 7: Former colonies, with average institutions for each country, instrumenting using settler mortality, only pre-1950 data

*Column 8: Former colonies, with average institutions for each country, instrumenting using settler mortality, only pre-1950

*Column 9: Former colonies, with average institutions for each country, instrumenting using settler mortality, only pre-1950 data

*Column 10: Former colonies, with average institutions for each country, instrumenting using settler mortality, only pre-1950 data


***---Panel A: Dependent variable is industrial production per capita

***---Panel B: Dependent variable is log GDP per capita	
