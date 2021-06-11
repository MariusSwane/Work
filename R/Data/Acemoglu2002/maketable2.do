***********************************************************
*Creates Table 2: Urbanization and Per Capita Income
***********************************************************
clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using maketable2, replace

/*Data Files Used
	maketable2
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
***----Column 1

use maketable2, clear

reg lgdppc13 uper19a

***----Column 2

reg lgdppc50 urbz1960

***----Column 3

reg logpgp95  urbz1995

***----Column 4

reg logpgp95  urbz1995 if ex2col==1

***----Column 5

reg logpgp95  urbz1995 if (f_nocolo==1 & ex2col!=1)

***----Column 6

reg logpgp95  urbz1995 africa asia america
