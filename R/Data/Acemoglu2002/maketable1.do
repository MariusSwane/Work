***********************************************************
*Creates Table 1: Summary Statistics
***********************************************************
clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using maketable1, replace

/*Data Files Used
	maketable1
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
***----Column 1

use maketable1, clear

summ logpgp95

summ urbz1995 

summ sjb1500

summ lpd1500s

summ pd1500s

summ lpd1000s

summ avexpr

summ cons90

summ cons1

summ  euro1900

***-----Column 2
	
use maketable1, clear

keep if baserf==1


summ logpgp95

summ urbz1995 

summ sjb1500

summ lpd1500s

summ pd1500s

summ lpd1000s

summ avexpr

summ cons90

summ cons1

summ  euro1900


***-----Column 3
	
use maketable1, clear

keep if (lpd1500s !=. & ex2col==1)


summ logpgp95

summ urbz1995 

summ sjb1500

summ lpd1500s

summ pd1500s

summ lpd1000s

summ avexpr

summ cons90

summ cons1

summ  euro1900


***----Column 4

use maketable1, clear

keep if baserf==1

egen med_urb=median(sjb1500)
g bmed=1 if sjb1500<=med_urb
g amed=1 if sjb1500>med_urb

summ bmed amed

keep if bmed==1

summ logpgp95

summ urbz1995 

summ sjb1500

summ lpd1500s

summ pd1500s

summ lpd1000s

summ avexpr

summ cons90

summ cons1

summ  euro1900




***----Column 5

use maketable1, clear

keep if baserf==1

egen med_urb=median(sjb1500)
g bmed=1 if sjb1500<=med_urb
g amed=1 if sjb1500>med_urb

summ bmed amed

keep if amed==1

summ logpgp95

summ urbz1995 

summ sjb1500

summ lpd1500s

summ pd1500s

summ lpd1000s

summ avexpr

summ cons90

summ cons1

summ  euro1900



***-----Column 6
	
use maketable1, clear

keep if (lpd1500s !=. & ex2col==1)


egen rank_dens=rank(lpd1500s), track
egen count_dens=count(lpd1500s)
gen p_dens=rank_dens/count_dens
g bmed=1 if p_dens<=.5
g amed=1 if p_dens>.5


keep if bmed==1

summ logpgp95

summ urbz1995 

summ sjb1500

summ lpd1500s

summ pd1500s

summ lpd1000s

summ avexpr

summ cons90

summ cons1

summ  euro1900

	


***-----Column 7
	
use maketable1, clear

keep if (lpd1500s !=. & ex2col==1)


egen rank_dens=rank(lpd1500s), track
egen count_dens=count(lpd1500s)
gen p_dens=rank_dens/count_dens
g bmed=1 if p_dens<=.5
g amed=1 if p_dens>.5


keep if amed==1

summ logpgp95

summ urbz1995 

summ sjb1500

summ lpd1500s

summ pd1500s

summ lpd1000s

summ avexpr

summ cons90

summ cons1

summ  euro1900

	