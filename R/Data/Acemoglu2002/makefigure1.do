***********************************************************
*Creates Figure 1: Urbanization in 1500
***********************************************************
clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using makefigure1, replace

/*Data Files Used
	makefigure1
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
	
use makefigure1, clear
keep if baserf==1

twoway (scatter logpgp95 sjb1500, msymbol(i) mlabel(shortnam) mlabsize(tiny)) (lfit logpgp95 sjb1500), yscale(r(7 10)) xscale(r(0 20)) ylab(7(1)10, nogrid labs(small)) xlab(0(5)20, labs(small)) legend(off) ytitle("Log GDP per capita, PPP, 1995", size(med)) xtitle(Urbanization in 1500, size(med)) saving(fig1.gph,replace)
graph export fig1.eps, replace



	
	