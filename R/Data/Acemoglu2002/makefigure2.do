***********************************************************
*Creates Figure 2: Population Density in 1500
***********************************************************
clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using makefigure2, replace

/*Data Files Used
	makefigure2
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
	
use makefigure2, clear
keep if ex2col==1

twoway (scatter logpgp95 lpd1500s, msymbol(i) mlabel(shortnam) mlabsize(tiny)) (lfit logpgp95 lpd1500s), yscale(r(6 10)) xscale(r(-5 5)) ylab(6(1)10, nogrid labs(small)) xlab(-5(1)5, labs(small)) legend(off)  ytitle("Log GDP per capita, PPP, 1995", size(med)) xtitle(Log Population Density in 1500, size(med)) saving(fig2.gph,replace)
graph export fig2.eps, replace



	
	