***********************************************************
*Creates Figure 3: Log GDP per Capita (PPP) in 1995 Against the Urbanization Rate in 1995
***********************************************************
clear
capture log close
cd "INSERT DIRECTORY PATH HERE"
log using makefigure3, replace

/*Data Files Used
	makefigure3
	
*Data Files Created as Final Product
	none
	
*Data Files Created as Intermediate Product
	none*/
	
	
use makefigure3, clear
keep if ex2col==1

twoway (scatter logpgp95 urbz1995, msymbol(i) mlabel(shortnam) mlabsize(tiny)) (lfit logpgp95 urbz1995), yscale(r(6 10)) xscale(r(0 100)) ylab(6(1)10, nogrid labs(small)) xlab(0(10)100, labs(small)) legend(off)  ytitle("Log GDP per capita, PPP, 1995", size(med)) xtitle(Urbanization in 1995, size(med)) saving(fig3.gph,replace)
graph export fig3.eps, replace



	
	