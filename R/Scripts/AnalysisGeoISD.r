#==============================================================================#
#	 █████╗ ███╗   ██╗ █████╗ ██╗  ██╗   ██╗███████╗██╗███████╗	       #
#	██╔══██╗████╗  ██║██╔══██╗██║  ╚██╗ ██╔╝██╔════╝██║██╔════╝	       #
#	███████║██╔██╗ ██║███████║██║   ╚████╔╝ ███████╗██║███████╗	       #
#	██╔══██║██║╚██╗██║██╔══██║██║    ╚██╔╝  ╚════██║██║╚════██║	       #
#	██║  ██║██║ ╚████║██║  ██║███████╗██║   ███████║██║███████║	       #
#	╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝   ╚══════╝╚═╝╚══════╝	       #
#==============================================================================#

# {{{ Header
#==============================================================================#
#	Loading Packages						       #
#==============================================================================#	

library(corrplot)
library(conflicted)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(Hmisc)
library(margins)
library(MASS)
library(pscl)
library(purrr)
library(rayshader)
library(RColorBrewer)
library(readr)
library(raster)
library(spdep)
library(sf)
library(sidedata)
library(summarytools)
library(texreg)
library(xtable)
library(notifier)

#==============================================================================#
#	Loading Data and functions					       #
#==============================================================================#	

load("../Data/GeoISDControls.Rdata")
source("goldenScatterCAtheme.r")

#==============================================================================#
#	Resolving conflicts 						       #
#==============================================================================#		

conflict_prefer("filter", "dplyr")  
conflict_prefer("select", "dplyr")
# }}}

# {{{ Pre-analysis
#==============================================================================#
#	Creating New Variables						       #
#==============================================================================#

prio_grid_isd <- prio_grid_isd %>% mutate(
				sqrtDeaths = sqrt(deaths),
				sqrtSBDeaths = sqrt(statebaseddeaths),
				sqrtState_based = sqrt(state_based),
				sqrtSpAll = sqrt(sp_os_i_sum),
				sqrtSpAny = sqrt(sp_i_sum_any),
				sqrtSpNoInt = sqrt(sp_os_sum),
				sqrtCapdist = sqrt(capdist),
				sqrtBDist = sqrt(bdist3),
				sqrtPopd = sqrt(popd),
				logDeaths = log(deaths + 1),
				logState_based = log(state_based + 1),
				logSpAll = log(sp_os_i_sum + 1),
				logSpAny = log(sp_i_sum_any + 1),
				logCapdist = log(capdist),
				logBDist = log(bdist3 + 1),
				logCDist = log(distcoast +1),
				logPopd = log(popd + 1),
				logOrg3 = log(org3 +1),
				sqrtorg3 = sqrt(org3),
				sqrtNonstate = sqrt(non_state),
				dumState = as.numeric(state_based > 0),    
				dumOrg3 = as.numeric(org3 > 0),
				SpAllXCapDist = logSpAll*capdist,
				SpAnyXCapDist = logSpAny*capdist,
				SpAll10 = sp_os_i_sum*10,
				region1 = as.numeric(factor(region)==1),
				region2 = as.numeric(factor(region)==2),
				region3 = as.numeric(factor(region)==3),
				region4 = as.numeric(factor(region)==4))

#==============================================================================#
#	Functions				    	             	       #
#==============================================================================#

# Function to carry out analysis
analysis <- function(dvs, ivs, controls, data, test_label){
  models_out_j <- NULL
  models_out <- NULL
  for(i in 1:length(dvs)) {
    dv <- dvs[i] 
    for(j in 1:length(ivs)) {
      m1 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls)),
               data=data)
      m2 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls,
				    extended_controls)), data=data)
      m3 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls,
				    extended_controls, even_more)), data=data)
      m4 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls,
				    extended_controls, even_more, rest)),
		   data=data)
      models_out_j[[j]] <- list(m1, m2, m3, m4)
    } 
    models_out[[i]] <- models_out_j %>% purrr::flatten()
  }
  models_out <- models_out
  final <- list('models'=models_out)
  return(final)
}

# Function to create pastel
lighten <- function (col, pct = 0.75, alpha = .8) 
{
    if (abs(pct) > 1) {
        print("Warning:  Error in Lighten; invalid pct")
        pcol <- col2rgb(col)/255
    }
    else {
        col <- col2rgb(col)/255
        if (pct > 0) {
            pcol <- col + pct * (1 - col)
        }
        else {
            pcol <- col * pct
        }
    }
    pcol <- rgb(pcol[1], pcol[2], pcol[3], alpha)
    pcol
}

cols <- c("red", "green", "blue")

pastels <- NULL
for (i in 1:length(cols)) {
	pastels[i] <- lighten(cols[i])
}

#==============================================================================#
#	Variables							       #
#==============================================================================#

controls <- c('+ mountains_mean + water_gc + barren_gc + logCDist')

extended_controls <- c('+ region3')

even_more <- c('+ logPopd')

rest <- c('+ logBDist')

climate_controls <- c('+ temp_sd + temp + prec_sd + prec_gpcc')

coefs <- list('sqrtSpAll' = 'Precolonial state presence (sqrt)', 
		'mountains_mean' = 'Mountainous terrain',
	     	'water_gc' = 'Water (%)', 
	 	'barren_gc' = 'Barren (%)', 
	      	'logCDist' = 'Distance to coast (log)',
	      	'logPopd' = 'Population density (log)', 
		'logBDist' = 'Distance to international boundary (log)', 
		'temp_sd' = 'Temperature (SD)',
	      	'temp' = 'Temperature (mean)', 
		'prec_sd' = 'Precipitation (SD)', 
		'prec_gpcc' = 'Precipitation (mean)', 
		'SpAll10' = 'Precolonial state presence (10)',
		'fra' = 'Former French colony',
		'gbr' = 'Former British colony',
		'region3' = 'North Africa',
		'sqrtSpAll' = 'Precolonial state presence (sqrt)')

# countcoefs <- list()
# 
# for( i in 1:length(coefs)) {
# 	names(countcoefs) <- names(coefs) 
# 	countcoefs[i] <- paste('Count model: ',coefs[i])
# }
# 
# logitcoefs <- list()
# 
# for( i in 1:length(coefs)) {
# 	logitcoefs[i] <- paste('Zero model: ',coefs[i])
# }
# 
# zinbcoefs <- c(logitcoefs, countcoefs)

control_names_int <-c('Baseline', 'Extended Controls')

control_names <-c('Geography', 'North Africa', 'Population densisty', 'Distance
		  to border')

dvs <- c('statebaseddeaths', 'state_based')

captions <- c('Fatalities', 'State based conflict events')

captions_int <- c('Fatalities * Distance to capital', 'Conflict events *
		  Distance to capital')

ivs <- c('sqrtSpAll')

interactions <- c('sqrtSpAll * logCapdist')
# }}}

# {{{ Descriptive Statistics
#==============================================================================#
#	Descriptive Statistics                                         	       #
#==============================================================================#

#==============================================================================#
# Correlation Matrix Plot

corrdata <- prio_grid_isd %>% select(bdist3, capdist, barren_gc, mountains_mean,
				     water_gc, sp_os_i_sum, distcoast, popd,
				     state_based, non_state, org3, deaths)

corrdata$geometry <- NULL

corrdata_mat <- cor(corrdata, method = c("spearman"), 
                     use="pairwise.complete.obs")

corrdata_mat <- rcorr(as.matrix(corrdata))

# TODO: Print to file
corrplot(corrdata_mat$r, method='color', diag=F, addCoef.col = "black", 
         p.mat = corrdata_mat$P, insig = "blank", tl.col = "black",
         tl.srt = 45)

# TODO: Print to file

# Summary Statistics

sumStats <- select(prio_grid_isd, statebaseddeaths, state_based, sp_os_i_sum,
		   bdist3, capdist, barren_gc, mountains_mean, water_gc,
		   distcoast, popd)

sumStats <- st_drop_geometry(sumStats)

stargazer(sumStats, median = FALSE, digits=1, title = "Summary Statistics",
	  column.sep.width = "1pt",
	  label = "summarystats",
	  covariate.labels = c("Fatalities", "Conflict events", "State presence", 
			       "Distance to boundary", "Distance to capital", 
			       "Barren", "Mountainous", "Water",
			       "Distance to coast"),
	  float.env = "sidewaystable",
	  out = "../Output/summaryStats.tex")
# }}}

# {{{ First analysis
#==============================================================================#
#	Analysis							       #
#==============================================================================#

main_models <- analysis(dvs = dvs, ivs = ivs, controls = controls, data =
			filter(prio_grid_isd, popd > 0), test_label = 'Linear
		Models')

interaction_models <- analysis(dvs = dvs, ivs = interactions, controls =
			       controls, data = filter(prio_grid_isd, popd > 0),
		       test_label = 'Interaction Models')

#==============================================================================#
# Controlling for French and British colonies

deathscol <- glm.nb(statebaseddeaths ~ sqrtSpAll * logCapdist + logPopd +
		    mountains_mean + water_gc + logCDist + logBDist + fra + gbr,
	    filter(prio_grid_isd, prio_grid_isd$popd > 0))
summary(deathscol)

#==============================================================================#
#	Marginal interacation plots 					       #
#==============================================================================#

deatsMainInt <- glm.nb(statebaseddeaths ~ sqrtSpAll * logCapdist +
		       mountains_mean + region3 + water_gc + barren_gc +
		       logCDist + logPopd + logBDist, data =
		       filter(prio_grid_isd, prio_grid_isd$popd > 0))


ggDeathsInt <- ggpredict(deathsMainInt, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 6.27, 7.817]"))

ggDeathsIntPlot <- ggplot(ggDeathsInt, aes(x^2, predicted, color = group)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels) +
		    name = 'Distance to capital', labels = c('Minimum', 'Mean',
							     'Maximum')) +
					xlab('State presence') +
					ylab('Predicted fatalities') +
					 geom_line() + goldenScatterCAtheme

pdf("../Output/deathsIntPlot.pdf",
    width = 10, height = 10/1.68)
ggDeathsIntPlot 
dev.off()

ggStateBased <- glm.nb(state_based ~ sqrtSpAll * logCapdist + mountains_mean +
		   region3 + water_gc + barren_gc + logCDist + logPopd + logBDist, data =
		    filter(prio_grid_isd, prio_grid_isd$popd > 0))

ggStateEffect <- ggeffect(ggStateBased, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 6.27, 7.817]"))

ggStatePlot <- ggplot(ggStateEffect, aes(x^2, predicted, color = group)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels,
		    name = 'Distance to capital', labels = c('Minimum', 'Mean',
							     'Maximum')) +
					 geom_line() + goldenScatterCAtheme

pdf("../Output/ggStatePlot.pdf",
    width = 10, height = 10/1.68)
ggStatePlot 
dev.off()


#==============================================================================#
# Regression tables
# TODO: Add coefficient maps for proper names in regression tables

for (i in 1:length(dvs)) { 
  name <- dvs[i]
  filename <- paste("../Output/",name,".tex",sep="") 
  texreg(main_models$models[[i]],
         file = filename,
         custom.model.names = control_names,
         custom.coef.map = coefs,
         stars = c(0.001, 0.01, 0.05, 0.1), 
         sideways = T, use.packages = F, scalebox = 1,
         caption = captions[i],
	 label = name,
	 booktabs = T,
         table = T)
}

#==============================================================================#
# Interaction tabels

for (i in 1:length(dvs)) { 
  name <- dvs[i]
  filename <- paste("../Output/interaction_",name,".tex",sep="") 
  intlabel <- paste("interaction_",name,sep="")
  texreg(interaction_models$models[[i]],
         file = filename,
         custom.model.names = control_names,
         custom.coef.map = coefs,
         stars = c(0.001, 0.01, 0.05, 0.1), 
         sideways = T, use.packages = F, scalebox = 1,
         caption = captions_int[i],
	 label = intlabel,
	 booktabs = T,
         table = T)
}
# }}}

# {{{ Spatial analysis
#==============================================================================#
#	Spatial analysis					       	       #
#==============================================================================#

prio_grid_shp <- st_read('../Data/PRIO-Grid/priogrid_cell.shp')

prio_grid <- read_csv('../Data/PRIO-Grid/priogridyv50-10.csv') %>% 
  as_tibble() %>% filter( (gwno >= 404 & gwno <= 626) | gwno == 651) %>% 
  group_by(gid) %>% 
  summarise(bdist3 = mean(bdist3),
  capdist = mean(capdist), excluded = mean(excluded), 
  	  temp_sd = sd(temp), gwno = last(gwno), temp = mean(temp), 
	  prec_sd = sd(prec_gpcc, na.rm = TRUE), prec_gpcc = mean(prec_gpcc))

prio_grid_static  <- read_csv('../Data/PRIO-Grid/PRIO-GRID Static Variables - 2021-06-04.csv') 

# Merging static and aggregated prio data
prio_grid  <- left_join(prio_grid, prio_grid_static, by = c("gid")) 

# Merging with the grid shape 
prio_grid <- left_join(prio_grid_shp, prio_grid, by = c("gid")) %>% 
  filter( (gwno >= 404 & gwno <= 626) | gwno == 651)

prio_sp <- as(prio_grid, Class = "Spatial")

nbQueen <- poly2nb(prio_sp, queen = T, row.names = prio_sp$gid)

lw = nb2listw(nbQueen, style = "W", zero.policy = TRUE)

moran.test(prio_grid_isd$deaths, listw = lw, zero.policy = TRUE)

prio_grid_isd$deaths_l <- lag.listw(lw, prio_grid_isd$deaths, zero.policy = T)

prio_grid_isd$state_based_l <- lag.listw(lw, prio_grid_isd$state_based, zero.policy = T)

prio_grid_isd$non_state_l <- lag.listw(lw, prio_grid_isd$non_state, zero.policy = T)

prio_grid_isd$org3_l <- lag.listw(lw, prio_grid_isd$org3, zero.policy = T)
# }}}

# {{{ ZINB analysis
#==============================================================================#
#	ZINB  								       #
#==============================================================================#		

zinb_deaths <- zeroinfl(deaths ~ sqrtSpAll * logCapdist + mountains_mean +
			region3 + water_gc + logCDist + logPopd + logBDist, data
		= filter(prio_grid_isd, popd > 0), dist = "negbin")

summary(zinb_deaths)


zinb_sb <- zeroinfl(state_based ~ sqrtSpAll * logCapdist + mountains_mean +
 			water_gc + logCDist + logPopd + logBDist, data =
			filter(prio_grid_isd, popd > 0), dist = "negbin")

summary(zinb_sb)

ggzinb<- ggpredict(zinb_deaths, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 6.27, 7.817]"))

ggzinbPlot <- ggplot(ggzinb, aes(x^2, predicted)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + 
		scale_fill_manual(values = pastels,
			name = 'Distance to capital', labels = c('Minimum',
				 'Mean', 'Maximum')) +
		xlab('State presence') +
		ylab('Predicted fatalities') +
		geom_line(aes(x^2, predicted, color = group),
			  show.legend = F) +
		goldenScatterCAtheme

pdf("../Output/zinbplot.pdf",
    width = 10, height = 10/1.68)
ggzinbPlot
dev.off()

ggzinbSB <- ggpredict(zinb_sb, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 6.27, 7.817]"))


ggzinbSBPlot <- ggplot(ggzinbSB, aes(x^2, predicted)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + 
		scale_fill_manual(values = pastels, 
				  name = 'Distance to capital', labels = 
					  c('Minimum', 'Mean', 'Maximum')) +
		xlab('State presence') +
		ylab('Predicted additional events') +
		geom_line(aes(x^2, predicted, color = group),
			  show.legend = F) +
		goldenScatterCAtheme

pdf("../Output/sbzinbplot.pdf",
    width = 10, height = 10/1.68)
ggzinbSBPlot
dev.off()

# Controlling for france and gbr
zinbgbr <- zeroinfl(deaths ~ sqrtSpAll * gbr + mountains_mean + water_gc +
		     region3 + logCDist + logPopd + logBDist, data =
		     filter(prio_grid_isd, popd > 0), dist = "negbin")

ggzinbgbr <- ggpredict(zinbgbr, terms = c("sqrtSpAll [0:15]", "gbr
						 [1, 0]"))

ggzinbgbrPlot <- ggplot(ggzinbgbr, aes(x^2, predicted, color = group)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels) +
					xlab('State presence') +
					ylab('Predicted fatalities') +
					 geom_line() + goldenScatterCAtheme

pdf("../Output/gbrzinbplot.pdf",
    width = 10, height = 10/1.68)
ggzinbgbrPlot
dev.off()

zinbfra <- zeroinfl(deaths ~ sqrtSpAll * fra + mountains_mean + water_gc +
		     region3 + logCDist + logPopd + logBDist, data =
		     filter(prio_grid_isd, popd > 0), dist = "negbin")

ggzinbfra <- ggpredict(zinbfra, terms = c("sqrtSpAll [0:15]", "fra
						 [1, 0]"))

ggzinbfraPlot <- ggplot(ggzinbfra, aes(x^2, predicted, color = group)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels) +
					xlab('State presence') +
					ylab('Predicted fatalities') +
					 geom_line() + goldenScatterCAtheme

pdf("../Output/frazinbplot.pdf",
    width = 10, height = 10/1.68)
ggzinbfraPlot
dev.off()

# Outputting regression tables

zinblist <- list(zinb_deaths, zinb_sb, zinbDcol)

texreg(zinblist, file = "../Output/zinb.tex", 
       stars = c(0.001, 0.01, 0.05, 0.1), 
       use.packages = F, 
       scalebox = .7,
       custom.coef.map = zinbcoefs,
       label = "zinb", 
       table = T,
       custom.model.names = c('Fatalities', 'Conflict events', 'Fatalities
			      (colonial controls)'),
       booktabs = T)
# }}}

# {{{ Histogram
#==============================================================================#
#	Histogram							       #
#==============================================================================#

histplot <- ggplot(filter(prio_grid_isd, statebaseddeaths < 20),
		   aes(statebaseddeaths)) + geom_histogram(binwidth = 1) +
					   goldenScatterCAtheme

histplot

pdf("../Output/histplot.pdf",
    width = 10, height = 10/1.68)
histplot
dev.off()
# }}}

# {{{ Table of atlas maps
#==============================================================================#
#	Table of atlas maps 						       #
#==============================================================================#

tb <- data.frame(x = c(
			  "\\citet{mcevedy1996penguin}",
			  "\\citet{Flint1976}",
			  "\\citet{Gailey1967}",
			  "\\citet{Oliver1985}",
			  "\\citet{Ajayi1985}",
			  "\\citet{Reid2012}",
			  "The Times atlas of world history
			  (\\citeyear{1978TTao})",
			  "\\citet{Kasule1998}"))

tb$x <- tb[order(tb$x), ]

names(tb) <- "Atlas maps"

xtb <- xtable(tb, caption = "List of maps from historical atlases
used in the Geo-ISD", label = "atlasmaps", auto = T)

print(xtb, tabular.environment = "tabularx", booktabs = T, width
      = "\\textwidth", include.rownames = F, file =
	      "../Output/atlasmaps.tex", sanitize.text.function = function(x) x)
# }}}

# {{{ Random sh*t
#==============================================================================#
#	Danger Zone!						      	       #
#==============================================================================#

testmodel <- glm.nb(dumState ~ sqrtSpAll * logCapdist + mountains_mean +
		    water_gc + barren_gc + distcoast + logPopd + bdist3, data =
		    prio_grid_isd) 
summary(testmodel)

# Just plotting the main independent variable

sqrtSpAllPlot <- ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes(fill = sqrtSpAll*logCapdist),
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()

pdf("../Output/sqrtSpAll.pdf",
	width = 10, height = 10/1.68)
  	sqrtSpAllPlot
dev.off()

logSpAllPlot <- ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes_string(fill = "logSpAll"),
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()

pdf("../Output/logSpAll.pdf",
	width = 10, height = 10/1.68)
  	logSpAllPlot
dev.off()

logOrg3 <- ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes_string(fill = "logOrg3"),
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()

pdf("../Output/logOrg3.pdf",
    width = 10, height = 10/1.68)
logOrg3
dev.off()

ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes_string(fill = "dumOrg3"),
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()

# Alternative measure of state presece


altmodelD <- glm.nb(statebaseddeaths ~ sqrtSpAny * logCapdist +
			mountains_mean + water_gc + barren_gc + logCDist +
			region3 + logPopd + logBDist, data =
			filter(prio_grid_isd, prio_grid_isd$popd > 0))

summary(altmodelD)

ggAltD <- ggeffect(altmodelD, terms = c("sqrtSpAny [0:15]", "logCapdist
						 [1.309, 6.27, 7.817]"))

ggAltDPlot <- ggplot(ggAltD, aes(x^2, predicted, color = group)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels) +
					xlab('State presence') +
					ylab('Predicted fatalities') +
					 geom_line() + goldenScatterCAtheme

ggAltDPlot 

ggStateBased <- glm.nb(state_based ~ sqrtSpAll * logCapdist + mountains_mean +
		    water_gc + barren_gc + logCDist + logPopd + logBDist, data =
		    filter(prio_grid_isd, prio_grid_isd$popd > 0))


# Colonial rulers

cols <- glm.nb(statebaseddeaths ~ sqrtSpAny * logCapdist +
			mountains_mean + water_gc + barren_gc + logCDist +
			region3 + logPopd + logBDist + gbr, data =
			filter(prio_grid_isd, prio_grid_isd$popd > 0))

summary(cols)

ggAltD <- ggeffect(altmodelD, terms = c("sqrtSpAny [0:15]", "logCapdist
						 [1.309, 6.27, 7.817]"))

ggAltDPlot <- ggplot(ggAltD, aes(x^2, predicted, color = group)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels) +
					xlab('State presence') +
					ylab('Predicted fatalities') +
					 geom_line() + goldenScatterCAtheme

ggAltDPlot 

# 3D plot experiment

geoisd <- st_read('../../QGIS/Geo-ISD.shp')

gisdnig <- filter(geoisd, COWID == 4798 | COWID == 4752 | COWID == 4521 | COWID
		  == 4327 | COWID == 4831 | COWID == 4776 | COWID == 4763 | COWID ==
434 | COWID == 4362 | COWID == 4768 | COWID == 4769 | COWID == 4751 | COWID == 4771 |
COWID == 4742 | COWID == 4773 | COWID == 4765 | COWID == 4775 | COWID == 4832)

gisdnig <- st_make_valid(gisdnig)

borders <- read_sf(dsn = "../Data/Shapes/Africa.shp") %>% 
	filter(ID == 626)

borders <- st_set_crs(borders, 4326)

ext <- extent(borders)

grid <- st_bbox(ext) %>% 
  st_make_grid(cellsize = (0.05), what = "polygons", flat_topped = T) %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))
notify(msg = c("Done!") # Comment this out if larger cell size


centroids <- st_centroid(grid)

grid$lon <- st_coordinates(centroids)[,1]
grid$lat <- st_coordinates(centroids)[,2]

#grid$sp1 <- lengths(st_intersects(grid, gisdnig))

#grid <- filter(grid, sp > 0)

grid <- grid %>% mutate(sp2 = lengths(st_within(grid, gisdnig)))
notify(msg = c("Done!") # Comment this out if larger cell size

grid <- na.omit(grid)
notify(msg = c("Done!") # Comment this out if larger cell size

# setNames(data.frame(coords[[1]], 
#                     matrix(unlist(coords[2]), ncol=2, byrow=TRUE)), 
#          c("ID", "lon", "lat"))

gridtest <- ggplot() +
		   geom_sf(data = grid,
			   linetype = 0,
			   aes(fill = sp2),
			   show.legend = F) +
    		   scale_fill_viridis_c() +
		   #geom_contour(data = filter(grid, sp2 > 1),
				#aes(x = round(lon, 1), y = round(lat, 1), z = sp2)) +
		   #geom_sf(data = borders, color = "gray", fill = NA) +
		   theme_minimal()

pdf("../Output/nigeriatest.pdf",
    width = 10, height = 10/1.683)
gridtest 
dev.off()

plot_gg(gridtest, pointcontract = 1, scale = 100, offset_edges = F, triangulate = F, verbose = T, zoom = 0.5, multicore = T)

render_snapshot("../Output/3DNigeria")

render_movie(filename = "../Output/Libyathemovie.gif")

#render_snapshot("../Output/3DLibya.html", webshot = T)
# }}}

