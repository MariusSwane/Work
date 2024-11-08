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
library(countrycode)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(haven)
library(Hmisc)
library(margins)
library(MASS)
#library(notifier)
library(pscl)
library(purrr)
library(rayshader)
library(RColorBrewer)
library(readr)
library(raster)
library(scales)
library(spdep)
library(sf)
library(sidedata)
library(stargazer)
library(summarytools)
library(texreg)
library(viridis)
library(xtable)

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

prio_grid_isd <- prio_grid_isd %>% filter( (gwno >= 404 & gwno <= 626) | gwno == 651 & gwno !=
			 581 & gwno != 590) %>% 
	mutate(
				sqrtDeaths = sqrt(deaths),
				sqrtSBDeaths = sqrt(statebaseddeaths),
				sqrtInterDeaths = sqrt(interdeaths),
				sqrtboth = sqrt(both),
				sqrtState_based = sqrt(state_based),
				sqrtSpAll = sqrt(sp_os_i_sum),
				sqrtSpAny = sqrt(sp_i_sum_any),
				sqrtSpNoInt = sqrt(sp_os_sum),
				sqrtCapdist = sqrt(capdist),
				sqrtBDist = sqrt(bdist3),
				sqrtPopd = sqrt(popd),
				logDeaths = log(deaths + 1),
				logInterDeaths = log(interdeaths + 1),
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
      m1 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls)), data=data)
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

virs <- viridis(3)

pastels <- NULL
for (i in 1:length(virs)) {
	pastels[i] <- lighten(virs[i])
}

#==============================================================================#
#	Variables							       #
#==============================================================================#

controls <- c('+ mountains_mean + water_gc + barren_gc + logCDist')

extended_controls <- c('+ region3')

even_more <- c('+ logPopd')

rest <- c('+ logBDist')

extra <- c('+ excluded')

#climate_controls <- c('+ temp_sd + temp + prec_sd + prec_gpcc')

coefs <- list('sqrtSpAll' = 'Precolonial state presence (sqrt)',
	      'sqrtSpAny' = 'Precolonial state presence (sqrt)(alt)',
	      	'logCapdist' = 'Distance to capital (log)',
		'sqrtSpAll:logCapdist' = 'Interaction term',
		'sqrtSpAny:logCapdist' = 'Interaction term (alt)',
		'mountains_mean' = 'Mountainous terrain',
	     	'water_gc' = 'Water (%)', 
	 	'barren_gc' = 'Barren (%)', 
	      	'logCDist' = 'Distance to coast (log)',
		'region3' = 'North Africa',
	      	'logPopd' = 'Population density (log)', 
		'logBDist' = 'Distance to international boundary (log)',
		'excluded' = 'Number of EPR excluded groups'
		)

control_names <-c('Geography', 'North Africa', 'Population densisty', 'Distance
		  to border')

dvs <- c('interdeaths', 'both')

captions <- c('Fatalities', 'Internal and internationalized conflict events')

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

#corrdata <- prio_grid_isd %>% select(bdist3, capdist, barren_gc, mountains_mean,
#				     water_gc, sp_os_i_sum, distcoast, popd,
#				     state_based, non_state, org3, deaths)
#
#corrdata$geometry <- NULL
#
#corrdata_mat <- cor(corrdata, method = c("spearman"), 
#                     use="pairwise.complete.obs")
#
#corrdata_mat <- rcorr(as.matrix(corrdata))

#corrplot(corrdata_mat$r, method='color', diag=F, addCoef.col = "black", 
#         p.mat = corrdata_mat$P, insig = "blank", tl.col = "black",
#         tl.srt = 45)


# Summary Statistics

sumStats <- select(prio_grid_isd, interdeaths, both, sp_os_i_sum,
		   bdist3, capdist, barren_gc, mountains_mean, water_gc,
		   distcoast, popd)

sumStats <- st_drop_geometry(sumStats)

stargazer(sumStats, median = FALSE, iqr = T, digits=1, title = "Summary Statistics",
	  column.sep.width = "1pt",
	  label = "summarystats",
	  covariate.labels = c("Fatalities", "Conflict events", "State presence", 
			       "Distance to boundary", "Distance to capital", 
			       "Barren", "Mountainous", "Water",
			       "Distance to coast", "Population density 1600"),
	  float.env = "sidewaystable",
	  out = "../Output/summaryStats.tex")

# Nigeriaplot

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

ext2 <- extent(2.17, 15.1, 3.7, 14.4)

grid <- st_bbox(ext2) %>% 
  st_make_grid(cellsize = (0.05), what = "polygons", flat_topped = T) %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))
notify(msg = c("Done!")) # Comment this out if larger cell size

centroids <- st_centroid(grid)

grid$lon <- st_coordinates(centroids)[,1]
grid$lat <- st_coordinates(centroids)[,2]

#grid$sp1 <- lengths(st_intersects(grid, gisdnig))

#grid <- filter(grid, sp > 0)

grid <- grid %>% mutate(sp2 = lengths(st_within(grid, gisdnig)))
notify(msg = c("Done!")) # Comment this out if larger cell size

#grid <- na.omit(grid)
#notify(msg = c("Done!")) # Comment this out if larger cell size

# setNames(data.frame(coords[[1]], 
#                     matrix(unlist(coords[2]), ncol=2, byrow=TRUE)), 
#          c("ID", "lon", "lat"))

gridtest <- ggplot() +
		   geom_sf(data = grid,
			   linetype = 0,
			   aes(fill = sp2)) +
    		   scale_fill_viridis_c("Pre-colonial
state presence") +
		   geom_contour(data = filter(grid, sp2 > 1),
				aes(x = round(lon, 1), y = round(lat, 1), z = sp2)) +
		   geom_sf(data = borders, color = "gray", size = 2, fill = NA) +
		   xlab("") + ylab("") +
		   theme_minimal()
gridtest 

pdf("../Output/nigeria.pdf",
    width = 10, height = 10/1.683)
gridtest 
dev.off()

# }}}

# {{{ First analysis
#==============================================================================#
#	Analysis							       #
#==============================================================================#

# TODO: does not converge atm
main_models <- analysis(dvs = dvs, ivs = ivs, controls = controls, data =
			filter(prio_grid_isd, popd > 0), test_label = 'Linear
		Models')

# Does converge
interaction_models <- analysis(dvs = dvs, ivs = interactions, controls =
			       controls, data = filter(prio_grid_isd, popd > 0),
		       test_label = 'Interaction Models')


#==============================================================================#
#	Marginal interacation plots 					       #
#==============================================================================#

deathsMainInt <- glm.nb(interdeaths ~ sqrtSpAll * logCapdist +
		       mountains_mean + region3 + water_gc + barren_gc +
		       logCDist + logPopd + logBDist, data =
		       filter(prio_grid_isd, prio_grid_isd$popd > 0))


ggDeathsInt <- ggpredict(deathsMainInt, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

ggDeathsIntPlot <- ggplot(ggDeathsInt, aes(x^2, predicted)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + 
			scale_fill_manual(values = pastels, 
					  name = 'Distance to capital', 
					  labels = c('Minimum', 'Maximum')) +
			xlab('State presence') +
			ylab('Predicted fatalities') +
			geom_line(aes(x^2, predicted, color = group), 
				      show.legend = F) +
			scale_color_manual(values = virs) +
			goldenScatterCAtheme

pdf("../Output/deathsInterPlot.pdf",
    width = 10, height = 10/1.68)
ggDeathsIntPlot 
dev.off()

ggboth <- glm.nb(both ~ sqrtSpAll * logCapdist + mountains_mean +
		   region3 + water_gc + barren_gc + logCDist + logPopd + logBDist, 
	   data = filter(prio_grid_isd, prio_grid_isd$popd > 0))

ggBothEffect <- ggeffect(ggboth, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

ggBothPlot <- ggplot(ggBothEffect, aes(x^2, predicted)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
				linetype = NA)) + 
		scale_fill_manual(values = pastels, 
				  name = 'Distance to capital', 
				  labels = c('Minimum', 'Maximum')) +
		xlab('State presence') +
		ylab('Predicted conflict events') +
		geom_line(aes(x^2, predicted, color = group), 
				show.legend = F) +
		scale_color_manual(values = virs) +
		goldenScatterCAtheme

pdf("../Output/ggBothPlot.pdf",
    width = 10, height = 10/1.68)
ggBothPlot 
dev.off()


#==============================================================================#
# Regression tables

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
         caption = captions[i],
	 label = intlabel,
	 booktabs = T,
         table = T)
}
# }}}

# {{{ Spatial analysis
#==============================================================================#
#	Spatial analysis					       	       #
#==============================================================================#

nbQueen <- poly2nb(prio_grid_isd, queen = T, row.names = prio_grid_isd$gid)

lw = nb2listw(nbQueen, style = "W", zero.policy = TRUE)

moran.test(prio_grid_isd$interdeaths, listw = lw, zero.policy = TRUE)

prio_grid_isd$deaths_l <- lag.listw(lw, prio_grid_isd$interdeaths, zero.policy = T)

prio_grid_isd$both_l <- lag.listw(lw, prio_grid_isd$both, zero.policy = T)

prio_grid_isd$non_state_l <- lag.listw(lw, prio_grid_isd$non_state, zero.policy = T)

prio_grid_isd$org3_l <- lag.listw(lw, prio_grid_isd$org3, zero.policy = T)

# Analysis

zinb_spatial_deaths <- zeroinfl(interdeaths ~ sqrtSpAll * logCapdist + mountains_mean +
			deaths_l + region3 + water_gc + logCDist + logPopd +
			logBDist, data = filter(prio_grid_isd, popd > 0), dist =
			"negbin")

summary(zinb_spatial_deaths)

gglagzinb <- ggpredict(zinb_spatial_deaths, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

gglagzinbPlot <- ggplot(gglagzinb, aes(x^2, predicted)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + 
		scale_fill_manual(values = pastels,
			name = 'Distance to capital', labels = c('Minimum',
				 'Maximum')) +
		xlab('State presence') +
		ylab('Predicted fatalities') +
		geom_line(aes(x^2, predicted, color = group),
			  show.legend = F) +
		scale_color_manual(values = virs) +
		goldenScatterCAtheme

pdf("../Output/spatialzinbdeathsplot.pdf",
    width = 10, height = 10/1.68)
gglagzinbPlot
dev.off()

zinb_spatial_both <- zeroinfl(both ~ sqrtSpAll * logCapdist + mountains_mean +
			both_l + region3 + water_gc + logCDist + logPopd +
			logBDist, data = filter(prio_grid_isd, popd > 0), dist =
			"negbin")

summary(zinb_spatial_both)

spatiallist <- list(zinb_spatial_deaths, zinb_spatial_both)

texreg(spatiallist,
        file = "../Output/spatial_count.tex",
        custom.model.names = c("Fatalities", "Events"),
        custom.coef.map = coefs,
	include.zero = F,
        stars = c(0.001, 0.01, 0.05, 0.1), 
        sideways = F, use.packages = F, scalebox = 1,
        caption = "Spatial lag models (count)",
	label = "spatialCount",
	booktabs = T,
        table = T)

texreg(spatiallist,
        file = "../Output/spatial_zero.tex",
        custom.model.names = c("Fatalities", "Events"),
        custom.coef.map = coefs,
	include.count = F,
        stars = c(0.001, 0.01, 0.05, 0.1), 
        sideways = F, use.packages = F, scalebox = 1,
        caption = "Spatial lag models (zero)",
	label = "spatialZero",
	booktabs = T,
        table = T)


bothlagzinb <- ggpredict(zinb_spatial_both, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

bothspatzinbPlot <- ggplot(bothlagzinb, aes(x^2, predicted)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + 
		scale_fill_manual(values = pastels,
			name = 'Distance to capital', labels = c('Minimum',
				 'Maximum')) +
		xlab('State presence') +
		ylab('Predicted fatalities') +
		geom_line(aes(x^2, predicted, color = group),
			  show.legend = F) +
		scale_color_manual(values = virs) +
		goldenScatterCAtheme

pdf("../Output/spatialzinbbothplot.pdf",
    width = 10, height = 10/1.68)
bothspatzinbPlot
dev.off()

# }}}

# {{{ ZINB analysis
#==============================================================================#
#	ZINB  								       #
#==============================================================================#		

zinb_deaths <- zeroinfl(interdeaths ~ sqrtSpAll * logCapdist + mountains_mean +
			region3 + water_gc + barren_gc + logCDist + logPopd + 
			logBDist, data
		= filter(prio_grid_isd, popd > 0), dist = "negbin")
summary(zinb_deaths)


zinb_both <- zeroinfl(both ~ sqrtSpAll * logCapdist + mountains_mean +
 			water_gc + barren_gc + logCDist + logPopd + logBDist,
		data = filter(prio_grid_isd, popd > 0), dist = "negbin")
summary(zinb_both)

ggzinb <- ggpredict(zinb_deaths, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

ggzinbPlot <- ggplot(ggzinb, aes(x^2, predicted)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + 
		scale_fill_manual(values = pastels,
			name = 'Distance to capital', labels = c('Minimum',
				 'Maximum')) +
		xlab('State presence') +
		ylab('Predicted additional fatalities') +
		geom_line(aes(x^2, predicted, color = group),
			  show.legend = F) +
		scale_color_manual(values = virs) +
		goldenScatterCAtheme

pdf("../Output/interdeathszinbplot.pdf",
    width = 10, height = 10/1.68)
ggzinbPlot
dev.off()

ggzinbboth <- ggpredict(zinb_both, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

ggzinbbothPlot <- ggplot(ggzinbboth, aes(x^2, predicted)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + 
		scale_fill_manual(values = pastels, 
				  name = 'Distance to capital', labels = 
					  c('Minimum', 'Maximum')) +
		xlab('State presence') +
		ylab('Predicted additional events') +
		geom_line(aes(x^2, predicted, color = group),
			  show.legend = F) +
		scale_color_manual(values = virs) +
		goldenScatterCAtheme

pdf("../Output/bothzinbplot.pdf",
    width = 10, height = 10/1.68)
ggzinbbothPlot
dev.off()

# Outputting regression tables

zinblist <- list(zinb_deaths, zinb_both)

texreg(zinblist, file = "../Output/zinbc.tex", 
       stars = c(0.001, 0.01, 0.05, 0.1), 
       use.packages = F, 
       scalebox = 1,
       custom.coef.map = coefs,
       include.zero = F,
       label = "zinbc", 
       table = T,
       caption = "ZINB models (count)",
       custom.model.names = c('Fatalities', 'Conflict events'),
       booktabs = T)

texreg(zinblist, file = "../Output/zinbz.tex", 
       stars = c(0.001, 0.01, 0.05, 0.1), 
       use.packages = F, 
       scalebox = 1,
       custom.coef.map = coefs,
       include.count = F,
       label = "zinbz", 
       table = T,
       caption = "ZINB models (zero)",
       custom.model.names = c('Fatalities', 'Conflict events'),
       booktabs = T)

# }}}

#==============================================================================#
# {{{ Robustness checks
#==============================================================================#

# EPR groups
#==============================================================================#

zinbepr <- zeroinfl(interdeaths ~ sqrtSpAll * logCapdist + mountains_mean + water_gc +
		     region3 + logCDist + logPopd + logBDist + excluded, data =
		     filter(prio_grid_isd, popd > 0), dist = "negbin")
summary(zinbepr)

# Alternative IV
#==============================================================================#

zinbaiv <- zeroinfl(interdeaths ~ sqrtSpAny * logCapdist + mountains_mean + water_gc +
		     region3 + logCDist + logPopd + logBDist, data =
		     filter(prio_grid_isd, popd > 0), dist = "negbin")
summary(zinbaiv)

ggaivzinb <- ggpredict(zinbaiv, terms = c("sqrtSpAny [0:10]", "logCapdist
						 [1.309, 7.817]"))

ggzinbaivPlot <- ggplot(ggaivzinb, aes(x^2, predicted)) +
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
			linetype = NA)) + 
		scale_fill_manual(values = pastels,
			name = 'Distance to capital', labels = c('Minimum',
				 'Maximum')) +
		xlab('State presence') +
		ylab('Predicted fatalities') +
		geom_line(aes(x^2, predicted, color = group),
			  show.legend = F) +
		scale_color_manual(values = virs) +
		goldenScatterCAtheme
ggzinbaivPlot

# Subsamples for france and gbr
#==============================================================================#

zinbgbr <- zeroinfl(interdeaths ~ sqrtSpAll * logCapdist + mountains_mean + water_gc +
		     region3 + logCDist + logPopd + logBDist, data =
		     filter(prio_grid_isd, popd > 0 & gbr == 0), dist = "negbin")
summary(zinbgbr)


zinbfra <- zeroinfl(interdeaths ~ sqrtSpAll * logCapdist + mountains_mean + water_gc +
		     region3 + logCDist + logPopd + logBDist, data =
		     filter(prio_grid_isd, popd > 0 & fra == 0), dist = "negbin")
summary(zinbfra)

# Prior conflict (🚧 WIP 🚧)
#==============================================================================#

conf <- read_tsv("../Data/dfo_historical_conflict_data_ehdr.tab")

conf <- subset(conf, conf$continent == "Sub-Saharan Africa" | conf$continent ==
	       "North Africa")

conf$gwno <- countrycode(conf$country, "country.name", "gwn") 

conf %>% group_by(gwno) %>% summarise(confls = length(gwno))

gisd_cntry %>% prio_grid_isd %>% group_by(gwno) %>% 
	summarise(interdeaths = sum(interdeaths),
		  spm = mean(sqrtSpAll),
		  sps = sum(sqrrtSpAll),


# Square root transformation
#==============================================================================#
sqrtmod <- lm(sqrt(interdeaths) ~ sqrtSpAll * logCapdist + mountains_mean + 
	      water_gc + barren_gc + region3 + logCDist + logPopd + logBDist,
      data = filter(prio_grid_isd, popd > 0))
summary(sqrtmod)

sqrtpred <- ggeffect(sqrtmod, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

sqrtplot <- ggplot(sqrtpred, aes(x^2, predicted^2, color = group)) +
	geom_ribbon(aes(ymin = conf.low^2, ymax = conf.high^2, fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels) +
					xlab('State presence') +
					ylab('Predicted fatalities') +
					 geom_line() + goldenScatterCAtheme
sqrtplot


# Sensitivity plots
#==============================================================================#

senseanalysis <- sensemakr(sqrtmod, "sqrtSpAll:logCapdist", benchmark_covariates
= "logPopd", kd = 1:2)

summary(senseanalysis)


pdf("../Output/senseplot.pdf",
    width = 10, height = 10/1.68)
plot(senseanalysis, sensitivity.of = "t-value")
dev.off()

# Country fixed effcts (not converging, "non finite value supplied by oprim")
#==============================================================================#

zinbfxe <- zeroinfl(interdeaths ~ sqrtSpAll * logCapdist + mountains_mean +
		     water_gc + region3 + logCDist + logPopd + logBDist +
		     factor(gwno),
	     data = filter(prio_grid_isd, popd > 0), dist = "negbin")
summary(zinbfxe)

nbfxe <- glm.nb(interdeaths ~ sqrtSpAll * logCapdist + mountains_mean +
		     water_gc + region3 + logCDist + logPopd + logBDist +
		     factor(gwno),
	     data = filter(prio_grid_isd, popd > 0))
summary(nbfxe)

# Log-transformed IV
zinblog <- zeroinfl(interdeaths ~ logSpAll * logCapdist + mountains_mean +
		     water_gc + region3 + logCDist + logPopd + logBDist + 
		     factor(gwno),
	     data = filter(prio_grid_isd, popd > 0), dist = "negbin")
summary(zinblog)

# Non-transformed IV (produces NaN's)
zinbnt <- zeroinfl(interdeaths ~ sp_os_i_sum * logCapdist + mountains_mean +
		     water_gc + region3 + logCDist + logPopd + logBDist +
		     factor(gwno),
	     data = filter(prio_grid_isd, popd > 0), dist = "negbin")
summary(zinbnt)

# Tse-tse fly
#==============================================================================#

tsetse <- read_dta("../Data/Alsan tse-tse/precolonial.dta", col_select =
		   c("TSI", "isocode", "malaria_index"))

tsetse$gwno <- countrycode(tsetse$isocode, "iso3c", "gwn")

prio_grid_isd <- left_join(prio_grid_isd, tsetse, by = "gwno")

#tsesubnat <- read_dta("../Data/Alsan tse-tse/subnational.dta")

zinbsick <- zeroinfl(interdeaths ~ sqrtSpAll * logCapdist + mountains_mean +
		     water_gc + region3 + logCDist + logPopd + logBDist +
		     TSI + malaria_index,
	     data = filter(prio_grid_isd, popd > 0), dist = "negbin")
summary(zinbsick)

sickpred <- ggpredict(zinbsick, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

sickplot <- ggplot(sickpred, aes(x^2, predicted^2, color = group)) +
	geom_ribbon(aes(ymin = conf.low^2, ymax = conf.high^2, fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels) +
					xlab('State presence') +
					ylab('Predicted fatalities') +
					 geom_line() + goldenScatterCAtheme
sickplot

nbsick <- glm.nb(interdeaths ~ sqrtSpAll * logCapdist + mountains_mean +
		     water_gc + region3 + logCDist + logPopd + logBDist +
		     TSI + malaria_index,
	     data = filter(prio_grid_isd, popd > 0))
summary(nbsick)

nbsickpred <- ggpredict(nbsick, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

nbsickplot <- ggplot(nbsickpred, aes(x^2, predicted^2, color = group)) +
	geom_ribbon(aes(ymin = conf.low^2, ymax = conf.high^2, fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels) +
					xlab('State presence') +
					ylab('Predicted fatalities') +
					 geom_line() + goldenScatterCAtheme
sickplot

# Nunn ruggedness
#==============================================================================#

ruggednunn <- read.csv("../Data/Nunn ruggedness/rugged_data.csv")

ruggednunn <- st_as_sf(ruggednunn, coords = c('lon', 'lat'), crs =
		       st_crs(prio_grid_shp))

# Robustness tables
#==============================================================================#

robustlist <- list(zinbepr, zinbaiv, zinbgbr, zinbfra, sqrtmod, zinbfxe, zinblog)

texreg(robustlist, file = "../Output/robust.tex", 
       stars = c(0.001, 0.01, 0.05, 0.1), 
       use.packages = F, 
       sideways = T,
       scalebox = .9,
       custom.coef.map = coefs,
       include.zero = F,
       label = "robustc", 
       table = T,
       caption = "Additional models (count)",
       custom.model.names = c('EPR groups', 'Alternativ IV', 
			      'Excluding former British colonies', 
			      'Excluding former French colonies',
       			'Squar root traansformed linear model',
			'Country fixed effects'),
       booktabs = T)

texreg(robustlist, file = "../Output/robustz.tex", 
       stars = c(0.001, 0.01, 0.05, 0.1), 
       use.packages = F, 
       sideways = T,
       scalebox = .9,
       custom.coef.map = coefs,
       include.count = F,
       label = "robustz", 
       table = T,
       caption = "Additional models (zero)",
       custom.model.names = c('EPR groups', 'Alternativ IV', 
			      'Excluding former British colonies', 
			      'Excluding former French colonies',
			'Country fixed effects'),
       booktabs = T)

# }}}

# {{{ Histogram
#==============================================================================#
#	Histogram							       #
#==============================================================================#

histplot <- ggplot(filter(prio_grid_isd, interdeaths < 20),
		   aes(interdeaths)) + geom_histogram(binwidth = 1) +
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

spPlot <- ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes(fill = sqrt(sp_os_i_sum))) +
    scale_fill_viridis_c("Pre-colonial
state presence",
			 labels = trans_format("identity", function(x) x^2)) +
    theme_minimal()
spPlot

pdf("../Output/spPlot.pdf",
	width = 10, height = 10/1.68)
  	spPlot
dev.off()

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

# Attempt at 3D plot
#plot_gg(gridtest, pointcontract = 1, scale = 100, offset_edges = F, triangulate = F, verbose = T, zoom = 0.5, multicore = T)
#notify(msg = c("Done!"))
#
#render_snapshot("../Output/3DNigeria")
#
#render_movie(filename = "../Output/Libyathemovie.gif")

#render_snapshot("../Output/3DLibya.html", webshot = T)


#==============================================================================#
#	Log and sqrt transformed DV					       #
#==============================================================================#		


logmod <- lm(log(interdeaths+1) ~ sqrtSpAll * logCapdist + mountains_mean + water_gc +
		     region3 + logCDist + logPopd + logBDist, data =
		     filter(prio_grid_isd, popd > 0))
summary(logmod)

logpred <- ggeffect(logmod, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 7.817]"))

logplot <- ggplot(logpred, aes(x^2, exp(predicted), color = group)) +
	geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high), fill = group,
			linetype = NA)) + scale_fill_manual(values = pastels) +
					xlab('State presence') +
					ylab('Predicted fatalities') +
					 geom_line() + goldenScatterCAtheme
logplot

# }}}

