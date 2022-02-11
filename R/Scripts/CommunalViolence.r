#==============================================================================#
#  ██████╗ ██████╗ ███╗   ███╗███╗   ███╗██╗   ██╗███╗   ██╗ █████╗ ██╗        # 
# ██╔════╝██╔═══██╗████╗ ████║████╗ ████║██║   ██║████╗  ██║██╔══██╗██║        # 
# ██║     ██║   ██║██╔████╔██║██╔████╔██║██║   ██║██╔██╗ ██║███████║██║        # 
# ██║     ██║   ██║██║╚██╔╝██║██║╚██╔╝██║██║   ██║██║╚██╗██║██╔══██║██║        # 
# ╚██████╗╚██████╔╝██║ ╚═╝ ██║██║ ╚═╝ ██║╚██████╔╝██║ ╚████║██║  ██║███████╗   # 
#  ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝   # 
#                                                                              # 
# ██╗   ██╗██╗ ██████╗ ██╗     ███████╗███╗   ██╗ ██████╗███████╗              # 
# ██║   ██║██║██╔═══██╗██║     ██╔════╝████╗  ██║██╔════╝██╔════╝              # 
# ██║   ██║██║██║   ██║██║     █████╗  ██╔██╗ ██║██║     █████╗                # 
# ╚██╗ ██╔╝██║██║   ██║██║     ██╔══╝  ██║╚██╗██║██║     ██╔══╝   	       #
#  ╚████╔╝ ██║╚██████╔╝███████╗███████╗██║ ╚████║╚██████╗███████╗              # 
#   ╚═══╝  ╚═╝ ╚═════╝ ╚══════╝╚══════╝╚═╝  ╚═══╝ ╚═════╝╚══════╝	       #
#==============================================================================#

# {{{ Header
#==============================================================================#
#	Loading Packages						       #
#==============================================================================#	

library(acled.api)
library(conflicted)
library(dplyr)
library(ggeffects)
library(ggplot2)
library(gridExtra)
library(MASS)
library(patchwork)
library(pscl)
library(purrr)
library(raster)
library(sidedata)
library(spdep)
library(sf)
library(stargazer)
library(terra)
library(texreg)

#==============================================================================#
#	Loading Data and functions					       #
#==============================================================================#	

load("../Data/GeoISDControls.Rdata")
source("goldenScatterCAtheme.r")
geoisd <- st_read('../../QGIS/Geo-ISD.shp')

#==============================================================================#
#	Resolving conflicts 						       #
#==============================================================================#		

conflict_prefer("filter", "dplyr")  
conflict_prefer("select", "dplyr")
conflict_prefer("extract", "raster")

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
				nopastor = as.factor(prec_gpcp > 400),
				region1 = as.numeric(factor(region)==1),
				region2 = as.numeric(factor(region)==2),
				region3 = as.numeric(factor(region)==3),
				region4 = as.numeric(factor(region)==4),
				region5 = as.numeric(factor(region)==5))

#==============================================================================#
#	Functions				    	             	       #
#==============================================================================#

#==============================================================================#
# Negative binomial models function

geoISDanalysis <- function(dvs, ivs, controls, data, test_label){
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
				    extended_controls, even_more)),
               data=data)
      m4 <-  glm.nb(as.formula(paste(dv, '~', ivs[j], controls,
				    extended_controls, even_more, rest)),
               data=data)
      m5 <-  glm.nb(as.formula(paste(dv, '~', ivs[j], controls,
				    extended_controls, even_more, rest,
				    pastoralism)),
               data=data)
      models_out_j[[j]] <- list(m1, m2, m3, m4, m5)
    } 
    models_out[[i]] <- models_out_j %>% flatten()
  }
  models_out <- models_out
  final <- list('models'=models_out)
  return(final)
}

# Creating function for lighter colors
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

# Creating colors, and pastels
cols <- c("red", "green", "blue")

pastels <- NULL
for (i in 1:length(cols)) {
	pastels[i] <- lighten(cols[i])
}

#==============================================================================#
#	Descriptive Statistics                                         	       #
#==============================================================================#

sumStats <- select(prio_grid_isd, org3, sp_os_i_sum, bdist3, capdist, barren_gc,
		   mountains_mean, water_gc, distcoast, popd)

sumStats <- st_drop_geometry(sumStats)

stargazer(sumStats, median = FALSE, digits=1, title = "Summary Statistics",
	  column.sep.width = "1pt", label = "summarystats", covariate.labels =
		  c("Conflict events", "State presence", "Distance to boundary",
		    "Distance to capital", "Barren", "Mountainous", "Water",
		    "Distance to coast"), out = "../Output/CVsummaryStats.tex")

#==============================================================================#
#	Variables							       #
#==============================================================================#

controls <- c('+ mountains_mean + water_gc + barren_gc + logCDist')

extended_controls <- c('+ region3')

even_more <- c('+ logPopd')

rest <- c('+ logBDist')

pastoralism <- c('+ nopastor')

climate_controls <- c('+ temp_sd + temp + prec_sd + prec_gpcc')

coefs_cv <- list('logSpAll' = 'Precolonial state presence (log)', 
		'mountains_mean' = 'Mountainous terrain',
	     	'water_gc' = 'Water (%)', 
	 	'barren_gc' = 'Barren (%)', 
	      	'logCDist' = 'Distance to coast (log)',
		'region3' = 'North Africa',
		'gbr' = 'Former British colony',
	      	'logPopd' = 'Population density (log)', 
		'logBDist' = 'Distance to border (log)', 
		'nopastorTRUE' = 'Land not suited for pastorial herding',
		'sqrtSpAll:gbr' = 'Interaction term',
		'sqrtSpAll' = 'Precolonial state presence (sqrt)')


control_names_cv <-c('Baseline', 'North Africa', 'Population density', 'Distance
		     to international boundary', 'Pastoralism')

cv_dvs <- c('non_state', 'org3')

cv_captions <- c('Non-state conflict events', 'Communal violence events')

ivs <- c('sqrtSpAll')

# }}}

# {{{ Main analysis
#==============================================================================#
#	Analysis							       #
#==============================================================================#

cv_nb_models <- geoISDanalysis(dvs = cv_dvs, ivs = ivs, controls = controls,
			       data = filter(prio_grid_isd, popd > 0),
			       test_label = 'Linear Models')

#==============================================================================#
#	Marginal interacation plots 					       #
#==============================================================================#

# Plotting
org3_NB <- glm.nb(org3 ~ sqrtSpAll + mountains_mean + water_gc + barren_gc +
		  logCDist + logBDist + logPopd + region3 + gbr, data =
		  filter(prio_grid_isd, popd > 0)) 
summary(org3_NB)

ggorg3 <- ggpredict(org3_NB, terms = "sqrtSpAll [0:15 by = .5] ")

org3plot <- ggplot(ggorg3, aes(x^2, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
	      fill = lighten("blue")) +
  geom_line(color = "blue") +
  labs(x = 'Precolonial state presence', y = 'Communial violence events') +
  goldenScatterCAtheme

# Printing to file
pdf("../Output/CommunalViolenceMargins.pdf",
    width = 10, height = 10/1.68)
org3plot
dev.off()

#==============================================================================#
#	Regression tables						       #
#==============================================================================#

for (i in 1:length(cv_dvs)) { 
  name <- cv_dvs[i]
  filename <- paste("../Output/",name,".tex",sep="") 
  texreg(cv_nb_models$models[[i]],
         file = filename,
         custom.model.names = control_names_cv,
         custom.coef.map = coefs_cv,
         stars = c(0.001, 0.01, 0.05, 0.1), 
         sideways = T, use.packages = F, scalebox = 1,
         #custom.note = "",
         caption = cv_captions[i],
	 label = name,
	 booktabs = T,
         table = T)
}
# }}}

# {{{ Spatial analysis
#==============================================================================#
#	Set up for spatial analysis				       	       #
#==============================================================================#

prio_grid_shp <- st_read('../Data/PRIO-Grid/priogrid_cell.shp')

prio_grid <- read_csv('../Data/PRIO-Grid/priogridyv50-10.csv') %>% 
  as_tibble() %>% filter( (gwno >= 404 & gwno <= 626) | gwno == 651) %>% 
  group_by(gid) %>% 
  summarise(bdist3 = mean(bdist3),
  capdist = mean(capdist), excluded = mean(excluded), 
  	  temp_sd = sd(temp), gwno = last(gwno), temp = mean(temp), 
	  prec_sd = sd(prec_gpcc, na.rm = TRUE), prec_gpcc = mean(prec_gpcc))

prio_grid_static  <- read_csv('../Data/PRIO-Grid/PRIO-GRID Static Variables -
			      2021-06-04.csv') 

# Merging static and aggregated prio data
prio_grid  <- left_join(prio_grid, prio_grid_static, by = c("gid")) 

# Merging with the grid shape 
prio_grid <- left_join(prio_grid_shp, prio_grid, by = c("gid")) %>% 
  filter( (gwno >= 404 & gwno <= 626) | gwno == 651)

prio_sp <- as(prio_grid, Class = "Spatial")

nbQueen <- poly2nb(prio_sp, queen = T, row.names = prio_sp$gid)

lw = nb2listw(nbQueen, style = "W", zero.policy = TRUE)

moran.test(prio_grid_isd$deaths, listw = lw, zero.policy = TRUE)

prio_grid_isd$non_state_l <- lag.listw(lw, prio_grid_isd$non_state, zero.policy
				       = T)

prio_grid_isd$org3_l <- lag.listw(lw, prio_grid_isd$org3, zero.policy = T)

# }}}

# {{{ ZINB
#==============================================================================#
#	ZINB  								       #
#==============================================================================#		

# Trying to be clever about it

# Part 1: models

zdvs <- c("org3", "non_state")

ziv <- "sqrtSpAll"

zcontrols <- c("+ mountains_mean + water_gc + barren_gc + logCDist + logBDist +
	       logPopd + region3")

exclusions <- c(475, 500, 501, 452)

#excluded <- c("Nigeria", "Uganda", "Kenya", "Ghana")

interactions <- c("gbr", "region1", "region5")

titles <- c('Main ZINB model', 'Excluding Nigeria', 'Excluding Uganda',
	    'Excluding Kenya', 'Excluding Ghana', 'Former British colony
	    interaction', 'East Africa interaction', 'West Africa interaction')

zinbanalysis <- function(zdvs, zivs, zcontrols, data, test_label){
  models_out <- NULL
  for(i in 1:length(zdvs)) {
    dv <- zdvs[i] 
# TODO: Do as in the next function to make prettier/more compact/clever
      m1 <- zeroinfl(as.formula(paste(dv, '~', ziv, zcontrols)),
               data=data)
      m2 <- zeroinfl(as.formula(paste(dv, '~', ziv, zcontrols)), 
		     data=filter(data, gwno != 475))
      m3 <- zeroinfl(as.formula(paste(dv, '~', ziv, zcontrols)),
               data=filter(data, gwno != 500))
      m4 <-  zeroinfl(as.formula(paste(dv, '~', ziv, zcontrols)),
               data=filter(data, gwno != 501))
      m5 <-  zeroinfl(as.formula(paste(dv, '~', ziv, zcontrols)),
               data=filter(data, gwno != 452))
      m6 <- zeroinfl(as.formula(paste(dv, '~', ziv, '*gbr',
				      zcontrols)), data=data)
      m7 <- zeroinfl(as.formula(paste(dv, '~', ziv, '*region1',
				      zcontrols)), data=data)
      m8 <- zeroinfl(as.formula(paste(dv, '~', ziv, '*region5',
				      zcontrols)), data=data)
      models_out[[i]] <- list(m1, m2, m3, m4, m5, m6, m7, m8)
  }
  #models_out <- models_out
  #final <- list('models'=models_out)
  return(models_out)
}

zmodels <- zinbanalysis(zdvs = zdvs, zivs = zivs, zcontrols = zcontrols,
			       data = filter(prio_grid_isd, popd > 0),
			       test_label = 'ZINB Models')

# Part 2: Figures

figs <- function(models) {
	figsout <- NULL
	for(i in 1:length(models)){
		mainfigs <- NULL
		mainplots <- NULL
		for(j in 1:5) {
			mainfigs[[j]] <- ggpredict(models[[i]][[j]], 
				terms = "sqrtSpAll [0:15 by = .5]") 
			mainplots[[j]] <- ggplot(mainfigs[[j]],
				aes(x^2, predicted)) +
  				geom_ribbon(aes(ymin = conf.low, 
					ymax = conf.high),
	      				fill = lighten("blue")) +
  				geom_line(color = "blue") +
  				labs(title = paste(titles[j]),
					x = 'Precolonial state presence', 
       					y = 'Communial violence events') +
  				goldenScatterCAtheme
		}
		intfigs <- NULL
		intplots <- NULL
		for(k in 1:3) {
			intfigs[[k]] <- ggpredict(models[[i]][[5+k]],
				terms = c("sqrtSpAll [0:15 by = .5]",
					  paste(interactions[k],"[0,1]")))
			intplots[[k]] <- ggplot(intfigs[[k]],
				aes(x^2, predicted, color = group)) +
				geom_ribbon(aes(ymin = conf.low, 
					ymax = conf.high, fill = group, 
					linetype = NA)) + 
				scale_fill_manual(values = pastels) +
				geom_line() + 
				labs(title = paste(titles[5+k]),
					x = 'Precolonial state presence', 
	     				y = 'Communial violence events') +
				goldenScatterCAtheme
		}
		figsout[[i]] <- flatten(list(mainplots, intplots))
	}
    #models_out[[i]] <- models_out_j %>% purrr::flatten()
	return(figsout)
}	

plots <- figs(models = zmodels)

org3plots <- do.call("grid.arrange", c(plots[[1]], ncol = 4))

ggsave("../Output/org3plots.pdf", org3plots, width = 15, height = 15/1.68)

nonstateplots <- do.call("grid.arrange", c(plots[[2]], ncol = 4))

ggsave("../Output/nonstateplots.pdf", nonstateplots, width = 15, height = 15/1.68)

# logOrg3 Plot
logOrg3 <- ggplot() +
	geom_sf(data = filter(prio_grid_isd),
            linetype = 0,
            aes(fill = log(org3)),
            show.legend = FALSE) + 
    labs(title = "Communal violence (log)") +
    scale_fill_viridis_c() +
    theme_minimal()

# }}}

# {{{ SIDE
#==============================================================================#
#	SIDE data					                       #
#==============================================================================#

#==============================================================================#
# Uganda

# # The below command only needs to be run once
# side_download(country = "Uganda", year = 2010, marker = "ethnic", dest.dir =
#	      "../Data/", conv.hull = T)

uga.ethnic <- side_load(country = "Uganda", year = 2010, marker = "ethnic",
			source.dir = "../Data")

uga.ethnic.meta.df <- sidemap2data(uga.ethnic)

names(uga.ethnic) <- uga.ethnic.meta.df$groupname

# Add Egypt?
gisdUga <- filter(geoisd, COWID == 5001 | COWID == 5003 | COWID == 4842 | COWID
		  == 517)

gisdUga <- st_make_valid(gisdUga)

crs(uga.ethnic) <- crs(geoisd)

ugaext <- extent(uga.ethnic)

ugaData <- extract(uga.ethnic, ugaext, df = T)
ugaData <- ugaData %>% mutate(id_cell = seq_len(nrow(.)))

grid <- st_bbox(ugaext) %>% 
  st_make_grid(cellsize = 0.00833334, what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

xy <- xyFromCell(uga.ethnic, as.integer(rownames(ugaData)))
result <- cbind(xy, ugaData)
colnames(result)[1:2] <- c("lon", "lat")

result <- st_as_sf(result, coords = c("lon","lat"))

result <- na.omit(result)

st_crs(result) <- crs(geoisd)

grid <- st_join(grid, result, join = st_contains)

grid <- na.omit(grid)

# Merging SIDE and Geo-ISD
ugaMerged <- grid %>% mutate(sp = lengths(st_within(grid, gisdUga)))

# Creating ethnic fractionalization index
for(i in 1:length(ugaMerged$id_cell.y)) {
		ugaMerged$ef[i] =
		(1 - (
		ugaMerged$acholi[i]^2 +
		ugaMerged$alur.jopahhola[i]^2 +
		ugaMerged$baganda[i]^2 + 
		ugaMerged$bagisu.sabiny[i]^2 +
		ugaMerged$bakiga[i]^2 +
		ugaMerged$banyankore[i]^2 +
		ugaMerged$banyoro[i]^2 +
		ugaMerged$basoga[i]^2 +
		ugaMerged$batoro[i]^2 +
		ugaMerged$iteso[i]^2 +
		ugaMerged$karimojong[i]^2 +
		ugaMerged$langi[i]^2 +
		ugaMerged$lugbara.madi[i]^2 +
		ugaMerged$other[i]^2))
}

# Adding water to UGA

ugaMerged <- st_join(ugaMerged, select(prio_grid_isd, geometry, water_gc), join
		     = st_within)

#ugaMerged <- grid %>% mutate(water = st_within(ugaMerged, prio_grid_isd$water_gc))

#==============================================================================#
# Burkina Faso

# # The below command only needs to be run once
#side_download(country = "Burkina Faso", year = 1999, marker = "ethnic", dest.dir =
      #"../Data/", conv.hull = T)

brk.ethnic <- side_load(country = "Burkina Faso", year = 1999, marker = "ethnic",
			source.dir = "../Data")

brk.ethnic.meta.df <- sidemap2data(brk.ethnic)

names(brk.ethnic) <- brk.ethnic.meta.df$groupname

gisdbrk <- filter(geoisd, COWID == 4327 | COWID == 4763 | COWID == 4321 | COWID ==
	      4751 | COWID == 4751 | COWID == 4325 | COWID == 4392 | COWID ==
	      4326 | COWID == 4395 | COWID == 4521 | COWID == 4328 | COWID ==
	      4395 | COWID == 4393 | COWID == 4399 | COWID == 4321 | COWID ==
	      4396)

gisdbrk <- st_make_valid(gisdbrk)

crs(brk.ethnic) <- crs(geoisd)

brkext <- extent(brk.ethnic)

brkData <- extract(brk.ethnic, brkext, df = T)
brkData <- brkData %>% mutate(id_cell = seq_len(nrow(.)))

grid <- st_bbox(brkext) %>% 
  st_make_grid(cellsize = 0.00833334, what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

xy <- xyFromCell(brk.ethnic, as.integer(rownames(brkData)))
result <- cbind(xy, brkData)
colnames(result)[1:2] <- c("lon", "lat")

result <- st_as_sf(result, coords = c("lon","lat"))

result <- na.omit(result)

st_crs(result) <- crs(geoisd)

grid <- st_join(grid, result, join = st_contains)

grid <- na.omit(grid)

# Merging SIDE and Geo-ISD
grid <- grid %>% mutate(sp = lengths(st_within(grid, gisdbrk)))

# Creating ethnic fractionalization index
for(i in 1:length(grid$id_cell.y)) {
		grid$ef[i] =
		(1 - (
		grid$bissa[i]^2 +
		grid$bobo[i]^2 +
		grid$dafing[i]^2 +
		grid$dagara[i]^2 +
		grid$dioula[i]^2 +
		grid$fulfulde.peul[i]^2 + 
		grid$gourmatche[i]^2 +
		grid$gourounsi[i]^2 +
		grid$lobi[i]^2 +
		grid$mossi[i]^2 +
		grid$other[i]^2 +
		grid$other.burkina[i]^2 +
		grid$samo[i]^2 +
		grid$senoufo[i]^2))
}

#==============================================================================#
# Nigeria

# # The below command only needs to be run once
# side_download(country = "Nigeria", year = 2013, marker = "ethnic", dest.dir =
#	      "../Data/", conv.hull = T)

nig.ethnic <- side_load(country = "Nigeria", year = 2013, marker = "ethnic",
			source.dir = "../Data")

nig.ethnic.meta.df <- sidemap2data(nig.ethnic)

names(nig.ethnic) <- nig.ethnic.meta.df$groupname

gisdnig <- filter(geoisd, COWID == 4798 | COWID == 4752 | COWID == 4521 | COWID
		  == 4327 | COWID == 4831 | COWID == 4776 | COWID == 4763 | COWID ==
434 | COWID == 4362 | COWID == 4768 | COWID == 4769 | COWID == 4751 | COWID == 4771 |
COWID == 4742 | COWID == 4773 | COWID == 4765 | COWID == 4775 | COWID == 4832)

gisdnig <- st_make_valid(gisdnig)

crs(nig.ethnic) <- crs(geoisd)

nigext <- extent(nig.ethnic)

nigData <- extract(nig.ethnic, nigext, df = T)
nigData <- nigData %>% mutate(id_cell = seq_len(nrow(.)))

grid <- st_bbox(nigext) %>% 
  st_make_grid(cellsize = 0.00833334, what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

xy <- xyFromCell(nig.ethnic, as.integer(rownames(nigData)))
result <- cbind(xy, nigData)
colnames(result)[1:2] <- c("lon", "lat")

result <- st_as_sf(result, coords = c("lon","lat"))

result <- na.omit(result)

st_crs(result) <- crs(geoisd)

grid <- st_join(grid, result, join = st_contains)

grid <- na.omit(grid)

# Merging SIDE and Geo-ISD
grid <- grid %>% mutate(sp = lengths(st_within(grid, gisdnig)))

# Creating ethnic fractionalization index
for(i in 1:length(grid$id_cell.y)) {
		grid$ef[i] =
		(1 - (
		grid$ebira.igbira[i]^2 +
		grid$fulfulde[i]^2 +
		grid$ibibio[i]^2 + 
		grid$igbo.ibo[i]^2 +
		grid$kanuri.beriberi[i]^2 +
		grid$ogoni[i]^2 +
		grid$tiv[i]^2 +
		grid$yoruba[i]^2 +
		grid$annang[i]^2 +
		grid$esan[i]^2 +
		grid$gbaju.gbagi[i]^2 +
		grid$idoma[i]^2 +
		grid$ijaw.izon[i]^2 +
		grid$mumuye[i]^2 +
		grid$urhobo[i]^2 +
		grid$bini.edo[i]^2 +
		grid$fulani[i]^2 +
		grid$hausa[i]^2 +
		grid$igala[i]^2 +
		grid$kambari[i]^2 +
		grid$nupe[i]^2 +
		grid$tarok[i]^2 +
		grid$wurkum[i]^2 +
		grid$other[i]^2))
}

#==============================================================================#
# Congo (DRC)

# # The below command only needs to be run once
# side_download(country = "Congo (DRC)", year = 2014, marker = "ethnic", dest.dir =
#	      "../Data/", conv.hull = T)

drc.ethnic <- side_load(country = "Congo (DRC)", year = 2014, marker = "ethnic",
			source.dir = "../Data")

drc.ethnic.meta.df <- sidemap2data(drc.ethnic)

names(drc.ethnic) <- gsub("\\W", "\\1", drc.ethnic.meta.df$groupname)

# Find relevant states
gisddrc <- filter(geoisd, COWID == 4776	| COWID == 4908 | COWID == 4842 | COWID
		  == 4910 | COWID == 4911 | COWID == 4904 | COWID == 4902 |
			  COWID == 5517 | COWID == 4703 | COWID == 4841 |
			  COWID == 4841)

gisddrc <- st_make_valid(gisddrc)

crs(drc.ethnic) <- crs(geoisd)

drcext <- extent(drc.ethnic)

drcData <- extract(drc.ethnic, drcext, df = T)
drcData <- drcData %>% mutate(id_cell = seq_len(nrow(.)))

# SOUTH WEST DRC
xy <- xyFromCell(drc.ethnic, as.integer(rownames(drcData)))
result <- cbind(xy, drcData)
colnames(result)[1:2] <- c("lon", "lat")

sw <- extent(min(result$lon), mean(result$lon), min(result$lat),
	     mean(result$lat))

result <- st_as_sf(result, coords = c("lon","lat"))

result <- na.omit(result)

st_crs(result) <- crs(geoisd)

grid <- st_bbox(sw) %>% 
  st_make_grid(cellsize = (0.00833334), what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

grid <- st_join(grid, result, join = st_contains)
	
grid <- na.omit(grid)

# Creating ethnic fractionalization index
for(i in 1:length(grid$id_cell.y)) {
		grid$ef[i] =
		(1 - (
		grid$bakongonordsud[i]^2 +
		grid$baselekmanetkivu[i]^2 +
		grid$baskasaietkwilukwngo[i]^2 + 
		grid$cuvettecentral[i]^2 +
		grid$kasaikatangatanganika[i]^2 +
		grid$lunda[i]^2 +
		grid$other[i]^2 +
		grid$ubangietitimbiri[i]^2 +
		grid$uelelacalbert[i]^2))
}

# NORTH WEST DRC
xy <- xyFromCell(drc.ethnic, as.integer(rownames(drcData)))
result <- cbind(xy, drcData)
colnames(result)[1:2] <- c("lon", "lat")

nw <- extent(mean(result$lon), max(result$lon), min(result$lat),
	     mean(result$lat))

result <- st_as_sf(result, coords = c("lon","lat"))

result <- na.omit(result)

st_crs(result) <- crs(geoisd)

grid <- st_bbox(nw) %>% 
  st_make_grid(cellsize = (0.00833334), what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

grid <- st_join(grid, result, join = st_contains)
	
grid <- na.omit(grid)

# Creating ethnic fractionalization index
for(i in 1:length(grid$id_cell.y)) {
		grid$ef[i] =
		(1 - (
		grid$bakongonordsud[i]^2 +
		grid$baselekmanetkivu[i]^2 +
		grid$baskasaietkwilukwngo[i]^2 + 
		grid$cuvettecentral[i]^2 +
		grid$kasaikatangatanganika[i]^2 +
		grid$lunda[i]^2 +
		grid$other[i]^2 +
		grid$ubangietitimbiri[i]^2 +
		grid$uelelacalbert[i]^2))
}

# NORTH EAST
xy <- xyFromCell(drc.ethnic, as.integer(rownames(drcData)))
result <- cbind(xy, drcData)
colnames(result)[1:2] <- c("lon", "lat")

se <- extent(mean(result$lon), max(result$lon), mean(result$lat),
	     max(result$lat))

result <- st_as_sf(result, coords = c("lon","lat"))

result <- na.omit(result)

st_crs(result) <- crs(geoisd)

grid <- st_bbox(se) %>% 
  st_make_grid(cellsize = (0.00833334), what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

grid <- st_join(grid, result, join = st_contains)
	
grid <- na.omit(grid)

# Creating ethnic fractionalization index
for(i in 1:length(grid$id_cell.y)) {
		grid$ef[i] =
		(1 - (
		grid$bakongonordsud[i]^2 +
		grid$baselekmanetkivu[i]^2 +
		grid$baskasaietkwilukwngo[i]^2 + 
		grid$cuvettecentral[i]^2 +
		grid$kasaikatangatanganika[i]^2 +
		grid$lunda[i]^2 +
		grid$other[i]^2 +
		grid$ubangietitimbiri[i]^2 +
		grid$uelelacalbert[i]^2))
}

# "SOUTH EAST"
xy <- xyFromCell(drc.ethnic, as.integer(rownames(drcData)))
result <- cbind(xy, drcData)
colnames(result)[1:2] <- c("lon", "lat")

ne <- extent(min(result$lon), mean(result$lon), mean(result$lat),
	     max(result$lat))

result <- st_as_sf(result, coords = c("lon","lat"))

result <- na.omit(result)

st_crs(result) <- crs(geoisd)

grid <- st_bbox(ne) %>% 
  st_make_grid(cellsize = (0.00833334), what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

grid <- st_join(grid, result, join = st_contains)
	
grid <- na.omit(grid)

# Creating ethnic fractionalization index
for(i in 1:length(grid$id_cell.y)) {
		grid$ef[i] =
		(1 - (
		grid$bakongonordsud[i]^2 +
		grid$baselekmanetkivu[i]^2 +
		grid$baskasaietkwilukwngo[i]^2 + 
		grid$cuvettecentral[i]^2 +
		grid$kasaikatangatanganika[i]^2 +
		grid$lunda[i]^2 +
		grid$other[i]^2 +
		grid$ubangietitimbiri[i]^2 +
		grid$uelelacalbert[i]^2))
}


# Merging SIDE and Geo-ISD
drcMerged <- grid %>% mutate(sp = lengths(st_within(grid, gisddrc)))


#==============================================================================#
# Test plot

gridtestwater <- ggplot() +
		   geom_sf(data = na.omit(ugaMerged),
			   linetype = 0,
			   aes(fill = water_gc),
			   show.legend = F) +
    		   scale_fill_viridis_c() +
		   theme_minimal() +
		   theme(plot.background = element_rect(fill = "white")) 


tiff("../Output/brkfasoSP.tiff",
    width = 10, height = 5, res = 300, units = 'in', compression = 'lzw')
gridtest
dev.off()

ugaplots <-  grid.arrange(gridtestef,gridtestsp, ncol = 2)

ggsave("../Output/ugaplots.tiff", ugaplots,
       device = "tiff", width = 10, height = 10/1.68, units = "in", dpi = 300,
       compression = "lzw")

# }}}

# {{{ Footer
#==============================================================================#
# Experimenting with tile package

lp <- lineplot(x = ggorg3zinbCol$x^2, y=ggorg3zinbCol$predicted, lower =
	 ggorg3zinbCol$conf.low, upper =ggorg3zinbCol$conf.high)

tile(lp)

#==============================================================================#
# Experimenting with 3D render

plot_gg(gridtest, multicore = T)


###################################
Org3 <- ggplot() +
	geom_sf(data = filter(prio_grid_isd, gwno != 452),
            linetype = 0,
            aes(fill = logOrg3),
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()

Org3

# Burkina Faso

geoisd <- st_read('../../QGIS/Geo-ISD.shp')

brk <- filter(geoisd, COWID == 4327 | COWID == 4763 | COWID == 4321 | COWID ==
	      4751 | COWID == 4751 | COWID == 4325 | COWID == 4392 | COWID ==
	      4326 | COWID == 4395 | COWID == 4521 | COWID == 4328 | COWID ==
	      4395 | COWID == 4393 | COWID == 4399 | COWID == 4321 | COWID ==
	      4396)

ext <- raster::extent(-5.7,2.5,9,15.5)

grid <- st_bbox(ext) %>% 
  st_make_grid(cellsize = (0.1), what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

brk <- st_make_valid(brk)

grid$sp <- lengths(st_intersects(grid, brk))


borders <- read_sf(dsn = "../Data/Shapes/Africa.shp") %>% 
	filter(ID == 22)

borders <- st_set_crs(borders, 4326)

gridtest <- ggplot() +
		   geom_sf(data = grid,
			   linetype = 0,
			   aes(fill = ef),
			   show.legend = F) +
    		   scale_fill_viridis_c() +
		   geom_sf(data = borders, color = "white", fill = NA) +
		   theme_minimal() +
		   theme(plot.background = element_rect(fill = "white")) 

pdf("../Output/burkinafasoEF.pdf",
    width = 10, height = 10/1.68)
gridtest
dev.off()

# }}}
