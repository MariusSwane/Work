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

#==============================================================================#
#	Loading Packages						       #
#==============================================================================#	

library(conflicted)
library(dplyr)
library(ggeffects)
library(ggplot2)
library(margins)
library(MASS)
library(pscl)
library(raster)
library(sidedata)
library(spdep)
library(sf)
library(stargazer)
library(texreg)

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
				region4 = as.numeric(factor(region)==4))

#==============================================================================#
#	Functions				    	             	       #
#==============================================================================#

#==============================================================================#
# Negative binomial models function

GeoISDanalysis <- function(dvs, ivs, controls, data, test_label){
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
    models_out[[i]] <- models_out_j %>% purrr::flatten()
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
	      	'distcoast' = 'Distance to coast',
		'region3' = 'North Africa',
	      	'logPopd' = 'Population density (log)', 
		'bdist3' = 'Distance to border', 
		'nopastorTRUE' = 'Land not suited for pastorial herding',
		'sqrtSpAll' = 'Precolonial state presence (sqrt)')


control_names_cv <-c('Baseline', 'North Africa', 'Population density', 'Distance
		     to international boundary', 'Pastoralism')

cv_dvs <- c('non_state', 'org3')

cv_captions <- c('Non-state conflict events', 'Communal violence events')

ivs <- c('sqrtSpAll')

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

#==============================================================================#
#	Analysis							       #
#==============================================================================#

cv_nb_models <- GeoISDanalysis(dvs = cv_dvs, ivs = ivs, controls = controls,
			       data = filter(prio_grid_isd, popd > 0),
			       test_label = 'Linear Models')

#==============================================================================#
#	Marginal interacation plots 					       #
#==============================================================================#

# Plotting
org3_NB <- glm.nb(org3 ~ sqrtSpAll + mountains_mean + water_gc + barren_gc +
		  logCDist + logBDist + logPopd + region3, data =
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
         table = T)
}

#==============================================================================#
#	ZINB  								       #
#==============================================================================#		

org3_zinb <- zeroinfl(org3 ~ sqrtSpAll + mountains_mean + water_gc + barren_gc +
		  logCDist + logBDist + logPopd + nopastor, data =
		  filter(prio_grid_isd, popd > 0)) 
summary(org3_zinb)

ggorg3zinb <- ggpredict(org3_zinb, terms = "sqrtSpAll [0:15 by = .5] ")

org3zinbplot <- ggplot(ggorg3zinb, aes(x^2, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
	      fill = lighten("blue")) +
  geom_line(color = "blue") +
  labs(x = 'Precolonial state presence', y = 'Communial violence events') +
  goldenScatterCAtheme

# Printing to file
pdf("../Output/CommunalViolenceZinbMargins.pdf",
    width = 10, height = 10/1.68)
org3zinbplot
dev.off()

# Controlling for French and British colonies
org3_zinbCol <- zeroinfl(org3 ~ sqrtSpAll + mountains_mean + water_gc + barren_gc +
		  logCDist + logBDist + logPopd + nopastor + fra + gbr, data =
		  filter(prio_grid_isd, popd > 0)) 
summary(org3_zinbCol)

ggorg3zinbCol <- ggpredict(org3_zinbCol, terms = "sqrtSpAll [0:15 by = .5] ")

org3zinbplotCol <- ggplot(ggorg3zinbCol, aes(x^2, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
	      fill = lighten("blue")) +
  geom_line(color = "blue") +
  labs(x = 'Precolonial state presence', y = 'Communial violence events') +
  goldenScatterCAtheme

org3zinbplotCol 
	
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

geoisd <- st_read('../../QGIS/Geo-ISD.shp')

# Add Egypt?
gisdUga <- filter(geoisd, COWID == 5001 | COWID == 5003 | COWID == 4842 | COWID
		  == 517)

gisdUga <- st_make_valid(gisdUga)

crs(uga.ethnic) <- crs(geoisd)

ext <- extent(uga.ethnic)

ugaData <- raster::extract(uga.ethnic, ext, df = T)
ugaData <- ugaData %>% mutate(id_cell = seq_len(nrow(.)))

grid <- st_bbox(ext) %>% 
  st_make_grid(cellsize = 0.00833334, what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

xy <- xyFromCell(uga.ethnic, as.integer(rownames(ugaData)))
result <- cbind(xy, ugaData)
colnames(result)[1:2] <- c("lon", "lat")

result <- st_as_sf(result, coords = c("lon","lat"))

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

# Creating population per grid cell
for(i in 1:length(ugaMerged$id_cell.y)) {
		ugaMerged$pop[i] =
		(ugaMerged$acholi[i] +
		ugaMerged$alur.jopahhola[i] +
		ugaMerged$baganda[i] + 
		ugaMerged$bagisu.sabiny[i] +
		ugaMerged$bakiga[i] +
		ugaMerged$banyankore[i] +
		ugaMerged$banyoro[i] +
		ugaMerged$basoga[i] +
		ugaMerged$batoro[i] +
		ugaMerged$iteso[i] +
		ugaMerged$karimojong[i] +
		ugaMerged$langi[i] +
		ugaMerged$lugbara.madi[i] +
		ugaMerged$other[i])
}

#==============================================================================#
# Congo (DRC)

# # The below command only needs to be run once
# side_download(country = "Congo (DRC)", year = 2014, marker = "ethnic", dest.dir =
#	      "../Data/", conv.hull = T)

drc.ethnic <- side_load(country = "Congo (DRC)", year = 2014, marker = "ethnic",
			source.dir = "../Data")

drc.ethnic.meta.df <- sidemap2data(drc.ethnic)

names(drc.ethnic) <- drc.ethnic.meta.df$groupname

# Find relevant states
gisddrc <- filter(geoisd, COWID == 4776	| COWID == 4908 | COWID == 4842 | COWID
		  == 4910 | COWID == 4911 | COWID == 4904 | COWID == 4902 |
			  COWID == 5517 | COWID == 4703 | COWID == 4841 |
			  COWID == 4841)

gisddrc <- st_make_valid(gisddrc)

crs(drc.ethnic) <- crs(geoisd)

ext <- extent(drc.ethnic)

drcData <- raster::extract(drc.ethnic, ext, df = T)
drcData <- drcData %>% mutate(id_cell = seq_len(nrow(.)))

grid <- st_bbox(ext) %>% 
  st_make_grid(cellsize = 0.00833334, what = "polygons") %>%
  st_set_crs(4326)
grid <- grid %>% st_sf() %>% mutate(id_cell = seq_len(nrow(.)))

xy <- xyFromCell(drc.ethnic, as.integer(rownames(drcData)))
result <- cbind(xy, drcData)
colnames(result)[1:2] <- c("lon", "lat")

result <- st_as_sf(result, coords = c("lon","lat"))

st_crs(result) <- crs(geoisd)

grid <- st_join(grid, result, join = st_contains)

grid <- na.omit(grid)

# Merging SIDE and Geo-ISD
drcMerged <- grid %>% mutate(sp = lengths(st_within(grid, gisddrc)))


#==============================================================================#
# Test plot

gridtest <- ggplot() +
		   geom_sf(data = prio_grid_isd,
			   linetype = 0,
			   aes(fill = pastor),
			   show.legend = F) +
    		   #scale_fill_viridis_c() +
		   theme_minimal()

gridtest
