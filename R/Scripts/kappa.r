# Loading packages
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(readxl)
library(sf)
library(tidyr)
library(viridis)

source("goldenScatterCAtheme.r")

#==============================================================================#
#	States per year							       #
#==============================================================================#
# Loading Data
isd <- read_excel("../Data/isd_2_1_1_7_2021.xlsx")

# Changing dates to years
isd$syear <- year(as.Date(isd$Start, '%d-%m-%Y'))
isd$eyear <- year(as.Date(isd$End, '%d-%m-%Y'))

# Selecting what is needed
isd <- select(isd, syear, eyear, Region, StateName)

# Omitting NAs
isd <- na.omit(isd)

# Changing to country-year format
isd <- isd %>% unnest(year = map2(syear, eyear, seq))

# Summing number of states per region-year 
isd <- isd %>% group_by(Region, year) %>% 
	summarise(Region = as.factor(Region),
		  total = length(unique(StateName)),
		  stateName = StateName) 

# Assigning names to regions
regions <- c("Africa", "Central & South Asia", "East Asia & Oceania",
	     "North & South America", "Europe", "Middle East")

levels(isd$Region) <- regions

# Plotting
plotty <- ggplot(isd, aes(x = year, y = total, colour = Region)) + 
	geom_line(size = 2) +
	scale_color_viridis_d() +
	ylab("Number of states") + xlab("Year") +
	goldenScatterCAtheme
# Preview plot
plotty

# Writing plot to file
pdf(file = "../Output/statesPerYear.pdf", 
        width=10, height=10/1.683)
plotty
dev.off()

#==============================================================================#
#	Civli conflict in Africa 					       #
#==============================================================================# 

# Loading data
load("../Data/GeoISDControls.Rdata")

# Plotting
ccplot <- ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes(fill = log(interdeaths+1)),
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()
# Preview plott
ccplot 

# Writing plot to file
pdf("../Output/civilconflictAfrica.pdf",
	width = 10, height = 10/1.68)
	ccplot
dev.off()

#==============================================================================#
#	Civil conflicts per year					       #
#==============================================================================#

# Loading data
#load("../Data/GEDEvent_v21_1.RData")
load("../Data/ucdp-prio-acd-211.rdata")

# Getting state based violence type from UCDP/PRIO
gedsb <- UcdpPrioConflict_v21_1 %>%  dplyr::select(conflict_id, type_of_conflict) %>%  
	mutate(conflict_id = as.integer(conflict_id)) %>% 
	mutate(sb = as.integer(type_of_conflict))

# Merging
ged <- left_join(GEDEvent_v21_1, gedsb, by = c("conflict_new_id" = 
					       "conflict_id"))

# Filtering out non relevant conflict and summing for region-year
ged <- filter(ged, sb == 3 | sb == 4) %>% 
	group_by(region, year) %>% 
	summarise(deaths = sum(best))

# Using UCDP/PRIO instead
# Converting form characters
UcdpPrioConflict_v21_1$type_of_conflict <- as.integer(UcdpPrioConflict_v21_1$type_of_conflict)
UcdpPrioConflict_v21_1$year <- as.integer(UcdpPrioConflict_v21_1$year)
UcdpPrioConflict_v21_1$region <- as.factor(as.integer(UcdpPrioConflict_v21_1$region))

# Assigning names to regions
ucdpregions <- c("Europe", "Middle East", "Asia", "Africa", "Americas")
levels(UcdpPrioConflict_v21_1$region) <- ucdpregions

# Summarising conflicts per year
ucdp <- filter(UcdpPrioConflict_v21_1, type_of_conflict == 3 | type_of_conflict == 4) %>% 
	group_by(region, year) %>% 
	summarise(ccount = length(type_of_conflict))

# Plotting
cclplot <- ggplot(ucdp, aes(x = year, y = ccount, color = region)) +
	geom_line(size = 2) +
	scale_color_viridis_d("Regoin") +
	ylab("Number of civil conflicts") + xlab("Year") +
	goldenScatterCAtheme
# Preview plot
cclplot

# Writing plot to file
pdf("../Output/civilConflictRegion.pdf",
	width = 10, height = 10/1.68)
	cclplot
dev.off()
