# Loading packages
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(readxl)
library(scales)
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
            show.legend = T) + 
    scale_fill_viridis_c("Fatalities",
			 labels=trans_format("identity", function(x)
					     round(exp(x),0)-1)) +
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
load("../Data/ucdp-prio-acd-211.rdata")

# Using UCDP/PRIO instead
# Converting from characters
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
	scale_color_viridis_d("Region") +
	ylab("Number of civil conflicts") + xlab("Year") +
	goldenScatterCAtheme
# Preview plot
cclplot

# Writing plot to file
pdf("../Output/civilConflictRegion.pdf",
	width = 10, height = 10/1.68)
	cclplot
dev.off()

#==============================================================================#
#	Communal violence fatalities per year				       #
#==============================================================================#

# Loading data
load("../Data/GEDEvent_v21_1.RData")
load("../Data/ucdp-nonstate-211.rdata")

# Getting non-state based violence
gedns <- Nonstate_v21_1 %>%  dplyr::select(conflict_id, org, year) %>%  
	mutate(conflict_id = as.integer(conflict_id)) %>% 
	mutate(org = as.integer(org))

# Merging
ged <- left_join(GEDEvent_v21_1, gedns, by = c("conflict_new_id" = 
					       "conflict_id"))

# Filtering out non relevant conflict and summing for region-year
ged <- filter(ged, org == 3) %>% 
	group_by(region, year.x) %>% 
	summarise(deaths = sum(best))


# Assigning names to regions
levels(ged$region) <- ucdpregions


# Plotting
cvplot <- ggplot(ged, aes(x = year.x, y = deaths, color = region)) +
	geom_line(size = 2) +
	scale_color_viridis_d("Region") +
	ylab("Fatalities") + xlab("Year") +
	goldenScatterCAtheme
# Preview plot
cvplot

# Writing plot to file
pdf("../Output/communalViolenceRegion.pdf",
	width = 10, height = 10/1.68)
	cvplot
dev.off()

#==============================================================================#
#	Borderness and frontierness	 				       #
#==============================================================================#

frontierplot <- ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes(fill = sqrt(sp_b_i_sum))) +
    scale_fill_viridis_c("Number of borders",
	labels=trans_format("identity", function(x) round(x^2,0))) +
    theme_minimal()
frontierplot

pdf("../Output/frontierplot.pdf",
    width = 10, height = 10/1.68)
	frontierplot
dev.off()

borderplot <- ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes(fill = sqrt(sp_o))) +
    scale_fill_viridis_c("Number of
overlapping borders",
	labels=trans_format("identity", function(x) round(x^2,0))) +
    theme_minimal()
borderplot

pdf("../Output/borderplot.pdf",
    width = 10, height = 10/1.68)
	borderplot
dev.off()

#==============================================================================#
#	Per state- state presence					       #
#==============================================================================#

gisd <- read.csv('../../QGIS/Geo-ISD.csv')

# SP Charles version

gisd <- gisd %>% 
  mutate(state_presence = 1) %>%
  group_by(gid, COWID, year) %>%
  summarise(state_presence_max = max(na.omit(state_presence)),
            no_unique_states = length(unique(COWID)),
            state_presence_all = sum(na.omit(state_presence))
  )

gisd <- gisd %>% 
  group_by(gid, COWID) %>%
  # Summing years with any state presence
  summarise(sp_sum_any_cb = sum(na.omit(state_presence_max)),
            # Summing state presence from all maps over years
            sp_sum_cb = sum(na.omit(state_presence_all)))

prio_grid_isd <- left_join(prio_grid_isd, gisd, by = c("gid")) %>% 
	mutate(state_presence = 1)

prio_grid_isd <- prio_grid_isd %>% group_by(COWID) %>% mutate(sp_sum_adj=sp_sum_cb/max(sp_sum_cb))

prio_grid_isd$sp_sum_adj[is.na(prio_grid_isd$sp_sum_adj) == T] <- 0

gisd_all <- ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes(fill = sp_sum_adj, alpha = sp_sum_adj)) +
    scale_fill_viridis_c("State presence
per pre-colonial state") +
    scale_alpha_continuous(guide = 'none') +
    theme_minimal()

pdf("../../Thesis/img/geo_isd_all.pdf",
    width = 10, height = 10/1.68)
	gisd_all
dev.off()

