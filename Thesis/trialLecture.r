library(countrycode)
library(dplyr)
library(ggplot2)
library(scales)
library(sf)
library(spData)
library(vdemdata)

source('../R/Scripts/goldenScatterCAtheme.r')

vd <- vdem

load("../R/Data/ucdp-prio-acd-211.rdata")

load("../R/Data/GEDEvent_v21_1.RData")

ged <- GEDEvent_v21_1 %>% group_by(country_id) %>% 
				   summarise(deaths = sum(best))

data(world)

world$country_id <- countrycode(world$iso_a2,"iso2c","gwn") 

ged <- left_join(world, ged, by = "country_id")

gedplot <- ggplot() +
	geom_sf(data = ged,
		aes(fill = sqrt(deaths))) +
	scale_fill_viridis_c("Fatalities",
    			 labels = trans_format("identity", function(x) x^2)) +
	theme_minimal() + 
	goldenScatterCAtheme
# Preview plot
gedplot

gdpplot <- ggplot() +
	geom_sf(data = ged,
		aes(fill = sqrt(gdpPercap))) +
	scale_fill_viridis_c("GDP per capita",
    			 labels = trans_format("identity", function(x) x^2)) +
	theme_minimal() + 
	goldenScatterCAtheme
# Preview plot
gdpplot

vdem <- select(vdem, v2x_polyarchy, country_name)

vdem$country_id <- countrycode(vdem$country_name,"country.name","gwn") 

vdem <- left_join(world, vdem, by = "country_id")

vdemplot <- ggplot() +
	geom_sf(data = vdem,
		aes(fill = v2x_polyarchy)) +
	scale_fill_viridis_c("Electoral Democracy Index")
	theme_minimal() + 
	goldenScatterCAtheme
# Preview plot
vdemplot
