library(countrycode)
library(dplyr)
library(ggplot2)
library(rworldmap)
library(rworldxtra)
library(scales)
library(sf)
library(spData)
library(vdemdata)

source('../R/Scripts/goldenScatterCAtheme.r')

load("../R/Data/ucdp-prio-acd-211.rdata")

load("../R/Data/GEDEvent_v21_1.RData")

ged <- filter(GEDEvent_v21_1, type_of_violence != 3) %>% group_by(country_id) %>% 
				   summarise(deaths = sum(best))

worldmap <- getMap(resolution = "high")

worldmap <- st_as_sf(worldmap)

worldmap$country_id <- countrycode(worldmap$ISO_N3,"iso3n","gwn") 

#st_crs(ged) <- st_crs(worldmap)

ged <- left_join(worldmap, ged, by = "country_id")

ggplot() +
	geom_sf(data = vd,
		aes(fill = country_id))

gedplot <- ggplot() +
	geom_sf(data = ged,
		aes(fill = sqrt(deaths))) +
	scale_fill_viridis_c("Fatalities",
    			 labels = trans_format("identity", function(x) x^2)) +
	theme_minimal() + 
	theme(text=element_text(size=16,
        family="Times New Roman")) +
	goldenScatterCAtheme
# Preview plot
gedplot

ggsave("img/gedplot.pdf", gedplot, device = cairo_pdf, width = 10, height = 10/1.68)

vd <- select(vdem, v2x_polyarchy, country_name, e_gdppc)

vd$country_id <- countrycode(vd$country_name,"country.name","gwn") 

vd <- left_join(worldmap, vd, by = "country_id")

vdemplot <- ggplot() +
	geom_sf(data = vd,
		aes(fill = v2x_polyarchy)) +
	scale_fill_viridis_c("Electoral Democracy Index") +
	theme_minimal() + 
	theme(text=element_text(size=16,
        family="Times New Roman")) +
	goldenScatterCAtheme

# Preview plot
vdemplot

ggsave("img/vdemplot2.pdf", vdemplot, device = cairo_pdf, width = 10, height = 10/1.68)

gdpplot <- ggplot() +
	geom_sf(data = ged,
		aes(fill = sqrt(gdpPercap))) +
	scale_fill_viridis_c("GDP per capita",
		labels = trans_format("identity", function(x) x^2)) +
	theme_minimal() + 
	theme(text=element_text(size=16,
        family="Times New Roman")) +
	goldenScatterCAtheme
# Preview plot
gdpplot

ggsave("img/gdpplot.pdf", gdpplot, device = cairo_pdf, width = 10, height = 10/1.68)

