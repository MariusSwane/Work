# Importing packages

library(tidyverse)
library(sf)
library(stargazer)
library(gridExtra)

ch <- st_read("christoffer/OneDrive_1_1-29-2024/icr_cwa.shp")

cr <- st_read("cristina/OneDrive_3_1-29-2024/ICRtest_cm.shp")

jh <- st_read("jan_helge/OneDrive_11_1-29-2024/icr_jhr.shp")

jh$cityerror[jh$cityerror == 0] <- NA
jh$coasterror[jh$coasterror == 0] <- NA

msw <- st_read("marius/MSW-ICR.shp")

oda <- st_read("oda/OneDrive_6_1-29-2024/icr_obb.shp")

stine <- st_read("stine/OneDrive_8_1-29-2024/icr.sr.shp")

icr <- list(ch,cr,jh,msw,oda,stine)

for(i in 1:length(icr)){
	icr[[i]] <- st_drop_geometry(icr[[i]])
}

icr <- Reduce(function(x, y) merge(x, y, all = T), icr)

icr <- icr %>% group_by(coder) %>% 
		summarise(
		 coastNAs = sum(is.na(coasterror)),
		 coasterror = mean(na.omit(coasterror)),
		 cityNAs = sum(is.na(cityerror)),
		 cityerror = mean(na.omit(cityerror)),
		 polygons = sum(!is.na(COWID)))
icr

gcpp <- read.csv("points.csv")

gcpp <- gcpp %>% group_by(coder) %>% 
		summarise(gcpp = mean(points))

icr <- merge(icr, gcpp)

icr <- icr %>% mutate(rate = (coasterror+cityerror)*gcpp,
		mean_error = (coasterror + cityerror)/2)

stargazer(select(icr, coder, mean_error, rate))

png("icr.png") 
grid.table(icr)
dev.off()

png("jh-city.png") 
hist(jh$cityerror)
dev.off()

png("jh-coast.png") 
hist(jh$coasterror)
dev.off()


