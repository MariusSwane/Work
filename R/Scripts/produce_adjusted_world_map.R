


############################################
######  Setting Working Directory  #########
############################################


library(rstudioapi)
setwd("C:/Users/charlesb/ARC_Project Dropbox/Charles Butcher/GEO-ISD/R")
setwd("~/ARC_Project Dropbox/Charles Butcher/GEO-ISD/R")

############################################
######  Loading packages  ##################
############################################

#library(raster)
library(spData)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(purrr)
library(RColorBrewer)

options(rgl.useNULL = FALSE)

############################################
######  Loading data  ######################
############################################

isd <- read.csv2('../ISDV2_Africa.csv')

geoisd <- st_read('../QGIS/Geo-ISD.shp')

geoisd_data <- read_csv('../QGIS/Geo-ISD.csv',
                        cols(
                          gid = col_double(),
                          xcoord = col_double(),
                          ycoord = col_double(),
                          col = col_double(),
                          row = col_double(),
                          ISDID = col_character(),
                          COWID = col_double(),
                          name = col_character(),
                          year = col_double(),
                          lyear = col_integer(),
                          hyear = col_integer(),
                          source = col_character(),
                          coder = col_character(),
                          note = col_character(),
                          error = col_double(),
                          layer = col_character(),
                          path = col_character()
                        ),
                        col_names = T,
)

prio_grid_shp <- st_read('../PRIO-Grid/priogrid_cell.shp')

prio_grid <- read_csv('../PRIO-Grid/prio_grid_2014.csv') %>% 
  filter( (gwno >= 404 & gwno <= 626) | gwno == 651)

prio_grid_merge <- left_join(prio_grid_shp, prio_grid, by = c("gid")) %>% 
  filter( (gwno >= 404 & gwno <= 626) | gwno == 651) %>% 
  dplyr::select(gid, xcoord, ycoord, col, row, year, gwno, geometry)


# State presence - Charles version
geoisd_data_gid_year <- geoisd_data %>% 
  mutate(state_presence = 1) %>%
  group_by(gid, COWID, year) %>%
  summarise(state_presence_max = max(na.omit(state_presence)),
            no_unique_states = length(unique(COWID)),
            state_presence_all = sum(na.omit(state_presence))
  )

geoisd_data_gid <- geoisd_data_gid_year %>% 
  group_by(gid, COWID) %>%
  # Summing years with any state presence
  summarise(sp_sum_any = sum(na.omit(state_presence_max)),
            # Summing state presence from all maps over years
            sp_sum = sum(na.omit(state_presence_all)))


prio_grid_isd <- left_join(prio_grid_merge, geoisd_data_gid, by = c("gid")) %>%
  mutate(state_presence=1)


### Produce a blank world  map ####

map_shape <- map_data(database = "world", regions = ".") 


ggplot() +
  geom_polygon(data = map_shape, aes(x=long, y = lat, group = group), colour="black", fill="white") +
  theme_minimal()


# Adjust the state presence measure so it is a proportion of the maximum level of state presence for that state. 

prio_grid_isd_t <- filter(prio_grid_isd) %>% group_by(COWID) %>% mutate(sp_sum_adj=sp_sum/max(sp_sum))


# World


world_plot <- ggplot() +
  geom_polygon(data = map_shape, aes(x=long, y = lat, group = group), colour="black", fill="white") +
  geom_sf(data = filter(prio_grid_isd_t, sp_sum>0),
          linetype = 0,
          aes(fill = sp_sum_adj, alpha=sp_sum_adj),
          show.legend = FALSE) +
  scale_fill_viridis_c(option="magma") +
  # geom_contour(data = filter(prio_grid_isd, sp_sum>0), aes(x=xcoord, y=ycoord, z=sp_sum), colour="black") +
  #  geom_point(data=cities_plot, aes(y=lat, x=lon),color="black", size=1) + geom_text(data=cities_plot, aes(y=lat, x=lon, label=city),hjust=0, vjust=2) +
  theme_minimal()


