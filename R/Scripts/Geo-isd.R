#==============================================================================#
#	 ______                 _________ ____   			       #
#	/ ____/__  ____        /  _/ ___// __ \  			       #
#      / / __/ _ \/ __ \______ / / \__ \/ / / /  			       #
#     / /_/ /  __/ /_/ /_____// / ___/ / /_/ /   			       #
#     \____/\___/\____/     /___//____/_____/    			       #
#                                                                              #
#==============================================================================#

rm(list = ls())

#==============================================================================#
#	Setting Working Directory  					       #
#==============================================================================#

# Only works in Windows
# Otherwise, consider using the RSTUDIO menu:
# Session -> Set Working Directory -> To Source File Location

# On linux it just works :)

#library(rstudioapi)
#setwd(dirname(getActiveDocumentContext()$path))

#==============================================================================#
#	Loading packages 						       # 
#==============================================================================#

library(conflicted)
library(dplyr)
library(ggplot2)
library(priogrid)
library(purrr)
library(raster)
library(RColorBrewer)
library(readr)
library(sf)
library(sidedata) # Some conflicts with raster (tail, stack, unstack, head)
library(sp)
library(spData)
library(tidyr)

#==============================================================================#
#	Setting preferences                                                    #
#==============================================================================#

conflict_prefer("filter", "dplyr")
conflict_prefer("last", "dplyr")  

#==============================================================================#
#	Loading data  							       #
#==============================================================================#

geoisd <- st_read('../../QGIS/Geo-ISD.shp')

geoisd_data <- read_csv('../../QGIS/Geo-ISD.csv',
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

geoisd_borders_data <- read_csv('../../QGIS/Geo-ISD_Borders2.csv',
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

# Cleaning
rm(prio_grid_static, prio_grid_shp)

# Interpolated years from atlas maps
geoisd_interpol <- geoisd_data %>% filter(is.na(year) == T & is.na(lyear) == F) %>%
  mutate(row_id = row.names(.),
         lyear_new = lyear, hyear_new=hyear)

geoisd_interpol <- geoisd_interpol %>% group_by(row_id) %>%
  nest(lyear_new, hyear_new) %>%
  mutate(data = map(data, ~seq(.x$lyear_new, .x$hyear_new, by = 1))) %>%
  unnest(data)

geoisd_interpol <- left_join(geoisd_data, geoisd_interpol)

geoisd_interpol <- geoisd_interpol %>% 
  mutate(year=ifelse(is.na(year) == T, data, year))


# Interpolating years for borders data

borders_interpol <- geoisd_borders_data %>% filter(is.na(year) == T & is.na(lyear) == F) %>%
  mutate(row_id = row.names(.),
         lyear_new = lyear, hyear_new=hyear)

borders_interpol <- borders_interpol %>% group_by(row_id) %>%
  nest(lyear_new, hyear_new) %>%
  mutate(data = map(data, ~seq(.x$lyear_new, .x$hyear_new, by = 1))) %>%
  unnest(data)

borders_interpol <- left_join(geoisd_borders_data, borders_interpol)

borders_interpol <- borders_interpol %>% 
  mutate(year=ifelse(is.na(year) == T, data, year))

#==============================================================================#
#  Checking for naming errors  						       #
#==============================================================================#

#isd <- read.csv2('../Data/ISDV2_Africa.csv')

#isd$isdcow <- paste(isd$COWNum, isd$COWID)
#
#geoisd$isdcow <- paste(geoisd$COWID, geoisd$ISDID)
#
#errorcheck <- left_join(geoisd, isd, by = c("isdcow"="isdcow"))
#errors <- filter(errorcheck, is.na(COWID.y))
                    

#==============================================================================#
#  Summarising the Data  						       #
#==============================================================================#

# Sate presence 
geoisd_data_gid_year <- geoisd_data %>% 
  mutate(state_presence = 1) %>%
  group_by(gid, year) %>%
  summarise(state_presence_max = max(na.omit(state_presence)),
            no_unique_states = length(unique(COWID)),
            state_presence_all = sum(na.omit(state_presence))
  )

geoisd_data_gid <- geoisd_data_gid_year %>% 
 group_by(gid) %>%
            # Summing years with any state presence
  summarise(sp_sum_any = sum(na.omit(state_presence_max)),
            # Summing state presence from all maps over years
            sp_sum = sum(na.omit(state_presence_all)))

# State presence with interpolated years from atlas maps

geoisd_data_gid_int <- geoisd_interpol %>% 
  mutate(state_presence = 1) %>%
  group_by(gid, year) %>%
  summarise(state_presence_max = max(na.omit(state_presence)),
            no_unique_states = length(unique(COWID)),
            state_presence_all = sum(na.omit(state_presence))
  )

geoisd_data_gid_int <- geoisd_data_gid_int %>% 
  group_by(gid) %>%
            # Summing years with any state presence
  summarise(sp_i_sum_any = sum(na.omit(state_presence_max)),
            # Summing state presence from all maps over years
            sp_i_sum = sum(na.omit(state_presence_all)))

# Borders sum 

borders <- geoisd_borders_data %>% 
  mutate(state_presence = 1) %>%
  group_by(gid, year) %>%
  summarise(state_presence_max = max(na.omit(state_presence)),
            no_unique_states = length(unique(COWID)),
            state_presence_all = sum(na.omit(state_presence))
  )

borders_summed <- borders %>% 
  group_by(gid) %>%
  # Summing years with any borders
  summarise(sp_b_sum_any = sum(na.omit(state_presence_max)),
            # Summing borders from all maps over years
            sp_b_sum = sum(na.omit(state_presence_all)))

# Borders sum of interpolated years from atlas maps

borders_int <- borders_interpol %>% 
  mutate(state_presence = 1) %>%
  group_by(gid, year) %>%
  summarise(state_presence_max = max(na.omit(state_presence)),
            no_unique_states = length(unique(COWID)),
            state_presence_all = sum(na.omit(state_presence))
  )

borders_int_summed <- borders_int %>% 
  group_by(gid) %>%
            # Summing years with any borders
  summarise(sp_b_i_sum_any = sum(na.omit(state_presence_max)),
            # Summing borders from all maps over years
            sp_b_i_sum = sum(na.omit(state_presence_all)))

# Overlapping sovereignty

overlap <- geoisd_data_gid_year %>% 
  group_by(gid, year) %>%
  summarise(overlap = max(na.omit(no_unique_states)))-1

overlap <- overlap %>% 
  group_by(gid) %>%
  summarise(sp_o = sum(na.omit(overlap)))

# Overlapping sovereignty (interpolated)

overlap_int <- geoisd_data_gid_int %>% 
  group_by(gid, year) %>%
  summarise(overlap_int = max(na.omit(no_unique_states)))-1

overlap_int <- overlap_int %>% 
  group_by(gid) %>%
  summarise(sp_o_i = sum(na.omit(overlap_int)))


# Max presence by any ONE state years

geoisd_one_state <- geoisd_data %>% 
  mutate(state_presence = 1) %>%
  group_by(gid, COWID, year) %>% 
  summarise(state_presence_max = max(na.omit(state_presence)),
            no_unique_states = length(unique(COWID)),
            state_presence_all = sum(na.omit(state_presence))) %>% 
  group_by(gid, COWID) %>%
            # Summing years with any state presence
  summarise(state_presence_sum_any = sum(na.omit(state_presence_max)),
            # Summing state presence from all maps over years
            state_presence_sum = sum(na.omit(state_presence_all))) %>% 
  group_by(gid) %>% 
            # Keeping only the states with the most presence
  summarise(sp_os_sum_any = max(na.omit(state_presence_sum_any)),
            sp_os_sum = max(na.omit(state_presence_sum)))

# Max presence by any ONE state with interpolated years

geoisd_one_state_int <- geoisd_interpol %>% 
  group_by(gid, COWID) %>%
  mutate(state_presence = 1) %>%
  group_by(gid, COWID, year) %>% 
  summarise(state_presence_max = max(na.omit(state_presence)),
            no_unique_states = length(unique(COWID)),
            state_presence_all = sum(na.omit(state_presence))) %>% 
  group_by(gid, COWID) %>%
  # Summing years with any state presence
  summarise(state_presence_sum_any = sum(na.omit(state_presence_max)),
            # Summing state presence from all maps over years
            state_presence_sum = sum(na.omit(state_presence_all))) %>% 
  group_by(gid) %>% 
  # Keeping only the states with the most presence
  summarise(sp_os_i_sum_any = max(na.omit(state_presence_sum_any)),
            sp_os_i_sum = max(na.omit(state_presence_sum)))

#==============================================================================#
#	Merge with PRIOGRID  						       #
#==============================================================================#

prio_grid_isd <- left_join(prio_grid, geoisd_data_gid, by = c("gid")) %>%
  mutate(sp_sum=ifelse(is.na(sp_sum) == T, 0, sp_sum),
         sp_sum_any=ifelse(is.na(sp_sum_any) == T, 0, sp_sum_any))

prio_grid_isd <- left_join(prio_grid_isd, geoisd_data_gid_int, by = c("gid")) %>%
  mutate(sp_i_sum=ifelse(is.na(sp_i_sum) == T, 0, sp_i_sum),
         sp_i_sum_any=ifelse(is.na(sp_i_sum_any) == T, 0, sp_i_sum_any))

prio_grid_isd <- left_join(prio_grid_isd, borders_summed, by = c("gid")) %>%
  mutate(sp_b_sum=ifelse(is.na(sp_b_sum) == T, 0, sp_b_sum),
         sp_b_sum_any=ifelse(is.na(sp_b_sum_any) == T, 0, sp_b_sum_any))

prio_grid_isd <- left_join(prio_grid_isd, borders_int_summed, by = c("gid")) %>%
  mutate(sp_b_i_sum=ifelse(is.na(sp_b_i_sum) == T, 0, sp_b_i_sum),
         sp_b_i_sum_any=ifelse(is.na(sp_b_i_sum_any) == T, 0, sp_b_i_sum_any))

prio_grid_isd <- left_join(prio_grid_isd, overlap, by = c("gid")) %>%
  mutate(sp_o=ifelse(is.na(sp_o) == T, 0, sp_o))

prio_grid_isd <- left_join(prio_grid_isd, overlap_int, by = c("gid")) %>%
  mutate(sp_o_i=ifelse(is.na(sp_o_i) == T, 0, sp_o_i))

prio_grid_isd <- left_join(prio_grid_isd, geoisd_one_state, by = c("gid")) %>%
  mutate(sp_os_sum=ifelse(is.na(sp_os_sum) == T, 0, sp_os_sum),
         sp_os_sum_any=ifelse(is.na(sp_os_sum_any) == T, 0, sp_os_sum_any))

prio_grid_isd <- left_join(prio_grid_isd, geoisd_one_state_int, by = c("gid")) %>%
  mutate(sp_os_i_sum=ifelse(is.na(sp_os_i_sum) == T, 0, sp_os_i_sum),
         sp_os_i_sum_any=ifelse(is.na(sp_os_i_sum_any) == T, 0, sp_os_i_sum_any))

#==============================================================================#
#	Cleaning							       #
#==============================================================================#		

rm(geoisd, geoisd_data, geoisd_interpol, prio_grid, geoisd_one_state_int,
   geoisd_one_state, overlap_int, overlap, borders_int_summed, borders_summed,
   geoisd_data_gid_int, geoisd_data_gid)

#==============================================================================#
#	Plotting state presence  					       #
#==============================================================================#

# Creating variables vector

vars <- c("sp_sum_any", "sp_sum", "sp_i_sum_any", "sp_i_sum", "sp_b_sum_any",
          "sp_b_sum", "sp_b_i_sum_any", "sp_b_i_sum", "sp_o",
          "sp_o_i", "sp_os_sum_any", "sp_os_sum", "sp_os_i_sum_any", 
          "sp_os_i_sum")

# Creating plots

plot_list = list()
for(i in 1:length(vars)) {
  p = ggplot() +
     geom_sf(data = prio_grid_isd,
              linetype = 0,
              aes_string(fill = vars[i]), 
              show.legend = FALSE) + 
      scale_fill_viridis_c() +
      theme_minimal()
  plot_list[[i]] = p 
}

# Printing plots to file

for (i in 1:length(plot_list)) {
  title <- vars[i]
  filename <- paste(title,".pdf",sep="") 
  pdf(filename, width=15, height=15/1.618)
      print(plot_list[[i]])
  dev.off()
}

# Creating logged plots | currently not working

for(i in 1:length(vars)) {
  ivar = as.name(vars[i])
  var = paste0("prio_grid_isd$",ivar)
  p = ggplot() +
    geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes_string(fill = log(var+1)), 
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()
  plot_list[[i]] = p 
}

# Printing logged plots to file

for (i in 1:length(plot_list)) {
  title <- vars[i]
  filename <- paste("ln_",title,".pdf",sep="") 
  pdf(filename, width=15, height=15/1.618)
  print(plot_list[[i]])
  dev.off()
}


#==============================================================================#
#	Loading coast data						       #
#==============================================================================#

coastline <- read_sf('../Data/GSHHS_shp/l/GSHHS_l_L1.shp') 

coastline <- sf::st_boundary(coastline)

#==============================================================================#
#	Calculating dist to coast and merging data			       #
#==============================================================================#

prio_grid_isd$distcoast <- get_closest_distance(prio_grid_isd, coastline)

# Cleaning
rm(coastline)

#==============================================================================#
#	Plotting to see it looks right					       #
#==============================================================================#

ggplot() +
    geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes_string(fill = "distcoast"),
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()

ggplot() +
	geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes_string(fill = "non_state"),
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()

#==============================================================================#
#	Rasterdata on pop-density					       #
#==============================================================================#

# Loading data
popdr <- raster("../Data/popd_1600AD.asc")

# Aggregating to Priogrid resolution and extent 
popdr  <- raster_to_pg(popdr, aggregation_function = "mean")

# Converting to tibble with PG-id
popdr  <- raster_to_tibble(popdr, add_pg_index = TRUE)

# Tidying -- not working
popdr  <- popdr %>% rename(popd = layer, gid = pgid) %>% dplyr::select(popd, gid)

# Merging
prio_grid_isd <- left_join(prio_grid_isd, popdr, by = c("gid"))

# Cleaning
rm(popdr)

#==============================================================================#
#	Afrobarometer							       #
#==============================================================================#

# Loading afrobarometer data
afroba  <- read.csv("../Data/afb_full_r3.csv")

afroba  <- st_as_sf(afroba, coords = c("longitude", "latitude"))

#afroba <- st_contains(afroba, prio_grid_shp)

vector_to_pg(afro, mean, need_aggregation = FALSE, missval = -99)

#==============================================================================#
#	UCDP Georeferenced Event Dataset (GED)				       #
#==============================================================================#

# Loading data
load("../Data/ged201.RData")

# Summarising
ged  <- ged201 %>% group_by(priogrid_gid) %>% 
	summarise(state_based = sum(type_of_violence == 1),
	non_state = sum(type_of_violence == 2),
	one_sided = sum(type_of_violence == 3))

# Merging
prio_grid_isd  <- left_join(prio_grid_isd, ged, by = c("gid" = "priogrid_gid"))

# NA's should be 0 for conflict events
prio_grid_isd  <- prio_grid_isd %>% 
	mutate(state_based = ifelse(is.na(state_based), 0, state_based)) %>% 
	mutate(non_state = ifelse(is.na(non_state), 0, non_state)) %>%  
	mutate(one_sided = ifelse(is.na(one_sided), 0, one_sided))

# Cleaning
rm(ged, ged201)
