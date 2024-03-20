# {{{ Loading Packages
library(dplyr)
library(sf)

# }}}

# {{{ Load Data

# Loading master shp file
master <- read_sf("master_shapefile.shp")

# Loading list of rows to be corrected
Bcorrected <- read.csv("2Bcorrected.csv")

# Loading list of corrected rows
corrections <- read.csv("corrections.csv")

# }}}

# {{{ Doing da ting

# Assigning polygon id from row number
master$id <- 1:nrow(master)

# Finding those to be corrected
merged <- merge(Bcorrected, master, all.x = T)

# Keeping the rest
#corrections$lyear <- merged$lyear
#corrections$hyear <- merged$hyear
#corrections$coder <- merged$coder
#corrections$note <- merged$note
#corrections$coastrerror <- merged$coasterror
#corrections$cityerror <- merged$cityerror
#corrections$sourcetype <- merged$sourcetype
#corrections$"Core/Great" <- merged$"Core/Great"
#corrections$core <- merged$core
corrections$id <- merged$id
#corrections$source_typ <- merged$source_typ
#corrections$layer <- merged$layer
#corrections$path <- merged$path
#corrections$geometry <- merged$geometry

# Applying corrections
#merged <- merge(corrections, merged, all.x = T)

# Fixing dumb column name
#colnames(master)[colnames(master) == "Core/Great"] <- "Core_Great"
colnames(merged)[colnames(merged) == "Core/Great"] <- "Core_Great"

# Changing to data frame-object to make next step work
masterdf <- data.frame(master)

# Merging back in with the master
masterdf <- rows_update(masterdf, corrections, by = "id", copy = T)

# }}}

# Testing to see if it worked
test <- filter(masterdf, id == 1136)
test

# {{{ Fixing Bagirmi/Borgu

masterdf$COWID[masterdf$COWNUM == 4831] <- "BAG"

# }}}

# {{{ Deleting misstakes

masterdf[masterdf$del
