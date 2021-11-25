library(foreign)
library(MASS)

groups <- read.dta("../Data/groupdata3.dta")

total <- length(unique(groups$cowgroupid))
total 
# 255 groups included

stategroups <- filter(groups, groups$JHordinal > 2) 
stategroups <- length(unique(stategroups$cowgroupid))
stategroups 
# 45 groups tied to states

# 73-79 groups in total make it into the final analysis
