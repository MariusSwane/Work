# Loading packages
library(dplyr)
library(lubridate)
library(readxl)

# Loading Data
isd <- read_excel("../Data/isd_2_1_1_7_2021.xlsx")

isd$Start_Year <- year(as.Date(isd$Start, '%d-%m-%Y'))
isd$End_Year <- year(as.Date(isd$End, '%d-%m-%Y'))

for (i in 1:length(minisd$StateName)) {
	minisd <- mutate(minisd, year = Start_Year[[i]]:End_Year[[i]])
}

minisd <- groupmutate(minisd, year = Start_Year[[i]]:minisd$End_Year[[i]])
